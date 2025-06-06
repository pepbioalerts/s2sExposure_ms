#SAMPLE WORKFLOW FOR A SPECIES

# load packages if missing and specific functions ====
load_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

load_if_missing ('terra')
load_if_missing ('tidyverse')
load_if_missing ('cNORM')
load_if_missing ('qs')
load_if_missing ('tools')
df_to_terra_stack <- function(df, nrow = NULL, ncol = NULL, ext = NULL, crs = NULL, model = NULL) {
  # Load terra
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Please install the 'terra' package.")
  }
  
  # If model is provided, extract raster parameters from it
  if (!is.null(model)) {
    nrow <- terra::nrow(model)
    ncol <- terra::ncol(model)
    ext  <- terra::ext(model)
    crs  <- terra::crs(model)
  } else {
    if (is.null(nrow) || is.null(ncol)) {
      stop("Either provide a model raster or specify nrow and ncol.")
    }
  }
  
  # Prepare template raster
  r_template <- terra::rast(nrow = nrow, ncol = ncol)
  if (!is.null(ext)) terra::ext(r_template) <- ext
  if (!is.null(crs)) terra::crs(r_template) <- crs
  
  # Create a list to hold layers
  layers <- vector("list", ncol(df) - 1)
  
  # Loop over columns (except first: cell numbers)
  for (i in 2:ncol(df)) {
    r <- r_template
    vals <- rep(NA, terra::ncell(r))
    vals[df[[1]]] <- df[[i]]
    terra::values(r) <- vals
    names(r) <- names(df)[i]
    layers[[i - 1]] <- r
  }
  # Combine into a SpatRaster (stack)
  terra::rast(layers)
}

#set working folder (CHOSE YOUR OWN FOLDER) -----
base_dir <- './script/'

#set parameters ==== 
vars <- c('T2M')
setwd(base_dir)
rasRef = terra::rast(xmin = -180,xmax=180,ymin=-90,ymax=90,res=0.5)
terra::values(rasRef)<- 1:ncell(rasRef)
climDB_dir = './s2s_forecasts_dB2/'
iniDate <- '198001'
endDate <- '202404'
quants<- c(0.01,0.99)
percs <- c(0,0.01,0.05,0.95,0.99,1)
monthforecast='may'
yearforecast='2024'
method <- c('wquant')

#  get a list of leadtimes of s2s  -----
# data access through 
s2s_climdb <- list.files (climDB_dir,pattern = "^lead.*\\.qs$",full.names = T)

# get species ranges as raster cells ====
# function to convert sp range polygons to cell ids and weights
# sp_range_ras <- terra::rasterize (terra::vect(sp_range),
#                                   y = climrefRas,
#                                   touches=T,cover=T)
# 
# 
# sp_range_cells <- tidyterra::as_tibble(sp_range_ras)
spf <-
  sp_range_cells <- 
  './sp_cells/Genus_species.qs'
spname <- tools::file_path_sans_ext(basename(spf))

#loop  leadtime ==== 
sp_forecast <- lapply (s2s_climdb,function (x){
  clim_lead <- qs::qread(x)
  leadtime <- tools::file_path_sans_ext(basename(x))
  print (paste0('====== Starting ',leadtime,' ======'))
  
  #get species cell ids
  spn <- tools::file_path_sans_ext(basename(spf))
  spn_extension <- tools::file_ext(basename(spf))
  sp_cell_info <- switch(spn_extension,
                         rds= readRDS(spf),
                         qs= qs::qread(spf)) 
  sp_cell <- sp_cell_info |> dplyr::pull(cell)
  sp_weights <- rep(1,times=length(sp_cell))
  
  #get species timeseries
  sp_clim <- clim_lead |> dplyr::filter (cell %in% sp_cell) |> 
    dplyr::mutate(year = substr(forecast_date,1,4),
                  month= substr(forecast_date,5,6)) |>  
    dplyr::select(tidyr::all_of(c('cell',vars,'forecast_date','year','month'))) 
  sp_clim_lead <- sp_clim |>  
    dplyr::filter (forecast_date >= iniDate & forecast_date < endDate )
  
  num_month_forecast <- match(monthforecast,tolower(month.abb)) |> as.character() 
  if (nchar(num_month_forecast)==1) num_month_forecast<- paste0('0',num_month_forecast)
  date_monthforecast <- lubridate::ym(paste0(yearforecast,'-',num_month_forecast))
  
  lt <- strsplit (leadtime,split ='lead')[[1]][2]
  target_leadDate <- date_monthforecast + months(as.numeric(lt))
  target_ym <- format(target_leadDate, "%Y%m") 
  sp_clim_leadForecast <- sp_clim |>  dplyr::filter (forecast_date == target_ym)
  
  
  #load variable characteristics and compute species thresholds
  var_char <- tibble::tibble (variable=c('T2M','PRECTOT'),year_stat=c('mean','sum'))
  sp_th_list  <- lapply (vars, function (vv){
    #get species timeseries YEARLY
    agg_stat <- var_char |> 
      dplyr::filter(variable==vv) |> 
      dplyr::pull(year_stat)
    if (length(agg_stat)==0) stop ('variable not considered')
    sp_clim_lead_yearly_vv <- sp_clim_lead |> 
      dplyr::filter (year< as.numeric(yearforecast)) |> #we xclude the year of the forecast to avoid partial years
      dplyr::group_by (cell,year) |>
      dplyr::summarize_at(vv,agg_stat,na.rm=T)
    sp_clim_lead_yearly_vv_ts <- tidyr::pivot_wider(sp_clim_lead_yearly_vv,
                                                    names_from='year',
                                                    values_from=all_of(vv))
    sp_m_ts_yr <- sp_clim_lead_yearly_vv_ts |> 
      dplyr::ungroup(cell)|> 
      dplyr::select (-cell) |> 
      as.matrix()
    
    combis <- expand.grid (quants,percs)
    names(combis) <- c('quants','percs')
    combis$value<-  sapply  (1:nrow(combis),function (i){
      p <-combis[i,'percs']
      q <-combis[i,'quants']
      th <- apply (sp_m_ts_yr,2,function (y){
        nonNA_cells <- !is.na(y)
        y <- y[nonNA_cells]
        w <- sp_weights[nonNA_cells]
        b = cNORM::weighted.quantile(x=y,weights = w,probs=q,type = 'Type7')  
      }) |> quantile(probs=p,na.rm=T)
    }) |> as.numeric()
    combis$timescale <- '00'
    #get species thresholds monthly
    months <- c('01','02','03','04','05','06','07','08','09','10','11','12')
    monthly_thresholds <- lapply (months, function (mo){
      sp_clim_lead_mo_vv <- sp_clim_lead |> 
        dplyr::filter (month==mo) |> 
        dplyr::select(all_of(c('cell','year',vv))) |> 
        dplyr::group_by (cell,year) 
      sp_clim_lead_mo_vv_ts <- tidyr::pivot_wider(sp_clim_lead_mo_vv ,
                                                  names_from='year',
                                                  values_from=all_of(vv))
      sp_m_ts_mo <- sp_clim_lead_mo_vv_ts |> 
        dplyr::ungroup(cell)|> 
        dplyr::select (-cell) |> 
        as.matrix()
      combis <- expand.grid (quants,percs)
      names(combis) <- c('quants','percs')
      combis$value<-  sapply  (1:nrow(combis),function (i){
        p <-combis[i,'percs']
        q <-combis[i,'quants']
        #quantile (matrixStats::rowQuantiles(sp_m_ts_mo,probs=q,na.rm=T),probs=p,na.rm=T)
        th <- apply (sp_m_ts_mo,2,function (y){
          nonNA_cells <- !is.na(y)
          y <- y[nonNA_cells]
          w <- sp_weights[nonNA_cells]
          b = cNORM::weighted.quantile(x=y,weights = w,probs=q,type = 'Type7')  
        }) |> quantile(probs=p,na.rm=T)
      }) |> as.numeric()
      combis$timescale <- mo
      combis
    }) |> dplyr::bind_rows()
    sp_thresholds <- dplyr::bind_rows(combis,monthly_thresholds)
    sp_thresholds$variable <- vv
    sp_thresholds$leadtime <- tools::file_path_sans_ext(basename(x))
    #get species MARGINAL CELLS
    marginal_cells <- lapply (percs, function (p){
      if (p>0.5) upper <- T else upper<-F
      if (upper) {
        th_val <- combis |> 
          dplyr::filter (quants==1,percs==p) |> 
          dplyr::pull(value)
        mc <-sp_clim_lead_yearly_vv_ts$cell [which (matrixStats::rowQuantiles(sp_m_ts_yr,probs=1,na.rm=T)  >= th_val)]
      } 
      if (!upper){
        th_val <- combis |> dplyr::filter (quants==0,percs==p) |> dplyr::pull(value)
        mc <-sp_clim_lead_yearly_vv_ts$cell [which (matrixStats::rowQuantiles(sp_m_ts_yr,probs=0,na.rm=T)  <= th_val)]
      } 
      return (mc)
    })
    names(marginal_cells) <- paste0(vv,'__q1-p',percs)
    #list results output
    return (list (sp_thresholds,marginal_cells))
  })
  
  sp_thresholds <- lapply(sp_th_list, `[[`, 1) |> 
    dplyr::bind_rows()
  
  #forecast 
  mm = lubridate::month (target_leadDate) |> as.character()
  if (nchar(mm)==1) mm<- paste0('0',mm)
  sp_forecast_monthExposure <- lapply (vars,function (vv){
    vals <- sp_clim_leadForecast |> dplyr::pull(c(vv)) 
    #monthly exposure
    target_th_month <- sp_thresholds |> 
      dplyr::filter (variable==vv & timescale==mm & leadtime==leadtime) |> 
      mutate (reasonable = ifelse (quants>0.5 & percs>0.5 |quants<0.5 & percs<0.5 ,TRUE,FALSE)) |> 
      dplyr::filter (reasonable==TRUE) 
    exposure_forecast_month <- lapply (1:nrow (target_th_month), function (rid){
      qq <- target_th_month$quants[rid]
      pp <- target_th_month$percs[rid]
      thval <- target_th_month$value[rid]
      exposure <- if (qq <0.5) as.numeric (vals<thval) else as.numeric (vals>=thval )
      out <- tibble::tibble (exposure)
      names (out)<-paste0(vv,'_month_q',qq,'-p',pp,'_',target_ym)
      out
    }) |> dplyr::bind_cols()
     #monthly exposure on most extreme threshold
    target_th_xtrmonths <- sp_thresholds |> 
      dplyr::filter (variable==vv & timescale!='00' & leadtime==leadtime) |>
      mutate (reasonable = ifelse (quants>0.5 & percs>0.5 |quants<0.5 & percs<0.5 ,TRUE,FALSE)) |> 
      dplyr::filter (reasonable==TRUE) |> 
      dplyr::group_by (quants,percs) |> 
      dplyr::summarize(min_value = min(value,na.rm=T),
                       max_value =  max (value,na.rm=T))
    exposure_forecast_xtrmonth <- lapply (1:nrow (target_th_xtrmonths), function (rid){
      qq <- target_th_xtrmonths$quants[rid]
      pp <- target_th_xtrmonths$percs[rid]
      #thval <- target_th_month$value[rid]
      thval <- ifelse (qq <0.5, 
                       as.numeric(target_th_xtrmonths$min_value[rid]), 
                       as.numeric (target_th_xtrmonths$max_value[rid]))
      exposure <- if (qq <0.5) as.numeric (vals<thval) else as.numeric (vals>=thval )
      out <- tibble::tibble (exposure)
      names (out)<-paste0(vv,'_xtrmonth_q',qq,'-p',pp,'_',target_ym)
      out
    }) |> dplyr::bind_cols()
    #merge all forecasts
    exposure_forecast_full<-dplyr::bind_cols(exposure_forecast_month,exposure_forecast_xtrmonth)
    exposure_forecast_full
  }) |> dplyr::bind_cols()
  sp_forecast_monthExposure$cells <- sp_clim_leadForecast$cell
  sp_forecast_monthExposure<- sp_forecast_monthExposure |> dplyr::relocate(cells)
  
  #create raster output
  outRaster <- df_to_terra_stack(sp_forecast_monthExposure,model = rasRef) 
  
  #return results
  return (outRaster)
}) |> terra::rast()

# save forecasts ====
forecast_dir <- paste0(base_dir,'/forecast/',monthforecast,'_',yearforecast)
dir.create(forecast_dir,showWarnings = F,recursive = T)
terra::writeRaster (sp_forecast, paste0(forecast_dir,'/',spname,'.tif'),overwrite=T)
