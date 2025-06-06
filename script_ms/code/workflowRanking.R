#### CODE TO CREATE REGIONAL (GID1 AND GID2) RANKS OF EXPOSURE
# LOAD DATA ====
library (ggplot2)
library (sf)
library (terra)
library (tidyverse)
library (tidyterra)
library (qs)
source ('./code/funcs_spExposure.R')
rmodel=terra::rast(xmin = -180,xmax=180,ymin=-90,ymax=90,res=0.5)
values(rmodel)<-1:ncell(rmodel)

# start by region (GID1) analysis----
outputDir <- './data_ms/analysis/regionMonitoring'
dir.create(outputDir)
region_codes <- qread ('./data_ms/attr/regionCodes.qs')
all_regions_output <- qs::qread('./data_ms/regionMonitoring/region_stats.qs')
## rank regions ----
## load exposure by regions ====
reg_exp <- qs::qread('./data_ms/regionMonitoring/gid1_regExpVars_withID.qs')
#load criteria
cb_vec <- readRDS('./data_ms/attr/criteria.rds')
dec_tib <- qs::qread('./data_ms/regionMonitoring/gid1_regExpVars.qs') #as reg_exp without id
varnames_dec_tib <-names(dec_tib)
dec_mat <- dec_tib|> as.matrix()

## build weights ====
#equal weights
eq_w <- rep(1/ncol(dec_mat),times=ncol(dec_mat))
#extreme x2 weight
xtr_w <- rep(0,times=ncol(dec_mat))
xtr_w [cb_vec$varname[-1] |> grep(pattern='extreme_')] <- 2
xtr_w <- xtr_w/sum(xtr_w)
#species x2 weight
N_target_cols = cb_vec$varname[-1] |> grep(pattern='_vulnerable') |> length()
N_nontarget_cols = (cb_vec$varname[-1] |> length()) - N_target_cols
sp_w <- rep(0,times=ncol(dec_mat))
sp_w [cb_vec$varname[-1] |> grep(pattern='_vulnerable')] <- 1/N_target_cols
sp_w <- sp_w/sum(sp_w)
#ecosys x2 weight
N_target_cols = cb_vec$varname[-1] |> grep(pattern='_widespread') |> length()
N_nontarget_cols = (cb_vec$varname[-1] |> length()) - N_target_cols
eco_w <- rep(0,times=ncol(dec_mat))
eco_w [cb_vec$varname[-1] |> grep(pattern='_widespread')] <- 1/N_target_cols
eco_w <- eco_w/sum(eco_w)

## compute ranks ====
library (MCDM)
library (caret)
source('./code/funcs_mcdm.R')
list_weights <- list(eqw=eq_w,spw=sp_w,ecow=eco_w,xtrw=xtr_w)
# add initial NA weight to variable attribute table
cb_vec$eq_w <- c(NA,eq_w)
cb_vec$xtr_w <- c(NA,xtr_w)
cb_vec$sp_w <- c(NA,sp_w)
cb_vec$eco_w <- c(NA,eco_w)

# get entropy weights
ll = lapply (c('eq_w','xtr_w','sp_w','eco_w'),function(wname){
  dimension_varialbes <- cb_vec$varname[which (cb_vec[,wname]>0)]
  dec_tib_subset <- dec_tib|> select(all_of(dimension_varialbes)) |> drop_na() 
  corr <- dec_tib_subset|> cor()
  hc =findCorrelation(corr,cutoff = 0.7)
  hc = sort(hc)
  reduced_Data = dec_tib_subset[,-c(hc)]
  var_direction <- cb_vec |> filter (varname %in% (reduced_Data |> names())) |> pull(cb)
  var_pos <- which(var_direction=='max')
  var_neg <- which(var_direction=='min')
  entWeight <- creditmodel:: entropy_weight(dat = reduced_Data,
                                            pos_vars = var_pos,
                                            neg_vars = var_neg)
  names(entWeight)<-c('varname',paste0('ew_',wname)) 
  other_vars <- tibble (varname = dimension_varialbes[!dimension_varialbes %in% entWeight$varname],
                        ew=0)
  names(other_vars)<-   names(entWeight)
  entWeight <-rbind(entWeight,other_vars)
  return(entWeight)
}) 
ew_all <- full_join(ll[[1]],ll[[2]]) |> full_join(ll[[3]]) |> full_join(ll[[4]])

#perform ranking
ranks_out <- lapply(names(ew_all)[-1],function(wname){
  print (wname)
  ww= ew_all |> pull(wname)
  ww[is.na(ww)]<-0
  topsis= TOPSISv_2(decision=dec_mat,weights =ww,cb = cb_vec$cb[-1] ) |> pull(Ranking)
  vikor = VIKOR_2(decision=dec_mat,weights =ww,cb = cb_vec$cb[-1] ,v=0.5) |> pull(Ranking)
  mmoora = MMOORA_2 (decision=dec_mat,weights =ww,cb = cb_vec$cb[-1]) |> pull(Ranking)
  waspas = WASPAS_2(decision=dec_mat,weights =ww,cb = cb_vec$cb[-1] ,lambda=0.5) |> pull(Ranking)
  df = data.frame (topsis,vikor,mmoora,waspas) 
  df$metaRank <-rowMeans(df,na.rm=T)
  df$metaRank <- dense_rank( df$metaRank )
  names(df) <- paste0(c('topsis','vikor','mmoora','waspas','metaRank'),'_',wname )
  return(df)
}) |> Reduce(f=cbind)
#merge ranks to regions_id
ranks_out <- cbind(reg_exp |> select(region_id),
                   ranks_out)
ranks_out <- left_join (ranks_out, 
                        all_regions_output |> select(region_id,region_gid1) |> distinct())
ranks_out <- ranks_out |> 
  left_join (region_codes,by=c('region_id'='reg_id','region_gid1'='GID_1')) |> 
  relocate(region_id,region_gid1,starts_with('metaRank')) |> rename(GID_1=region_gid1) |> 
  left_join(gadmCodesCells |> select (UID,GID_0,GID_1) |> drop_na() |> distinct())
# mapping of ranks and metarank (e.g. fields starting with metaRank_ew_*) can be done by downloading GADM database and joining

# start by subregion (GID2) analysis----
rm (list=ls())
region_codes <- qread ('./data_ms/attr/subregion_codes.qs')
all_regions_output <- qs::qread('./data_ms/subregionMonitoring/subregion_stats.qs')

## rank regions ----
## load exposure by reagions ====
reg_exp <- qs::qread('./data_ms/subregionMonitoring/gid2_regExpVars_withID.qs')

#load criteria
cb_vec <- readRDS('./data_ms/attr/criteria.rds')
dec_tib <- qs::qread('./data_ms/subregionMonitoring/gid2_regExpVars.qs') #as reg_exp without id
varnames_dec_tib <-names(dec_tib)
dec_mat <- dec_tib|> as.matrix()

## build weights ====
#equal weights
eq_w <- rep(1/ncol(dec_mat),times=ncol(dec_mat))
#extreme x2 weight
xtr_w <- rep(0,times=ncol(dec_mat))
xtr_w [cb_vec$varname[-1] |> grep(pattern='extreme_')] <- 2
xtr_w <- xtr_w/sum(xtr_w)
#species x2 weight
N_target_cols = cb_vec$varname[-1] |> grep(pattern='_vulnerable') |> length()
N_nontarget_cols = (cb_vec$varname[-1] |> length()) - N_target_cols
sp_w <- rep(0,times=ncol(dec_mat))
sp_w [cb_vec$varname[-1] |> grep(pattern='_vulnerable')] <- 1/N_target_cols
sp_w <- sp_w/sum(sp_w)
#ecosys x2 weight
N_target_cols = cb_vec$varname[-1] |> grep(pattern='_widespread') |> length()
N_nontarget_cols = (cb_vec$varname[-1] |> length()) - N_target_cols
eco_w <- rep(0,times=ncol(dec_mat))
eco_w [cb_vec$varname[-1] |> grep(pattern='_widespread')] <- 1/N_target_cols
eco_w <- eco_w/sum(eco_w)

## compute ranks ====
library (MCDM)
library (caret)
source('./code/funcs_mcdm.R')
list_weights <- list(eqw=eq_w,spw=sp_w,ecow=eco_w,xtrw=xtr_w)
# add initial NA weight to variable attribute table
cb_vec$eq_w <- c(NA,eq_w)
cb_vec$xtr_w <- c(NA,xtr_w)
cb_vec$sp_w <- c(NA,sp_w)
cb_vec$eco_w <- c(NA,eco_w)

# get entropy weights
ll = lapply (c('eq_w','xtr_w','sp_w','eco_w'),function(wname){
  dimension_varialbes <- cb_vec$varname[which (cb_vec[,wname]>0)]
  dec_tib_subset <- dec_tib|> select(all_of(dimension_varialbes)) |> drop_na() 
  corr <- dec_tib_subset|> cor()
  hc =findCorrelation(corr,cutoff = 0.7)
  hc = sort(hc)
  reduced_Data = dec_tib_subset[,-c(hc)]
  var_direction <- cb_vec |> filter (varname %in% (reduced_Data |> names())) |> pull(cb)
  var_pos <- which(var_direction=='max')
  var_neg <- which(var_direction=='min')
  entWeight <- creditmodel:: entropy_weight(dat = reduced_Data,
                                            pos_vars = var_pos,
                                            neg_vars = var_neg)
  names(entWeight)<-c('varname',paste0('ew_',wname)) 
  other_vars <- tibble (varname = dimension_varialbes[!dimension_varialbes %in% entWeight$varname],
                        ew=0)
  names(other_vars)<-   names(entWeight)
  entWeight <-rbind(entWeight,other_vars)
  return(entWeight)
}) 
ew_all <- full_join(ll[[1]],ll[[2]]) |> full_join(ll[[3]]) |> full_join(ll[[4]])
#perform ranking
ranks_out <- lapply(names(ew_all)[-1],function(wname){
  print (wname)
  ww= ew_all |> pull(wname)
  ww[is.na(ww)]<-0
  topsis= TOPSISv_2(decision=dec_mat,weights =ww,cb = cb_vec$cb[-1] ) |> pull(Ranking)
  vikor = VIKOR_2(decision=dec_mat,weights =ww,cb = cb_vec$cb[-1] ,v=0.5) |> pull(Ranking)
  mmoora = MMOORA_2 (decision=dec_mat,weights =ww,cb = cb_vec$cb[-1]) |> pull(Ranking)
  waspas = WASPAS_2(decision=dec_mat,weights =ww,cb = cb_vec$cb[-1] ,lambda=0.5) |> pull(Ranking)
  df = data.frame (topsis,vikor,mmoora,waspas) 
  df$metaRank <-rowMeans(df,na.rm=T)
  df$metaRank <- dense_rank( df$metaRank )
  names(df) <- paste0(c('topsis','vikor','mmoora','waspas','metaRank'),'_',wname )
  return(df)
}) |> Reduce(f=cbind)
#merge ranks to regions_id
ranks_out <- cbind(reg_exp |> select(subregion_id),
                   ranks_out)
ranks_out <- left_join (ranks_out, 
                        all_regions_output |> select(subregion_id,subregion_gid2) |> distinct())
gadmCodesCells <- qread('./data_ms/attr/gadmCodesCells.qs')
ranks_out <- ranks_out |> 
  left_join (region_codes,by=c('subregion_id'='reg_id','subregion_gid2'='GID_2')) |> 
  relocate(subregion_id,subregion_gid2,starts_with('metaRank')) |> rename(GID_2=subregion_gid2) |> 
  left_join(gadmCodesCells |> select (UID,GID_0,GID_1,GID_2,HASC_2) |> drop_na() |> distinct())
# mapping of ranks and metarank (e.g. fields starting with metaRank_ew_*) can be done by downloading GADM database and joining
