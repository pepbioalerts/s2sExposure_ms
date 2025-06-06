# make map
make_map <- function (data,col_cells='cells',col_attr='exposure',layer_name=NULL,
                      rmodel=terra::rast(xmin = -180,xmax=180,ymin=-90,ymax=90,res=0.5)){
  #require(terra)
  #require (dplyr)
  rmodel [data |> dplyr::pull(get(col_cells))] <- data |> dplyr::pull(get(col_attr))
  if (is.null(layer_name)) names(rmodel) <- col_attr
  if (!is.null(layer_name)) names(rmodel) <- layer_name
  return(rmodel)
}

#pheno stats
duration_impact = function (x,value){
  x <-unlist(x)
  seq_exp_duration <-rle (abs(x))
  duration_tib <- tibble (value=seq_exp_duration$values,lengths= seq_exp_duration$lengths)
  if (all(x==0))  {maxdur_exp <- NA} else {
    maxdur_exp <-   duration_tib |> filter (value==1) |> pull(lengths) |> max(na.rm=T)
  }
  if(value=='maxdur_exp') return (maxdur_exp)
  if (all(x==0)) {
    time_to_impact <- NA
  } else {
    rid_maxtarget <- duration_tib |> mutate (target = value==1 & lengths==maxdur_exp ) |> pull(target) |> which()
    time_to_impact <- if (any(rid_maxtarget==1)) 0 else sum (duration_tib$lengths[1:(min(rid_maxtarget)-1)])
  }
  if (value=='time_to_impact') return (time_to_impact)
}
first_exp <- function (x) {
  x <- unlist(x)
  id <- which (abs(x)==1)
  if (length(id)==0) NA else
    min(id)
}

