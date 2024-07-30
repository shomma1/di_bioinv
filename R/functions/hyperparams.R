# Script for transforming hyper parameters

trans.hyper <- function (model_info) {

  source("./1_script/functions/transform_hyperparams.R")
  
  df.hyper <- data.frame()
  margs <- list()
  for (i in 1:nrow(model_info)) {
    model <- model_info[i, ]
    res <- res.list[[model$model.name]]
    
    if (model$model != "GLM") {
      if (model$model %in% c("MICAR-dt","MICAR-sp")) {
        res <- transform.hyper.micar(res)
        summary.hyper <- res[["zmargs"]]
        margs[[model$model.name]] <- res[["margs"]]
      } else if (model$model %in% c("MPCAR-dt","MPCAR-sp")) {
        res <- transform.hyper.mlerouxcar(res)
        summary.hyper <- res[["zmargs"]]
        margs[[model$model.name]] <- res[["margs"]]
        
      } else if (model$model == "ICAR") {
        res <- transform.hyper.icar(res)
        summary.hyper <- res[["zmargs"]]
        margs[[model$model.name]] <- res[["margs"]]
        
      } else if (model$model == "PCAR") {
        res <- transform.hyper.lerouxcar(res)
        summary.hyper <- res[["zmargs"]]
        margs[[model$model.name]] <- res[["margs"]]
        
      } else {
        res <- res$summary.hyperpar
        summary.hyper <- res[["zmargs"]]
        margs[[model$model.name]] <- res[["margs"]]
        
      } 
      
      summary.hyper <- cbind(params = rownames(summary.hyper), summary.hyper)
      rownames(summary.hyper) <- NULL
      df.hyper <- rbind(df.hyper, 
                        cbind( model[rep(1, nrow(summary.hyper)), ],
                               summary.hyper))
    } 
  }
  
  return(list(df.hyper = df.hyper, marges = margs))

}
