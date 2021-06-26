df_to_layer <- function(x, lons, lats){
    mat <- matrix(data = x, nrow = uniqueN(lats), ncol = uniqueN(lons))
    layer <- raster(
        mat[nrow(mat):1,],
        xmn=min(lons), xmx=max(lons), 
        ymn=min(lats), ymx=max(lats)
     )
    return(layer)
}

summary_inner <- function(object){
    # Fit
    if(class(object)=='bart') {
        fitobj <- object$fit
    } else if(class(object)=='rbart') {
        fitobj <- object$fit[[1]] 
    }

    # Call
    (objcall <- object$call)

    # Predictors
    (predictors <- paste(attr(fitobj$data@x, "term.labels"), sep=' '))
    
    # Labels
    (true.vector <- fitobj$data@y)

    # Predictions
    (predictions <- colMeans(pnorm(object$yhat.train)))
    
    # Prediction instance
    (pred <- prediction(predictions, true.vector))
    # Performance instance
    (perf.tss <- performance(pred, "sens", "spec"))
    # TSS values list
    tss.list <- (perf.tss@x.values[[1]] + perf.tss@y.values[[1]] - 1)
    # TSS values ~ threshold dataframe
    (tss.df <- tibble(alpha=perf.tss@alpha.values[[1]],tss=tss.list))

    # AUC
    (auc <- performance(pred,"auc")@y.values[[1]])
    
    # Threshold
    (thresh <- min(tss.df$alpha[which(tss.df$tss==max(tss.df$tss))]))
    
    # TSS
    (tss <- tss.df[which(tss.df$alpha==thresh),'tss'][[1]])
    
    # Type I error rate (false positive)
    (type_I <-  1 - perf.tss@y.values[[1]][which(perf.tss@alpha.values[[1]] == thresh)])
    # Type II error rate (false negative)
    (type_II <- 1 - perf.tss@x.values[[1]][which(perf.tss@alpha.values[[1]] == thresh)])
    
    diagnostics <- list(
        fit = fitobj,
        call = objcall,
        predictors = predictors,
        labels = true.vector,
        predictions = predictions,
        auc = auc,
        threshold = thresh,
        tss = tss,
        type_I = type_I,
        type_II = type_II
    )
    return(diagnostics)
}