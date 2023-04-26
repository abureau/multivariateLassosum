plot.validate.lassosum <- function(obj, legend.x="bottomright", log.scale=TRUE, base.plot = FALSE) {
  #' @title Plot function for \code{validate.lassosum} objects
  #' @param obj A validate.lassosum object
  #' @param legend.x A character specifying where to place the legend when using R base plot
  #' @param log.scale Should a log x scale should be used?
  #' @param base.plot Should we force the use of a R base plot? If \code{FALSE}, ggplot2 is used.
  #' @export
  t <- obj$validation.table
  if(!any(is.finite(t$value))) return(invisible(NULL)) # If nothing to plot
  installed <- installed.packages()[,1]
  ggplot2.installed <- ("ggplot2" %in% installed)
  if(ggplot2.installed & !base.plot){
    ggplot(data=t, aes(x=lambda, y=value, group=factor(s), color=factor(s))) +
      geom_line(size = 1)+
      geom_point(size = 1.5) +
      {if(log.scale) scale_x_continuous(trans="log10")}+
      labs(x = "Lambda", y = obj$validation.type, color = "s")
  }else{
    if(!ggplot2.installed & !base.plot)cat("Please install ggplot2...\n")
    plot(unlist(t$lambda), 
         unlist(t$value), 
         type="n", log=ifelse(log.scale, yes = "x", no = ""),
         xlab="lambda", ylab=obj$validation.type)
    ss <- unique(t$s)
    for(i in 1:length(ss)) {
      subs <- subset(t, s == ss[i])
      points(subs$lambda, subs$value, col=i, type="o")
    }
    legend(x=legend.x, col=1:length(ss), pch=1, 
           legend=paste0("s=", ss))
  }
}