mybarplot = function(dat, device="X11", OutputFile=NULL) {
  margin = max(nchar(as.character(dat[,1])))
  w = 4*log((max(dat[,2], na.rm=T)-min(dat[,2], na.rm=T))+ 0.2*margin)
  h = 3*log(nrow(dat))
  #cexR = max(0.3,6*hhm/nrow(dataset))
#rowmargin = 1+cexR*0.4*max(nchar(rownames(dataset)))
  switch(device,
         "pdf"=pdf(OutputFile, width=w, height=h),
         "png"=png(OutputFile, width=96*w, height=96*h),
         "jpg"=jpeg(OutputFile, width=96*w, height=96*h),
         X11(width=w, height=h))
  par(las=2)
  par(mar=c(4, 0.6*margin, 2, 2))
  barplot2(dat[,2], names.arg=dat[,1], horiz=T)
  switch(device,
         "pdf"=dev.off(),
         "png"=dev.off(),
         "jpg"=dev.off())

}

mybarplot2 = function(dat, device="X11", OutputFile=NULL, wf=3, hf=4, legendText, legendPos="topright", ...) {
  ### some useful parameters to pass on to barplot2 are
  ### legend.text
  ### inside
  margin = max(nchar(as.character(dat[,1])))
  w = wf*log((max(dat[,-1], na.rm=T)-min(dat[,-1], na.rm=T))+ 0.2*margin)
  h = hf*log(nrow(dat))
  #cexR = max(0.3,6*hhm/nrow(dataset))
#rowmargin = 1+cexR*0.4*max(nchar(rownames(dataset)))
  switch(device,
         "pdf"=pdf(OutputFile, width=w, height=h),
         "png"=png(OutputFile, width=96*w, height=96*h),
         "jpg"=jpeg(OutputFile, width=96*w, height=96*h),
         X11(width=w, height=h))
  par(las=2)
  par(mar=c(4, 0.55*margin, 0, 2))
  barHeights <- apply(dat[,-1], 2, as.numeric)
  barplot2(t(barHeights), names.arg=dat[,1], horiz=T, beside=T,
     col=as.character(1:ncol(barHeights)), ...)
  abline(v=0, col="grey")
  legend(legendPos, legendText, fill=as.character(1:ncol(barHeights)))
  switch(device,
         "pdf"=dev.off(),
         "png"=dev.off(),
         "jpg"=dev.off())

}
