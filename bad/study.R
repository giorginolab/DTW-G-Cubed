load("bad/sw_matrix.RData")


image(1:nrow(sw_matrix), 1:ncol(sw_matrix), sw_matrix, useRaster = T, asp=1)

c<-"grey"
for(i in seq(0,5000,by=200)) {
  abline(a=i, b=0.5, col=c)
  abline(a=i, b=2, col=c)
  abline(a=-i, b=0.5, col=c)
  abline(a=-i, b=2, col=c)
}


load("bad/gcm.RData")
cm <- gcm$costMatrix
image(1:nrow(sw_matrix), 1:ncol(sw_matrix), sw_matrix, useRaster = T, asp=1)

