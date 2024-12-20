# Import packages

#library(dtwMultiAlign)
library(dtw)
library(DescTools)
library(astrochron)

# Import Picard1 and U1464 datasets

Picard1 <- read.csv("bad/PICARD_1.csv", header=TRUE, stringsAsFactors=FALSE)
Picard1=Picard1[c(83:7950),]# Eocene-Miocene Unconformity
head(Picard1)
plot(Picard1[,c(1,2)], type="l", xlim = c(150, 1400), ylim = c(0, 50))

U1464 <- read.csv("bad/U1464-HSGR.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
M1 = Gmean(U1464[c(1:553),2])
M2 = Gmean(U1464[c(553:600),2])
SD1 = Gsd(U1464[c(1:553),2])
SD2 = Gsd(U1464[c(553:600),2])
U1464[c(1:553),2]=(U1464[c(1:553),2]+(M2-M1))*(SD1/SD2)

head(U1464)
plot(U1464[,c(1,2)], type="l", xlim = c(0, 800), ylim = c(0, 60))

#### Resampling of the data ####

# Normalization of the data

# Picard1

Picard1_interpolated <- linterp(Picard1[,c(1,2)], dt = 0.2, genplot = F)
Picard1_interpolated1 <- linterp(Picard1[,c(1,3)], dt = 0.2, genplot = F)
Picard1_interpolated[,3] = Picard1_interpolated1[,2]

Pmean = Gmean(Picard1_interpolated$GR)
Pstd = Gsd(Picard1_interpolated$GR)

Picard1_scaled = (Picard1_interpolated$GR - Pmean)/Pstd
Picard1_rescaled = data.frame(Picard1_interpolated$DEPT, Picard1_scaled, Picard1_interpolated$V3)

Picard1_scaled = mwStats(Picard1_rescaled, cols = 2, win=3, conv = 1, ends = T)
Picard1_standardized = data.frame(Picard1_scaled$Center_win, Picard1_scaled$Average, Picard1_interpolated$V3)

plot(Picard1_standardized[,c(1,2)], type="l", xlim = c(150, 1400), ylim = c(-20, 20), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")

# U1464

U1464_interpolated <- linterp(U1464[,c(1,2)], dt = 0.2, genplot = F)
U1464_interpolated1 <- linterp(U1464[,c(1,3)], dt = 0.2, genplot = F)
U1464_interpolated[,3] = U1464_interpolated1[,2]

Mmean = Gmean(U1464_interpolated$HSGR)
Mstd = Gsd(U1464_interpolated$HSGR)

U1464_scaled = (U1464_interpolated$HSGR - Mmean)/Mstd
U1464_rescaled = data.frame(U1464_interpolated$DEPTH_WMSF, U1464_scaled, U1464_interpolated$V3)

U1464_scaled = mwStats(U1464_rescaled, cols = 2, win=3, conv = 1, ends = T)
U1464_standardized = data.frame(U1464_scaled$Center_win, U1464_scaled$Average, U1464_interpolated$V3)

plot(U1464_standardized[,c(1,2)], type="l", xlim = c(0, 800), ylim = c(-20, 20), xlab = "U1464 Resampled Depth", ylab = "Normalized GR")

# DTW with asymmetricP01 and open-close ends

system.time(al_U1464_p1_ap1 <- dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
#plot(al_U1464_p1_ap1, "density")

# Tuning the standardized data on reference depth scale
U1464_on_Picard1_depth = tune(U1464_standardized[,c(1,2)], cbind(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap1$index2s]), extrapolate = F)

#dev.off()

# Plotting the data

plot(Picard1_standardized[,c(1,2)], type = "l", ylim = c(-20, 20), xlim = c(150, 1400), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR")
lines(U1464_on_Picard1_depth, col = "red")

# DTW Distance and RMSE

al_U1464_p1_ap1$normalizedDistance
al_U1464_p1_ap1$distance


#### Custom Windowing with Slack ####

window_custom <- function (iw, jw, window_iw, window_jw,query.size,reference.size,...)
{
  lean_window = outer(window_iw, window_jw,FUN = "==")
  slack_window = lean_window
  slack = 500
  
  # TONI this looks like a 2d convolution with a square kernel to me
  for (i in 1:length(window_iw)){
    for (j in 1:length(window_jw)){
      if (lean_window[i,j]==TRUE) {
        if(i+slack <= length(window_iw)) {slack_window[(i+1):(i+slack),j]=TRUE}
        if(i-slack >= 1) {slack_window[(i-slack):(i-1),j]=TRUE}
        if(j+slack <= length(window_jw)) {slack_window[i,(j+1):(j+slack)]=TRUE}
        if(j-slack >= 1) {slack_window[i,(j-slack):(j-1)]=TRUE}
      }
    }
  }
  return(slack_window)
}

sw_matrix <- window_custom(NA, NA, U1464_standardized$U1464_interpolated.V3, Picard1_standardized$Picard1_interpolated.V3, NA, NA)
str(sw_matrix)

image(1:nrow(sw_matrix), 1:ncol(sw_matrix), sw_matrix, useRaster = T, asp=1)

sw_fun <- function(iw, jw, query.size, reference.size, ...) sw_matrix

al_U1464_p1_ap1 <- dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1,
                   window.type = sw_fun, open.end = T, open.begin = T)


#### Custom-Windowing ####

compare.window <- matrix(data=TRUE,nrow=nrow(U1464_standardized),ncol=nrow(Picard1_standardized))
image(x=Picard1_standardized[,1],y=U1464_standardized[,1],z=t(compare.window),useRaster=TRUE)

#install.packages("DescTools")
library(DescTools)

base_1_x <- Closest(260, Picard1_standardized[,1],which=TRUE)
base_1_y <- Closest(46, U1464_standardized[,1],which=TRUE)

base_2_x <- Closest(370, Picard1_standardized[,1],which=TRUE)
base_2_y <- Closest(74, U1464_standardized[,1],which=TRUE)

base_3_x <- Closest(430, Picard1_standardized[,1],which=TRUE)
base_3_y <- Closest(163, U1464_standardized[,1],which=TRUE)

base_4_x <- Closest(545, Picard1_standardized[,1],which=TRUE)
base_4_y <- Closest(242, U1464_standardized[,1],which=TRUE)

base_5_x <- Closest(1010, Picard1_standardized[,1],which=TRUE)
base_5_y <- Closest(570, U1464_standardized[,1],which=TRUE)

base_6_x <- Closest(1190, Picard1_standardized[,1],which=TRUE)
base_6_y <- Closest(700, U1464_standardized[,1],which=TRUE)

compare.window <- matrix(data = TRUE, nrow = nrow(U1464_standardized), ncol = nrow(Picard1_standardized))

compare.window[(base_1_y+300):nrow(U1464_standardized),1:(base_1_x-300)] <- 0
compare.window[1:(base_1_y-200),(base_1_x+200):ncol(compare.window)] <- 0

compare.window[(base_2_y+400):nrow(U1464_standardized),1:(base_2_x-400)] <- 0
compare.window[1:(base_2_y-200),(base_2_x+200):ncol(compare.window)] <- 0

compare.window[(base_3_y+800):nrow(U1464_standardized),1:(base_3_x-800)] <- 0
compare.window[1:(base_3_y-800),(base_3_x+800):ncol(compare.window)] <- 0

compare.window[(base_4_y+800):nrow(U1464_standardized),1:(base_4_x-800)] <- 0
compare.window[1:(base_4_y-800),(base_4_x+800):ncol(compare.window)] <- 0

compare.window[(base_5_y+800):nrow(U1464_standardized),1:(base_5_x-800)] <- 0
compare.window[1:(base_5_y-800),(base_5_x+800):ncol(compare.window)] <- 0

compare.window[(base_6_y+200):nrow(U1464_standardized),1:(base_6_x-200)] <- 0
compare.window[1:(base_6_y-200),(base_6_x+200):ncol(compare.window)] <- 0

image(x=Picard1_standardized[,1],y=U1464_standardized[,1],z=t(compare.window),useRaster=TRUE)

compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))
image(x=Picard1_standardized[,1],y=U1464_standardized[,1],z=t(compare.window),useRaster=TRUE)

win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

system.time(al_U1464_p1_ap1 <- dtw(U1464_standardized$U1464_scaled.Average, Picard1_standardized$Picard1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = T, open.begin = T))

al_U1464_p1_ap1$normalizedDistance
al_U1464_p1_ap1$distance

#plot(al_U1464_p1_ap1, type = "density")

image(y = Picard1_standardized[,1], x = U1464_standardized[,1], z = compare.window, useRaster = T)
lines(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap1$index1], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap1$index2], col = "white", lwd = 2)

U1464_on_Picard1_depth = tune(U1464_standardized[,c(1,2)], cbind(U1464_standardized$U1464_scaled.Center_win[al_U1464_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_U1464_p1_ap1$index2s]), extrapolate = F)

#dev.off()

plot(Picard1_standardized[,c(1,2)], type = "l", ylim = c(-20, 20), xlim = c(150, 1400), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (U1464)")
lines(U1464_on_Picard1_depth, col = "red")

plot(U1464_on_Picard1_depth, type = "l", ylim = c(-20, 20), xlim = c(150, 1400), xlab = "Picard1 Resampled Depth", ylab = "U1464")
