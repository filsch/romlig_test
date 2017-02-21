library(geoR); library(akima); library(fields); library(MASS)
data = read.table("https://www.math.ntnu.no/emner/TMA4250/2017v/Exercise1/topo.dat") 
set.seed(100)
#Setting plot settings and interpolating a colorpalette
par(mar=c(5,5,5,5), mgp=c(1.5,.5,0), par(mfrow=c(1,1))); n_col = 100
signature = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
                               "yellow", "#FF7F00", "red", "#7F0000"))

#Plotting the datavalues at positions
plot(x='',y='', xlim=c(0,315), ylim=c(0,315), main="Datapoints",
     xlab="x coords", ylab="y coords")
text(data$x, data$y, labels=data$z, cex= 0.7)

#Plotting the data, interpolated middlepoints
int_data = interp(data$x,data$y,data$z, extrap = TRUE); 
filled.contour(int_data, xlim=c(0,315), ylim=c(0,315), color.palette = signature,
               nlevels = n_col, plot.axes = {axis(1)
                 axis(2)
                 contour(int_data,add=T,lwd=2,labcex=1)   }, main="Dataset Y",xlab="x coords", ylab="y coords")
par(mfrow=c(1,2))
#Making grid and performing universal kriging
grid <- expand.grid(seq(0,315), seq(0,315))
kriger = krige.conv(coords=cbind(data$x,data$y), data=data$z, locations = grid,
                    krige=krige.control(type.krige= "OK", trend.d = "2nd", trend.l = "2nd", cov.model="exponential", cov.pars=c(2500, 100)));

#Plotting predictions and std.deviations
plot_values = list(z=matrix(kriger$predict,nrow=316,ncol=316,byrow=FALSE), x = seq(0,315), y = seq(0,315))
image.plot(plot_values, main="Kriging predictions, no noise on data Y", zlim=c(670,980))

noisefree_z_pred = plot_values$z

plot_values$z = matrix(sqrt(kriger$krige.var),nrow=316,ncol=316,byrow=FALSE) 
image.plot(plot_values, main="Kriging std.dev., no noise on data Y", zlim=c(0,45))

noisefree_z_std = plot_values$z
#Kommentar: We see that the kriging operation interpolates so that the predictions fit nicely with the contours, with the same points for max and min.
#Also observe that the standard deviation of the predicted data are biggest in the areas devoid of datapoints, and zero in the datapoints itself. This
#corresponds nicely with theoretic results. 

#We have that the kriging predictions for the actual values Y(s) are Gaussian with an expected value and variance estimated by kriging.
expected_so = kriger$predict[316*100 + 101]; variance_so = kriger$krige.var[316*100+101] #A bit worried about the point
probability_over_700 = 1 - dnorm(700, expected_so, sqrt(variance_so))
below_90_percent = qnorm(0.9, expected_so, sqrt(variance_so))

#New kriging predictions for noisy data with plots, sigma^2 = 5
kriger = krige.conv(coords=cbind(data$x,data$y), data=data$z + mvrnorm(52,0,5), locations = grid,
                    krige=krige.control(type.krige= "OK", trend.d = "2nd", trend.l = "2nd", cov.model="exponential", cov.pars=c(2500, 100)));

plot_values$z = matrix(kriger$predict,nrow=316,ncol=316,byrow=FALSE) 
image.plot(plot_values, main="Kriging predictions, noise sigma^2 = 5 on data Y", zlim=c(670,980))
image.plot(x = seq(0,315),y = seq(0,315),z = plot_values$z - noisefree_z_pred, main="Kriging pred. difference no noise vs noise sigma^2 = 5", zlim=c(-7,7), xlab="x coords",ylab="y coords")

#New kriging predictions for noisy data with plots, sigma^2 = 15 with comparison
kriger = krige.conv(coords=cbind(data$x,data$y), data=data$z + mvrnorm(52,0,5), locations = grid,
                    krige=krige.control(type.krige= "OK", trend.d = "2nd", trend.l = "2nd", cov.model="exponential", cov.pars=c(2500, 100)));

plot_values$z = matrix(kriger$predict,nrow=316,ncol=316,byrow=FALSE) 
image.plot(plot_values, main="Kriging predictions, noise sigma^2 = 15 on data Y", zlim=c(670,980))
image.plot(x = seq(0,315),y = seq(0,315),z = plot_values$z - noisefree_z_pred, main="Kriging pred. difference no noise vs noise sigma^2 = 15", zlim=c(-7,7), xlab="x coords",ylab="y coords")


#Kommentar: No difference in variance as variance is independent of measurements
