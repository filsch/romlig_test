library(MASS); library(fields); library(akima); library(geoR); library(spatial)

#Task 2a)
#Setting some initial parameters
lambda_m = 5
lambda_c = 10
sigma = 0.1
wrapping = FALSE
lambda_2x2 = lambda_m*2*2
title = 'Neymann-Scott PP with wrapping'

#Generating number of mother points and their "placement" in D_s
lambda = lambda_m
if (!wrapping){
  lambda = lambda_2x2
  title = 'Neymann-Scott PP with actual sampling over (-0.5, 1.5) x (-0.5, 1.5)'
}
num_mother_points = rpois(1,lambda = lambda);
place_mothers = runif(num_mother_points*2) #matrix(runif(num_mother_points*2), ncol = 2, byrow=TRUE)

if (!wrapping){
  place_mothers = place_mothers*2 - 0.5
}

#Generating actual points by daughter process per mother points,
# expected placement bivariate gaussian with mean mother position
#Note: Be vary of dimensions and meaning of rows. Even = x pos, odd = y pos

num_daughter_points = rpois(num_mother_points,lambda = lambda_c)
total_daughter_points = sum(num_daughter_points)
if (length(num_daughter_points) > 1){
  mu = rep(place_mothers[1:2], times = num_daughter_points[1])
} else {
  mu = rep(place_mothers[1:2], times = num_daughter_points)
}

for (i in 2:length(num_daughter_points)){
  mu = c(mu, rep(place_mothers[(i*2 - 1):(i*2)], times = num_daughter_points[i]))
}

placement = mu + sigma*rnorm(n = total_daughter_points)
plot(x = '', y = '', xlim=c(0,1),ylim=c(0,1), main=title, xlab ="X-pos", ylab ="Y-pos")
#plot(place_mothers[seq(1,length(place_mothers),2)],place_mothers[seq(2,length(place_mothers),2)],
#     col='red', pch=4, xlab ='Y-pos', ylab='X-pos', main=title, xlim=c(0,1), ylim=c(0,1))

if (wrapping){
  #Wrapping the datapoints
  if (sum(placement > 1) > 1){
  placement[placement > 1] = placement[placement > 1]%%1
  }
  if (sum(placement < 0) > 1){
    placement[placement < 0] = placement[placement < 0]%%1
  }
  points(placement[seq(1,length(placement),2)],placement[seq(2,length(placement),2)],
         xlab ='Y-pos', ylab='X-pos', main=title, xlim=c(0,1), ylim=c(0,1))
} else {
  #Removing points outside of D_s
  placement = matrix(placement, ncol=2, nrow=total_daughter_points, byrow=TRUE)
  placement = placement[placement[,1] > 0,]
  placement = placement[placement[,1] < 1,]
  placement = placement[placement[,2] > 0,]
  placement = placement[placement[,2] < 1,]
  
  points(placement[,1],placement[,2], xlab ='Y-pos',
         ylab='X-pos', main=title, xlim=c(0,1), ylim=c(0,1))
}

#Task 2b)



