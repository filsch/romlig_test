library(geoR); library(akima); library(fields)

set.seed(100)
#For task b
pointsOnGRF = 1:50
variances = c(1,2,5,10); expectancy = 20; ranges = c(1,10,25,50,100)
distances = matrix(numeric(50*50),50,50)
for (i in 1:50){
  for (j in i:50){
  distances[i,j] = abs(pointsOnGRF[i]-pointsOnGRF[j])
  }
}
distances = distances + t(distances)
color = c("red","green","blue","orange","yellow")
plot(x = '', y = '', xlim=c(0,50),ylim=c(15,25), main="Exponential")
for (i in 1:length(ranges)){
    correlationExponential = cov.spatial(distances,cov.model="exponential", cov.pars=c(1, 50));
    SigmaE = correlationExponential * variances[4]
    realizationsE = mvrnorm(1,numeric(50)+20,SigmaE)
    lines(1:50,realizationsE, col=color[i])
}
plot(x = '', y = '', xlim=c(0,50),ylim=c(15,25), main="Matern")
for (i in 1:length(ranges)){
  correlationMatern = cov.spatial(distances,cov.model = "matern", cov.pars=c(1, 50), kappa=1);
  SigmaM = correlationMatern * variances[1]
  realizationsM = mvrnorm(1,numeric(50) + 20,SigmaM)
  lines(1:50,realizationsM, col=color[i])
}


  #For task d:
correlationMatern = cov.spatial(distances,cov.model = "matern", cov.pars=c(1, 50), kappa=1);
correlationExponential = cov.spatial(distances,cov.model="exponential", cov.pars=c(1, 50));
SigmaM = correlationMatern * variances[1]
SigmaE = correlationExponential * variances[1]
realizationsM = mvrnorm(1,numeric(50) + 20,SigmaM)
realizationsE = mvrnorm(1,numeric(50) + 20,SigmaE)

obsZ = mvrnorm(1,c(realizationsE[10],realizationsE[25],realizationsE[30]), diag(3))
pointsZ = c(10,25,30); 
distances12 = matrix(numeric(50*3),50,3); distances12[,1] = abs(1:50 - pointsZ[1]); distances12[,2] = abs(1:50 - pointsZ[2]); distances12[,3] = abs(1:50 - pointsZ[3])
distances22 = matrix(numeric(9),3,3); distances22[,1] = abs(pointsZ - pointsZ[1]); distances22[,2] = abs(pointsZ - pointsZ[2]); distances22[,3] = abs(pointsZ - pointsZ[3])

Sigma11 = cov.spatial(distances,cov.model="exponential", cov.pars=c(1, 50));
Sigma12 = cov.spatial(distances12,cov.model="exponential", cov.pars=c(1, 50));
Sigma22 = cov.spatial(distances22,cov.model="exponential", cov.pars=c(1, 50)) + diag(3);

conditional_mu = numeric(50) + 20 - Sigma12%*%solve(Sigma22)%*%(obsZ-20);
conditional_Sigma = Sigma11 - Sigma12%*%solve(Sigma22)%*%t(Sigma12)

polygon_y = cbind(conditional_mu - 2*sqrt(diag(conditional_Sigma)),conditional_mu + 2*sqrt(diag(conditional_Sigma)));
polygon_x = cbind(1:50,50:1)

par(mfrow=c(2,1))
plot(x = '', y = '', xlim=c(0,50),ylim=c(17,23), main="Exponential conditional") 
polygon(polygon_x,polygon_y, col="red", lty="dashed"); lines(1:50,conditional_mu) 


#For task e:

sim = 50
simulations = mvrnorm(sim,conditional_mu, conditional_Sigma)
mean_simulations = colMeans(simulations); var_simulations = colMeans(simulations^2) - mean_simulations^2
polygon_y = cbind(mean_simulations - 2*sqrt(var_simulations), mean_simulations + 2*sqrt(var_simulations));
polygon_x = cbind(1:50,50:1)

plot(x = '', y = '', xlim=c(0,50),ylim=c(17,23), main="Exponential, simulated") 
polygon(polygon_x,polygon_y, col="red", lty="dashed"); lines(1:50, mean_simulations) 

