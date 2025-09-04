library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(arrow)
library(gstat)
library(sp)
library(raster)
library(MASS)
library(geoR)
library(mapview)
library(gridExtra)

#set larger font:
par(cex.lab = 1.5,    # axis labels
    cex.axis = 1.3,   # tick labels
    cex.main = 1.6)   #title

#plot semivariance:
library(sp)
library(sf)
library(lattice)
data(meuse)
meuse.df <- meuse
meuse.sf <- st_as_sf(meuse.df, 
                     coords = c("x", "y"),     
                     crs = 28992  # CRS code
)
ccol = 'black' #grey(.5)
v <- variogram(log(zinc) ~ 1, meuse.sf)
v.fit <- fit.variogram(v, vgm(1, "Sph", 800, 1))
plot(v, v.fit, pch = 16, panel = function(x,y,subscripts,...) {
  larrows(0,v.fit$psill[1], v.fit$range[2], v.fit$psill[1], 
          col=ccol, ends = 'both', length=.1, angle=15)
  larrows(v.fit$range[2],0, v.fit$range[2], v.fit$psill[1], 
          col=ccol, ends = 'both', length=.1, angle=15)
  larrows(v.fit$range[2],v.fit$psill[1], v.fit$range[2], 
          sum(v.fit$psill), 
          col=ccol, ends = 'both', length=.1, angle=15)
  ltext(v.fit$rang[2]/2, 1.2*v.fit$psill[1], "range", col=ccol,
        adj = c(.5, 0), cex=.9)
  ltext(1.02 * v.fit$rang[2], 0.5 *v.fit$psill[1], "nugget", col=ccol,
        adj = c(0, 0.5), cex=.9)
  ltext(1.02 * v.fit$rang[2], v.fit$psill[1] + 0.5 * v.fit$psill[2], 
        "partial sill", col=ccol, adj = c(0, 0.5), cex=.9)
  vgm.panel.xyplot(x,y,subscripts,...)
}
)


covmodel <- "matern"
sigma2 <- 1

# Plot covariance function
curve(cov.spatial(x, cov.pars = c(sigma2, 1),
                  cov.model = covmodel), lty = 1,
      from = 0, to = 1, ylim = c(0, 1), main = "Matérn",
      xlab = "distance", ylab = "Covariance(distance)")
curve(cov.spatial(x, cov.pars = c(sigma2, 0.5),
                  cov.model = covmodel), lty = 2, add = TRUE)
curve(cov.spatial(x, cov.pars = c(sigma2, 0.2),
                  cov.model = covmodel), lty = 3, add = TRUE)
curve(cov.spatial(x, cov.pars = c(sigma2, 0.01),
                  cov.model = covmodel), lty = 4, add = TRUE)

legend(cex = 1.5, "topright", lty = c(1, 1, 2, 3),
       col = c("white", "black", "black", "black"),
       lwd = 2, bty = "n", inset = .01,
       c(expression(paste(sigma^2, " = 1 ")),
         expression(paste(phi, " = 1")),
         expression(paste(phi, " = 0.5")),
         expression(paste(phi, " = 0.2")),
         expression(paste(phi, " = 0.01"))))

# Simulate Gaussian random field in a regular grid 32 X 32
sim1 <- grf(1024, grid = "reg",
            cov.model = covmodel, cov.pars = c(sigma2, 1))
sim2 <- grf(1024, grid = "reg",
            cov.model = covmodel, cov.pars = c(sigma2, 0.2))
sim3 <- grf(1024, grid = "reg",
            cov.model = covmodel, cov.pars = c(sigma2, 0.01))

sim1_sp <- data.frame(x = sim1$coords[,1],
                      y = sim1$coords[,2],
                      value = sim1$data)
coordinates(sim1_sp) <- ~x + y

sim2_sp <- data.frame(x = sim2$coords[,1],
                      y = sim2$coords[,2],
                      value = sim2$data)
coordinates(sim2_sp) <- ~x + y

sim3_sp <- data.frame(x = sim3$coords[,1],
                      y = sim3$coords[,2],
                      value = sim3$data)
coordinates(sim3_sp) <- ~x + y
#variograms
v1<-variogram(value~ 1, sim1_sp)
v2<-variogram(value~ 1, sim2_sp)
v3<-variogram(value~ 1, sim3_sp)
library(lattice)


#plot with large font:
old_theme <- trellis.par.get()

trellis.par.set(
  axis.text = list(cex = 1.3),      # tick labels
  par.main.text = list(cex = 1.6),  # title
  par.xlab.text = list(cex = 1.5),  # x-axis label
  par.ylab.text = list(cex = 1.5)   # y-axis label
)
plot(v1,
     main = expression("Empirical variogram with " * phi == 1),
     xlab = "Distance",
     ylab = "Semivariance")
plot(v2,
     main = expression("Empirical variogram with " * phi == 0.2),
     xlab = "Distance",
     ylab = "Semivariance")
plot(v3,
     main = expression("Empirical variogram with " * phi == 0.01),
     xlab = "Distance",
     ylab = "Semivariance")

# Restore the original theme 
trellis.par.set(old_theme)

plot(v1, main = expression("Empirical variogram with " * phi == 1))
plot(v2, main = expression("Empirical variogram with " * phi == 0.2))
plot(v3, main = expression("Empirical variogram with " * phi == 0.01))




#get multiple realisations for variogram cloud
#first for phi=1
sims1 <- vector("list", 30)
for (i in 1:30) {
  sims1[[i]] <- grf(1024, grid = "reg", cov.model = covmodel, cov.pars = c(sigma2, 1))
}
cloud_points1 <- NULL

for (i in 1:30) {
  v <- variog(sims1[[i]], cloud = TRUE, messages = FALSE)
  cloud_points1 <- rbind(cloud_points1, data.frame(dist = v$u, gamma = v$v))
}
par(mfrow = c(1,1))
par(mar = c(5,5,4,2))   #resets size
plot(cloud_points1$dist, cloud_points1$gamma,
     pch = 21, bg = "white", col = "black",
     xlab = "Distance", ylab = "Semivariance",
     main = expression("Variogram cloud for " * phi == 1))

# for phi=0.2
sims2 <- vector("list", 30)
for (i in 1:30) {
  sims2[[i]] <- grf(1024, grid = "reg", cov.model = covmodel, cov.pars = c(sigma2, 0.2))
}
cloud_points2 <- NULL

for (i in 1:30) {
  v <- variog(sims2[[i]], cloud = TRUE, messages = FALSE)
  cloud_points2 <- rbind(cloud_points2, data.frame(dist = v$u, gamma = v$v))
}
par(mfrow = c(1,1))
par(mar = c(5,5,4,2))   #resets size
plot(cloud_points2$dist, cloud_points2$gamma,
     pch = 21, bg = "white", col = "black",
     xlab = "Distance", ylab = "Semivariance",
     main = expression("Variogram cloud for " * phi == 0.2))


# for phi=0.01
sims3 <- vector("list", 30)
for (i in 1:30) {
  sims3[[i]] <- grf(1024, grid = "reg", cov.model = covmodel, cov.pars = c(sigma2, 0.01))
}
cloud_points3 <- NULL

for (i in 1:30) {
  v <- variog(sims3[[i]], cloud = TRUE, messages = FALSE)
  cloud_points3 <- rbind(cloud_points3, data.frame(dist = v$u, gamma = v$v))
}
par(mfrow = c(1,1))
par(mar = c(5,5,4,2))   #resets size
plot(cloud_points3$dist, cloud_points3$gamma,
     pch = 21, bg = "white", col = "black",
     xlab = "Distance", ylab = "Semivariance",
     main = expression("Variogram cloud for " * phi == 0,01))


# Number of bins to keep consistent lags
nbins <- 30

# Initialize matrix to store gamma values (rows = bins, cols = sims)
gamma_mat <- matrix(NA, nrow = nbins, ncol = length(sims3))
lags <- NULL

for (i in 1:length(sims3)) {
  v <- variog(sims3[[i]], max.dist = max(sims3[[i]]$coords[,1]), uvec = NULL, messages = FALSE)
  # Save lags and gamma
  if (is.null(lags)) lags <- v$u[1:nbins]
  gamma_mat[, i] <- v$v[1:nbins]
}

# Average semivariances across simulations
avg_gamma <- rowMeans(gamma_mat, na.rm = TRUE)

# Plot average variogram
plot(lags, avg_gamma, type = "b", pch = 19,
     xlab = "Distance", ylab = "Semivariance",
     main = "Average Empirical Variogram (30 realizations)")




#isotropy vs anisotropy
data(meuse)
coordinates(meuse) <- ~x + y
v.dir <- variogram(log(zinc) ~ 1, meuse, alpha = (0:3) * 45)
v.anis <- vgm(psill = 0.6, model = "Exp", range = 1000, nugget = 0.05,
              anis = c(45, 0.3))
plot(v.dir, v.anis,
     xlab = "Distance (m)",
     ylab = "Semivariance",
     pch = 16,
     col = "darkblue")

