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
library(RColorBrewer)

#prepare data: read, omit NA's, keep only needed columns
data_January <-read.csv(file="C:\\Users\\Sofija\\Desktop\\Durham\\dissertation\\R\\data_January.csv")
data_Jan_multivar<-data_January[
  , c("latitude", "longitude", "avg_temp_c", "precipitation_mm", "avg_wind_speed_kmh")
]
data_multivar<-na.omit(data_Jan_multivar)
data_multivar <- st_as_sf(data_multivar, coords = c("longitude", "latitude"), crs = 4326)
europe <- ne_countries(continent = "Europe", returnclass = "sf")


#plots avg temp:
ggplot() +
  geom_sf(data = europe, fill = "gray90", color = "white") +
  geom_sf(data = data_multivar, aes(color = avg_temp_c), size = 2) +
  scale_color_viridis_c(option = "plasma", name = "Avg Temp (°C)") +
  coord_sf(xlim = c(-10, 35), ylim = c(35, 70), expand = FALSE) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  labs(
    #title = "Average January Temperature",
    x = "Longitude", y = "Latitude")

#plot precipitation
ggplot() +
  geom_sf(data = europe, fill = "gray90", color = "white") +
  geom_sf(data = data_multivar, aes(color = precipitation_mm), size = 2) +
  scale_color_viridis_c(option = "plasma", name = "Precip (mm)") +
  coord_sf(xlim = c(-10, 35), ylim = c(35, 70), expand = FALSE) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  labs(
    #title = "Average daily precipitation in January",
    x = "Longitude", y = "Latitude")
#plot avg wind speed
ggplot() +
  geom_sf(data = europe, fill = "gray90", color = "white") +
  geom_sf(data = data_multivar, aes(color = avg_wind_speed_kmh), size = 2) +
  scale_color_viridis_c(option = "plasma", name = "wind speed (kmh)") +
  coord_sf(xlim = c(-10, 35), ylim = c(35, 70), expand = FALSE) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  labs(
    #title = "Average wind speed in January",
    x = "Longitude", y = "Latitude")



#look at histograms to inspect if transformation is needed
par(mfrow = c(3, 1))
hist(
  data_multivar$avg_temp_c,
  main = "Average Temperature",
  xlab = "Temperature (°C)",
  cex.main = 1.6,
  cex.lab = 1.6,
  cex.axis = 1.3
)
hist(
  data_multivar$precipitation_mm,
  main = "Precipitation",
  xlab = "Precipitation (mm)",
  cex.main = 1.6,
  cex.lab = 1.6,
  cex.axis = 1.3
)
hist(
  data_multivar$avg_wind_speed_kmh,
  main = "Average Wind Speed",
  xlab = "Speed (km/h)",
  cex.main = 1.6,
  cex.lab = 1.6,
  cex.axis = 1.3
)

par(mfrow = c(1, 1)) 


#apply box-cox transformation to be able to do predictions
transform_columns <- c("avg_temp_c", "precipitation_mm", "avg_wind_speed_kmh")
data_boxcox <- data_multivar
boxcox_params <- list()

for (col in transform_columns) {
  shift <- abs(min(data_boxcox[[col]], na.rm = TRUE)) + 1
  shifted <- data_boxcox[[col]] + shift
  bc <- boxcox(lm(shifted ~ 1), lambda = seq(-2, 2, 0.1))
  lambda <- bc$x[which.max(bc$y)]
  if (abs(lambda) < 1e-6) {
    transformed <- log(shifted)
  } else {
    transformed <- (shifted^lambda - 1) / lambda
  }
  data_boxcox[[col]] <- transformed
  boxcox_params[[col]] <- list(shift = shift, lambda = lambda) #save parameters to unshift later
}
# look at histograms after transformation
par(mfrow = c(3,1))
hist(
  data_boxcox$avg_temp_c,
  main = "Average temperature",
  xlab = "temperature",
  cex.main = 1.6,
  cex.lab = 1.6,
  cex.axis = 1.3
)
hist(data_boxcox$precipitation_mm, 
     main = "Precipitation", 
     xlab = "precip (mm)",
     cex.main = 1.6,
     cex.lab = 1.6,
     cex.axis = 1.3)
hist(data_boxcox$avg_wind_speed_kmh, 
     main = "Avg wind speed", 
     xlab = "speed (kmh)",
     cex.main = 1.6,
     cex.lab = 1.6,
     cex.axis = 1.3)
par(mfrow = c(1, 1))

hist(data_boxcox$avg_temp_c, main = "Average temperature", xlab = "temperature")
hist(data_boxcox$precipitation_mm, main = "Precipitation", xlab = "precip (mm)")
hist(data_boxcox$avg_wind_speed_kmh, main = "Avg wind speed", xlab = "speed (kmh)")
par(mfrow = c(1, 1))


#create grid of europe
europe<-ne_countries(continent = "Europe", returnclass = "sf")
europe <- europe %>%
  filter(!admin %in% c("Russia", "Turkey", "Kazakhstan", "Azerbaijan", "Georgia", "Armenia"))

box<-st_bbox(europe)

r <- raster(xmn = box["xmin"], xmx = box["xmax"],
            ymn = box["ymin"], ymx = box["ymax"],
            res = 0.7)   #gives point density 

grid_data <- rasterToPoints(r, spatial = TRUE)
europe_sp <- as(europe, "Spatial")
grid_europe <- grid_data[europe_sp, ]
grid_europe <- grid_europe[coordinates(grid_europe)[,1] > -10 & coordinates(grid_europe)[,2] > 36]
grid_europe_sf <- st_as_sf(grid_europe)
ggplot() +
      geom_sf(data = europe, fill = "gray90", color = "white") +
     geom_sf(data = grid_europe_sf, color = "blue", size = 0.5, alpha = 0.6) +
     coord_sf(xlim = c(-10, 35), ylim = c(35, 70), expand = FALSE) +
     theme_minimal() +
     labs(x = "Longitude", y = "Latitude")
grid_europe_sf$x <- st_coordinates(grid_europe_sf)[,1]
grid_europe_sf$y <- st_coordinates(grid_europe_sf)[,2]

# direct variogram modeling
data_boxcox$x <- st_coordinates(data_boxcox)[, 1]
data_boxcox$y <- st_coordinates(data_boxcox)[, 2]
trellis.par.set(
  list(
    par.xlab.text = list(cex = 1.6),  
    par.ylab.text = list(cex = 1.6),  
    axis.text = list(cex = 1.3),      
    add.text = list(cex = 1.4)        
  )
)
svg_temp <- variogram(avg_temp_c ~ 1+x+y, data = data_boxcox, cutoff = 7000)
vgm_model_temp <- fit.variogram(svg_temp, model = vgm("Sph"))
plot(svg_temp, model=vgm_model_temp)

svg_prec <- variogram(precipitation_mm ~ 1, data = data_boxcox, cutoff = 10000)
vgm_model_prec <- fit.variogram(svg_prec, model = vgm("Exp"))
plot(svg_prec, model=vgm_model_prec)

svg_wind <- variogram(avg_wind_speed_kmh ~ 1, data = data_boxcox, cutoff = 8000)
vgm_model_wind <- fit.variogram(svg_wind, model = vgm("Exp"))
plot(svg_wind, model=vgm_model_wind)

vgm_model_wind1 <- vgm(psill = 0.035, model = "Nug")
vgm_model_wind1 <- vgm(psill = 0.075, model = "Sph", range = 3000, add.to = vgm_model_wind1)

vgm_model_wind_1 <- fit.variogram(svg_wind, model = vgm_model_wind1)
plot(svg_wind, model = vgm_model_wind1)


hscat(avg_temp_c~1,data_boxcox,(0:9)*100)

#Check for anisotropy:
svgm_anis <- variogram(object = avg_temp_c ~ 1,
                       data = data_multivar,
                       alpha= c(45, 60, 90, 135),
                       cutoff = 1000,
                       width = 10
)
plot(svgm_anis)
svgm_anis <- variogram(object = precipitation_mm ~ 1,
                       data = data_multivar,
                       alpha= c(45, 60, 90, 135),
                       cutoff = 1000,
                       width = 10
)
plot(svgm_anis)
svgm_anis <- variogram(object = avg_wind_speed_kmh ~ 1,
                       data = data_multivar,
                       alpha= c(45, 60, 90, 135),
                       cutoff = 1000,
                       width = 10
)
plot(svgm_anis)
#no anisotropy




#krige
data_boxcox_sp <- as(data_boxcox, "Spatial") #first convert to Spatial

uk_model_temp <- krige(avg_temp_c ~ x + y, data_boxcox_sp, newdata = grid_europe_sf, model = vgm_model_temp)
k_model_prec <- krige(precipitation_mm ~ 1,locations = data_boxcox_sp,newdata = grid_europe_sf,model = vgm_model_prec)
k_model_wind <- krige(avg_wind_speed_kmh ~ 1,locations = data_boxcox_sp,newdata = grid_europe_sf,model = vgm_model_wind)

#untransform every variable:
shift_temp <- boxcox_params[["avg_temp_c"]]$shift
lambda_temp <- boxcox_params[["avg_temp_c"]]$lambda
pred_temp <- uk_model_temp$var1.pred
if (abs(lambda_temp) < 1e-6) {
    untrans_temp <- exp(pred_temp) - shift_temp
    } else {
      untrans_temp <- (pred_temp * lambda_temp + 1)^(1 / lambda_temp) - shift_temp
    }
uk_model_temp$var1.pred_untrans <- untrans_temp
uk_model_temp_sp <- as(uk_model_temp, "Spatial")
#spplot(uk_model_temp_sp, "var1.pred_untrans", main = "Prediction for average January temperature")

uk_model_temp_sf <- st_as_sf(uk_model_temp_sp)
ggplot(uk_model_temp_sf) +
  geom_sf(aes(fill = var1.pred_untrans), color = "black", shape = 21, size = 2) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey50") +
  theme_minimal(base_size = 14) + 
  labs(
    title = "Predicted average temperature",
    x = "Longitude",
    y = "Latitude",
    fill = "Temperature (°C)"
  )


shift_prec <- boxcox_params[["precipitation_mm"]]$shift
lambda_prec <- boxcox_params[["precipitation_mm"]]$lambda
pred_prec <- k_model_prec$var1.pred
if (abs(lambda_prec) < 1e-6) {
  untrans_prec <- exp(pred_prec) - shift_prec
} else {
  untrans_prec <- (pred_prec * lambda_prec + 1)^(1 / lambda_prec) - shift_prec
}
k_model_prec$var1.pred_untrans <- untrans_prec
#k_model_prec_sp <- as(k_model_prec, "Spatial")
#spplot(k_model_prec_sp, "var1.pred_untrans", main = "Prediction for daily precipitation")
k_model_prec_sf <- st_as_sf(k_model_prec)
ggplot(k_model_prec_sf) +
  geom_sf(aes(fill = var1.pred_untrans), color = "black", shape = 21, size = 2) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey50") +
  theme_minimal(base_size = 14) + 
  labs(
    title = "Predicted daily precipitaion",
    x = "Longitude",
    y = "Latitude",
    fill = "Precipitation (mm)"
  )



shift_wind <- boxcox_params[["avg_wind_speed_kmh"]]$shift
lambda_wind <- boxcox_params[["avg_wind_speed_kmh"]]$lambda
pred_wind <- k_model_wind$var1.pred
if (abs(lambda_wind) < 1e-6) {
  untrans_wind <- exp(pred_wind) - shift_wind
} else {
  untrans_wind <- (pred_wind * lambda_wind + 1)^(1 / lambda_wind) - shift_wind
}
k_model_wind$var1.pred_untrans <- untrans_wind
#k_model_wind_sp <- as(k_model_wind, "Spatial")
#spplot(k_model_wind_sp, "var1.pred_untrans", main = "Prediction for average wind speed")
k_model_wind_sf <- st_as_sf(k_model_wind)
ggplot(k_model_wind_sf) +
  geom_sf(aes(fill = var1.pred_untrans), color = "black", shape = 21, size = 2) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey50") +
  theme_minimal(base_size = 14) + 
  labs(
    title = "Predicted wind speed",
    x = "Longitude",
    y = "Latitude",
    fill = "Wind speed (kmh)"
  )

#univariate analysis:
#cross validation(cv, untransform, plot):
cv_temp <- krige.cv(avg_temp_c ~ x + y, data_boxcox_sp, model = vgm_model_temp, nfold = nrow(data_boxcox_sp))

pred_trans <- cv_temp$var1.pred
obs_trans <- cv_temp$observed
pred_untrans <- if (abs(lambda_temp) < 1e-6) {
  exp(pred_trans) - shift_temp
} else {
  (pred_trans * lambda_temp + 1)^(1 / lambda_temp) - shift_temp
}

obs_untrans <- if (abs(lambda_temp) < 1e-6) {
  exp(obs_trans) - shift_temp
} else {
  (obs_trans * lambda_temp + 1)^(1 / lambda_temp) - shift_temp
}

residuals_untrans <- obs_untrans - pred_untrans

cv_temp$pred_untrans <- pred_untrans
cv_temp$obs_untrans <- obs_untrans
cv_temp$residual_untrans <- residuals_untrans
cv_temp_clean <- cv_temp[!is.na(cv_temp$residual_untrans), ]
bubble(cv_temp_clean, "residual_untrans",
       main = "Temperature residuals",
       maxsize = 2.5)

plot(cv_temp_clean$obs_untrans, cv_temp_clean$residual_untrans,
     xlab = "Observed Temperature (°C)",
     ylab = "Residual (°C)",
     main = "Residuals vs Observed")
abline(h = 0, col = "red", lty = 2)


cv_prec <- krige.cv(precipitation_mm ~ 1, data_boxcox_sp, model = vgm_model_prec, nfold = nrow(data_boxcox_sp))

pred_trans <- cv_prec$var1.pred
obs_trans <- cv_prec$observed
pred_untrans <- if (abs(lambda_prec) < 1e-6) {
  exp(pred_trans) - shift_prec
} else {
  (pred_trans * lambda_prec + 1)^(1 / lambda_prec) - shift_prec
}

obs_untrans <- if (abs(lambda_prec) < 1e-6) {
  exp(obs_trans) - shift_prec
} else {
  (obs_trans * lambda_prec + 1)^(1 / lambda_prec) - shift_prec
}
residuals_untrans <- obs_untrans - pred_untrans
cv_prec$pred_untrans <- pred_untrans
cv_prec$obs_untrans <- obs_untrans
cv_prec$residual_untrans <- residuals_untrans
cv_prec_clean <- cv_prec[!is.na(cv_prec$residual_untrans), ]
bubble(cv_prec_clean, "residual_untrans",
       main = "Precipitation residuals",
       maxsize = 2.5)



cv_wind <- krige.cv(avg_wind_speed_kmh ~ 1, data_boxcox_sp, model = vgm_model_wind, nfold = nrow(data_boxcox_sp))

pred_trans <- cv_wind$var1.pred
obs_trans <- cv_wind$observed
pred_untrans <- if (abs(lambda_wind) < 1e-6) {
  exp(pred_trans) - shift_wind
} else {
  (pred_trans * lambda_wind + 1)^(1 / lambda_wind) - shift_wind
}

obs_untrans <- if (abs(lambda_wind) < 1e-6) {
  exp(obs_trans) - shift_wind
} else {
  (obs_trans * lambda_wind + 1)^(1 / lambda_wind) - shift_wind
}
residuals_untrans <- obs_untrans - pred_untrans
cv_wind$pred_untrans <- pred_untrans
cv_wind$obs_untrans <- obs_untrans
cv_wind$residual_untrans <- residuals_untrans
cv_wind_clean <- cv_wind[!is.na(cv_wind$residual_untrans), ]
bubble(cv_wind_clean, "residual_untrans",
       main = "Wind speed residuals",
       maxsize = 2.5)











#multivariate analysis
#first check correlation
cor(as.data.frame(data_multivar)[c("avg_temp_c", "precipitation_mm", "avg_wind_speed_kmh")])

#first standartiize down to the same scale so we can plot cross variograms
data_boxcox$precipitation_mm_z <- scale(data_boxcox$precipitation_mm)
data_boxcox$avg_wind_speed_kmh_z   <- scale(data_boxcox$avg_wind_speed_kmh)
data_boxcox$avg_temp_c_z   <- scale(data_boxcox$avg_temp_c)


g<-gstat(NULL, "avg_temp", avg_temp_c_z~ 1+y, data=data_boxcox)
g<-gstat(g, "precipitation", precipitation_mm_z~ 1+y, data=data_boxcox)
g<-gstat(g, "wind_speed", avg_wind_speed_kmh_z~ 1+y, data=data_boxcox)
#g
vm<-variogram(g, cutoff = 4000)
#plot(vm)
vm.fit <- fit.lmc(vm, g, vgm(1, "Sph", 1500, 1)) 
plot(vm, vm.fit)
grid_europe_sf$x <- st_coordinates(grid_europe_sf)[,1]
grid_europe_sf$y <- st_coordinates(grid_europe_sf)[,2]


vgm1 <- vgm(psill = 0.4,
            model = "Exp", 
            range = 1000
)
vgm1 <- vgm(psill = 0.3,
            model = "Gau", 
            range = 5000,
            add.to = vgm1
)
vgm1_fit<-fit.lmc(vm, g, vgm1)
plot(vm,vgm1_fit)     #this
cok.maps <- predict(vgm1_fit, grid_europe_sf)  #WORKS


#Unstandartize and untransfrom:
# Extract predicted standardized Box-Cox values for temperature
pred_temp_z <- cok.maps$avg_temp.pred
# Unstandardize function
unstandardize <- function(z, mean_val, sd_val) {
    return(z * sd_val + mean_val)}
# Inverse Box-Cox function
inv_boxcox <- function(y, lambda, shift) {
     if (abs(lambda) < 1e-6) {
     return(exp(y) - shift)
     } else {
         return((lambda * y + 1)^(1 / lambda) - shift)}
      }
# Get mean and sd from original scaling of temperature
mean_temp <- attr(scale(data_boxcox$avg_temp_c), "scaled:center")
sd_temp <- attr(scale(data_boxcox$avg_temp_c), "scaled:scale")
# Box-Cox parameters saved earlier
lambda_temp <- boxcox_params[["avg_temp_c"]]$lambda
shift_temp <- boxcox_params[["avg_temp_c"]]$shift
 
#1: Unstandardize
pred_temp_y <- unstandardize(pred_temp_z, mean_temp, sd_temp)

#2: Inverse Box-Cox to original units
pred_temp_original <- inv_boxcox(pred_temp_y, lambda_temp, shift_temp)
# Precipitation
pred_precip_z <- cok.maps$precipitation.pred
mean_precip <- attr(scale(data_boxcox$precipitation_mm), "scaled:center")
sd_precip <- attr(scale(data_boxcox$precipitation_mm), "scaled:scale")
lambda_precip <- boxcox_params[["precipitation_mm"]]$lambda
shift_precip <- boxcox_params[["precipitation_mm"]]$shift
pred_precip_y <- unstandardize(pred_precip_z, mean_precip, sd_precip)
pred_precip_original <- inv_boxcox(pred_precip_y, lambda_precip, shift_precip)
 
# Wind speed
pred_wind_z <- cok.maps$wind_speed.pred
mean_wind <- attr(scale(data_boxcox$avg_wind_speed_kmh), "scaled:center")
sd_wind <- attr(scale(data_boxcox$avg_wind_speed_kmh), "scaled:scale")
lambda_wind <- boxcox_params[["avg_wind_speed_kmh"]]$lambda
shift_wind <- boxcox_params[["avg_wind_speed_kmh"]]$shift
 
pred_wind_y <- unstandardize(pred_wind_z, mean_wind, sd_wind)
pred_wind_original <- inv_boxcox(pred_wind_y, lambda_wind, shift_wind)

cok.maps$avg_temp_original <- pred_temp_original
cok.maps$precipitation_original <- pred_precip_original
cok.maps$wind_speed_original <- pred_wind_original
  
#plot: 
ggplot(cok.maps) +
      geom_point(aes(x = st_coordinates(cok.maps)[,1], y = st_coordinates(cok.maps)[,2], fill = avg_temp_original), 
      shape = 21, size = 3, color = "black") +
      scale_fill_viridis_c(option = "plasma", na.value = "grey50") +
  theme_minimal(base_size = 14) + 
  labs(title = "Predicted Average Temperature",
       x = "Longitude",
       y = "Latitude",
       fill = "Temperature (°C)")
  
ggplot(cok.maps) +
  geom_point(aes(x = st_coordinates(cok.maps)[,1], y = st_coordinates(cok.maps)[,2], fill = precipitation_original), 
             shape = 21, size = 3, color = "black") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey50") +
  theme_minimal(base_size = 14) + 
  labs(title = "Predicted daily precipitation",
       x = "Longitude",
       y = "Latitude",
       fill = "Precipitation (mm)")

ggplot(cok.maps) +
  geom_point(aes(x = st_coordinates(cok.maps)[,1], y = st_coordinates(cok.maps)[,2], fill = wind_speed_original), 
             shape = 21, size = 3, color = "black") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey50") +
  theme_minimal(base_size = 14) + 
  labs(title = "Predicted wind speed",
    x = "Longitude",
    y = "Latitude",
    fill = "Wind speed (km/h)")
  
  

#cross validation: create gstat model using variogram;
#run cross validatio; evaluate results
#temp:
g_temp<-gstat(NULL, "avg_temp", avg_temp_c_z~ 1+y, data=data_boxcox)
g_temp<-gstat(g_temp, "wind_speed", avg_wind_speed_kmh_z~ 1+y, data=data_boxcox)
g_temp<-gstat(g_temp, "precipitation", precipitation_mm_z~ 1+y, data=data_boxcox)
v_temp<-variogram(g_temp, cutoff = 4000)
vgm1 <- vgm(psill = 0.1,
            model = "Exp", 
            range = 1000
)
v_temp_fit<-fit.lmc(v_temp, g_temp, vgm1)
plot(v_temp,v_temp_fit)     #this
temp_cv <- gstat.cv(v_temp_fit, id = "avg_temp", nfold = 5, nmax = 40)

g_wind<-gstat(NULL, "wind_speed", avg_wind_speed_kmh_z~ 1+y, data=data_boxcox)
g_wind<-gstat(g_wind, "precipitation", precipitation_mm_z~ 1+y, data=data_boxcox)
g_wind<-gstat(g_wind, "avg_temp", avg_temp_c_z~ 1+y, data=data_boxcox)
v_wind<-variogram(g_wind, cutoff = 4000)
vgm1 <- vgm(psill = 0.4,
            model = "Exp", 
            range = 1000
)
vgm1 <- vgm(psill = 0.3,
            model = "Gau", 
            range = 5000,
            add.to = vgm1
)
v_wind_fit<-fit.lmc(v_wind, g_wind, vgm1)
plot(v_wind,v_wind_fit)     #this
wind_cv <- gstat.cv(v_wind_fit, id = "wind_speed", nfold = 5, nmax = 40)


g_prec<-gstat(NULL, "precipitation", precipitation_mm_z~ 1+y, data=data_boxcox)
g_prec<-gstat(g_prec, "avg_temp", avg_temp_c_z~ 1+y, data=data_boxcox)
g_prec<-gstat(g_prec, "wind_speed", avg_wind_speed_kmh_z~ 1+y, data=data_boxcox)
v_prec<-variogram(g_prec, cutoff = 4000)
vgm1 <- vgm(psill = 0.1,
            model = "Exp", 
            range = 1000
)
v_prec_fit<-fit.lmc(v_prec, g_prec, vgm1)
plot(v_prec,v_prec_fit)     #this
prec_cv <- gstat.cv(v_prec_fit, id = "precipitation", nfold = 5, nmax = 40)

#untransform:
temp_cv$pred_unscaled <- unstandardize(temp_cv$avg_temp.pred, mean_temp, sd_temp)
temp_cv$obs_unscaled <- unstandardize(temp_cv$observed, mean_temp, sd_temp)
temp_cv$pred_original <- inv_boxcox(temp_cv$pred_unscaled, lambda_temp, shift_temp)
temp_cv$obs_original <- inv_boxcox(temp_cv$obs_unscaled, lambda_temp, shift_temp)
temp_cv$resid_original <- temp_cv$obs_original - temp_cv$pred_original

prec_cv$pred_unscaled <- unstandardize(prec_cv$precipitation.pred, mean_precip, sd_precip)
prec_cv$obs_unscaled <- unstandardize(prec_cv$observed, mean_precip, sd_precip)
prec_cv$pred_original <- inv_boxcox(prec_cv$pred_unscaled, lambda_precip, shift_precip)
prec_cv$obs_original <- inv_boxcox(prec_cv$obs_unscaled, lambda_precip, shift_precip)
prec_cv$resid_original <- prec_cv$obs_original - prec_cv$pred_original

wind_cv$pred_unscaled <- unstandardize(wind_cv$wind_speed.pred, mean_wind, sd_wind)
wind_cv$obs_unscaled <- unstandardize(wind_cv$observed, mean_wind, sd_wind)
wind_cv$pred_original <- inv_boxcox(wind_cv$pred_unscaled, lambda_wind, shift_wind)
wind_cv$obs_original <- inv_boxcox(wind_cv$obs_unscaled, lambda_wind, shift_wind)
wind_cv$resid_original <- wind_cv$obs_original - wind_cv$pred_original

#plot residuals:
temp_cv_clean <- temp_cv[!is.na(temp_cv$resid_original), ]
bubble(temp_cv_clean, "resid_original",
     main = "Temperature residuals",
    maxsize = 2.5) 
bubble(prec_cv, "resid_original",
       main = "Precipitation residuals",
       maxsize = 2.5) 
bubble(wind_cv, "resid_original",
       main = "Wind speed residuals",
       maxsize = 2.5) 


#compare
# temp
rmse <- function(resid) sqrt(mean(resid^2))
rmse_univ_temp <- rmse(cv_temp_clean$residual_untrans)
rmse_multi_temp <- rmse(temp_cv_clean$resid_original)

# MAE
mae <- function(resid) mean(abs(resid))

mae_uni_temp <- mae(cv_temp_clean$residual_untrans)
mae_multi_temp <- mae(temp_cv_clean$resid_original)

# Print results
cat("Univariate Kriging - RMSE:", rmse_univ_temp, "MAE:", mae_uni_temp, "\n")
cat("Cokriging - RMSE:", rmse_multi_temp, "MAE:", mae_multi_temp, "\n")

#precip
rmse_univ_precip<- rmse(cv_prec_clean$residual_untrans)
rmse_multi_precip<- rmse(prec_cv$resid_original)

# MAE
mae_uni_precip <- mae(cv_prec_clean$residual_untrans)
mae_multi_precip <- mae(prec_cv$resid_original)

# Print results
cat("Univariate Kriging - RMSE:", rmse_univ_precip, "MAE:", mae_uni_precip, "\n")
cat("Cokriging - RMSE:", rmse_multi_precip, "MAE:", mae_multi_precip, "\n")


#wind
rmse_univ_wind<- rmse(cv_wind_clean$residual_untrans)
rmse_multi_wind<- rmse(wind_cv$resid_original)

# MAE
mae_uni_wind <- mae(cv_wind_clean$residual_untrans)
mae_multi_wind <- mae(wind_cv$resid_original)

# Print results
cat("Univariate Kriging - RMSE:", rmse_univ_wind, "MAE:", mae_uni_wind, "\n")
cat("Cokriging - RMSE:", rmse_multi_wind, "MAE:", mae_multi_wind, "\n")

#print results:
results <- data.frame(
  Variable = c("Temperature", "Temperature", "Precipitation", "Precipitation", "Wind Speed", "Wind Speed"),
  Method = c("Univariate Kriging", "Cokriging", "Univariate Kriging", "Cokriging", "Univariate Kriging", "Cokriging"),
  RMSE = c(2.090809, 2.732024, 1.19024, 1.256789, 4.717023, 4.732893),
  MAE = c(1.26939, 1.748593, 0.6229189, 0.6861993, 3.445825, 3.468553)
)
