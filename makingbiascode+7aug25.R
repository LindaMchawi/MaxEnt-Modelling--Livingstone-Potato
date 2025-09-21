library(raster)

#Model 1
# Load the bias raster
bias <- raster("set1bias.asc")

# Replace all zero or negative values with a small positive number (e.g., 0.01)
bias[bias <= 0] <- 0.01

# Optionally normalize (optional but helpful)
bias <- bias / cellStats(bias, stat = 'max')

# Save corrected bias raster
writeRaster(bias, "set1bias_fixed.asc", format="ascii", overwrite=TRUE)


#Model 2
# Load the bias raster
bias <- raster("set2bias.asc")

# Replace all zero or negative values with a small positive number (e.g., 0.01)
bias[bias <= 0] <- 0.01

# Optionally normalize (optional but helpful)
bias <- bias / cellStats(bias, stat = 'max')

# Save corrected bias raster
writeRaster(bias, "set2bias_fixed.asc", format="ascii", overwrite=TRUE)


#Model 3
# Load the bias raster
bias <- raster("set3bias.asc")

# Replace all zero or negative values with a small positive number (e.g., 0.01)
bias[bias <= 0] <- 0.01

# Optionally normalize (optional but helpful)
bias <- bias / cellStats(bias, stat = 'max')

# Save corrected bias raster
writeRaster(bias, "set3bias_fixed.asc", format="ascii", overwrite=TRUE)


#Model 4
# Load the bias raster
bias <- raster("set4bias.asc")

# Replace all zero or negative values with a small positive number (e.g., 0.01)
bias[bias <= 0] <- 0.01

# Optionally normalize (optional but helpful)
bias <- bias / cellStats(bias, stat = 'max')

# Save corrected bias raster
writeRaster(bias, "set4bias_fixed.asc", format="ascii", overwrite=TRUE)


# For pointsample alltime
#Model 1
# Load the bias raster
bias <- raster("allset1bias.asc")

# Replace all zero or negative values with a small positive number (e.g., 0.01)
bias[bias <= 0] <- 0.01

# Optionally normalize (optional but helpful)
bias <- bias / cellStats(bias, stat = 'max')

# Save corrected bias raster
writeRaster(bias, "allset1bias_fixed.asc", format="ascii", overwrite=TRUE)


#Model 2
# Load the bias raster
bias <- raster("allset2bias.asc")

# Replace all zero or negative values with a small positive number (e.g., 0.01)
bias[bias <= 0] <- 0.01

# Optionally normalize (optional but helpful)
bias <- bias / cellStats(bias, stat = 'max')

# Save corrected bias raster
writeRaster(bias, "allset2bias_fixed.asc", format="ascii", overwrite=TRUE)


#Model 3
# Load the bias raster
bias <- raster("allset3bias.asc")

# Replace all zero or negative values with a small positive number (e.g., 0.01)
bias[bias <= 0] <- 0.01

# Optionally normalize (optional but helpful)
bias <- bias / cellStats(bias, stat = 'max')

# Save corrected bias raster
writeRaster(bias, "allset3bias_fixed.asc", format="ascii", overwrite=TRUE)


#Model 4
# Load the bias raster
bias <- raster("allset4bias.asc")

# Replace all zero or negative values with a small positive number (e.g., 0.01)
bias[bias <= 0] <- 0.01

# Optionally normalize (optional but helpful)
bias <- bias / cellStats(bias, stat = 'max')

# Save corrected bias raster
writeRaster(bias, "allset4bias_fixed.asc", format="ascii", overwrite=TRUE)