# Load libraries
library(raster)
library(sp)
library(MASS)
library(ENMeval)

# ==== 1. Load and stack bioclimatic variables ====
env <- stack(
  raster("2.5mbio2.asc"),
  raster("2.5mbio3.asc"),
  raster("2.5mbio4.asc"),
  raster("2.5mbio6.asc"),
  raster("2.5mbio13.asc"),
  raster("2.5mbio14.asc"),
  raster("2.5mbio15.asc")
)

# ==== 2. Load and prepare occurrence points ====
occ_raw <- read.csv("Alltimes.csv")
occ <- na.omit(occ_raw[, c("x", "y")])  # Adjust column names if needed
coordinates(occ) <- ~x + y
crs(occ) <- crs(env)

# ==== 3. Create KDE bias raster (with fallback if < 2 points) ====
message("Creating bias file...")
occur.ras <- rasterize(occ, env, 1)

# Extract presence coordinates
presences <- which(values(occur.ras) == 1)
pres.locs <- xyFromCell(occur.ras, presences)

if (!is.null(pres.locs) && nrow(pres.locs) >= 2) {
  dens <- kde2d(
    pres.locs[,1], pres.locs[,2],
    n = c(ncol(env), nrow(env)),
    lims = c(extent(env)[1:2], extent(env)[3:4])
  )
  bias.ras <- raster(list(x = dens$x, y = dens$y, z = dens$z))
  crs(bias.ras) <- crs(env)
  extent(bias.ras) <- extent(env)
  bias.ras <- resample(bias.ras, env)
  writeRaster(bias.ras, "allset1bias.asc", format = "ascii", overwrite = TRUE)
  message("Bias file created: allset1bias.asc")
} else {
  warning("Fewer than 2 presence points; using uniform bias raster.")
  bias.ras <- env[[1]]
  values(bias.ras) <- 1
}

# ==== 4. Sample background points using bias ====
message("Sampling background points...")
bg.n <- min(10000, length(which(!is.na(values(env[[1]])))))
bg.cells <- sample(
  which(!is.na(values(env[[1]]))),
  size = bg.n,
  prob = values(bias.ras)[!is.na(values(env[[1]]))]
)
bg <- xyFromCell(bias.ras, bg.cells)
colnames(bg) <- c("x", "y")

# ==== 5. Run ENMevaluate with block partitioning ====
message("Running ENMevaluate...")
occ.df <- as.data.frame(coordinates(occ))

enmeval_results <- ENMevaluate(
  occ = occ.df,
  env = env,
  bg.coords = bg,
  tune.args = list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), rm = 1:5),
  partitions = "jackknife",
  algorithm = "maxnet"
)

# ==== 6. Save results ====
write.csv(enmeval_results@results, "allset1enmeval_results_block.csv")
message("Model tuning complete. Results saved to allset1enmeval_results_block.csv")

# ==== 7. (Optional) Plot model evaluation ====
plot(enmeval_results)


#set 2
# Load libraries
library(raster)
library(sp)
library(MASS)
library(ENMeval)

# ==== 1. Load and stack bioclimatic variables ====
env <- stack(
  raster("2.5mbio2.asc"),
  raster("2.5mbio3.asc"),
  raster("2.5mbio4.asc"),
  raster("2.5mbio6.asc")
)

# ==== 2. Load and prepare occurrence points ====
occ_raw <- read.csv("Alltimes.csv")
occ <- na.omit(occ_raw[, c("x", "y")])  # Adjust column names if needed
coordinates(occ) <- ~x + y
crs(occ) <- crs(env)

# ==== 3. Create KDE bias raster (with fallback if < 2 points) ====
message("Creating bias file...")
occur.ras <- rasterize(occ, env, 1)

# Extract presence coordinates
presences <- which(values(occur.ras) == 1)
pres.locs <- xyFromCell(occur.ras, presences)

if (!is.null(pres.locs) && nrow(pres.locs) >= 2) {
  dens <- kde2d(
    pres.locs[,1], pres.locs[,2],
    n = c(ncol(env), nrow(env)),
    lims = c(extent(env)[1:2], extent(env)[3:4])
  )
  bias.ras <- raster(list(x = dens$x, y = dens$y, z = dens$z))
  crs(bias.ras) <- crs(env)
  extent(bias.ras) <- extent(env)
  bias.ras <- resample(bias.ras, env)
  writeRaster(bias.ras, "allset2bias.asc", format = "ascii", overwrite = TRUE)
  message("Bias file created: allset2bias.asc")
} else {
  warning("Fewer than 2 presence points; using uniform bias raster.")
  bias.ras <- env[[1]]
  values(bias.ras) <- 1
}

# ==== 4. Sample background points using bias ====
message("Sampling background points...")
bg.n <- min(10000, length(which(!is.na(values(env[[1]])))))
bg.cells <- sample(
  which(!is.na(values(env[[1]]))),
  size = bg.n,
  prob = values(bias.ras)[!is.na(values(env[[1]]))]
)
bg <- xyFromCell(bias.ras, bg.cells)
colnames(bg) <- c("x", "y")

# ==== 5. Run ENMevaluate with block partitioning ====
message("Running ENMevaluate...")
occ.df <- as.data.frame(coordinates(occ))

enmeval_results <- ENMevaluate(
  occ = occ.df,
  env = env,
  bg.coords = bg,
  tune.args = list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), rm = 1:5),
  partitions = "jackknife",
  algorithm = "maxnet"
)

# ==== 6. Save results ====
write.csv(enmeval_results@results, "allset2enmeval_results_block.csv")
message("Model tuning complete. Results saved to allset2enmeval_results_block.csv")

# ==== 7. (Optional) Plot model evaluation ====
plot(enmeval_results)


#set3
# Load libraries
library(raster)
library(sp)
library(MASS)
library(ENMeval)

# ==== 1. Load and stack bioclimatic variables ====
env <- stack(
  raster("2.5mbio13.asc"),
  raster("2.5mbio14.asc"),
  raster("2.5mbio15.asc")
)

# ==== 2. Load and prepare occurrence points ====
occ_raw <- read.csv("Alltimes.csv")
occ <- na.omit(occ_raw[, c("x", "y")])  # Adjust column names if needed
coordinates(occ) <- ~x + y
crs(occ) <- crs(env)

# ==== 3. Create KDE bias raster (with fallback if < 2 points) ====
message("Creating bias file...")
occur.ras <- rasterize(occ, env, 1)

# Extract presence coordinates
presences <- which(values(occur.ras) == 1)
pres.locs <- xyFromCell(occur.ras, presences)

if (!is.null(pres.locs) && nrow(pres.locs) >= 2) {
  dens <- kde2d(
    pres.locs[,1], pres.locs[,2],
    n = c(ncol(env), nrow(env)),
    lims = c(extent(env)[1:2], extent(env)[3:4])
  )
  bias.ras <- raster(list(x = dens$x, y = dens$y, z = dens$z))
  crs(bias.ras) <- crs(env)
  extent(bias.ras) <- extent(env)
  bias.ras <- resample(bias.ras, env)
  writeRaster(bias.ras, "allset3bias.asc", format = "ascii", overwrite = TRUE)
  message("Bias file created: allset3bias.asc")
} else {
  warning("Fewer than 2 presence points; using uniform bias raster.")
  bias.ras <- env[[1]]
  values(bias.ras) <- 1
}

# ==== 4. Sample background points using bias ====
message("Sampling background points...")
bg.n <- min(10000, length(which(!is.na(values(env[[1]])))))
bg.cells <- sample(
  which(!is.na(values(env[[1]]))),
  size = bg.n,
  prob = values(bias.ras)[!is.na(values(env[[1]]))]
)
bg <- xyFromCell(bias.ras, bg.cells)
colnames(bg) <- c("x", "y")

# ==== 5. Run ENMevaluate with block partitioning ====
message("Running ENMevaluate...")
occ.df <- as.data.frame(coordinates(occ))

enmeval_results <- ENMevaluate(
  occ = occ.df,
  env = env,
  bg.coords = bg,
  tune.args = list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), rm = 1:5),
  partitions = "jackknife",
  algorithm = "maxnet"
)

# ==== 6. Save results ====
write.csv(enmeval_results@results, "allset3enmeval_results_block.csv")
message("Model tuning complete. Results saved to allset3enmeval_results_block.csv")

# ==== 7. (Optional) Plot model evaluation ====
plot(enmeval_results)


#set4
# Load libraries
library(raster)
library(sp)
library(MASS)
library(ENMeval)

# ==== 1. Load and stack bioclimatic variables ====
env <- stack(
  raster("2.5mbio3.asc"),
  raster("2.5mbio6.asc"),
  raster("2.5mbio13.asc"),
  raster("2.5mbio14.asc")
)

# ==== 2. Load and prepare occurrence points ====
occ_raw <- read.csv("Alltimes.csv")
occ <- na.omit(occ_raw[, c("x", "y")])  # Adjust column names if needed
coordinates(occ) <- ~x + y
crs(occ) <- crs(env)

# ==== 3. Create KDE bias raster (with fallback if < 2 points) ====
message("Creating bias file...")
occur.ras <- rasterize(occ, env, 1)

# Extract presence coordinates
presences <- which(values(occur.ras) == 1)
pres.locs <- xyFromCell(occur.ras, presences)

if (!is.null(pres.locs) && nrow(pres.locs) >= 2) {
  dens <- kde2d(
    pres.locs[,1], pres.locs[,2],
    n = c(ncol(env), nrow(env)),
    lims = c(extent(env)[1:2], extent(env)[3:4])
  )
  bias.ras <- raster(list(x = dens$x, y = dens$y, z = dens$z))
  crs(bias.ras) <- crs(env)
  extent(bias.ras) <- extent(env)
  bias.ras <- resample(bias.ras, env)
  writeRaster(bias.ras, "allset4bias.asc", format = "ascii", overwrite = TRUE)
  message("Bias file created: allset4bias.asc")
} else {
  warning("Fewer than 2 presence points; using uniform bias raster.")
  bias.ras <- env[[1]]
  values(bias.ras) <- 1
}

# ==== 4. Sample background points using bias ====
message("Sampling background points...")
bg.n <- min(10000, length(which(!is.na(values(env[[1]])))))
bg.cells <- sample(
  which(!is.na(values(env[[1]]))),
  size = bg.n,
  prob = values(bias.ras)[!is.na(values(env[[1]]))]
)
bg <- xyFromCell(bias.ras, bg.cells)
colnames(bg) <- c("x", "y")

# ==== 5. Run ENMevaluate with block partitioning ====
message("Running ENMevaluate...")
occ.df <- as.data.frame(coordinates(occ))

enmeval_results <- ENMevaluate(
  occ = occ.df,
  env = env,
  bg.coords = bg,
  tune.args = list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), rm = 1:5),
  partitions = "jackknife",
  algorithm = "maxnet"
)

# ==== 6. Save results ====
write.csv(enmeval_results@results, "allset4enmeval_results_block.csv")
message("Model tuning complete. Results saved to allset4enmeval_results_block.csv")

# ==== 7. (Optional) Plot model evaluation ====
plot(enmeval_results)
