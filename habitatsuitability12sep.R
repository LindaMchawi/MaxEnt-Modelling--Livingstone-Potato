# ================================
# Malawi-only Habitat Suitability (P10) & Change
# ================================
# Outputs:
#  - current_P10_binary_Malawi.asc / .png
#  - future_P10_binary_Malawi.asc  / .png
#  - change_LSG_Malawi.asc         / .png   (Loss=-1, Stable=0, Gain=1)
#  - area_summary_P10_Malawi.csv
#
# Requirements: terra, sf, rnaturalearth, rnaturalearthdata, (optional) readr
# ================================

##50yrs ssp245 model set 1
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._4.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._4_Set 1.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.677  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##50yrs ssp245 model set 2
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._6.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._6_Set 2.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.578  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##50yrs ssp245 model set 3
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._0.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._0_Set 3.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.606  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##50yrs ssp245 model set 4
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._8.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._8_Set 4.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.493  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##50yrs ssp585 model set 1
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._4.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._4_Set 1.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.603  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")



##50yrs ssp585 model set 2
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._2.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._2_Set 2.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.593  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##50yrs ssp585 model set 3
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._7.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._7_Set 3.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.517  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##50yrs ssp585 model set 4
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._1.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._1_Set 4.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.671  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##Alltime ssp245 model set 1
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._4.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._4_set1.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.667  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##Alltime ssp245 model set 2
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._8.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._8_set2.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.541  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##Alltime ssp245 model set 3
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._5.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._5_set3.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.510  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##Alltime ssp245 model set 4
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._0.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._0_set4.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.371  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##Alltime ssp585 model set 1
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._8.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._8_set1.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.572  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##50yrs ssp585 model set 2
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._1.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._1_set2.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.563  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##50yrs ssp585 model set 3
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._2.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._2_set3.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.455  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")


##50yrs ssp585 model set 4
# --- 0) Libraries ---
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(readr)  # only needed if you want to auto-compute P10 from presence CSV
})

# --- 1) USER INPUTS ---
# Continuous MaxEnt rasters (0–1) in ASCII (.asc)
current_path <- "Plectranthus_esculentus_N.E.Br._5.asc"
future_path  <- "Plectranthus_esculentus_N.E.Br._5_set4.asc"

# If your ASCII has NoData flagged as -9999, set here; else set to NA
NA_value <- -9999

# (A) Auto-compute P10 from presence predictions CSV (training presences)
presence_csv <- "Plectranthus_esculentus_N.E.Br._4_samplePredictions.csv"  # or set to NULL to skip
prediction_col <- NULL  # set if you know it (e.g., "logistic"). Else script auto-picks.

# (B) OR provide a manual P10 (set to NA to ignore and use CSV instead)
p10_manual <- 0.355  # e.g., 0.421  or NA_real_ to compute from CSV

# --- 2) Helper: get Malawi boundary (admin-0) ---
mw <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")
mw <- st_make_valid(mw)
mw_vect <- vect(st_transform(mw, crs = "+proj=longlat +datum=WGS84"))

# --- 3) Load rasters & harmonize ---
r_cur <- rast(current_path)
r_fut <- rast(future_path)

# Handle NoData
if (!is.na(NA_value)) {
  NAflag(r_cur) <- NA_value
  NAflag(r_fut) <- NA_value
}

# Ensure same grid
if (!compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_cur, method = "bilinear")
}

# Clip to Malawi
r_cur_mw <- mask(crop(r_cur, mw_vect), mw_vect)
r_fut_mw <- mask(crop(r_fut, mw_vect), mw_vect)

# --- 4) Compute P10 (if not supplied) ---
get_p10_from_csv <- function(csv_path, pred_col = NULL) {
  if (is.null(csv_path) || !file.exists(csv_path)) return(NA_real_)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!is.null(pred_col) && pred_col %in% names(df)) {
    sc <- df[[pred_col]]
  } else {
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) == 0) stop("No numeric columns in presence CSV to compute P10.")
    sc <- df[[tail(num_cols, 1)]]
  }
  as.numeric(quantile(sc, probs = 0.10, na.rm = TRUE))
}

if (is.na(p10_manual)) {
  p10 <- get_p10_from_csv(presence_csv, prediction_col)
  if (is.na(p10)) stop("P10 not provided and could not be computed from CSV. Set p10_manual.")
} else {
  p10 <- p10_manual
}
message(sprintf("Using P10 threshold = %.6f", p10))

# --- 5) Binary maps (0=unsuitable, 1=suitable) ---
r_cur_bin <- ifel(r_cur_mw >= p10, 1, 0)
r_fut_bin <- ifel(r_fut_mw >= p10, 1, 0)
names(r_cur_bin) <- "current_P10"
names(r_fut_bin) <- "future_P10"

# --- 6) Change map: Loss/Stable/Gain ---
chg <- r_fut_bin - r_cur_bin
# chg = -1 (Loss: 1→0), 0 (Stable: 0→0 or 1→1), +1 (Gain: 0→1)
names(chg) <- "change_LSG"

# --- 7) Area calculations (km²) ---
aea <- "+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=0 +lon_0=25 +datum=WGS84 +units=m"
cur_bin_aea <- project(r_cur_bin, aea, method = "near")
fut_bin_aea <- project(r_fut_bin, aea, method = "near")
chg_aea     <- project(chg,      aea, method = "near")

cell_km2 <- cellSize(cur_bin_aea, unit = "km")

area_suitable_current <- global(cur_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]
area_suitable_future  <- global(fut_bin_aea * cell_km2, "sum", na.rm = TRUE)[[1]]

area_loss   <- global((chg_aea == -1) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_stable <- global((chg_aea ==  0) * cell_km2, "sum", na.rm = TRUE)[[1]]
area_gain   <- global((chg_aea ==  1) * cell_km2, "sum", na.rm = TRUE)[[1]]

percent_change <- 100 * (area_suitable_future - area_suitable_current) / area_suitable_current

# Ensure integer types for ASCII grids
r_cur_bin <- round(r_cur_bin)        # 0/1
r_fut_bin <- round(r_fut_bin)        # 0/1
chg       <- round(chg)              # -1/0/1

NA_value <- -9999  # keep whatever you set earlier

# Write ESRI ASCII by using .asc extension (no filetype argument)
writeRaster(r_cur_bin, "current_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(r_fut_bin, "future_P10_binary_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

writeRaster(chg, "change_LSG_Malawi.asc",
            overwrite = TRUE, NAflag = NA_value)

# --- 9) PNG plots ---
png("current_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_cur_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Current Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("future_P10_binary_Malawi.png", width = 1600, height = 1200, res = 200)
plot(r_fut_bin, col = c("grey85", "darkgreen"), legend = FALSE, main = "Future Suitability (P10) — Malawi")
legend("topright", fill = c("grey85","darkgreen"), legend = c("Unsuitable","Suitable"), bty = "n")
dev.off()

png("change_LSG_Malawi.png", width = 1600, height = 1200, res = 200)
plot(chg, col = c("red","grey70","forestgreen"), legend = FALSE, main = "Habitat Change (P10) — Malawi")
legend("topright", fill = c("red","grey70","forestgreen"),
       legend = c("Loss (1→0)","Stable (0→0,1→1)","Gain (0→1)"), bty = "n")
dev.off()

# --- 10) CSV summary ---
summary_df <- data.frame(
  metric = c("Suitable area (current)", "Suitable area (future)", "Loss", "Stable", "Gain", "Percent change"),
  value  = c(area_suitable_current, area_suitable_future, area_loss, area_stable, area_gain, percent_change),
  unit   = c("km^2","km^2","km^2","km^2","km^2","%")
)
write.csv(summary_df, "area_summary_P10_Malawi.csv", row.names = FALSE)

# --- 11) Console summary ---
cat("==== Malawi P10 Summary ====\n")
cat(sprintf("P10 threshold ...........: %.6f\n", p10))
cat(sprintf("Current suitable area ...: %.2f km²\n", area_suitable_current))
cat(sprintf("Future suitable area ....: %.2f km²\n", area_suitable_future))
cat(sprintf("Loss  ...................: %.2f km²\n", area_loss))
cat(sprintf("Stable ..................: %.2f km²\n", area_stable))
cat(sprintf("Gain  ...................: %.2f km²\n", area_gain))
cat(sprintf("Percent change ..........: %.2f %%\n", percent_change))
cat("ASCII outputs: current_P10_binary_Malawi.asc, future_P10_binary_Malawi.asc, change_LSG_Malawi.asc\n")
cat("PNGs: current_P10_binary_Malawi.png, future_P10_binary_Malawi.png, change_LSG_Malawi.png\n")
cat("CSV: area_summary_P10_Malawi.csv\n")
