
# ------Calcul en boucle de NDVI, MNDWI, NDBI sur 32 images Landsat------
# Prétraitement : reprojection + rééchantillonnage + découpe
#                 sur un raster de référence
# Visualisation : composite saison sèche (mai-oct) vs
#                 saison humide (nov-avril)

library(terra)

# ---- 1. Paramètres à adapter ----

# Dossier contenant toutes les bandes des 32 images (B3, B4, B5, B6)
dossier_images <- "C:/Users/Spagnoli/Desktop/DATA/LANDSAT/bandes"

# Raster de référence : définit résolution, emprise et projection cible
ref <- rast("C:/Users/Spagnoli/Desktop/R/predicteurs/lst_160225_1k.tif")

# Dossiers de sortie
dossier_sortie <- "C:/Users/Spagnoli/Desktop/R"
dir.create(file.path(dossier_sortie, "NDVI"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dossier_sortie, "MNDWI"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dossier_sortie, "NDBI"),  recursive = TRUE, showWarnings = FALSE)

# ---- 2. Fonction de prétraitement (reprojection + résolution + découpe) ----
# On reprojette d'abord sur la CRS de ref, puis on rééchantillonne
# directement sur la grille de ref : cela aligne en une seule étape
# la résolution, l'emprise (découpe/extension) et la projection.
aligner_raster <- function(r, ref, methode = "bilinear") {
  if (crs(r) != crs(ref)) {
    r <- project(r, crs(ref))
  }
  r <- resample(r, ref, method = methode)
  return(r)
}

# ---- 3. Repérer les 32 images et regrouper les bandes par date ----
# Récupère la date au format YYYYMMDD quelque part dans le nom de fichier
# (convention Landsat standard, ex: LC08_L2SP_..._20210615_..._B3.TIF)

fichiers <- list.files(dossier_images, pattern = "\\.tif$|\\.TIF$", full.names = TRUE)

extraire_date <- function(nom_fichier) {
  m <- regmatches(nom_fichier, regexpr("[0-9]{8}", nom_fichier))
  if (length(m) == 0) return(NA)
  as.Date(m, format = "%Y%m%d")
}

extraire_bande <- function(nom_fichier) {
  m <- regmatches(nom_fichier, regexpr("B[3-6](?!\\d)", nom_fichier, perl = TRUE))
  if (length(m) == 0) return(NA)
  m
}

infos <- data.frame(
  fichier = fichiers,
  date    = sapply(fichiers, extraire_date),
  bande   = sapply(fichiers, extraire_bande)
)
infos$date <- as.Date(infos$date, origin = "1970-01-01")
infos <- infos[!is.na(infos$date) & !is.na(infos$bande), ]

dates_uniques <- sort(unique(infos$date))
cat("Nombre d'images (dates) détectées :", length(dates_uniques), "\n")

# ---- 4. Boucle principale : prétraitement + calcul des indices ----

resultats <- data.frame(date = as.Date(character()), saison = character(),
                        stringsAsFactors = FALSE)

# Empilements pour stocker les rasters NDVI en mémoire (utile pour les composites)
pile_ndvi <- list()

for (i in seq_along(dates_uniques)) {
  
  d <- dates_uniques[i]
  cat("Traitement image", i, "/", length(dates_uniques), "-", as.character(d), "\n")
  
  sous <- infos[infos$date == d, ]
  
  chemin_b3 <- sous$fichier[sous$bande == "B3"]
  chemin_b4 <- sous$fichier[sous$bande == "B4"]
  chemin_b5 <- sous$fichier[sous$bande == "B5"]
  chemin_b6 <- sous$fichier[sous$bande == "B6"]
  
  if (length(chemin_b3) == 0 || length(chemin_b4) == 0 ||
      length(chemin_b5) == 0 || length(chemin_b6) == 0) {
    warning("Bande manquante pour la date ", d, " - image ignorée")
    next
  }
  
  # Chargement
  b3 <- rast(chemin_b3[1])  # Vert
  b4 <- rast(chemin_b4[1])  # Rouge
  b5 <- rast(chemin_b5[1])  # PIR (NIR)
  b6 <- rast(chemin_b6[1])  # SWIR1
  
  # Prétraitement : reprojection + résolution + découpe sur ref
  b3 <- aligner_raster(b3, ref)
  b4 <- aligner_raster(b4, ref)
  b5 <- aligner_raster(b5, ref)
  b6 <- aligner_raster(b6, ref)
  
  # ---- Calcul des indices ----
  ndvi  <- (b5 - b4) / (b5 + b4)          # (NIR - Rouge) / (NIR + Rouge)
  mndwi <- (b3 - b6) / (b3 + b6)          # (Vert - SWIR1) / (Vert + SWIR1)
  ndbi  <- (b6 - b5) / (b6 + b5)          # (SWIR1 - NIR) / (SWIR1 + NIR)
  
  names(ndvi)  <- paste0("NDVI_",  format(d, "%Y%m%d"))
  names(mndwi) <- paste0("MNDWI_", format(d, "%Y%m%d"))
  names(ndbi)  <- paste0("NDBI_",  format(d, "%Y%m%d"))
  
  # Sauvegarde sur disque
  writeRaster(ndvi,  file.path(dossier_sortie, "NDVI",  paste0(names(ndvi), ".tif")),
              overwrite = TRUE)
  writeRaster(mndwi, file.path(dossier_sortie, "MNDWI", paste0(names(mndwi), ".tif")),
              overwrite = TRUE)
  writeRaster(ndbi,  file.path(dossier_sortie, "NDBI",  paste0(names(ndbi), ".tif")),
              overwrite = TRUE)
  
  # Détermination de la saison : mai-octobre = sèche, nov-avril = humide
  mois <- as.integer(format(d, "%m"))
  saison <- ifelse(mois >= 5 & mois <= 10, "seche", "humide")
  
  pile_ndvi[[as.character(d)]] <- ndvi
  resultats <- rbind(resultats, data.frame(date = d, saison = saison))
}

cat("\nTraitement terminé.\n")
print(table(resultats$saison))
