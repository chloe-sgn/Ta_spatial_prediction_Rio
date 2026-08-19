# Calcul LST sur Landsat 8 et 9 lvl 1 la nuit

library(terra)
library(tidyverse)


DOSSIER_B10  <- "C:/Users/Spagnoli/Desktop/DATA/LANDSAT/landsat_nuit"
DOSSIER_LST  <- "C:/Users/Spagnoli/Desktop/DATA/LANDSAT/LST_n"
RASTER_EMPRISE <- "C:/Users/Spagnoli/Desktop/DATA/LANDSAT/LST_n/LST_060225.tif"
dir.create(DOSSIER_LST, showWarnings = FALSE, recursive = TRUE)

# Constantes LST
LAMBDA  <- 10.895  # longueur d'onde moy bande 10 (Sholihah & Shibata, 2019)
RHO     <- 14388   # en µK
EPSILON <- 0.97    # émissivité constante (Amorim & Dubreuil, 2016)
KELVIN  <- 273.15  # conversion K -> °C

# ----FONCTION LECTURE MTL----
#Images Landsat Level-1 sont livrées avec un fichier *_MTL.txt
#contenant tous les paramètres radiométriques. On lit ce fichier
#texte pour en extraire les valeurs dont on a besoin.

lire_param_MTL <- function(chemin_mtl, cle) {
  # Lit toutes les lignes du fichier MTL
  lignes <- readLines(chemin_mtl)
  # Cherche la ligne contenant exactement la clé (ex: "RADIANCE_MULT_BAND_10")
  ligne  <- lignes[grepl(paste0("\\b", cle, "\\b"), lignes)]
  if (length(ligne) == 0) stop(paste("Clé MTL introuvable :", cle))
  # Extrait la valeur numérique après le "="
  as.numeric(trimws(sub(".*=\\s*", "", ligne[1])))
}

# ----CHARGEMENT DE L'EMPRISE----

cat("=== Chargement de l'emprise ===\n")
raster_ref <- rast(RASTER_EMPRISE)
cat("    CRS     :", crs(raster_ref, describe = TRUE)$name, "\n")
cat("    Étendue :", paste(round(as.vector(ext(raster_ref)), 4), collapse = " | "), "\n")

# ----LISTE DES FICHIERS B10----

fichiers <- list.files(DOSSIER_B10,
                       pattern = "LC0[89].*B10\\.TIF$",  # Level-1 : _B10.TIF (pas ST_B10)
                       full.names = TRUE)

cat("=== Fichiers trouvés :", length(fichiers), "===\n")
for (f in fichiers) cat(" •", basename(f), "\n")

# ----FONCTION PRINCIPALE : B10 -> Radiance TOA -> BT -> LST----

calculer_LST_nuit <- function(chemin_fichier) {
  
  nom     <- basename(chemin_fichier)
  capteur <- ifelse(grepl("LC08", nom), "L8", "L9")
  date    <- substr(nom, 18, 25)
  date_f  <- format(as.Date(date, "%Y%m%d"), "%d/%m/%Y")
  
  cat("\n--- Traitement :", nom, "---\n")
  cat("    Capteur :", capteur, "| Date :", date_f, "\n")
  
  # ÉTAPE 1 — Localisation du fichier MTL associé
  # Le fichier MTL est dans le même dossier que B10 et porte
  # le même préfixe scène (17 premiers caractères = LXSS_LLLL_PPPRRR_YYYYMMDD)
  prefixe_scene <- substr(nom, 1, 40)  # conserve assez de caractères pour être unique
  prefixe_scene <- sub("_B10\\.TIF$", "", nom)  # supprime le suffixe _B10.TIF
  
  chemin_mtl <- list.files(
    dirname(chemin_fichier),
    pattern    = paste0(prefixe_scene, ".*MTL\\.txt$"),
    full.names = TRUE
  )
  
  if (length(chemin_mtl) == 0) {
    # Certaines livraisons utilisent l'extension _MTL.txt ou _MTL.TXT
    chemin_mtl <- list.files(
      dirname(chemin_fichier),
      pattern    = paste0(prefixe_scene, ".*MTL.*\\.txt$"),
      full.names = TRUE,
      ignore.case = TRUE
    )
  }
  
  if (length(chemin_mtl) == 0) stop(paste("Fichier MTL introuvable pour :", nom))
  chemin_mtl <- chemin_mtl[1]
  cat("    MTL     :", basename(chemin_mtl), "\n")
  
  # ÉTAPE 2 — Lecture des paramètres radiométriques dans le MTL
  # Paramètres Radiance TOA (formule USGS : Lλ = Ml * Qcal + Al)
  Ml <- lire_param_MTL(chemin_mtl, "RADIANCE_MULT_BAND_10")
  Al <- lire_param_MTL(chemin_mtl, "RADIANCE_ADD_BAND_10")
  cat("    Ml (RADIANCE_MULT_BAND_10) :", Ml, "\n")
  cat("    Al (RADIANCE_ADD_BAND_10)  :", Al, "\n")
  
  # Constantes thermales pour la température de brillance TOA
  # BT = K2 / ln((K1 / Lλ) + 1)  — en Kelvins
  K1 <- lire_param_MTL(chemin_mtl, "K1_CONSTANT_BAND_10")
  K2 <- lire_param_MTL(chemin_mtl, "K2_CONSTANT_BAND_10")
  cat("    K1 (K1_CONSTANT_BAND_10)   :", K1, "\n")
  cat("    K2 (K2_CONSTANT_BAND_10)   :", K2, "\n")
  
  # LECTURE RASTER B10 + DÉCOUPE SUR EMPRISE
  
  dn <- rast(chemin_fichier)  # Qcal : comptes numériques bruts
  
  ref_reproj <- if (crs(raster_ref) != crs(dn)) {
    cat("    Reprojection de l'emprise...\n")
    project(raster_ref, crs(dn))
  } else {
    raster_ref
  }
  
  dn_crop <- crop(dn, ext(ref_reproj))
  cat("    Découpe effectuée\n")
  
  # ÉTAPE 3 — Radiance TOA  (Albut 2024 / USGS)
  # Lλ = Ml * Qcal + Al
  # Lλ en W/(m² · sr · µm)
  
  L_toa <- (Ml * dn_crop) + Al
  
  # ÉTAPE 4 — Température de brillance TOA (en Kelvins)
  # BT = K2 / ln((K1 / Lλ) + 1)
  
  BT_K <- K2 / log((K1 / L_toa) + 1)
  
  # ÉTAPE 5 — LST avec émissivité constante (Albut 2024)
  # St = BT / (1 + (λ * BT / ρ) * ln(ε))
  # Résultat en Kelvins, puis conversion en °C
  
  LST_K <- BT_K / (1 + (LAMBDA * BT_K / RHO) * log(EPSILON))
  LST_C <- LST_K - KELVIN
  
  # CONTRÔLE QUALITÉ
  
  vals <- values(LST_C, na.rm = TRUE)
  cat("    Min :", round(min(vals),  1), "°C |",
      "Max :", round(max(vals),  1), "°C |",
      "Moyenne :", round(mean(vals), 1), "°C\n")
  cat("    Pixels valides :", length(vals), "\n")
  
  # EXPORT
  
  nom_sortie    <- gsub("_B10\\.TIF$", "_LST_celsius.TIF", nom)
  chemin_sortie <- file.path(DOSSIER_LST, nom_sortie)
  chemin_temp   <- file.path(DOSSIER_LST, gsub("_B10\\.TIF$", "_LST_temp.TIF", nom))
  
  # Écriture initiale
  writeRaster(LST_C, chemin_sortie,
              overwrite = TRUE,
              datatype  = "FLT4S",
              NAflag    = -9999,
              gdal      = c("COMPRESS=LZW",
                            "PREDICTOR=2",
                            "NODATA=-9999"))
  
  # Reconstruction avec gdal_translate (pour avoir les stats internes dans le raster)
  system(paste(
    "gdal_translate",
    "-of GTiff",
    "-a_nodata -9999",
    "-co COMPRESS=LZW",
    "-co PREDICTOR=2",
    "-stats",
    chemin_sortie,
    chemin_temp
  ))
  
  # Remplacement
  file.remove(chemin_sortie)
  file.rename(chemin_temp, chemin_sortie)
  
  cat("    Exporté :", nom_sortie, "\n")
  
  # Vérification finale
  lst_verif <- rast(chemin_sortie)
  cat("    Vérif fichier — Mean :",
      round(global(lst_verif, "mean", na.rm = TRUE)[[1]], 2), "°C |",
      "NAflag :", NAflag(lst_verif), "\n")
  
  tibble(
    fichier  = nom,
    date     = date_f,
    Ml       = Ml,
    Al       = Al,
    K1       = K1,
    K2       = K2,
    min_LST  = round(min(vals),  2),
    max_LST  = round(max(vals),  2),
    mean_LST = round(mean(vals), 2)
  )
}

# ----TRAITEMENT DE TOUS LES FICHIERS---

resultats <- map_dfr(fichiers, possibly(calculer_LST_nuit, otherwise = tibble()))

cat("\n=== Résultats finaux ===\n")
print(resultats)
