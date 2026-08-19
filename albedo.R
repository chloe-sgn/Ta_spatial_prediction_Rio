# Calcul de l'albédo de surface à partir d'images Sentinel-2 (Niveau 2A)
# Méthode : Bonafoni & Sekertekin (2020)

library(terra)
library(sf)
library(ggplot2)
library(viridis)
library(tidyterra)

# Dossier racine contenant tous les fichiers .jp2 Sentinel-2
# (le script cherche récursivement dans les sous-dossiers)
dossier_sentinel <- "C:/Users/Spagnoli/Desktop/QGIS/albedo/sec"

# Raster LST de référence définit l'emprise, la résolution et le CRS cibles
chemin_lst <- "C:/Users/Spagnoli/Desktop/DATA/LANDSAT/landsat_nuit/LST_060225.tif"

# Répertoire de sortie
dir_sortie <- "C:/Users/Spagnoli/Desktop/QGIS/albedo/sec/"
dir.create(dir_sortie, showWarnings = FALSE, recursive = TRUE)

# Facteur d'échelle L2A : DN / 10000 = réflectance de surface
facteur_echelle <- 10000

# ------DÉTECTION DES FICHIERS PAR PATTERN------
# Format attendu : B02/B03/B11/B12 → _20m.jp2 | B04/B08 → _10m.jp2
# Les 4 tuiles de chaque bande doivent être dans dossier_sentinel
# (peu importe la profondeur de l'arborescence, recursive = TRUE)

message("=== Détection des fichiers Sentinel-2 ===")

bandes_B02 <- list.files(
  path       = dossier_sentinel,
  pattern    = "^T[0-9]{2}[A-Z]{3}_[0-9]{8}T[0-9]{6}_B02_20m\\.jp2$",
  full.names = TRUE,
  recursive  = TRUE
)
bandes_B03 <- list.files(
  path       = dossier_sentinel,
  pattern    = "^T[0-9]{2}[A-Z]{3}_[0-9]{8}T[0-9]{6}_B03_20m\\.jp2$",
  full.names = TRUE,
  recursive  = TRUE
)
bandes_B04 <- list.files(
  path       = dossier_sentinel,
  pattern    = "^T[0-9]{2}[A-Z]{3}_[0-9]{8}T[0-9]{6}_B04_10m\\.jp2$",
  full.names = TRUE,
  recursive  = TRUE
)
bandes_B08 <- list.files(
  path       = dossier_sentinel,
  pattern    = "^T[0-9]{2}[A-Z]{3}_[0-9]{8}T[0-9]{6}_B08_10m\\.jp2$",
  full.names = TRUE,
  recursive  = TRUE
)
bandes_B11 <- list.files(
  path       = dossier_sentinel,
  pattern    = "^T[0-9]{2}[A-Z]{3}_[0-9]{8}T[0-9]{6}_B11_20m\\.jp2$",
  full.names = TRUE,
  recursive  = TRUE
)
bandes_B12 <- list.files(
  path       = dossier_sentinel,
  pattern    = "^T[0-9]{2}[A-Z]{3}_[0-9]{8}T[0-9]{6}_B12_20m\\.jp2$",
  full.names = TRUE,
  recursive  = TRUE
)

#Affichage des fichiers trouvés
cat("\n  Bandes B02 trouvées :", length(bandes_B02), "\n")
for (f in bandes_B02) cat("    -", basename(f), "\n")

cat("\n  Bandes B03 trouvées :", length(bandes_B03), "\n")
for (f in bandes_B03) cat("    -", basename(f), "\n")

cat("\n  Bandes B04 trouvées :", length(bandes_B04), "\n")
for (f in bandes_B04) cat("    -", basename(f), "\n")

cat("\n  Bandes B08 trouvées :", length(bandes_B08), "\n")
for (f in bandes_B08) cat("    -", basename(f), "\n")

cat("\n  Bandes B11 trouvées :", length(bandes_B11), "\n")
for (f in bandes_B11) cat("    -", basename(f), "\n")

cat("\n  Bandes B12 trouvées :", length(bandes_B12), "\n")
for (f in bandes_B12) cat("    -", basename(f), "\n")

# ------CHARGEMENT DU RASTER LST DE RÉFÉRENCE------

message("\n=== Raster LST de référence ===")
lst_ref <- rast(chemin_lst)

cat(sprintf("  CRS        : %s\n", crs(lst_ref, describe = TRUE)$name))
cat(sprintf("  Résolution : %.1f m\n", res(lst_ref)[1]))
cat(sprintf("  Dimensions : %d lignes x %d colonnes\n", nrow(lst_ref), ncol(lst_ref)))
cat(sprintf("  Emprise    : xmin=%.4f  xmax=%.4f  ymin=%.4f  ymax=%.4f\n",
            ext(lst_ref)[1], ext(lst_ref)[2], ext(lst_ref)[3], ext(lst_ref)[4]))

# ------FONCTION : MOSAÏQUE + REPROJECTION + ALIGNEMENT SUR LST------

construire_mosaique <- function(nom_bande, liste_fichiers, lst_reference,
                                est_scl   = FALSE,
                                natif_10m = FALSE,   # TRUE pour B04 et B08
                                facteur   = 10000) {
  
  message(sprintf("\n--- Bande %s (%d tuile(s)) ---", nom_bande, length(liste_fichiers)))
  methode <- ifelse(est_scl, "near", "bilinear")
  
  # Charger et reprojeter chaque tuile vers le CRS du LST
  tuiles <- lapply(seq_along(liste_fichiers), function(i) {
    message(sprintf("  Chargement tuile %d : %s", i, basename(liste_fichiers[i])))
    r <- rast(liste_fichiers[i])
    if (!same.crs(r, lst_reference)) {
      r <- project(r, crs(lst_reference), method = methode)
    }
    return(r)
  })
  
  # Assembler en mosaique
  # fun = "first" : sur les zones de chevauchement entre tuiles,
  # on prend la valeur de la premiere tuile (evite les artefacts de moyenne)
  if (length(tuiles) == 1) {
    mosaique <- tuiles[[1]]
  } else {
    message("  Assemblage mosaique...")
    mosaique <- do.call(mosaic, c(tuiles, list(fun = "first")))
  }
  
  # Agregation 10m -> 20m (uniquement pour B04 et B08)
  # aggregate() avec fun = "mean" : chaque pixel 20m = moyenne des 4 pixels
  # 10m qui le composent. Methode correcte car elle conserve l energie spectrale
  # (pas d interpolation, pas de perte d info).
  if (natif_10m) {
    message("  Agregation 10m -> 20m (moyenne des pixels 2x2)...")
    mosaique <- aggregate(mosaique, fact = 2, fun = "mean", na.rm = TRUE)
    cat(sprintf("    Resolution apres agregation : %.1f m\n", res(mosaique)[1]))
  }
  
  # Recadrer sur l emprise du LST
  mosaique <- crop(mosaique, ext(lst_reference))
  
  # Reechantillonner sur la grille exacte du LST
  # -> meme resolution, meme origine, meme nombre de pixels que le LST
  message("  Alignement sur la grille LST...")
  mosaique <- resample(mosaique, lst_reference, method = methode)
  
  # Conversion DN -> reflectance de surface [0, 1]
  if (!est_scl) {
    mosaique <- mosaique / facteur
    mosaique <- clamp(mosaique, 0, 1)
  }
  
  names(mosaique) <- nom_bande
  return(mosaique)
}

# ------CONSTRUCTION DES MOSAÏQUES------

message("\n=== Construction des mosaïques ===")

rho2  <- construire_mosaique("B02", bandes_B02, lst_ref)
rho3  <- construire_mosaique("B03", bandes_B03, lst_ref)
rho4  <- construire_mosaique("B04", bandes_B04, lst_ref, natif_10m = TRUE)
rho8  <- construire_mosaique("B08", bandes_B08, lst_ref, natif_10m = TRUE)
rho11 <- construire_mosaique("B11", bandes_B11, lst_ref)
rho12 <- construire_mosaique("B12", bandes_B12, lst_ref)

# Vérification de l'alignement
cat("\n=== Vérification de l'alignement spatial ===\n")
for (nm in c("lst_ref", "rho2", "rho3", "rho4", "rho8", "rho11", "rho12")) {
  r <- get(nm)
  cat(sprintf("  %-8s : %d x %d px | res=%.1f m\n",
              nm, nrow(r), ncol(r), res(r)[1]))
}

# ------CALCUL DE L'ALBÉDO — Bonafoni & Sekertekin (2020)------
# α = 0.2266·ρ₂ + 0.1236·ρ₃ + 0.1573·ρ₄ + 0.3417·ρ₈ + 0.1170·ρ₁₁ + 0.0338·ρ₁₂

message("\n=== Calcul de l'albédo (Bonafoni & Sekertekin 2020) ===")

albedo <- 0.2266 * rho2 +
  0.1236 * rho3 +
  0.1573 * rho4 +
  0.3417 * rho8 +
  0.1170 * rho11 +
  0.0338 * rho12

albedo <- clamp(albedo, 0, 1)
names(albedo) <- "albedo_bonafoni"

# Statistiques descriptives
vals <- values(albedo, na.rm = TRUE)
cat(sprintf("  Min     : %.4f\n", min(vals)))
cat(sprintf("  Max     : %.4f\n", max(vals)))
cat(sprintf("  Moyenne : %.4f\n", mean(vals)))
cat(sprintf("  Médiane : %.4f\n", median(vals)))
cat(sprintf("  Écart-t : %.4f\n", sd(vals)))
cat(sprintf("  N pixels: %d\n",   length(vals)))

# ------SAUVEGARDE------

message("\n=== Sauvegarde ===")

chemin_sortie <- file.path(dir_sortie, "albedo_sec.tif")

writeRaster(albedo,
            chemin_sortie,
            overwrite = TRUE,
            datatype  = "FLT4S",
            gdal      = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER"))

message("  Fichier : ", chemin_sortie)
