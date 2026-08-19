# ----Calcul LST sur Landsat 8 et 9 lvl 2 le jour---

library(terra)      # Pour manipuler les rasters
library(tidyverse)  # Pour la manipulation des données

DOSSIER_B10 <- "C:/Users/Spagnoli/Desktop/DATA/LANDSAT/landsat_nuit"
DOSSIER_LST <- "C:/Users/Spagnoli/Desktop/DATA/LANDSAT/LST_n"
RASTER_EMPRISE <- "C:/Users/Spagnoli/Desktop/DATA/LANDSAT/LST_n/LST_060225.tif"
dir.create(DOSSIER_LST, showWarnings = FALSE, recursive = TRUE)

#constantes pour calculs 
SCALE_FACTOR <- 0.00341802 # issu du calcul de ST de https://www.usgs.gov/landsat-missions/landsat-collection-2-surface-temperature
OFFSET <- 149.0 # issu du calcul de ST de https://www.usgs.gov/landsat-missions/landsat-collection-2-surface-temperature
LAMBDA <- 10.895 # longueur d'onde moy de la bande 10   
RHO <- 14388 # en uK
EPSILON <- 0.97 # constante d'emissivité (Amorim & Dubreuil, 2016) 
KELVIN <- 273.15 # pour convertir en °C

# ----CHARGEMENT DE L'EMPRISE---

cat("=== Chargement de l'emprise ===\n") # fonction baseR : affiche texte dans console comme print(),
                                        #"\n" = retour à la ligne

raster_ref <- rast(RASTER_EMPRISE) # fonction terra : lit chemin raster et charge en mémoire comme objet SpatRaster  

cat("    CRS     :", crs(raster_ref, describe = TRUE)$name, "\n") # fonction terra : retourne syst coord du raster, 
                                                                # describe = TRUE retourne un df avec nom du CRS 

cat("    Étendue :", paste(round(as.vector(ext(raster_ref)), 4), collapse = " | "), "\n") # ext() fonction terra : 
                # retourne étendue spatiale du raster (xmin, xmax, ymin, ymax) utilisé pour définir zone de découpe et retourne un objet SpatRaster 
              # as.vector() fonction BaseR : convertit le SpatRatser en vecteur num de 4 valeurs (xmin, xmax, ymin, ymax)
            # round(,4) arrondit chaque valeur à 4 décimales 
          # paste(, collapse = " | ") colle les 4 valeurs séparées par | dans une chaîne de texte 

# ----LISTE DES FICHIERS----

fichiers <- list.files(DOSSIER_B10, 
                       pattern = "LC0[89].*ST_B10\\.TIF$", 
                       full.names = TRUE) # fonction baseR : liste tous les fichiers d'un dossier qui correspondent à un pattern

cat("=== Fichiers trouvés :", length(fichiers), "===\n") # fonction baseR : compte le nb d'elements dans la liste de fichiers
                                                    
for (f in fichiers) cat(" •", basename(f), "\n") # boucle for parcourt chaque fichier f de la liste de fichiers un par un
                                                # basename (chemin) fonction baseR: extrait juste nom du fichier 

# ----FONCTION CALCUL LST----

calculer_LST_crop <- function(chemin_fichier) { # function permet de créer sa propre fonction, ici : 1 fonction avec 1 argument (chemin_fichier)
  
  nom     <- basename(chemin_fichier) 
  capteur <- ifelse(grepl("LC08", nom), "L8", "L9") # grepl() cherche si le texte "LC08" est présent dans nom du fichier, retourne TRUE ou FALSE
                                                   # ifelse () : le TRUE ou FALSE retourné de grepl() est la condition,
                                                  # si ça retourne TRUE la valeur est L8, si FALSE la valeur est L9 
  
  date    <- substr(nom, 18, 25)   # extrait la date du nom de fichier (de la 18 à la 25 position des caractères du nom du fichier)
  
  date_f  <- format(as.Date(date, "%Y%m%d"), "%d/%m/%Y") # as.Date fonction baseR : convertit date en vrai objet date R 
                                                        # format () passe de "20230118" à "18/01/2023"
  
  cat("\n--- Traitement :", nom, "---\n")
  cat("    Capteur :", capteur, "| Date :", date_f, "\n")
  
  # Lecture du raster
  dn <- rast(chemin_fichier) #dn = bande 10 en comptes num 
  
  # Reprojection de l'emprise dans le CRS du raster Landsat si nécessaire
  ref_reproj <- if (crs(raster_ref) != crs(dn)) { # fonction terra : retourne syst coord du raster ref 
    cat("    Reprojection de l'emprise...\n")
    project(raster_ref, crs(dn)) # fonction terra : reprojette dn dans syst coord raster_ref 
  } else {
    raster_ref
  } # bloc if {} else {} retourne une valeur dans ref_reproj : si syst coord sont ≠ -> ref_reproj reçoit raster reproj 
  # si syst coord sont identiques -> ref_reproj reçoit raster original sans modif   
  
  # Découpe sur l'emprise
  dn_crop <- crop(dn, ext(ref_reproj)) # ext() extrait l'étendue spatiale du raster de ref et crop ( crop () : fonction terra) le dn à partir du ref  
  cat("    Découpe effectuée\n")
  
  # Étape 1 — Température de brillance (Kelvins)
  BT <- (dn_crop * SCALE_FACTOR) + OFFSET # BT = (B10 * 0.00341802) + 149
  
  # Étape 2 — LST (Kelvins) avec émissivité constante  
  LST_K <- BT / (1 + (LAMBDA * BT / RHO) * log(EPSILON)) # LSTk = BT / (1 + (10.895 * BT / 14388) * ln (0.97))
  
  # Étape 3 — Conversion en °C
  LST_C <- LST_K - KELVIN # LSTc = LSTk - 273.15
  
  # Vérification plage de valeurs
  vals      <- values(LST_C, na.rm = TRUE) # values(raster, na.rm) fonction terra : extrait toutes valeurs raster sous vecteur num R 
                                          # permet d'utiliser toutes les fonctions stats de base 
                                        # na.rm = TRUE exclut les pixels NA du vecteur 
  cat("    Min :", round(min(vals), 1), "°C |",
      "Max :", round(max(vals), 1), "°C |",
      "Moyenne :", round(mean(vals), 1), "°C\n") # calculs des valeurs min, max, et mean, arrondies avec round()
  cat("    Pixels valides :", length(vals), "\n")  # longueur des valeurs, étape contrôle
  
# -----EXPORT — avec modif pour obtenir stats du raster----
  
  nom_sortie    <- gsub("ST_B10\\.TIF$", "LST_celsius.TIF", nom) # fonction baseR : remplace le pattern avec le nouveau nom (pattern, remplacement, texte)  
  chemin_sortie <- file.path(DOSSIER_LST, nom_sortie) # fonction baseR : construit un chemin avec dossier et nom de fichier  
  chemin_temp   <- file.path(DOSSIER_LST, gsub("ST_B10\\.TIF$", "LST_temp.TIF", nom)) # fichier temporaire le temps que gdal_translate construise le fichier final avec les stats
  # writeRaster va écrire le fichier pixel par pixel en flux continu. À la fin de l'écriture des pixels, au lieu de refaire une passe complète pour connaitre les stats, 
  # il écrit -9999 dans les stats de la bande interne (il ne peut faire qu'une passe).
  # Il n'existe pas de fonction capable de lire et modifier le même fichier en même temps. 
  # Donc on va utiliser gdal_translate qui lit tous les pixels écrits par writeRaster depuis chemin_sortie puis recalcule les stats dans chemin_temp    
  
  # Étape 1 — Écriture initiale
  writeRaster(LST_C, chemin_sortie,  
              overwrite = TRUE,
              datatype  = "FLT4S", #32bit 
              NAflag    = -9999,
              gdal      = c("COMPRESS=LZW",
                            "PREDICTOR=2", #predictor pour 32bit
                            "NODATA=-9999"))
  
  # Étape 2 — Reconstruction avec gdal_translate 
  system(paste(
    "gdal_translate", #gdal_translate reconstruit le fichier 
    "-of GTiff",
    "-a_nodata -9999", #ignore les NA -9999
    "-co COMPRESS=LZW",
    "-co PREDICTOR=2",
    "-stats", #recalcule les stats internes de la bande 
    chemin_sortie,
    chemin_temp
  ))
  
  # Étape 3 — Remplacer l'original par le fichier corrigé
  file.remove(chemin_sortie) # gdal_translate ne peut pas écraser le fichier de writeRaster, donc il le supprime
  file.rename(chemin_temp, chemin_sortie) # renomme le fichier de chemin_sortie avec celui de chemin_temp
  
  cat("    Exporté :", nom_sortie, "\n")
  
  # Vérification finale du fichier exporté
  lst_verif <- rast(chemin_sortie)
  cat("    Vérif fichier — Mean :",
      round(global(lst_verif, "mean", na.rm = TRUE)[[1]], 2), "°C |", # global() fonction terra : calcule moy pixels en ignorant NA. Retourne df (1l / 1c)
                # [[1]] extrait le 1er element du DF retourné par global() pour obtenir que une valeur de moy et pas un DF complet
              # round(..., 2) arrondit à 2 décimales
      "NAflag :", NAflag(lst_verif), "\n") # fonction terra : lit valeur NA du raster. Retourne -9999 si bien défini, NaN si absent 
  
  # Retourne un résumé
  tibble(
    fichier  = nom,
    date     = date_f,
    min_LST  = round(min(vals),  2),
    max_LST  = round(max(vals),  2),
    mean_LST = round(mean(vals), 2)
  ) # fonction tidyverse : créé DF moderne. Retourne résumé par image qui sera collé avec map_dfr() 
}

# ----TRAITEMENT DE TOUS LES FICHIERS----
# map_dfr(liste, fonction) fonction tidyverse : applique une fonction à chaque élément d'une liste et colle les résultats en 1 DF
# équivalent boucle for 
# possibly(fonction, otherwise) fonction tidyverse : évite de faire planter script si calcul_LST ne fonctionne pas pour 1 raster avec otherwise = tibble vide
resultats <- map_dfr(fichiers, possibly(calculer_LST_crop, otherwise = tibble())) 
