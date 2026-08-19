# DIAGRAMME OMBROTHERMIQUE — Rio de Janeiro
# Normale climatique sur PERIODE PERSONNALISABLE (ex : 1990-2020, 2000-2020)
# Source : "Historical monthly weather data" WorldClim
#          (CRU-TS-4.09 recalé avec WorldClim 2.1, couverture 1950-2024)
# Pas compte ni API mais téléchargement automatique

library(terra)
library(ggplot2)

annee_debut <- 1990   # ex : 1990, 2000...
annee_fin   <- 2020   # ex : 2020
resolution  <- "10"   # "2.5", "5" ou "10" (minutes d'arc). 10 = le plus léger/rapide
rio_lon     <- -43.17
rio_lat     <- -22.91
dossier_dl  <- "./worldclim_hist_data"

dir.create(dossier_dl, showWarnings = FALSE)

# Augmenter le délai de téléchargement (les zips font 100-300+ Mo selon la résolution)
options(timeout = 1200)

# ---- 1. Déterminer les décennies nécessaires ----
# Les données sont livrées par blocs de 10 ans : 1950-1959, ..., 2010-2019, 2020-2024
# (le dernier bloc n'a que 5 ans -> il faut le lister séparément, pas fusionné avec 2010-2019)
debuts_decennies <- c(seq(1950, 2010, 10), 2020)
fins_decennies    <- debuts_decennies + 9
fins_decennies[length(fins_decennies)] <- 2024  # dernier bloc : 2020-2024

bornes_decennies <- data.frame(
  debut = debuts_decennies,
  fin   = fins_decennies
)
bornes_decennies$label <- paste0(bornes_decennies$debut, "-", bornes_decennies$fin)

decennies_utiles <- bornes_decennies[
  bornes_decennies$fin >= annee_debut & bornes_decennies$debut <= annee_fin,
]

cat("Décennies à télécharger pour la période", annee_debut, "-", annee_fin, ":\n")
print(decennies_utiles$label)

# ---- 2. Téléchargement + dézippage des décennies nécessaires ---
base_url <- "https://geodata.ucdavis.edu/climate/worldclim/2_1/hist/cts4.09/"

telecharger_decennie <- function(var, label_decennie) {
  fname <- paste0("wc2.1_cruts4.09_", resolution, "m_", var, "_", label_decennie, ".zip")
  url   <- paste0(base_url, fname)
  dest  <- file.path(dossier_dl, fname)
  
  # motif des tifs déjà extraits pour cette décennie (évite un re-téléchargement)
  prefixe_annee <- substr(label_decennie, 1, 3)
  motif <- paste0("wc2.1_cruts4.09_", resolution, "m_", var, "_", prefixe_annee, ".*\\.tif$")
  deja_la <- list.files(dossier_dl, pattern = motif, full.names = TRUE)
  
  if (length(deja_la) == 0) {
    message("Téléchargement : ", fname, " (peut prendre plusieurs minutes selon la résolution)")
    ok <- tryCatch({
      download.file(url, dest, mode = "wb", quiet = FALSE)
      TRUE
    }, error = function(e) {
      message("  -> échec : ", conditionMessage(e))
      if (file.exists(dest)) file.remove(dest)  # nettoyer un fichier partiel
      FALSE
    })
    if (ok) {
      utils::unzip(dest, exdir = dossier_dl)
      file.remove(dest)
    }
  }
  invisible(NULL)
}

for (var in c("tmin", "tmax", "prec")) {
  for (lbl in decennies_utiles$label) {
    telecharger_decennie(var, lbl)
  }
}

# ---- 3. Extraction des valeurs mensuelles pour Rio, année par année-----
rio_pt   <- matrix(c(rio_lon, rio_lat), nrow = 1)
rio_vect <- terra::vect(rio_pt, crs = "EPSG:4326")

extraire_serie <- function(var) {
  fichiers <- list.files(
    dossier_dl,
    pattern = paste0("^wc2\\.1_cruts4\\.09_", resolution, "m_", var, "_\\d{4}-\\d{2}\\.tif$"),
    full.names = TRUE
  )
  
  infos <- data.frame(
    fichier = fichiers,
    annee   = as.integer(sub(".*_(\\d{4})-(\\d{2})\\.tif$", "\\1", fichiers)),
    mois    = as.integer(sub(".*_(\\d{4})-(\\d{2})\\.tif$", "\\2", fichiers))
  )
  
  # ne garder que la période demandée
  infos <- infos[infos$annee >= annee_debut & infos$annee <= annee_fin, ]
  
  r <- terra::rast(infos$fichier)
  val <- terra::extract(r, rio_vect)[, -1]
  val <- as.numeric(val)
  
  # si NA (point pile en bordure de pixel) -> repli sur un buffer de 50 km
  if (any(is.na(val))) {
    message("Correction NA (", var, ") — extraction par buffer 50 km...")
    val <- as.numeric(terra::extract(r, rio_vect, buffer = 50000,
                                     fun = mean, na.rm = TRUE)[, -1])
  }
  
  infos$valeur <- val
  infos
}

serie_tmin <- extraire_serie("tmin")
serie_tmax <- extraire_serie("tmax")
serie_prec <- extraire_serie("prec")

# ---- 4. Moyenne mensuelle sur la période (= la "normale" perso)-----
moyenne_mensuelle <- function(serie) {
  tapply(serie$valeur, serie$mois, mean, na.rm = TRUE)
}

tmin_moy <- moyenne_mensuelle(serie_tmin)
tmax_moy <- moyenne_mensuelle(serie_tmax)
prec_moy <- moyenne_mensuelle(serie_prec)   # moyenne des totaux mensuels

tavg_moy <- (tmin_moy + tmax_moy) / 2   # tavg n'existe pas dans ce jeu -> recalculé

normale <- data.frame(
  mois       = 1:12,
  mois_label = factor(1:12, levels = 1:12,
                      labels = c("Jan","Fév","Mar","Avr","Mai","Jun",
                                 "Jul","Aoû","Sep","Oct","Nov","Déc")),
  temp_moy   = as.numeric(tavg_moy[as.character(1:12)]),
  temp_max   = as.numeric(tmax_moy[as.character(1:12)]),
  temp_min   = as.numeric(tmin_moy[as.character(1:12)]),
  precip_tot = as.numeric(prec_moy[as.character(1:12)])
)

cat("\nNormale climatique Rio de Janeiro (", annee_debut, "-", annee_fin, ") :\n", sep = "")
print(normale[, c("mois_label", "temp_min", "temp_moy", "temp_max", "precip_tot")])

write.csv(normale,
          paste0("./normale_Rio_", annee_debut, "_", annee_fin, ".csv"),
          row.names = FALSE)

# ---- 5. Diagramme ombrothermique ----
# Convention de Gaussen : un mois est considéré "sec" quand P <= 2T
# -> l'axe des précipitations est fixé à 2x l'échelle de température
coeff <- 2

# (ancienne version : coefficient dynamique pour que les 2 courbes tiennent
#  bien dans le graphique, sans respecter la convention de Gaussen)
# coeff <- max(normale$precip_tot, na.rm = TRUE) / max(normale$temp_moy, na.rm = TRUE)

diagramme <- ggplot(normale, aes(x = mois_label)) +
  
  geom_bar(
    aes(y = precip_tot / coeff, fill = "Précipitations (mm)"),
    stat = "identity", width = 0.6, alpha = 0.75
  ) +
  geom_ribbon(
    aes(ymin = temp_min, ymax = temp_max, group = 1,
        fill = "Amplitude Tn-Tx"),
    alpha = 0.25
  ) +
  geom_line(
    aes(y = temp_moy, colour = "Température moy. (°C)", group = 1),
    linewidth = 1.2
  ) +
  geom_point(
    aes(y = temp_moy, colour = "Température moy. (°C)"),
    size = 3, shape = 21, fill = "white", stroke = 1.5
  ) +
  scale_y_continuous(
    name     = "Température (°C)",
    sec.axis = sec_axis(~ . * coeff, name = "Précipitations (mm)")
  ) +
  scale_fill_manual(values = c(
    "Précipitations (mm)" = "#4A90D9",
    "Amplitude Tn-Tx"      = "#E8A0B4"
  )) +
  scale_colour_manual(values = c(
    "Température moy. (°C)" = "#E84040"
  )) +
  labs(
    title    = "Diagramme ombrothermique — Rio de Janeiro",
    subtitle = paste0("Normale climatique ", annee_debut, "-", annee_fin,
                      " (CRU-TS 4.09 recalé WorldClim 2.1)"),
    x        = "Mois",
    caption  = "Source : CRU-TS 4.09 (Harris et al., 2020) ; downscaling WorldClim 2.1 (Fick & Hijmans, 2017)",
    fill     = NULL, colour = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(colour = "grey40", size = 11),
    axis.title.y       = element_text(colour = "#E84040", face = "bold"),
    axis.title.y.right = element_text(colour = "#4A90D9", face = "bold"),
    axis.text.y        = element_text(colour = "#E84040"),
    axis.text.y.right  = element_text(colour = "#4A90D9"),
    legend.position    = "bottom",
    panel.grid.minor   = element_blank()
  )

print(diagramme)
