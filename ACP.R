library(terra)
library(ggplot2)
library(tidyverse)
library(readxl)
library(FactoMineR) # pour ACP
library(factoextra)

# ------PRÉPARATION DU TABLEAU POUR L'ACP ------

df <- read.csv("C:/Users/Spagnoli/Desktop/R/DF_jour1k")

conflicted::conflicts_prefer(dplyr::select)
# Variables continues pour l'ACP
vars_continues <- df %>%
  select(LST_pixel, LST_buffer, Ta_era5, RH975_era5, NDVI, NDBI, MNDWI, albedo, MNT, TPI, 
         dist_ocean, slope, aspect, SVF, longitude) %>%
  na.omit()

# Variables illustratives qualitatives (saison + LCZ + station)
# Elles doivent être dans le même tableau, après les variables actives
df_acp <- df %>%
  select(LST_pixel, LST_buffer, Ta_era5, RH975_era5, NDVI, NDBI, MNDWI, albedo, MNT, TPI, 
         dist_ocean, slope, aspect, SVF, longitude, saison, date, nom) %>%
  na.omit() %>%
  mutate(
    date = as.factor(date),
    saison = as.factor(saison),
    nom = as.factor(nom)
  )

# ==============================================================================
# DISTRIBUTION DES VARIABLES PAR STATION
# ==============================================================================

# Moyennes par station pour chaque variable
df_stations <- df %>%
  group_by(nom) %>%
  summarise(
    LST_pixel = mean(LST_pixel, na.rm = TRUE),
    Ta_era5   = mean(Ta_era5, na.rm = TRUE),
    RH975_era5 = mean(RH975_era5, na.rm = TRUE),
    NDVI      = mean(NDVI,      na.rm = TRUE),
    NDBI      = mean(NDBI,      na.rm = TRUE),
    MNDWI     = mean(MNDWI,     na.rm = TRUE),
    albedo    = mean(albedo,    na.rm = TRUE),
    MNT       = mean(MNT,       na.rm = TRUE),
    TPI       = mean(TPI,       na.rm = TRUE),
    dist_ocean  = mean(dist_ocean,  na.rm = TRUE),
    slope     = mean(slope,  na.rm = TRUE),
    aspect    = mean(aspect,  na.rm = TRUE),
    SVF       = mean(SVF,  na.rm = TRUE),
    longitude = mean(longitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(nom)

# GRAPHIQUES PAR VARIABLE

variables <- c("LST_pixel", "Ta_era5", "RH975_era5", "NDVI", "NDBI", "MNDWI", "albedo", "MNT", 
               "TPI", "dist_ocean", "slope", "aspect", "SVF", "longitude")

for (var in variables) {
  
  p <- df_stations %>%
    arrange(.data[[var]]) %>%
    mutate(nom = factor(nom, levels = nom)) %>%  # ordre croissant
    ggplot(aes(x = nom, y = .data[[var]], fill = .data[[var]])) +
    geom_col() +
    geom_text(aes(label = round(.data[[var]], 2)),
              hjust = -0.1, size = 3) +
    scale_fill_gradient(low = "#99CCFF", high = "#FF9966") +
    coord_flip() +
    labs(
      title = paste("Distribution de", var, "par station"),
      x     = NULL,
      y     = var
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p)
}

# TABLEAU RÉCAPITULATIF TRIÉ PAR VARIABLE

for (var in variables) {
  cat("\n=== Classement des stations par", var, "===\n")
  df_stations %>%
    arrange(desc(.data[[var]])) %>%
    select(nom, all_of(var)) %>%
    print()
}

# ==============================================================================
# 4. ACP
# ==============================================================================
# L'ACP : les colonnes actives sont utilisées via le df_acp, les dernières sont illustratives qualitatives
acp <- PCA(
  df_acp,
  scale.unit  = TRUE,  # centrage-réduction indispensable (unités différentes)
  ncp         = 5,   # nombre de composantes conservées
  quali.sup   = c(16, 17, 18),   # indices des colonnes saison, date et station 
  graph       = FALSE
)

# ==============================================================================
# 5. VISUALISATIONS
# ==============================================================================

# --- Variance expliquée par composante
fviz_eig(acp, addlabels = TRUE, ylim = c(0, 60),
         main = "Variance expliquée par composante",
         barfill = "#99CCFF",
         barcolor = "#99CCFF")

# --- Cercle des corrélations (les prédicteurs)
fviz_pca_var(
  acp,
  col.var = "contrib",        # couleur selon contribution
  gradient.cols = c("#99CCFF", "#FF9966"),
  repel = TRUE,
  title = "Cercle des corrélations"
)

# Seulement les individus colorés par saison (sans les variables)
fviz_pca_ind(
  acp,
  habillage     = df_acp$saison,
  addEllipses   = TRUE,
  ellipse.level = 0.95,
  geom.ind      = "point",    # points seulement, sans labels
  title         = "Individus par saison"
)

# Seulement les individus colorés par date (sans les variables)
fviz_pca_ind(
  acp,
  habillage     = df_acp$date,
  addEllipses   = TRUE,
  ellipse.level = 0.95,
  geom.ind      = "point",    # points seulement, sans labels
  title         = "Individus par date"
)

# Seulement les individus colorés par LCZ
fviz_pca_ind(
  acp,
  habillage     = df_acp$LCZ,
  addEllipses   = TRUE,
  ellipse.level = 0.95,
  geom.ind      = "point",
  title         = "Individus par LCZ"
)

# Seulement les individus colorés par station
fviz_pca_ind(
  acp,
  habillage     = df_acp$nom,
  addEllipses   = TRUE,
  ellipse.level = 0.95,
  geom.ind      = "point",
  title         = "Individus par station"
)

# ==============================================================================
# 6. RÉSULTATS CHIFFRÉS
# ==============================================================================

# Corrélations variables / composantes
print(acp$var$cor)

# Contributions des variables aux axes
print(acp$var$contrib)

# Cos2 : qualité de représentation des variables
print(acp$var$cos2)

# Matrice de corrélation entre prédicteurs -> ce qui guide si on garde ou pas les variables
cor_matrix <- cor(df_acp[, 1:12], use = "complete.obs")
print(round(cor_matrix, 2)) # si 2 variables ont correlation > 0.85 on n'en garde qu'une 

library(corrplot)

corrplot(cor_matrix,
         method = "color",          # cases colorées
         type = "upper",            # triangle supérieur seulement (évite la redondance)
         order = "hclust",          # regroupe les variables corrélées ensemble
         addCoef.col = "black",     # affiche les valeurs numériques
         number.cex = 0.7,          # taille des chiffres
         tl.col = "black",          # couleur des labels
         tl.srt = 45,               # angle des labels
         col = colorRampPalette(c("#99CCFF", "#FF9966"))(200),
         diag = FALSE)              # masque la diagonale (toujours 1)
