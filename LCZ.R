# ------ Créer un raster LCZ rapidement sur toute l'image de référence------
# Ce raster va nous aider notamment à faire le masque urbain|rural pour analyser l'ICU

library(LCZ4r)
library(terra)
library(sf)

# ------CHARGER LE RASTER LST DE RÉFÉRENCE------
lst_ref <- rast("C:/Users/Spagnoli/Desktop/DATA/LANDSAT/LST_n/LST_060225.tif")

# Vérifications
print(lst_ref)
cat("CRS     :", crs(lst_ref, describe = TRUE)$name, "\n")
cat("Résolution :", res(lst_ref), "\n")
cat("Étendue :", as.character(ext(lst_ref)), "\n")

# ------CONVERTIR L'EMPRISE DU LST EN SHAPEFILE (ROI pour LCZ4r)------
# LCZ4r attend un objet sf en WGS84 (EPSG:4326)
roi_sf <- st_as_sf(as.polygons(ext(lst_ref), crs = crs(lst_ref))) |>
  st_transform(4326)

# ------TÉLÉCHARGER LA CARTE LCZ SUR L'EMPRISE DU LST------
lcz_roi <- lcz_get_map(roi = roi_sf)

lcz_plot_map(lcz_roi, title = "LCZ – Emprise LST")

# ------RÉCUPÉRER LES PARAMÈTRES------
lcz_params <- lcz_get_parameters(lcz_roi)

# ------REPROJETER ET RÉÉCHANTILLONNER TOUT SUR LA GRILLE LST------
# Reprojection du raster LCZ vers le CRS du LST
lcz_reproj <- project(lcz_roi, crs(lst_ref), method = "near")

# Rééchantillonnage sur la grille exacte du LST
lcz_final <- resample(lcz_reproj, lst_ref, method = "near")

# Même chose pour les paramètres continus
params_continus <- c("SVFmean", "ISFmean", "TSFmean", "AHmean",
                     "BSFmean", "HREmean", "SALmean", "ARmean",
                     "SADmean", "z0")

params_reproj <- project(lcz_params[[params_continus]],
                         crs(lst_ref),
                         method = "bilinear")

params_final <- resample(params_reproj, lst_ref, method = "bilinear")

# Vérification alignement
cat("\n--- Vérification alignement ---\n")
cat("LST      :", res(lst_ref),   "| nrow:", nrow(lst_ref),   "ncol:", ncol(lst_ref),   "\n")
cat("LCZ      :", res(lcz_final), "| nrow:", nrow(lcz_final), "ncol:", ncol(lcz_final), "\n")
cat("Params   :", res(params_final), "| nrow:", nrow(params_final), "\n")

# ------EXPORT RASTER LCZ + QML------
output_dir <- "C:/Users/Spagnoli/Desktop/QGIS/LCZ_Rio_State"
params_dir <- file.path(output_dir, "params_lst_grid")
dir.create(params_dir, recursive = TRUE, showWarnings = TRUE)
stopifnot(dir.exists(params_dir))

# Export LCZ
tif_lcz <- file.path(output_dir, "lcz_lst_grid.tif")
qml_lcz <- file.path(output_dir, "lcz_lst_grid.qml") # pour légende
writeRaster(lcz_final, tif_lcz, overwrite = TRUE)

lcz_colors <- data.frame(
  value = 1:17,
  label = c("1 - Compact high-rise","2 - Compact midrise","3 - Compact low-rise",
            "4 - Open high-rise","5 - Open midrise","6 - Open low-rise",
            "7 - Lightweight low-rise","8 - Large low-rise","9 - Sparsely built",
            "10 - Heavy industry","A - Dense trees","B - Scattered trees",
            "C - Bush/scrub","D - Low plants","E - Bare rock/paved",
            "F - Bare soil/sand","G - Water"),
  color = c("#8c0000","#d10000","#ff0000","#bf4d00","#ff6600","#ff9955",
            "#faee05","#bcbcbc","#ffccaa","#555555","#006a00","#00aa00",
            "#648525","#b9db79","#000000","#fbf7ae","#6a6aff")
)

qml_items <- paste(sapply(1:nrow(lcz_colors), function(i) {
  sprintf('        <paletteEntry alpha="255" color="%s" label="%s" value="%d"/>',
          lcz_colors$color[i], lcz_colors$label[i], lcz_colors$value[i])
}), collapse = "\n")

writeLines(sprintf(
  '<!DOCTYPE qgis PUBLIC "http://mrcc.com/qgis.dtd" "SYSTEM">
<qgis version="3.0" styleCategories="AllStyleCategories">
  <pipe>
    <rasterrenderer alphaBand="-1" band="1" opacity="1" type="paletted">
      <rasterTransparency/>
      <colorPalette>
%s
      </colorPalette>
    </rasterrenderer>
  </pipe>
  <blendMode>0</blendMode>
</qgis>', qml_items), qml_lcz)

cat("LCZ exporté :", tif_lcz, "\n")

# Export paramètres
for (param in params_continus) {
  
  raster_layer <- params_final[[param]]
  tif_out <- file.path(params_dir, paste0("lcz_", param, ".tif"))
  qml_out <- file.path(params_dir, paste0("lcz_", param, ".qml"))
  
  writeRaster(raster_layer, tif_out, overwrite = TRUE)
  
  vals    <- values(raster_layer, na.rm = TRUE)
  min_val <- round(min(vals), 4)
  max_val <- round(max(vals), 4)
  
  writeLines(sprintf(
    '<!DOCTYPE qgis PUBLIC "http://mrcc.com/qgis.dtd" "SYSTEM">
<qgis version="3.0" styleCategories="AllStyleCategories">
  <pipe>
    <rasterrenderer band="1" classificationMax="%s" classificationMin="%s"
                    opacity="1" type="singlebandpseudocolor">
      <rasterTransparency/>
      <rastershader>
        <colorrampshader classificationMode="2" clip="0" colorRampType="INTERPOLATED">
          <item alpha="255" color="#440154" label="%s" value="%s"/>
          <item alpha="255" color="#31688e" label=""    value="%s"/>
          <item alpha="255" color="#35b779" label=""    value="%s"/>
          <item alpha="255" color="#fde725" label="%s" value="%s"/>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
  </pipe>
  <blendMode>0</blendMode>
</qgis>',
    max_val, min_val,
    min_val, min_val,
    round(min_val + (max_val - min_val) * 0.33, 4),
    round(min_val + (max_val - min_val) * 0.66, 4),
    max_val, max_val
  ), qml_out)
  
  cat("Exporté :", param, "\n")
}

cat("\nTout est aligné sur la grille LST et exporté dans :", output_dir, "\n")