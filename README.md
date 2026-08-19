# Cartographie de la température de l'air à Rio de Janeiro

Ce projet vise à prédire spatialement la température de l'air (Ta) à Rio de 
Janeiro à partir de la température de surface (LST) issue de l'imagerie 
satellite, combinée à des variables environnementales et topographiques 
(NDVI, NDBI, MNDWI, albédo, MNT, pente, distance à la mer, Sky View Factor, etc.)
et des données de modèles climatiques globaux désagrégées (ERa5).

Dans un premier temps, l'étude explore l'évolution spatio-temporelle de la relation Ta–LST selon 
la typologie urbaine, le cycle jour/nuit et la saisonnalité (saison humide/sèche). Dans un second temps,
une méthode adaptée à notre jeu de données à été appliquée pour identifier les prédicteurs de la Ta
la nuit et le jour à Rio, puis pour prédire spatialement la Ta validée par une méthode LOSO.
Ensuite, l'analyse de l'îlot de chaleur urbain (ICU) a permis de caractériser ses variations 
spatio-temporelles à Rio.

## Contenu à exécuter dans l'ordre
- Script du calcul de la normale climatique de Rio 
- Scripts du calcul des LST de jour et de nuit
- Scripts de calcul des indices spectraux et autres variables (NDVI, NDBI, MNDWI, albédo; LCZ)
- Script de création du DataFrame à 1km
- Script de l'analyse de la relation Ta~LST selon la typologie
- Script de l'ACP pour visualiser les prédicteurs
- Script du modèle RF de prédiction spatiale de la Ta
- Script de l'analyse de l'ICU (tests de Wilcoxon appariés)
- Cartes et visualisations des résultats
