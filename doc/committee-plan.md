## Plan comité-conseil - Paragraphe par paragraphe

### Introduction

1. Le calcul des LCBD permet l’identification des sites qui recèlent une diversité
   exceptionnelle, mais l’utilisation de la méthode est actuellement limitée à échelle
   spatiale discontinue.

1. Les modèles de distribution d’espèces permettent de faire des prédictions à échelle
   continue, sur lesquelles il serait possible de calculer les LCBD.

1. Les bases de données BIOCLIM et eBird sont fréquemment utilisées dans des SDM, et
   pourraient donc servir au calcul des LCBD à échelle spatiale continue.

1. L’approche prédictive suggérée permettrait de pallier au manque d’informations et
   d’échantillonnages dans certaines régions.

1. La méthode pourrait être appliquée aux scénarios de changements climatiques du GIEC pour
   identifier les sites où la diversité bêta risque de changer de façon importante.

### Méthodes

1. Récolte des données : Les données disponibles ont été récoltées pour l’Amérique du Nord et
   pour la famille des Parulines.

1. Traitement des données : Les données ont été converties en présence-absence par pixel,
   pour une résolution de 10 arc-minutes (puis la transformation de Hellinger a été appliquée
   temporairement, mais la transformation appropriée reste à déterminer).

1. SDM : La méthode BIOCLIM se base sur les quantiles des variables environnementales pour
   effectuer les prédictions.

1. LCBD : Les LCBD ont été calculées par la variance totale de la matrice Y, et non par une
   matrice de dissimilarité.

1. Validité des résultats : La façon exacte pour évaluer la validité des prédictions et de la
   méthode reste à déterminer.

1. Méthodes alternatives : Certaines méthodes alternatives sont envisagées pour les SDM, ex.
   Random Forests, réseaux de neurones, etc.

1. Scénarios changements climatiques & diversité bêta temporelle : Nous souhaitons appliquer
   la méthode aux conditions environnementales prévues par les scénarios de changements
   climatiques, puis utiliser les indices temporels de diversité bêta (TBI) pour identifier
   les sites avec les changements les plus importants

1. Résultats préliminaires : Les résultats préliminaires permettent de comparer les données
   brutes aux prédictions faites par les SDM. Les figures présentées porteront sur :
   - Distribution des espèces individuelles
   - Richesse
   - LCBD
   - Relation richesse-LCBD
