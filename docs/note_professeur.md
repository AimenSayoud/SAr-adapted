# Note d'avancement — mesure du déplacement vertical de la tourbière de Rzecin

**Objet :** état d'avancement, résultat principal, et demande d'accès aux
données de terrain archivées (laser/WTD, DoD drone).

---

## 1. Objectif

Mesurer le déplacement vertical (subsidence / gonflement) de la surface de la
tourbière de Rzecin (fen transitionnel à tapis de Sphagnum flottant, ~90 ha)
sur 2022–2024, par télédétection, en s'inspirant d'une méthode publiée sur un
site polonais comparable (Ghezelayagh et al. 2024, *Ecological Indicators*,
Biebrza) et de la méthode de référence sur tourbières (Hrysiewicz et al. 2024,
*Remote Sensing of Environment*).

## 2. Ce qui a été réalisé (pipeline complet, reproductible, sur GitHub)

- Chaîne InSAR Sentinel-1 (bande C) complète : acquisition, interférogrammes
  HyP3, réseau, masques d'eau Sentinel-2, classification, inversions.
- **Trois approches d'inversion testées :** (a) SBAS dense classique (MintPy) ;
  (b) ISBAS sur mesure (pixels intermittents) ; (c) paires saisonnières-
  annuelles optimisées par données atmosphériques ERA5 ; (d) **réseau hybride
  court + annuel** (méthode Hrysiewicz), inversé conjointement.
- Correction troposphérique ERA5 ; demande GACOS en cours.
- Validations internes sans terrain : cohérence par saison, fermeture de
  chaîne annuelle, corrélation à un proxy hydrologique, contrôle Sentinel-2.

## 3. Résultat principal (honnête)

**L'InSAR Sentinel-1 en bande C, seul, ne fournit pas de déplacement vertical
fiable sur le cœur de la tourbière**, quelle que soit la méthode :

- Le réseau hybride « résout » 14 238 pixels, mais après filtre de qualité
  (résidu d'inversion + nombre de paires), **0 pixel sur les 564 de la
  tourbière est fiable** : le résidu d'inversion y est de 2.0–2.9 radians
  (≈ 18–26 mm de dispersion par interférogramme) — les phases ne forment pas
  une série temporelle cohérente.
- Les paires annuelles se contredisent (taux 2022→2023 = −24 mm/an,
  2023→2024 = +32 mm/an ; défaut de fermeture 11 mm/an) → erreurs de
  déroulement de phase.
- La cohérence reste modérée mais insuffisante (été : 0.30 ; hiver : 0.37) et
  n'est que faiblement liée à l'eau de surface (r = −0.17) : la décorrélation
  vient surtout de la **diffusion de volume par la végétation** sur tout le
  tapis, pas seulement du lac.

**Ce résultat n'est pas un échec de traitement : c'est une limite physique
documentée** (bande C + végétation humide + petit site), cohérente avec la
littérature. Il constitue en soi un résultat publiable.

## 4. Interprétation et question ouverte

Deux causes possibles au résidu élevé, **que seule une donnée de terrain peut
départager** :

1. **décorrélation / erreurs de déroulement réelles** (le signal est perdu) ;
2. **mouvement réversible du tapis flottant** (la surface monte/descend avec la
   nappe) : dans ce cas la phase peut être correcte, mais une inversion en
   déplacement cumulé donne un résidu élevé *par construction*.

Distinguer les deux change tout : dans le cas 2, l'InSAR contiendrait un vrai
signal de « respiration » exploitable ; dans le cas 1, non.

## 5. Ce dont j'ai besoin pour trancher et publier

Pour lever cette ambiguïté et transformer ce travail en étude publiable, il me
faut les données de terrain archivées du site :

- **Laser / capteur d'élévation de surface (et niveau de nappe, WTD)** — même
  un seul point suffit pour le test décisif : la série InSAR au pixel du laser
  reproduit-elle la respiration mesurée ? (méthode validée par Hrysiewicz avec
  des caméras in situ). Le WTD sert aussi de proxy physique de la nappe.
- **Modèles numériques de surface drone (DTM/DSM) 2022 et 2024** — pour un
  changement net d'élévation (DoD) directement mesuré, indépendant de l'InSAR ;
  ce sera la mesure verticale principale, l'InSAR devenant une covariable.

**Questions concrètes :**
- Ces archives (laser/WTD continu, vols drone 2022–2024) sont-elles
  disponibles ? Sous quel format, quelles dates exactes, quelles coordonnées
  du point laser ?
- Les vols drone sont-ils du **LiDAR** (pénètre la végétation → sol) ou de la
  **photogrammétrie** (voit la canopée) ? C'est déterminant pour l'exactitude.
- Existe-t-il des repères géodésiques (GNSS/nivellement) sur le site ?

## 6. Le plan satellite d'ici là (validé)

En attendant les données de terrain, le volet satellite se poursuit de façon
utile et autonome :

- appliquer la correction atmosphérique **GACOS** (demandée) et quantifier son
  effet ;
- utiliser Sentinel-1/2 comme **capteurs de l'état hydrologique** (cohérence,
  rétrodiffusion, humidité/inondation) plutôt que de déformation — un rôle
  bien étayé par la littérature (modélisation de la nappe par S1/S2) ;
- améliorer le proxy de nappe avec **ERA5-Land (humidité du sol)** via Google
  Earth Engine, meilleur que la précipitation cumulée ;
- produire les cartes de caractérisation (cohérence saisonnière, dynamique
  d'inondation) qui serviront de contexte à la fusion finale.

**Architecture cible** (dès réception des données terrain) : *laser = modèle
temporel de la respiration ; drone = vérité spatiale du changement net ;
Sentinel-1/2 = covariable hydrologique spatio-temporelle pour corriger et
étendre.* C'est cette fusion qui donnera le déplacement vertical net défendable.

---

*Dépôt du projet (code, notebooks, rapports, revues) : AimenSayoud/SAr-adapted.*
