# Vérification des faits & challenge scientifique complet

**Objet.** Revenir aux faits, vérifier chaque affirmation dans le dépôt et la
littérature, puis soumettre notre conclusion actuelle (« Rzecin est un site
extrêmement difficile ») à réfutation active. Étiquetage : [VÉRIFIÉ] contrôlé
dans le dépôt/littérature · [NON VÉRIFIÉ] non accessible · [INT] interprétation
· [HYP] hypothèse.

---

## 1. Superficie de l'AOI — vérification

- [VÉRIFIÉ] La valeur « 89.7 ha » n'est **pas une donnée officielle** : elle est
  **calculée** par `aoi.area_ha()` à partir du polygone du fichier
  `data/aoi/Rzecin_corrected.geojson` (exporté d'un KML ArcGIS).
- [VÉRIFIÉ] Deux méthodes indépendantes donnent la même aire pour CE polygone :
  équirectangulaire locale = 89.74 ha ; géodésique ENU = 89.69 ha. Donc
  **89.7 ha est l'aire réelle du tracé que nous utilisons**, pas une erreur de
  calcul.
- [VÉRIFIÉ] L'attribut `Shape_Area = 0.90` du geojson est en unités incohérentes
  (degrés², non convertibles ici) — **non fiable**, à ignorer.
- [VÉRIFIÉ] Le polygone **n'a aucun trou intérieur** : le **lac résiduel est
  INCLUS** dans l'AOI. Idem pour d'éventuelles marges d'eau libre.

**Conséquence sur « 84.7 ha vs 89.7 ha » :** l'écart (~5 ha, ~6 %) ne se règle
PAS en remplaçant le texte — le geojson lui-même calcule 89.7. Deux
possibilités :
1. [HYP] le tracé KML est légèrement trop large (inclut lac + marges) par
   rapport à la délimitation officielle de la tourbière (84.7 ha) ;
2. [HYP] « 84.7 ha » désigne une autre entité (réserve, tourbe active,
   empreinte eddy-covariance) que ce contour.

**Action recommandée (avant toute correction chiffrée) :** obtenir la source
autoritative du 84.7 ha (quel document/désignation ?) et, si possible, le
**shapefile officiel de la tourbière**. Alors : remplacer le geojson, relancer
le masquage, et **toutes les valeurs se mettront à jour automatiquement**
(elles dérivent toutes de `area_ha`/`aoi_mask`). Corriger « 89.7 → 84.7 » dans
les textes sans changer le polygone créerait une incohérence code/rapport.
→ **En attente de la délimitation officielle. Ne PAS figer 84.7 tant que le
polygone n'est pas remplacé.**

*(Note d'hygiène : le dossier `insar-wetlands-pipeline-mqd0qv/notebooks/`
contient 3 notebooks dupliqués obsolètes suivis par git — à supprimer pour
éviter d'éditer la mauvaise copie.)*

## 2. Étude Hrysiewicz et al. (2024, RSE 291) — fiche vérifiée

**[VÉRIFIÉ] (résultats de recherche, résumé + notices) :**
- Sites : **Cors Fochno (264 ha)** et **Cors Caron** (plus vaste), mid-Wales,
  Royaume-Uni.
- Type : **bogs ombrotrophes de plaine tempérée (raised bogs)**, à Sphagnum,
  **en cours de restauration** (donc historiquement drainés/perturbés).
- Régime hydrologique : nappe sous la surface, alimentation par les pluies
  (ombrotrophe) ; PAS de grand plan d'eau libre central.
- Période Sentinel-1 : **mi-2015 → début-2023 (~8 ans)**.
- Réseau temporel : **combinaison de lignes de base LONGUES ET courtes**
  (point-clé explicite de l'article).
- Validation : **caméras in situ** sub-mm, mi-2019 → mi-2022.
- Performance : Pearson **0.8–0.9** ; 76 % des écarts < ±5 mm ; 93 % < ±10 mm.
- Vitesses : Cors Caron subsidence max **−7 mm/an** ; Cors Fochno **−9 à
  +5 mm/an** (subsidence au centre, soulèvement aux marges).

**[NON VÉRIFIÉ] (article payant, preprint/ESA renvoient 403 via le proxy) :**
logiciel exact (SNAP/ISCE/GAMMA/LiCSAR/StaMPS/DS-InSAR ?), **polarisation**
(VV ? VV+VH ?), **orbite(s)** (asc/desc, numéros), **multilooking/résolution**,
seuils de cohérence, **méthode de correction atmosphérique** (GACOS ?),
**stratégie de point de référence**. Ces points sont déterminants et doivent
être obtenus (accès institutionnel à l'article, ou contact des auteurs).

## 3. Comparaison rigoureuse Rzecin / Cors Fochno-Caron / Biebrza

| Facteur | Rzecin (nous) | Cors Fochno/Caron (Hrysiewicz) | Biebrza (Ghezelayagh) |
|---|---|---|---|
| Type | Fen transitionnel/pauvre, tapis Sphagnum flottant | Raised bogs ombrotrophes, Sphagnum | Fen minerotrophe drainé |
| Végétation | Sphagnum + laîches | **Sphagnum** (similaire !) | herbacée/roselière |
| Eau libre | Lac résiduel INCLUS dans l'AOI | pas de grand plan d'eau | canaux de drainage |
| Superficie | ~85–90 ha | **264 ha (Fochno)**, Caron + grand | 96 610 ha |
| État | quasi naturel (+ drainage partiel) | perturbés, en restauration | drainé actif |
| Signal net attendu | quasi nul | modéré (−7 à +5 mm/an) | fort (−14 mm/an) |
| Période S1 | **3 ans (2022–24)** | **~8 ans (2015–23)** | 7 ans (2015–22) |
| Satellite | **S1A seul, 12 j** (S1B en panne) | S1A+S1B sur ~6 ans, 6 j | S1A+S1B |
| Validation terrain | aucune (à ce jour) | caméras in situ | piézomètres |
| Réseau | court ≤48 j, puis hybride +3 annuelles | **long+court sur 8 ans** | annuel avril→avril |

**[INT] Lecture :** la végétation (Sphagnum) est **SIMILAIRE** entre Rzecin et
les bogs de Hrysiewicz — cela **affaiblit** l'idée d'un site « uniquement
difficile par la végétation ». Les différences réellement défavorables chez
nous sont : **période 2,7× plus courte** (3 vs 8 ans), **échantillonnage 2×
plus pauvre** (S1A seul 12 j vs S1A+B 6 j), **AOI 3× plus petite** (moins de
moyennage), **eau libre incluse**, et un **réseau long-baseline testé
tardivement** (3 paires annuelles seulement). Ce sont majoritairement des
différences de **données/méthode**, pas de physique intrinsèque du site.

## 4. Challenge : toutes les explications possibles, classées

Objectif : essayer de PROUVER que « le site est extrêmement difficile » est
faux. On liste les explications de l'observation-clé (**résidu d'inversion
2–2.9 rad, 0 pixel fiable sur l'AOI**), rangées par probabilité.

**H1 — Notre pipeline n'a pas testé les méthodes adaptées aux faibles
cohérences (le plus probable, le plus sous-estimé).**
Nous avons utilisé HyP3 (interférogrammes GAMMA, looks 10×2 figés) + MintPy +
ISBAS maison. **Jamais testé : DS-InSAR / phase-linking (SqueeSAR-like),
LiCSBAS, StaMPS, ni un traitement SNAP/ISCE sur mesure (filtrage Goldstein,
multilooking optimisé).** Le phase-linking est PRÉCISÉMENT conçu pour les
cibles distribuées à cohérence faible (végétation, tourbe) — c'est la classe
de méthode que les études modernes utilisent. [INT] Une partie du résidu
élevé peut être des erreurs de déroulement / un traitement sous-optimal, pas
une perte physique.
*Pour/contre :* pour — méthode designée pour ça, non testée ; contre — la
cohérence brute reste basse (0.30 été), qu'aucun logiciel n'augmente.
*Expérience discriminante :* **rejouer le stack en DS-InSAR/phase-linking
(ex. via ISCE2+FRInGE, ou LiCSBAS)** ; si le résidu chute et le signal
apparaît, c'était la méthode.

**H2 — L'AOI inclut de l'eau libre qui dégrade les statistiques (partiel).**
[VÉRIFIÉ] Le lac résiduel est dans l'AOI. [INT] Mais le contrôle qualité est
PAR PIXEL : un pixel-lac mauvais ne rend pas un pixel-tourbe voisin mauvais.
La corrélation cohérence–eau S2 était faible (r=−0.17) → l'eau n'est pas le
décorrélateur principal ; la décorrélation est étalée sur tout le tapis.
*Expérience :* masquer précisément le lac + marges, recalculer la qualité sur
le **cœur Sphagnum seul** (lié au point 84.7 ha : le vrai cœur est peut-être
plus petit et plus cohérent).

**H3 — Le site est génuinement difficile RELATIVEMENT à notre traitement
(plausible, mais non isolé).**
Cohérence 0.30–0.41 étalée sur le tapis, faiblement liée à l'eau → décorrélation
de volume par la végétation partout. [INT] Mais « difficile » est relatif :
Hrysiewicz réussit sur Sphagnum à des cohérences comparables **avec un meilleur
traitement et 8 ans de données**. Donc H3 n'est établi que si H1 est écarté.

**H4 — Le mouvement est réel mais réversible/non-linéaire → résidu élevé par
construction (non départageable de H3 sans terrain).**
Le tapis flottant bouge avec la nappe ; une inversion en déplacement CUMULÉ a
un résidu élevé même avec une phase correcte. *Expérience :* ajuster un modèle
**saisonnier+tendance** au lieu d'incréments cumulés — si le résidu chute, le
mouvement est réversible, pas du bruit ; **le laser tranche définitivement**.

**H5 — Correction atmosphérique insuffisante (partiel).** Rampes jusqu'à
179 mm/an, seulement deramp planaire, **GACOS pas encore appliqué**. [INT]
Mais l'atmosphère est spatialement lisse : elle explique le bruit de tendance
des paires annuelles, **pas** un résidu inter-paires de 2.5 rad par pixel.
*Expérience :* appliquer GACOS (en cours), quantifier.

**H6 — Point de référence instable sur site flottant (partiel).** Référence
jamais validée. [INT] Mais un mauvais point de référence décale TOUT
uniformément — il n'inflige pas un résidu de fit par pixel. Affecte
l'exactitude de la tendance, pas le résidu. *Expérience :* DS-InSAR/PSI sur
diffuseurs stables + ancrage GNSS.

### Conclusion du challenge (révision honnête)

[INT] **Nous ne pouvons PAS conclure « Rzecin est extrêmement difficile ».**
Trois explications restent vivantes et non départagées :
- **H1** (pipeline incomplet — surtout l'absence de DS-InSAR/phase-linking) ;
- **H3** (site difficile *relativement à notre traitement et nos 3 ans*) ;
- **H4** (mouvement réversible → mésappariement de modèle).

Les différences avec les études qui réussissent sont **majoritairement de
données et de méthode** (8 ans vs 3, S1A+B vs S1A, 264 ha vs 90, DS-InSAR vs
non, long-baseline systématique vs tardif) — **pas** une preuve de supériorité
de difficulté physique intrinsèque. La végétation Sphagnum, elle, est
similaire à un site qui réussit.

## 5. Vérification de l'AOI (notre découpage)

- [VÉRIFIÉ] GeoJSON unique `Rzecin_corrected.geojson`, tracé KML, 89.7 ha, **lac
  inclus**, sommets simplifiés pour les requêtes (aoi_wkt).
- [HYP] Le contour est probablement **trop large et trop inclusif** : il englobe
  le lac et des marges. Un découpage sur le **cœur Sphagnum** (excluant eau
  libre + laîches de bordure) donnerait une AOI plus petite (possiblement
  proche de 84.7 ha) et **plus cohérente** — ce qui pourrait changer le verdict
  qualité.
- [INT] Notre AOI influence donc DIRECTEMENT la performance observée : on
  évalue la qualité sur une zone qui inclut des pixels structurellement
  incohérents (eau). *Action :* produire une **AOI « cœur »** (masque
  eau-permanente + classe Sphagnum) et recomparer.

## 6. Comparaison écologique fondée sur la littérature

- [VÉRIFIÉ, littérature site] Rzecin : poor/transitional fen, tapis Sphagnum
  flottant par terrestrialisation d'un lac, nappe 0–30 cm sous la surface,
  **stabilité hydrologique exceptionnelle** (Juszczak et al. 2013), quasi
  naturel mais influence de drainage régionale.
- [VÉRIFIÉ] Cors Fochno/Caron : raised bogs ombrotrophes à Sphagnum, en
  restauration (donc perturbés). Cors Fochno montre subsidence-centre /
  soulèvement-marges — un patron de « respiration » spatialement structuré.
- **Ressemblances réellement importantes pour l'InSAR :** dominance **Sphagnum**
  (même régime de décorrélation de volume en bande C) ; climat tempéré humide ;
  échelle de quelques centaines d'ha ; mouvement de surface de type respiration.
- **Différences réellement importantes pour l'InSAR :** (i) **tapis flottant
  libre** à Rzecin (mouvement potentiellement plus ample et plus réversible →
  H4) vs bog « respirant » mais non flottant ; (ii) **lac résiduel** (eau libre
  centrale) absent des raised bogs ; (iii) **taille** (90 vs 264+ ha) ;
  (iv) **durée/échantillonnage** (3 ans S1A vs 8 ans S1A+B).
- **Différences probablement PEU importantes pour l'InSAR :** trophie
  (minérotrophe vs ombrotrophe) en soi ; statut « naturel vs restauration »
  (les deux ont de la subsidence/respiration mesurable).
- [INT] Biebrza est le MOINS bon analogue (fen drainé, signal fort, immense) ;
  **Cors Fochno/Caron sont les meilleurs analogues** — et l'InSAR y marche.
  Cela déplace le poids de la preuve vers « données/méthode » plutôt que
  « site intrinsèquement impossible ».

## 7. Matrice d'expériences pour départager (priorisées)

| # | Expérience | Départage | Coût |
|---|---|---|---|
| 1 | **DS-InSAR / phase-linking** (ISCE2+FRInGE ou LiCSBAS) sur le stack existant | H1 (méthode) vs H3 (site) | Moyen, données en main |
| 2 | **EGMS** (European Ground Motion Service) : vérifier la couverture Rzecin — produit S1 européen déjà optimisé, indépendant | H1 (notre pipeline) : un tiers voit-il un signal ? | **Faible, immédiat** |
| 3 | **AOI cœur** (masque eau + Sphagnum) + recalcul qualité | H2/AOI + lien 84.7 ha | Faible, code en main |
| 4 | Modèle **saisonnier+tendance** vs incréments cumulés | H4 (réversible) vs H3 (bruit) | Faible, code en main |
| 5 | **GACOS** (demandé) | H5 (atmosphère) | En cours |
| 6 | **Laser au pixel** (dès réception) | H3 vs H4, test décisif | dépend du terrain |
| 7 | Obtenir les **détails de traitement de Hrysiewicz** (accès article) | cadre la comparaison H1 | Faible (biblio) |

**[INT] Les deux prochaines actions les plus rentables, sans terrain :**
(2) **vérifier EGMS sur Rzecin** — si le service européen, avec son traitement
PSI/SBAS optimisé, ne voit rien non plus, cela renforce fortement « site
difficile » de façon indépendante de notre pipeline ; s'il voit un signal,
c'était notre méthode ; et (1) **tester le phase-linking/DS-InSAR**, la classe
de méthode que nous n'avons jamais essayée et qui est conçue pour ce cas.

---

### Synthèse en une phrase
La vérification factuelle **affaiblit** notre conclusion « site extrêmement
difficile » : l'analogue le plus proche (Sphagnum, Cors Fochno) réussit en
bande C, et nos désavantages sont surtout de **données et de méthode** (3 ans
vs 8, S1A vs S1A+B, pas de DS-InSAR, AOI trop large incluant le lac) — trois
explications (méthode incomplète, site-relatif, mouvement réversible) restent
ouvertes et se départagent par EGMS + phase-linking + AOI-cœur + laser, pas par
davantage d'affirmations.
