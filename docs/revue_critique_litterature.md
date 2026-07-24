# Revue critique et confrontation à la littérature

**Objet.** Soumettre nos propres conclusions à réfutation, confronter chaque
affirmation importante à la littérature internationale (succès *et* échecs),
expliciter les hypothèses implicites et les biais, et assigner à chaque
conclusion un **niveau de confiance** avec les expériences nécessaires pour
la confirmer ou l'infirmer. Standard visé : RSE / ISPRS J. / IEEE TGRS.

Convention d'étiquetage utilisée partout :
**[FAIT]** observé/mesuré · **[HYP]** hypothèse · **[INT]** interprétation ·
**[SPEC]** spéculation.

---

## 0. Résumé des révisions imposées par la littérature

| Conclusion de notre rapport | Verdict après revue | Confiance révisée |
|---|---|---|
| « L'InSAR C-band ne PEUT PAS mesurer le vertical sur tourbière » | **RÉFUTÉE comme énoncé général** | Faible → à retirer |
| « Notre échec est d'origine physique (bande C) » | **CONFONDUE** avec nos choix de méthode | Faible-modérée |
| « Le centre de phase sur végétation inondée n'est pas le sol » | **SOUTENUE** (wetland-InSAR) mais nuancée | Élevée |
| « Le taux net à Rzecin est ≈ 0 » | Plausible mais **non démontré** (SNR) | Modérée |
| « S1/S2 = covariable hydrologique/nappe » | **FORTEMENT SOUTENUE** | Élevée |
| « La bande L est la solution physique » | Soutenue, mais bénéfice à quantifier | Modérée-élevée |

Le message : notre récit « échec = limite physique de la bande C » est
**partiellement un biais de confirmation**. Des équipes réussissent l'InSAR
C-band sur tourbière. Nous n'avons pas isolé si notre échec vient du **site**
ou de notre **méthode**.

---

## 1. La conclusion la plus dangereuse : « C-band ne peut pas »

### Ce que dit notre rapport
Que l'incompatibilité bande C + fen végétalisé rend le déplacement vertical
inobservable par Sentinel-1.

### Contre-preuves directes dans la littérature
- **[FAIT] Hrysiewicz et al. (2024, *Remote Sensing of Environment* 291,
  S0034425724002505)** : InSAR Sentinel-1 **C-band** sur deux tourbières
  hautes tempérées (Cors Fochno, Cors Caron, Pays de Galles), validé par
  caméras in situ sub-millimétriques. Corrélation de Pearson **0.8–0.9**,
  **76 %** des écarts < ±5 mm, **93 %** < ±10 mm ; RMSE multi-annuel
  **~7 mm/an** (~3.5 mm sur mesures individuelles). Surtout : la méthode
  **récupère les oscillations annuelles de surface de 10–40 mm** — la
  « respiration » que NOUS avions qualifiée de simple bruit confondant.
- **[FAIT]** Études tourbières tropicales (SBAS Sentinel-1, seuil de
  cohérence 0.25) couvrant les 2/3 de la zone ; performance **comparable à
  l'ALOS L-band** (0.3 → 62 %) sur le même site.
- **[FAIT]** Nombreuses études SBAS C-band réussies sur tourbières
  dégradées/restaurées (Frontiers Env. Sci. 2023 ; IOP 2025).

### Différence méthodologique cruciale que nous avons ignorée
- **[FAIT]** Hrysiewicz combine **des réseaux à lignes de base temporelles
  LONGUES ET courtes**. **[FAIT]** Notre configuration n'a utilisé que des
  paires **≤ 48 jours** (`max_temporal_baseline_days: 48`). Nous n'avons
  **jamais** testé le réseau long-baseline qui fait le succès des autres.
- **[INT]** Une partie de notre « limite physique » pourrait n'être qu'une
  **limite de conception de réseau** — un choix, pas une fatalité.

### Différence de site (la nuance qui pourrait nous sauver)
- **[FAIT]** Les sites où l'InSAR C-band réussit sont des **bogs
  ombrotrophes** (nappe SOUS la surface, tapis non inondé). Rzecin est un
  **fen transitionnel minerotrophe** à **nappe affleurante/émergente**.
- **[HYP]** Le régime hydrologique (inondé vs non inondé), pas la bande,
  serait le vrai discriminant. Cette hypothèse est **testable** (§6) et
  n'est PAS démontrée à ce stade.

### Confiance
Notre affirmation initiale « C-band impossible » : **confiance faible, à
retirer**. Reformulation défendable (voir §3).

---

## 2. Le centre de phase et le « wetland InSAR » : notre point le plus solide, mais à recadrer

### Ce que dit notre rapport
Sur végétation inondée, le centre de phase est le double-bounce eau-tige ; il
suit le niveau d'eau et la biomasse, pas le sol tourbeux.

### Appui de la littérature
- **[FAIT]** Champ mature du **wetland InSAR** (Wdowinski & Hong ; Hong et
  al. 2014 TGRS ; Alsdorf ; Kim) : sur végétation émergente inondée, le
  double-bounce domine et l'InSAR mesure le **changement de niveau d'eau**,
  exploité comme signal utile.
- **[FAIT]** « Avec la végétation herbacée émergente, le double-bounce
  végétation-eau **augmente** la cohérence au fil de la saison de
  croissance » — ce qui **contredit** notre idée simple « inondation =
  décorrélation totale ».

### Nuance qui affaiblit notre version
- **[FAIT]** Le **C-band en VV** (notre polarisation) est **MOINS sensible
  au double-bounce** que le HH et les longues longueurs d'onde. Donc :
  - soit on veut le niveau d'eau → il aurait fallu **HH** (indisponible sur
    notre trace S1, VV/VH seulement) → **limitation instrumentale réelle** ;
  - soit on veut le sol → en VV on a surtout de la **diffusion de volume**
    (décorrélante), pas un double-bounce stable.
- **[INT]** Donc à Rzecin, en VV, sur zones inondées : ni bon signal de
  niveau d'eau (VV faible en double-bounce) ni bon signal de sol (volume
  décorrélant). C'est un argument **plus fin et plus défendable** que
  « le centre de phase n'est pas le sol » énoncé sèchement.

### Confiance
Le mécanisme (phase ≠ sol sur zones inondées) : **confiance élevée**. Mais
« donc bruit pur » est **trop fort** : c'est potentiellement un signal de
niveau d'eau mal capté en VV, pas du bruit. **Hypothèses implicites
corrigées** : (i) inondation ≠ décorrélation systématique ; (ii) tout le
fen n'est pas inondé en permanence.

---

## 3. Reformulation défendable de notre « échec »

**[INT] Énoncé prudent, publiable :** *Sur ce fen transitionnel à nappe
affleurante, avec une pile Sentinel-1 VV mono-satellite à cycle 12 j
(2022–2024), un réseau small-baseline ≤ 48 j et un traitement MintPy
standard, nous n'avons pas extrait de taux vertical net fiable sur le cœur
de la tourbière. Deux causes NON encore séparées : (a) régime hydrologique
inondé/émergent déplaçant le mécanisme de diffusion ; (b) choix
méthodologiques (réseau uniquement court, VV, petite AOI ~564 px limitant le
moyennage spatial). Le taux net attendu (~mm/an pour un fen préservé) est de
plus **à la limite du plancher de bruit InSAR documenté (~3.5–7 mm/an)**.*

Cette version distingue **fait**, **cause site**, **cause méthode** et
**limite SNR** — au lieu de tout imputer à « la physique de la bande C ».

---

## 4. Le taux net ≈ 0 : plausible, pas démontré

### Ce que nous avons
- **[FAIT]** Pipeline B : AOI médiane −4.5 mm/an, **IQR ±18.9 mm/an**.

### Critique
- **[FAIT]** L'incertitude (18.9) écrase le signal (4.5). Statistiquement,
  **on ne peut pas distinguer −4.5 mm/an de 0, ni de −14 mm/an** (Biebrza).
  Notre « ≈ 0 » et notre « cohérent avec un fen préservé » sont des
  **[INT]**, pas des **[FAIT]**.
- **[HYP]** Biais possible : le deramp planaire par paire peut **retirer une
  partie d'un vrai signal** à grande longueur d'onde s'il existe. Nous
  supposons que le signal est concentré sur l'AOI et la rampe hors-AOI —
  non vérifié.
- **[SPEC]** Un fen minerotrophe alimenté par la nappe pourrait accumuler
  (gonfler) plutôt que s'affaisser — signe opposé à Biebrza. Nos données ne
  tranchent pas.

### Confiance
« Taux net non détectable avec cette méthode » : **confiance élevée**.
« Taux net physiquement ≈ 0 » : **confiance modérée à faible** — c'est une
absence de preuve, pas une preuve d'absence.

---

## 5. La reformulation S1/S2 = covariable hydrologique : notre point le mieux étayé

- **[FAIT]** « Modelling water table depth at rewetted peatlands with
  Sentinel-1 and Sentinel-2 » (2025) : la profondeur de nappe est
  modélisable directement depuis S1/S2. Appui **direct** à notre réorientation.
- **[FAIT]** Cohérence InSAR utilisée comme indice de végétation/humidité
  (Sentinel-1 coherence as a vegetation index, RSE 2022) et pour la
  cartographie d'inondation (nombreux travaux).
- **[FAIT]** La motion de surface tourbeuse est un **proxy reconnu du niveau
  de nappe** (Hrysiewicz ; Alshammari) — ce qui valide l'idée que la
  quantité pertinente couple surface et hydrologie.
- **[FAIT]** Rzecin lui-même a déjà été étudié en optique (PlanetScope
  2017–2023 ; EVI/phénologie vs nappe, MDPI *Remote Sensing* 2026),
  confirmant l'intérêt du suivi multi-capteurs sur ce site précis.

### Confiance
**Élevée.** C'est la partie la plus solide de notre stratégie — et,
ironiquement, celle que nous avions présentée comme un repli.

---

## 6. Expériences décisives pour lever les ambiguïtés

Chaque expérience vise à **séparer site vs méthode** et à **falsifier** une
de nos conclusions.

1. **Test de réplication Hrysiewicz (le plus important).** Reconstruire un
   réseau combinant lignes de base **longues (annuelles) ET courtes**, comme
   RSE 2024, au lieu de ≤ 48 j seul. → Si l'InSAR récupère alors la
   respiration, notre « limite physique » était une **limite de méthode**.
   *Falsifie ou confirme §1.*
2. **Test de la respiration au point laser (test décisif).** Le laser in
   situ donne l'oscillation vraie (façon caméras de Hrysiewicz). Comparer la
   série InSAR au pixel du laser : si l'InSAR récupère la respiration
   (corrélation élevée) même sans trend net fiable, alors **C-band
   fonctionne** ici et seul le trend est limité par le SNR. Sinon, le cœur
   du fen est **génuinement décorrélé**. *Discrimine §1 vs §3 vs §4.*
3. **Stratification de la cohérence.** Cartographier cohérence vs
   niveau de nappe / classe d'inondation / phénologie S2. → Teste
   l'hypothèse « inondé = décorrélé » (§2) et localise les zones exploitables.
4. **DS-InSAR / PSI sur diffuseurs stables.** Chercher des cibles stables
   (passerelles, instruments, marges plus sèches) pour ancrer une référence
   fiable — ce qui manquait à MintPy.
5. **Budget SNR quantitatif.** Estimer le plancher de bruit réel sur notre
   pile (dispersion sur zones stables) et le comparer au signal net attendu.
   → Établit si « ≈ 0 » est une détection ou une non-détection (§4).
6. **Validation croisée drone.** Le DoD 2022→2024, dé-aliasé par la nappe,
   fournit une vérité indépendante du signe et de l'ordre de grandeur du
   trend net — arbitre externe de §4.

---

## 7. Hypothèses implicites restantes (au-delà du rapport précédent)

- **[HYP]** Que « paires annuelles avril→avril » implique un état de nappe
  comparable. **Non vérifié** ; le laser/ERA5 doivent le contrôler.
- **[HYP]** Que le deramp planaire n'ampute pas de vrai signal long-onde.
- **[HYP]** Que la petite AOI (~564 px) suffit statistiquement ; les succès
  publiés portent sur 10³–10⁶ ha (moyennage bien supérieur).
- **[HYP]** Que VV est neutre. Il ne l'est pas pour un régime double-bounce ;
  HH aurait été préférable mais est indisponible sur cette trace.
- **[HYP]** Que 3 ans (2022–2024) suffisent pour un trend ; Hrysiewicz
  utilise ~2015–2023 (~8 ans), Biebrza 7 ans. Notre fenêtre est **courte**.

---

## 8. Biais méthodologiques identifiés dans notre propre démarche

- **Biais de confirmation** : après l'échec MintPy, nous avons cherché des
  raisons « physiques » de l'échec et sous-pondéré les succès publiés en
  C-band. La présente revue corrige ce déséquilibre.
- **Biais de disponibilité** : nous avons calqué la comparaison sur
  Biebrza (drainé) parce que c'était l'article en main, alors que les
  analogues pertinents sont les **bogs suivis par caméra** (Hrysiewicz).
- **Biais de méthode unique** : un seul design de réseau (court) testé sous
  MintPy ; pas d'ISCE/GAMMA/StaMPS/DS-InSAR en comparaison.
- **Sur-généralisation** : d'un échec sur *un* site/une config vers « la
  bande C ne peut pas ».

---

## 9. Conclusion de la revue

Notre travail expérimental et notre réorientation stratégique (drone/laser
comme mesure principale, S1/S2 comme covariable hydrologique) restent
**valides et bien étayés**. En revanche, l'**explication** que nous en
donnions — « limite physique de la bande C » — ne résiste pas à la
littérature : l'InSAR C-band mesure la respiration des tourbières avec une
précision de quelques mm. La formulation correcte, publiable, est que **nous
n'avons pas encore séparé la cause-site de la cause-méthode**, et que les
six expériences du §6 — au premier chef la réplication long+short baseline et
la comparaison au laser — sont nécessaires avant toute affirmation forte.

**Ce que cela change pour une soumission :** la contribution n'est pas
« l'InSAR échoue sur les fens » (invérifié et probablement faux), mais soit
(a) « conditions d'observabilité de l'InSAR C-band selon le régime
hydrologique de la tourbière (bog non inondé vs fen inondé) », étayée par les
expériences §6, soit (b) « architecture de fusion laser–drone–S1/S2 pour le
déplacement net corrigé de la nappe ». La (b) est la plus solide aujourd'hui.

---

### Références clés mobilisées

- Hrysiewicz A. et al. (2024). Estimation and validation of InSAR-derived
  surface displacements at temperate raised peatlands. *Remote Sensing of
  Environment* 291. (C-band, r=0.8–0.9, RMSE ~7 mm/an, respiration 10–40 mm.)
- Hong S.-H., Wdowinski S. (2014). Double-bounce component in
  cross-polarimetric SAR. *IEEE TGRS*. (Mécanisme double-bounce en zone humide.)
- Wdowinski S., Hong S.-H. Wetland InSAR: a review of the technique and
  applications. (Phase = niveau d'eau sur végétation inondée.)
- Ghezelayagh P. et al. (2024). *Ecological Indicators* 166:112305. (Biebrza,
  seasonal-annual, −1.44 cm/an — fen drainé, référence initiale.)
- « Modelling water table depth at rewetted peatlands with Sentinel-1 and
  Sentinel-2 » (2025). (Appui direct à la reformulation covariable.)
- « Sentinel-1 interferometric coherence as a vegetation index » (2022),
  *RSE*. (Cohérence comme covariable de surface.)
- Étude PlanetScope sur Rzecin (2017–2023), MDPI *Remote Sensing* 18(4):593,
  2026. (Site étudié en optique : EVI/phénologie vs nappe.)
- Revue « Remote sensing of peatland degradation in temperate and boreal
  zones » (2024), *Ecological Indicators*. (État de l'art succès/lacunes.)
