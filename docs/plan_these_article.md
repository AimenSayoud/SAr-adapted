# Plan détaillé — Thèse / Article publiable

**Titre de travail (article) :**
*Vertical displacement of a near-natural temperate fen from C-band InSAR: why dense
SBAS time-series fail on small peatlands, and a robust seasonal-annual alternative
validated by UAV and in-situ laser measurements (Rzecin, Poland)*

**Titre de travail (thèse) :**
*Mesure du déplacement vertical de la tourbière de Rzecin par télédétection
multi-capteurs : limites physiques de l'InSAR en bande C et stratégie
d'observation robuste*

---

## Fil narratif en une phrase

> Nous avons tenté la voie standard (SBAS/MintPy) sur une petite tourbière
> quasi naturelle, documenté précisément **pourquoi et où elle échoue**
> (limites physiques de la bande C + hypothèses logicielles inadaptées aux
> zones humides), construit une alternative qui fonctionne (paires
> saisonnières-annuelles optimisées par ERA5, améliorant Ghezelayagh et al.
> 2024), et validé le tout par deux mesures indépendantes (DoD drone +
> laser in situ) — produisant à la fois un résultat de site et un guide
> méthodologique transférable.

C'est une structure « échec instructif → solution → validation » : l'échec
n'est pas caché, il **est** une partie du résultat. C'est ce qui rend l'étude
publiable même si le taux de déplacement final est proche de zéro (un fen
quasi naturel qui ne s'affaisse pas est un résultat cohérent et utile —
c'est le point bas du gradient hydrologique de Ghezelayagh et al., Fig. 6).

---

## 1. Introduction

1.1. Tourbières = stock de carbone ; le tassement (subsidence) comme
     indicateur intégrateur de dégradation (oxydation + retrait) ;
     émissions de CO₂ associées. [refs : Wösten 1997 ; Hooijer 2012 ;
     Couwenberg & Hooijer 2013 ; Hoyt 2020 ; Tiemeyer 2020]

1.2. L'InSAR comme solution au problème d'échelle des mesures terrain ;
     succès récents sur tourbières **grandes et drainées** :
     Ghezelayagh et al. 2024 (Biebrza, 96 610 ha, −1.44 cm/an, R²=0.7
     contre piézomètres) ; Hrysiewicz 2023 ; Hoyt 2020 ; Tampuu et al.

1.3. **Le trou dans la littérature** que nous attaquons :
     - les sites étudiés avec succès sont vastes (moyennage spatial massif)
       et fortement drainés (signal fort, 1–3 cm/an) ;
     - qu'en est-il d'un site **petit (~90 ha)**, **quasi naturel** (signal
       attendu faible, 0–1 cm/an), à végétation dense et surface
       intermittemment inondée ? C'est le cas de la majorité des tourbières
       protégées d'Europe — précisément celles que les programmes de
       restauration doivent surveiller.
     - les chaînes de traitement standard (SBAS/MintPy) ont été conçues
       pour des cibles cohérentes (roche, urbain) : leurs hypothèses
       implicites n'ont jamais été confrontées systématiquement à ce type
       de terrain.

1.4. Questions de recherche :
     - **QR1** : la chaîne SBAS dense standard (HyP3 + MintPy) peut-elle
       produire une série temporelle exploitable sur un petit fen quasi
       naturel en bande C ? Si non, à quel niveau exact se situe le blocage
       (données, algorithme, physique) ?
     - **QR2** : une stratégie « seasonal-annual » à petites paires
       optimisées (type Ghezelayagh 2024) peut-elle produire un taux
       vertical robuste sur ce même site, et peut-on l'améliorer en
       remplaçant la sélection météo manuelle par un critère ERA5
       quantitatif ?
     - **QR3** : que disent les mesures indépendantes (DoD drone
       2022–2024, laser in situ continu) du déplacement réel, et
       valident-elles l'InSAR ?
     - **QR4** : quelles recommandations générales pour le suivi InSAR des
       petites tourbières, aujourd'hui (bande C) et demain (bande L) ?

## 2. Site et données

2.1. **Rzecin** (Pologne, 52.76°N 16.31°E) : fen quasi naturel de ~89.7 ha,
     site éddy-covariance historique ; contexte hydrologique ; pourquoi ce
     site est un « cas difficile idéal ».

2.2. Données (tableau récapitulatif) :
     | Source | Période | Rôle |
     |---|---|---|
     | Sentinel-1 A, orbite 175 asc., burst 175_374052_IW1, VV | 2022–2024, 91 dates (cycle 12 j) | InSAR |
     | HyP3 `INSAR_ISCE_BURST`, looks 10×2, ~80 m | 346 paires ≤48 j + paires annuelles | interférogrammes prêts à l'emploi |
     | Sentinel-2 L2A | 2022–2024 | masque eau/inondation, NDVI, phénologie |
     | ERA5 (tcwv, t2m, tp) | 2022–2024 | score atmosphérique, gel, hydrologie |
     | **UAV DTM/DSM** | vols 2022 et 2024 | changement net d'élévation (DoD) |
     | **Laser in situ** (1 point) | continu | vérité terrain temporelle, respiration |

2.3. Le choix « produits cloud ASF HyP3 plutôt que traitement SLC local » :
     reproductibilité, coût, cohérence avec l'article de référence — et sa
     conséquence (pas de contrôle sur le déroulement SNAPHU).

## 3. Méthodes I — la chaîne SBAS/ISBAS dense (la voie standard)

3.1. Réseau small-baseline : 91 dates → 354 paires ≤48 j ; contrôle de
     connexité ; statistiques de cohérence par paire.

3.2. Pré-traitement contextuel (l'apport « zones humides » du pipeline) :
     - masque eau/inondation fusionné S2 (NDWI/MNDWI/SCL) + S1 RTC
       (rétrodiffusion, double-bounce) ;
     - masque de gel par T2m ERA5 ;
     - **classification comportementale A–E** des pixels (stable / …/ eau
       permanente) — grille de lecture utilisée dans tout le reste ;
     - indice de qualité W_i par paire.

3.3. SBAS via MintPy (config, point de référence, correction de
     déroulement) ; **ISBAS maison** (WLS pixel par pixel tolérant aux
     paires manquantes, pondéré par cohérence) ; fermeture de phase par
     triplets comme diagnostic de déroulement.

## 4. Résultats I — anatomie d'un échec instructif (QR1)

*Le chapitre « négatif » assumé — c'est lui qui différencie l'étude.*

4.1. **Chronologie des 8 blocages MintPy**, chacun classé
     {bug de code / limite de l'outil / limite physique} (tableau déjà
     construit dans le rapport). Points saillants :
     - la correction `bridging` exige le pixel de référence dans la
       composante connexe de **chaque** interférogramme (limite documentée,
       discussion GitHub #819) — impossible ici ;
     - `phase_closure` : « 0 common regions » — ce n'est pas un bug, c'est
       le terrain : **aucune région du site ne reste cohérente sur
       l'ensemble du réseau**.

4.2. **Le résultat chiffré de l'échec** : 182 pixels résolus seulement
     (SBAS), concentrés hors tourbière ; gain de couverture ISBAS
     (Phase 10) réel mais insuffisant pour une carte ; fermeture de phase
     confirmant des erreurs de déroulement massives en été.

4.3. **La physique derrière (le cœur de l'analyse)** :
     - *Bande C (λ=5.55 cm) vs végétation de fen* : diffusion de volume
       dans le couvert herbacé/mousses, décorrélation temporelle en
       quelques jours dès la croissance printanière ; l'onde ne « voit »
       le sol qu'en fin d'hiver / début de printemps.
     - *Rapport signal/bruit intrinsèque des paires courtes* : à 12–48 j,
       le déplacement attendu (0.05–0.2 mm) est 100× sous le bruit de
       phase ; le réseau dense ne moyenne le bruit que si la cohérence
       survit — ici elle ne survit pas l'été.
     - *Surface intermittente* : inondation/gel rendent des sous-ensembles
       de dates inutilisables par pixel → le « réseau connexe unique »
       supposé par SBAS n'existe pas physiquement ; l'ISBAS corrige
       l'algorithme mais pas le manque d'information.
     - *Hypothèses MintPy héritées des cibles cohérentes* : référence
       visible partout, composantes connexes stables, correction tropo en
       aval — toutes inadaptées à une zone humide petite et intermittente.
     - Synthèse : **l'échec n'est pas un défaut d'implémentation mais une
       incompatibilité physique** entre bande C + cycle 12 j + fen végétalisé
       et l'architecture SBAS dense. (Réponse complète à QR1.)

## 5. Méthodes II — la voie qui marche : paires saisonnières-annuelles optimisées (QR2)

5.1. Principe (d'après Ghezelayagh et al. 2024, §2.3.1) : comparer avril à
     avril (phénologie identique, végétation minimale, signal accumulé sur
     1 an ≫ bruit), sans réseau ni inversion globale.

5.2. **Nos quatre améliorations** (tableau comparatif) :
     1. sélection atmosphérique **quantitative** par ERA5 (|Δtcwv| entre
        les deux dates + pénalité précipitation 24 h) au lieu de relevés
        météo manuels ;
     2. **top-k=3 paires par transition** (au lieu d'une seule) → médiane
        d'ensemble pixel à pixel qui rejette une paire corrompue (démontré
        synthétiquement : paire polluée à +40 mm/an, médiane exacte à
        ±0.15 mm/an) ;
     3. **deramp explicite** par paire sur pixels stables (rampe
        orbitale/atmosphérique résiduelle, non traitée dans l'article) ;
     4. **incertitude par pixel** (IQR inter-paires) au lieu d'un RMSE
        global.

5.3. Estimateur de contrôle : ajustement WLS conjoint de toutes les paires
     (réutilise le cœur ISBAS) — deux estimateurs indépendants qui doivent
     converger.

## 6. Méthodes III — validation multi-capteurs (QR3)

6.1. **DoD drone** : co-registration des DTM 2022/2024 sur les surfaces
     stables (classes A/B), retrait de plan résiduel, différence, masquage
     S2 (eau, changement de végétation), agrégation par classe ; budget
     d'erreur (l'erreur relative après co-registration ≪ erreur absolue).

6.2. **Laser in situ** : décomposition tendance / cycle saisonnier
     (« respiration ») ; test de l'hypothèse avril→avril (même phase du
     cycle chaque année ?) ; calibration absolue de l'InSAR au pixel
     correspondant.

6.3. Schéma de fusion pyramidal : laser (mm, 1 point, continu) → DoD (cm,
     toute l'AOI, 2 époques) → InSAR (relatif, régional) ; critères de
     concordance.

## 7. Résultats II — le déplacement vertical de Rzecin

7.1. Paires retenues (tableau façon Table 1 de l'article : dates, Δt,
     |Δtcwv|, score) ; cartes par paire brut/deramp.

7.2. Carte finale : médiane d'ensemble ± IQR ; statistiques par classe
     comportementale (boxplots, pendant de leur Fig. 6) ; concordance
     ensemble vs WLS.

7.3. DoD drone 2022→2024 : carte, stats par classe, incertitude.

7.4. Laser : tendance nette, amplitude de respiration, comparaison au
     pixel InSAR ; R²/RMSE de la validation croisée (l'équivalent de leur
     Fig. 3/4 — notre R² à nous).

7.5. Tableau de synthèse trois capteurs par zone : « ça donne ça » —
     le résultat central de la thèse. Interprétation attendue : taux
     faible/nul sur le fen quasi naturel (cohérent avec le bas du gradient
     hydrologique de Biebrza), respiration saisonnière dominante sur la
     tendance.

## 8. Discussion

8.1. **Rzecin dans le gradient régional** : nos valeurs vs Biebrza
     (−14.4 mm/an moyen, fens = classe la moins subsidente) ; un fen
     protégé et hydraté ne s'affaisse pas → l'indicateur InSAR détecte
     aussi la *stabilité*, ce qui compte pour le suivi de restauration.

8.2. **Guide méthodologique transférable** (la contribution générique) :
     arbre de décision « quelle stratégie InSAR pour quelle tourbière » :
     - grande + drainée → seasonal-annual simple (article original) ;
     - petite + naturelle → notre variante (pool ERA5 + top-k + deramp +
       validation locale) ; SBAS dense : seulement si cohérence hivernale
       continue et signal > bruit par paire ;
     - le bruit SBAS/MintPy sur tourbière n'est pas gaussien-réductible :
       il est structurel (décorrélation saisonnière + déroulement) — le
       densifier ne le moyenne pas.

8.3. **Perspective bande L** (réponse à QR4, et justification honnête de
     nos limites) :
     - la bande L (λ≈24 cm) pénètre le couvert herbacé et maintient la
       cohérence sur des mois — la solution *physique* au problème
       documenté au §4.3 ;
     - **NISAR** (NASA/ISRO, lancé 2025, L+S bande, libre) et **ROSE-L**
       (ESA, ~2028) fourniront exactement cela ; ALOS-2 existe mais
       archive clairsemée sur l'Europe et accès restreint ;
     - **mais aucune archive L n'existe pour 2022–2024 sur notre site** :
       notre étude documente donc la meilleure stratégie possible avec les
       données réellement disponibles pour la période, et fournit la
       *baseline* (site caractérisé, classes, validation terrain) sur
       laquelle brancher NISAR dès que ses séries couvriront quelques
       années.

8.4. Limites : 3 ans (vs 7 à Biebrza) ; ~90 ha (moyennage spatial limité) ;
     un seul point laser ; DoD sensible à la végétation ; pas de
     piézomètres distribués.

## 9. Conclusions

- QR1 : non — et la cause est physique (bande C / fen), pas logicielle ;
  8 blocages classés, reproductibles, utiles à la communauté.
- QR2 : oui — pipeline B, améliorations mesurables vs l'état de l'art.
- QR3 : validation trois capteurs, R² local, concordance par zone.
- QR4 : arbre de décision + perspective NISAR/ROSE-L.

---

## Annexes prévues

A. Configuration MintPy complète + journal des 8 blocages (reproductibilité).
B. Détail du score ERA5 et sensibilité au choix de la fenêtre ±15 j.
C. Budget d'erreur DoD (co-registration, végétation).
D. Code : dépôt GitHub `SAr-adapted` (pipeline complet, tests synthétiques).

## Liste des figures clés

1. Site + AOI + classes comportementales (carte).
2. Réseau SBAS 354 paires + cohérence moyenne par paire (chute estivale).
3. « Anatomie de l'échec » : composantes connexes par saison + fermeture de
   phase + les 182 pixels SBAS.
4. Sélection pipeline B : tcwv aux dates candidates + paires retenues.
5. Cartes par paire brut vs deramp (le sauvetage du signal).
6. Carte finale médiane ± IQR + boxplots par classe.
7. DoD drone 2022→2024.
8. Série laser : respiration + tendance + point InSAR superposé.
9. Synthèse trois capteurs (le graphe de validation, façon Fig. 3–4 de
   Ghezelayagh).
10. Arbre de décision méthodologique (la figure « à citer »).

## Cibles de publication (par ordre de pertinence)

1. *Ecological Indicators* — même revue que l'article de référence ;
   continuité directe (« companion study » sur le bas du gradient).
2. *Remote Sensing of Environment* — si l'angle méthodologique (échec
   SBAS documenté + guide) est mis en avant.
3. *Remote Sensing* (MDPI) / *International Journal of Applied Earth
   Observation and Geoinformation* — alternatives plus rapides.

## Correspondance chapitres de thèse ↔ phases du dépôt

| Chapitre | Phases/notebooks |
|---|---|
| 2 (données) | phase01 |
| 3 (méthodes SBAS) | phase02–07 |
| 4 (échec instructif) | phase08–10 + rapport |
| 5 (pipeline B) | phase01b–04b (+15) |
| 6–7 (validation) | phase05b–07b (à coder : drone + laser + fusion) |
| 8 (discussion) | phase13 recyclée (hydro/ERA5) |
