# Revue par comité scientifique — évaluation type RSE / ISPRS J. / Nature Communications

**Mandat.** Évaluer le rapport général du projet Rzecin comme un panel de
reviewers spécialisés (physique radar/InSAR, télédétection des zones humides,
hydrologie des tourbières, géodésie/validation, statistiques, écologie des
tourbières). Objectif : identifier toutes les faiblesses, chercher activement
les erreurs, hypothèses non démontrées, biais et explications alternatives ;
puis proposer un plan précis pour atteindre le niveau publiable.

**Recommandation éditoriale globale (simulée) : *Major revision / reject &
resubmit*.** Le travail est reproductible et honnête, mais dans son état
actuel il ne constitue pas encore une contribution publiable dans une revue
de premier plan : le résultat InSAR est une non-détection non diagnostiquée,
la caractérisation du site est partiellement erronée, les données in situ
disponibles ne sont pas exploitées, et la validation (N=1) est insuffisante.
Le potentiel existe pour un excellent article de fusion multi-capteurs
(cadre A ci-dessous), sous conditions.

---

## Reviewer 1 — Physique radar / InSAR

**R1.1 [Majeur] Le « résultat » InSAR est une non-détection non diagnostiquée.**
182 pixels résolus hors zone d'intérêt et −4.5 ± 18.9 mm/an ne sont pas des
mesures mais l'absence de mesure. Aucune courbe de décorrélation cohérence
vs ligne de base temporelle par saison n'est fournie — l'affirmation
« décorrélation en quelques jours » est asservie, non mesurée. RSE exige la
statistique de cohérence quantitative du site.

**R1.2 [Majeur] Réseau interférométrique sous-optimal et non comparé.**
Seules des paires ≤ 48 j ont été utilisées. Les succès publiés (Hrysiewicz
et al. 2024) combinent lignes de base LONGUES et courtes. L'échec ne peut
être attribué à la physique tant que le design de réseau des études qui
réussissent n'a pas été répliqué. C'est le test manquant central.

**R1.3 [Majeur] Point de référence non validé sur un site flottant.**
Un fen flottant n'offre AUCUN point stable en son sein ; la référence a été
prise hors AOI mais sa stabilité n'a jamais été vérifiée (ni GNSS, ni
nivellement, ni laser). Toute mesure « relative à un pixel stable » est
suspecte tant que ce pixel n'est pas prouvé stable.

**R1.4 [Majeur] Correction atmosphérique insuffisante.**
Les paires annuelles montrent des rampes jusqu'à 179 mm/an. Un deramp
planaire ne capture pas l'atmosphère turbulente. GACOS ou une correction
troposphérique ERA5 (pyaps) doit être appliquée et le résidu quantifié avant
toute interprétation.

**R1.5 [Majeur] Erreurs de déroulement non contrôlées.**
Le changement de signe entre paires de 2 ans prouve des sauts de déroulement.
La fermeture de phase sur les triplets annuels n'est pas rapportée. Sans
diagnostic de qualité de déroulement, les taux par paire ne sont pas fiables.

**R1.6 [Mineur] Erreur topographique résiduelle.**
Le terme phase topographique (erreur du DEM Copernicus × ligne de base
perpendiculaire) n'est pas discuté ; pertinent pour les paires annuelles à
plus grande ligne de base.

**R1.7 [Contrainte] Polarisation VV.**
VV est peu sensible au double-bounce ; HH serait préférable en régime
inondé mais est indisponible sur cette trace S1 (VV/VH). À énoncer comme
limitation instrumentale motivant la bande L.

## Reviewer 2 — Télédétection des zones humides

**R2.1 [Majeur — erreur factuelle] Mauvaise caractérisation hydrologique du site.**
Le rapport affirme une « nappe au-dessus de la végétation » à l'échelle du
site. La littérature spécifique à Rzecin indique une nappe typiquement à
**0–30 cm SOUS la surface** durant la saison de croissance, sur un tapis de
Sphagnum flottant, avec une **stabilité hydrologique exceptionnelle**
(Juszczak et al. 2013 ; « Exceptional hydrological stability of a
Sphagnum-dominated peatland », Quaternary/Holocene). L'eau libre concerne le
lac résiduel et les dépressions, pas l'ensemble du fen. Conséquence lourde :
un tapis de Sphagnum émergent ressemble DAVANTAGE aux bogs où l'InSAR C-band
réussit qu'à un marais inondé — ce qui affaiblit l'argument « régime inondé =
inobservable » et renforce l'hypothèse « échec = méthode ».

**R2.2 [Majeur] Classes comportementales non validées et circulaires.**
Les classes A–E sont dérivées de la cohérence/eau puis utilisées pour
interpréter l'InSAR — raisonnement circulaire. Elles doivent être validées
contre une carte de végétation/microforme indépendante (hummock/hollow,
Sphagnum vs cariçaie vs eau libre).

**R2.3 [Modéré] Complétude du masque d'eau optique non évaluée.**
Sur un site très humide et souvent nuageux, la couverture S2 exploitable et
les lacunes temporelles ne sont pas quantifiées.

## Reviewer 3 — Hydrologie des tourbières

**R3.1 [Majeur] Données de nappe in situ disponibles et NON utilisées.**
Rzecin est un site d'eddy-covariance instrumenté (Juszczak, Chojnicki,
Lamentowicz) avec suivi de nappe (WTD) et expériences de manipulation du
niveau d'eau. Le rapport propose ERA5 comme proxy hydrologique alors que des
séries WTD in situ existent presque certainement. Ne pas utiliser la donnée
de terrain disponible est rédhibitoire pour ce volet.

**R3.2 [Majeur] Contexte climatique absent.**
2022 fut une année de sécheresse en Europe centrale. Un « taux net ≈ 0 » sur
2022–2024 doit être lu au regard de la trajectoire WTD et des anomalies
climatiques : la stabilité hydrologique documentée du site (résistance aux
phases sèches) pourrait EXPLIQUER l'absence de subsidence nette — hypothèse
positive et testable, non un simple « rien détecté ».

**R3.3 [Modéré] Aucun bilan de masse tourbe.**
Densité apparente, profondeur de tourbe, partition oxydation/retrait :
absents. Si l'ambition est le carbone, le lien déplacement → perte de masse
n'est pas établi.

## Reviewer 4 — Géodésie / validation

**R4.1 [Majeur] Validation N=1.**
Un seul point laser ne valide qu'une série temporelle en un lieu. Aucune
statistique spatiale (R²/RMSE) n'est possible, contrairement à Ghezelayagh
(11 sites, 337 relevés). Pour un produit spatial, c'est insuffisant en l'état.

**R4.2 [Majeur] Budget d'erreur UAV non chiffré et LiDAR vs photogrammétrie non tranché.**
Un DSM photogrammétrique voit la canopée, pas le sol, et a une précision
verticale de ~5–15 cm — supérieure au signal (~cm). Faut-il un LiDAR UAV
(pénètre) ou du SfM (canopée) ? Aucune stratégie GCP/points de contrôle,
aucun RMSE de co-registration, aucun masque de changement de hauteur de
végétation n'est fourni.

**R4.3 [Majeur] Absence d'ancrage géodésique absolu.**
Ni RTK-GNSS ni nivellement de précision pour lier les mesures à un référentiel
et valider la stabilité de la référence InSAR (cf. R1.3).

## Reviewer 5 — Statistiques

**R5.1 [Majeur] Confusion dispersion / incertitude.**
« IQR ±18.9 mm/an » est présenté comme une incertitude. L'IQR est la
dispersion inter-paires, PAS l'erreur-type de la médiane (~IQR/1.35/√n_eff).
La conclusion de non-détection tient, mais la quantité rapportée est mal
définie.

**R5.2 [Majeur] Degrés de liberté effectifs surestimés.**
Les 9 paires partagent des dates → fortement corrélées ; n_eff ≪ 9. Aucune
propagation de cette corrélation.

**R5.3 [Majeur] Absence de test d'hypothèse formel.**
Aucun intervalle de confiance ni test « −4.5 diffère-t-il de 0 ? de −14 ? ».
Le rapport conclut qualitativement.

**R5.4 [Modéré] Deramp = risque de sur-ajustement / fuite de signal.**
Ajuster 3 paramètres par paire puis les appliquer à l'AOI peut absorber du
vrai signal ; aucune validation croisée du modèle de rampe.

**R5.5 [Modéré] Biais de sélection top-k.**
Sélectionner les paires sur un score atmosphérique ERA5 (un proxy) peut
biaiser les estimations de déplacement conservées ; non examiné.

## Reviewer 6 — Écologie des tourbières

**R6.1 [Majeur] Interprétation écologique indigente.**
Les classes sont radar, pas écologiques. Aucun lien aux communautés
végétales, aux microformes, ni aux flux de GES pourtant mesurés sur site
(eddy-covariance). La signification écologique du mouvement de surface d'un
poor fen à Sphagnum n'est pas discutée.

**R6.2 [Modéré] Statut de conservation/gestion non décrit.**

## Faiblesses transversales / conceptuelles

- **C1 [Majeur] Confusion de portée.** Le rapport est à la fois un papier de
  méthode InSAR, une étude de site et un cadre de fusion. Une revue de haut
  niveau exige UNE contribution focalisée.
- **C2 [Majeur] La « narration d'échec » n'est pas encore une contribution.**
  Documenter un échec InSAR de plus sur un site humide n'est pas nouveau tant
  que la cause (site vs méthode) n'est pas isolée expérimentalement.
- **C3 [Positif] Reproductibilité forte** (code, tests synthétiques, git) —
  un atout à mettre en avant (dépôt public, données ouvertes).
- **C4 [Majeur] Littérature spécifique au site absente** (Juszczak,
  Lamentowicz, Chojnicki, Milecka, groupe PEATBOG/WETMAN) — signale une
  connaissance incomplète de son propre terrain.

---

## Plan de transformation en étude publiable de très haut niveau

Deux cadres viables. Le cadre A est le plus solide et faisable ; le cadre B
est conditionnel aux expériences.

### Cadre A (recommandé) — « Water-table-corrected surface motion of a floating poor fen by laser–UAV–Sentinel fusion »
Contribution : un cadre transférable séparant la respiration réversible
(pilotée par la nappe) du changement net, sur une classe de tourbières
(fens minerotrophes flottants) sous-étudiée par rapport aux bogs.
Cible : *Remote Sensing of Environment*, *ISPRS J.*, ou *JGR-Biogeosciences*.

### Cadre B (conditionnel) — « Observability of C-band InSAR across peatland hydrological regimes »
Contribution : un cadre de décision « quelle stratégie InSAR pour quelle
hydromorphologie de tourbière », étayé par l'expérience site-vs-méthode.
Cible : *ISPRS J.* ou *IEEE TGRS*. Publiable seulement si les expériences
InSAR sont réellement conduites (et même un résultat négatif rigoureux
devient alors une contribution).

### Actions requises (checklist priorisée)

**Bloc 1 — Corriger et documenter le site (indispensable, rapide)**
1. Réécrire la caractérisation hydrologique avec la littérature Rzecin
   (nappe 0–30 cm sous la surface, tapis Sphagnum flottant, stabilité
   hydrologique — Juszczak et al. 2013 ; « Exceptional hydrological
   stability… »). Corriger l'affirmation « eau au-dessus de la végétation ».
2. Intégrer une carte de végétation/microforme et le statut de gestion.

**Bloc 2 — Exploiter les données in situ existantes (indispensable)**
3. Récupérer et utiliser les séries WTD in situ 2022–2024 (réseau
   eddy-covariance) ; établir la fonction de transfert élévation = f(WTD).
4. Ajouter le contexte climatique (sécheresse 2022) et relier au WTD.
5. Relier le mouvement net aux flux de GES (eddy-covariance) — le lien
   carbone qui hausse l'impact.

**Bloc 3 — Le test InSAR décisif (tranche site vs méthode ; cœur du cadre B)**
6. Répliquer le réseau long+court baseline (façon Hrysiewicz), hors MintPy si
   besoin (SNAPHU/ISCE, DS-InSAR/PSI), avec correction atmosphérique
   GACOS/ERA5.
7. Produire les courbes de décorrélation cohérence vs baseline par saison et
   les cartes de cohérence ; quantifier la fraction d'AOI jamais cohérente.
8. Test décisif : corréler la série InSAR au pixel du laser (récupère-t-on la
   respiration, comme les caméras de Hrysiewicz ?). Rapporter r et RMSE.
9. Diagnostic de déroulement : fermeture de phase sur les triplets annuels.

**Bloc 4 — Géodésie et validation (indispensable pour un produit spatial)**
10. Clarifier LiDAR vs SfM pour l'UAV ; budget d'erreur vertical complet
    (GCP, points de contrôle, RMSE de co-registration) ; masque de changement
    de végétation (NDVI multitemporel) ; DoD avec incertitude par pixel ;
    dé-aliasing par le WTD aux dates de vol.
11. Ancrage géodésique : RTK-GNSS / nivellement d'au moins quelques repères ;
    valider la stabilité de la référence InSAR.
12. Élargir la validation au-delà de N=1 : utiliser le DoD UAV comme
    validation SPATIALE de l'InSAR là où cohérent ; ajouter tout nivellement
    disponible.

**Bloc 5 — Rigueur statistique (indispensable)**
13. Remplacer « IQR = incertitude » par l'erreur-type de la médiane avec
    n_eff tenant compte des dates partagées ; fournir IC et tests formels
    (−4.5 vs 0 ; vs −14).
14. Valider par validation croisée le modèle de deramp ; tester la fuite de
    signal ; quantifier le biais de sélection top-k.
15. Propagation d'incertitude de bout en bout (InSAR + UAV + laser + fusion).

**Bloc 6 — Cadrage, comparaison et perspective**
16. Focaliser sur UNE contribution (cadre A recommandé).
17. Tableau comparatif vs études InSAR tourbières publiées (bog vs fen,
    longueur d'onde, réseau, validation, RMSE).
18. Plan d'acquisition NISAR (bande L) : positionner Rzecin comme site
    supersite de cal/val multi-capteurs.
19. Disponibilité des données/code (atout reproductibilité).

### Barre pour *Nature Communications* spécifiquement
Un single-site, même bien fait, y passe rarement. Il faudrait soit
(i) plusieurs sites couvrant le gradient hydromorphologique (bog → fen
flottant → marais) pour établir une loi d'observabilité générale, soit
(ii) une avancée conceptuelle forte (p. ex. séparation quantitative
respiration/subsidence validée et reliée au bilan carbone à l'échelle
régionale). Sinon, viser RSE/ISPRS J., plus adaptés et déjà ambitieux.

---

### Synthèse en une phrase pour l'auteur
Le projet a le squelette d'un très bon article de fusion multi-capteurs, mais
doit : (1) corriger la description du site avec sa propre littérature,
(2) utiliser les données in situ (WTD, GES) existantes, (3) conduire le test
InSAR site-vs-méthode au lieu de conclure sur la physique, (4) chiffrer la
géodésie UAV et l'incertitude, (5) durcir la statistique. Ces cinq blocs
transforment une non-détection non diagnostiquée en une étude de référence.
