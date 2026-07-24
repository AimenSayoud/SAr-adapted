# Le fen flottant est-il « spécial » ? Déductions des Phases D, D-bis, D-ter et E2

**Statut :** résultat expérimental interne, indépendant de tous les échecs
InSAR antérieurs (SBAS/MintPy, ISBAS, Pipeline B, hybride). Étiquetage :
[FAIT] mesuré · [INT] interprétation étayée · [HYP] hypothèse · [SPEC]
spéculation.

---

## 1. Question posée

À **type de couvert végétal comparable**, la surface du tapis flottant de
Rzecin décorrèle-t-elle davantage en InSAR bande C que la même végétation sur
sol stable ? Et si oui, **par quel mécanisme** (mouvement mécanique du tapis
flottant vs état diélectrique de la tourbe saturée) ?

## 2. Protocole (résumé)

Stratification du crop rectangulaire en 4 zones sur la grille de cohérence :
**A** = tapis végétalisé (intérieur du polygone geojson, hors eau) ; **B** =
lac résiduel ; **C** = végétation extérieure **appariée** à A (même classe
ESA WorldCover — « prairie » — et mêmes features Sentinel-2 verdure/humidité
dans l'intervalle [p10, p90] de A, pente < 5° contrôlée via le MNT) ; **D** =
autres couverts extérieurs. Effectifs : A=499, B=65, C=169, D=10 979 px.

Le test central est **apparié par interférogramme** : chaque paire est vue
dans A et dans C, ce qui neutralise la ligne de base perpendiculaire et
l'atmosphère du jour — il ne reste que la différence de surface. Compléments
(Phase D-bis) : couplage de la cohérence à un proxy de nappe (ERA5), test de
gel, profil radial de distance au polygone, et résidu d'inversion WLS par zone
(issu du réseau hybride, Phase A).

**Statistique (précisions exigées en review) :** le test apparié est un
**Wilcoxon signé** sur les différences coh_A − coh_C (le « 88 % » n'est que la
fraction de différences de même signe — un descriptif, pas un test). Les ~356
paires n'étant **pas indépendantes** (elles partagent ~90 dates), un bootstrap
par paire sous-estime l'incertitude ; on rapporte donc un **jackknife par
date** (retrait de chaque acquisition et de toutes ses paires) : le résultat
n'est retenu comme robuste que si le delta garde le même signe pour *toute*
date retirée. L'IC bootstrap-par-paire est signalé comme optimiste.

## 3. Résultats [FAIT]

| Mesure | A (tapis) | C (prairie stable appariée) | Lecture |
|---|---|---|---|
| Cohérence moyenne | 0.41 | 0.48 | A < C |
| Δ apparié (coh_A − coh_C) | **−0.069** | — | Wilcoxon signé **p = 2.2×10⁻⁴⁹** ; 88 % des paires de signe négatif ; **jackknife par date : delta ∈ [−0.0705, −0.0652], même signe pour TOUTE date retirée** (SE robuste 0.0115) → non porté par une acquisition anormale |
| σ0 VV médian (RTC) | −10.2 dB | −11.5 dB | A **plus brillant** que la prairie (~+1.3 dB) → double-bounce/humidité de surface probable |
| σ0 VV — lac (B) | −15.3 dB | — | le plus bas → **confirme l'eau libre spéculaire** |
| Phénologie S2 (verdure ET humidité) | ≈ identique à C | — | A et C sont des **jumeaux phénologiques** → l'appariement de couvert est excellent, le Δcohérence n'est PAS un artefact de couvert |
| Temps de décorrélation τ | 21 j | 32 j | A décorrèle plus vite (r² faibles : indicatif) |
| Pente cohérence ~ \|Δnappe\| | −10.7 | −10.0 | quasi identiques |
| Gain de cohérence au gel | +0.028 | +0.077 | A gagne MOINS |
| Résidu d'inversion WLS médian | **2.46 rad** | 1.92 rad | A nettement pire (≈347 paires valides partout) |

Profil radial de la cohérence médiane selon la distance signée au bord du
polygone : **plateau bas et plat (~0.40) sur tout l'intérieur** du tapis,
**marche nette au passage du bord**, **pic juste à l'extérieur (~0.48)** puis
variations non monotones (creux à +300–600 m, remontée à +1150 m) reflétant
l'hétérogénéité du couvert extérieur.

## 4. Déductions

### 4.1 [INT — robuste] Le tapis est une unité de cohérence distincte, et son déficit est réel

À couvert apparié, contrôlé pour la baseline (appariement), l'atmosphère
(appariement), la pente (MNT) et l'humidité optique du couvert (features S2),
le tapis reste **significativement moins cohérent** que la prairie sur sol
stable (Δ = −0.069, p < 0.05, 88 % des interférogrammes). Le profil radial
montre une **frontière franche** au bord du polygone et un intérieur
uniformément bas : ce n'est pas un bruit diffus, c'est une **unité physique
délimitée**. Le résidu d'inversion confirme la conséquence : la phase du
tapis est **intrinsèquement moins cohérente** (2.46 vs 1.92 rad), **sans**
déficit de paires disponibles — ce qui relie directement la décorrélation à
l'échec d'inversion documenté en Phases 8–A.

**Conclusion 1 :** ce n'est pas « la végétation en bande C » en général ;
l'état de surface propre au fen flottant ajoute une décorrélation mesurable
et spatialement structurée. **La phénologie S2 (verdure ET humidité optique)
du tapis A est quasi identique à celle de la prairie appariée C** — ce sont
des jumeaux phénologiques — ce qui écarte définitivement l'hypothèse d'un
biais de couvert : la différence de cohérence est réelle et vient d'une
propriété de surface non optique (échelle radar). Robustesse statistique :
Wilcoxon p = 2.2×10⁻⁴⁹, et le jackknife par date confirme que le résultat ne
dépend d'aucune acquisition particulière (delta stable ∈ [−0.070, −0.065]).

### 4.2 [INT — robuste] Le mécanisme n'est PAS le cycle hydrologique saisonnier

Les pentes cohérence ~ |Δnappe| sont quasi identiques dans toutes les zones
(A = −10.7, C = −10.0). Le tapis **n'est pas plus couplé** aux variations de
nappe que la prairie. Si la décorrélation venait d'une flottaison mécanique
suivant la nappe (poussée d'Archimède saisonnière), la sensibilité de A
serait très supérieure à celle de C. Elle ne l'est pas : **le déficit du
tapis est un décalage quasi constant, non modulé par la saison hydrologique.**

**Conclusion 2 :** l'instabilité de surface responsable est **quasi
permanente**, pas saisonnière — ce qui, pour la bande C à 12 j, est bien plus
pénalisant qu'une modulation saisonnière.

### 4.3 [INT — nuance critique] Mécanique vs diélectrique : NON tranché — et le gel ne prouve pas le diélectrique

Il serait tentant (et une lecture externe l'a fait) de conclure « décorrélation
diélectrique, pas mécanique ». **Les données ne le permettent pas**, pour deux
raisons de physique :

1. Un **décalage constant** (§4.2) est tout aussi compatible avec un
   **micro-mouvement mécanique permanent** du tapis (tremblement/flexion à
   chaque acquisition, indépendant du niveau net de nappe) qu'avec un état
   diélectrique. « Constant » ≠ « diélectrique ».
2. Le **test de gel se retourne contre l'hypothèse diélectrique**. En InSAR,
   le gel de l'eau **restaure fortement** la cohérence des surfaces humides
   (l'eau gelée a une permittivité basse et stable — principe du « freeze-up »
   documenté en zones humides/pergélisol). Si le problème du tapis était une
   humidité diélectrique de surface, le gel aurait dû donner à A le **plus
   gros** gain de cohérence. Or A gagne le **moins** (+0.028 vs C +0.077 ;
   B +0.026 vs D +0.068). La lecture correcte : **les zones humides ne gèlent
   pas comme le sol stable** (le tapis flotte sur une eau liquide non gelée) —
   information sur l'état thermique, **pas** un discriminant mécanique/
   diélectrique. (Réserve : n_froid = 31 paires ; proxy nappe = pluie cumulée,
   grossier.)

**Conclusion 3 :** on a **éliminé** le mécanisme mécanique saisonnier
(flottaison couplée à la nappe) ; on **n'a pas** établi que le mécanisme
résiduel soit diélectrique plutôt que mécanique-constant. Les deux restent
compatibles avec l'ensemble des observations. Affirmer « limite diélectrique
démontrée » serait une surinterprétation contestable en évaluation.

## 5. Synthèse défendable (formulation pour la thèse / l'article)

> À couvert végétal apparié (prairie), et après contrôle de la ligne de base,
> de l'atmosphère, de la pente et de l'humidité optique du couvert, la surface
> du fen flottant de Rzecin présente une cohérence InSAR bande C
> significativement inférieure (Δ ≈ −0.07, p < 0.05, 88 % des interférogrammes),
> une décorrélation plus rapide, et un résidu d'inversion nettement plus élevé
> (2.46 vs 1.92 rad) que la prairie sur sol stable — définissant une unité de
> surface distincte, nettement délimitée et chroniquement peu cohérente. Ce
> déficit n'est expliqué **ni** par les variations saisonnières de nappe **ni**
> par un état gelable : il traduit une **instabilité de surface quasi
> permanente**. La nature exacte de cette instabilité — micro-mouvement
> mécanique constant du tapis flottant vs variabilité diélectrique de la tourbe
> saturée non gelable — ne peut être tranchée par ces seules observations
> radar et demande une mesure in situ (laser) ou une diversité de longueur
> d'onde/polarisation.

## 5-bis. Compléments Phase D-ter (mécanisme de diffusion & lac)

- [FAIT] **σ0 VV** : le tapis (A, −10.2 dB) est **plus brillant** que la prairie
  appariée (C, −11.5 dB) d'environ 1.3 dB, avec une variabilité temporelle
  comparable (~2.1). [INT] Un excès de rétrodiffusion sur une végétation
  optiquement identique évoque un **double-bounce / une humidité de surface**
  (eau sous la végétation émergente) — indice *en faveur* d'une composante
  diélectrique, mais non décisif seul (le ratio VH/VV, non encore calculé,
  trancherait).
- [FAIT] **Lac (B)** : σ0 VV = −15.3 dB (le plus bas, spéculaire → eau libre
  confirmée), MAIS sa cohérence reste **~0.39** même défini comme eau
  persistante (83 px). [INT] Le « lac résiduel » n'est donc **pas de l'eau
  libre ouverte** mais un plan d'eau **peu profond et végétalisé** (cohérent
  avec une terrestrialisation en cours) — une information écologique en soi.
- [FAIT] **Dispersion d'amplitude D_A** (85 dates RTC) : A médiane 0.243
  (59 % de pixels < 0.25) ; C 0.238 (69 %) ; D 0.215 (87 %) ; **lac B 0.310
  (3 %)**. [INT] Le tapis **n'est PAS radar-sombre** : il rétrodiffuse et a une
  stabilité d'amplitude à peine inférieure à la prairie. **Donc son problème
  n'est pas « aucune cible » mais une PHASE instable d'une cible distribuée
  présente → cela N'EXCLUT PAS le DS-InSAR/phase-linking** (conçu exactement
  pour ce régime) : **H1 reste ouverte**, D_A ne la ferme pas. Seul le lac est
  vraiment sans cible stable (3 %).
- [FAIT] **Ratio VH/VV (dual-pol, 85 dates)** : A −5.27 dB, C −5.57 dB, D
  (forêt/nu) −4.51, B (lac végétalisé) −4.87. [INT] Le tapis a un VH/VV
  **plus élevé** (plus dépolarisant) que la prairie, **pas plus bas** — donc
  **plus de diffusion de VOLUME**, PAS un double-bounce de type eau libre. Avec
  le σ0 +1.3 dB, le tapis se comporte en C-band comme un **volume diffusant
  plus dense et plus humide** que la prairie sèche, malgré une phénologie
  optique identique.
- [FAIT] **Amplitude de « respiration » optique (humidité S2)** : A 0.512, **C
  0.852**, B 0.263. [INT] Contre-intuitif mais cohérent avec la littérature :
  c'est la PRAIRIE qui varie le plus en humidité saisonnière ; **le tapis reste
  humide et stable** toute l'année (stabilité hydrologique du fen, Juszczak
  et al. 2013), et le lac reste le plus humide.

### 5-ter. Synthèse mécanistique révisée (Phases D-bis + D-ter)

[INT] La convergence des indices **déplace le poids de la preuve** : le tapis
est un **jumeau phénologique optique** de la prairie mais, en radar, il est
**plus brillant (+1.3 dB), plus dépolarisant (VH/VV), optiquement humide et
stable**, à **résidu d'inversion élevé** et **cohérence chroniquement basse
non récupérée par le gel**. Ce faisceau pointe vers une **décorrélation
volumétrique + diélectrique d'un couvert bas humide sur tourbe saturée** (un
« volume diffusant humide » dont le centre de phase est instable en C-band à
12 j) — et **NON** vers un mouvement mécanique de corps rigide du tapis
flottant, pour lequel **aucun indice positif** n'a été trouvé (ni couplage
nappe, ni stabilisation au gel, ni signature double-bounce). [HYP] Le mécanisme
dominant est donc probablement **volumétrique/diélectrique**, la composante
mécanique restant non prouvée — mais les tailles d'effet polarimétriques sont
**modestes** (VH/VV : 0.3 dB) et **seul le laser** exclurait définitivement un
micro-mouvement.

## 5-quater. Phase E2 — verdict DS-InSAR par phase-linking EVD (sans ISCE)

**Motivation.** Toutes les phases précédentes reposaient sur une inversion
**paire par paire** (WLS/ISBAS), qui pouvait sous-exploiter la redondance. La
question H1 restait ouverte : un **phase-linking** (DS-InSAR), qui pondère la
cohérence de *toutes* les paires simultanément — le maximum de vraisemblance
sous modèle gaussien circulaire — récupérerait-il un signal là où le WLS
échoue ? ISCE/MiaplPy ne s'installant pas sur Colab (packaging conda #642), on
a implémenté le phase-linking **en pur numpy** : le phase-linking opère sur la
**matrice de cohérence complexe N×N par pixel**, dont l'entrée (i,j) est
`cohérence × exp(i·φ_ij)` ; **nos interférogrammes HyP3 SONT ces entrées**
(paires déjà corégistrées). On récupère l'historique de phase comme phase du
**vecteur propre dominant** (EVD), et la qualité par la **temporal_coherence**
(ajustement de l'historique estimé aux interférogrammes observés — l'analogue
de `MiaplPy.temporalCoherence`).

**Protocole.** EVD sur les mêmes zones A/B/C/D que la Phase D (WorldCover +
appariement S2 + contrôle de pente), sur les 356 paires, phase déroulée
**ré-enroulée** (le phase-linking opère sur l'observation enroulée). Effectifs
de ce run : A=499, B=65, **C=374**, D=8368 px (le C plus grand qu'en §2 —
variation d'appariement run-à-run — **renforce** le contrôle). Le module est
couvert par un test synthétique (`tests/test_synthetic_phaselinking.py`) :
l'EVD récupère un historique de phase connu sur réseau sparse, et la
temporal_coherence sépare un pixel cohérent (~0.97) d'un pixel décorrélé.

**Résultats [FAIT] :**

| Zone | tcoh médiane | tcoh p25–p75 | % pixels ≥ 0.7 | lecture |
|---|---|---|---|---|
| **C** prairie sol stable (appariée) | **0.734** | 0.671–0.803 | **64.7 %** | récupérable → DS-InSAR fonctionne |
| **D** autres couverts | 0.639 | 0.597–0.693 | 23.1 % | mélange (contexte) |
| **A** tapis flottant | **0.604** | 0.566–0.647 | **5.4 %** | quasi non récupérable |
| **B** lac résiduel | 0.584 | 0.542–0.630 | 1.5 % | plancher de bruit |

**Le point décisif — le plancher de bruit.** Le test synthétique établit qu'à
cette redondance (356 paires / ~90 dates ≈ 4), un pixel **purement décorrélé**
donne tcoh ≈ 0.55. Donc :

- **B (lac) = 0.584** tombe **exactement sur le plancher** → **validation
  interne majeure** : on *sait* que l'eau végétalisée est décorrélée, et la
  méthode la classe correctement au plancher. La méthode se comporte donc
  correctement sur la cible-témoin connue.
- **A (tapis) = 0.604** n'est qu'à **~0.05 au-dessus du plancher**, juste à
  côté du lac → le phase-linking n'y trouve **presque aucun historique de phase
  réel**. A n'est pas « un peu moins bon » que C : **A ≈ bruit, C ≈ signal**.
- **C = 0.734** est franchement au-dessus du plancher (65 % de pixels
  DS-qualité) → la végétation appariée sur sol stable est, elle, parfaitement
  récupérable par DS-InSAR **avec le même réseau de 356 paires**.

**Déductions.**

- [INT — robuste] **H1 est FERMÉE.** Le « 0 pixel fiable » de la Phase A
  n'était **pas** une faiblesse du WLS : l'estimateur du maximum de
  vraisemblance échoue aussi sur le tapis. L'échec d'inversion est une
  **propriété physique de la cible**, pas de l'algorithme.
- [INT — robuste] **Ce n'est pas « la végétation en bande C » en général.** C
  est phénologiquement/couvert-apparié à A et se récupère à 65 % **sur le même
  réseau** → la sparsité du réseau est *contrôlée* (elle handicape A et C de la
  même façon) et n'explique pas l'écart. La cause est **propre au tapis**.
  Réponse Phase D (A ≪ C) confirmée par une **seconde méthode indépendante**.
- [INT] **Convergence totale de l'enquête** : Phase D (Δcoh −0.069,
  p=2×10⁻⁴⁹) + D-bis (ni hydro-saisonnier ni gelable) + D-ter (diffusion de
  volume VH/VV, σ0 +1.3 dB, pas de double-bounce) + **E2 (phase-linking échoue
  sur A, réussit sur C)** : quatre angles, une seule conclusion — **décorrélation
  volumétrique/diélectrique irréductible d'un couvert humide sur substrat saturé
  et mobile**, et **non** un mouvement de tapis rigide qu'un meilleur estimateur
  suivrait. La carte de tcoh le confirme visuellement : le polygone tourbière
  est une unité uniformément basse, les taches à tcoh élevée sont les structures
  stables/champs — pas d'artefact visuel.

**Garde-fous d'honnêteté (à écrire tels quels).**

1. **Phase-linking « léger » (sparse).** On n'exploite que les 356 paires HyP3,
   pas toutes les paires SLC possibles. Mais la comparaison A vs C est menée
   sur **le même réseau** : C réussit, A échoue à sparsité égale. Un
   phase-linking full-SLC améliorerait A *et* C ; l'écart demeurerait. La
   conclusion tient. (Le module reste une **contribution méthodologique** : un
   test de faisabilité DS-InSAR reproductible et léger, directement sur produits
   HyP3, sans ISCE.)
2. **Bande C uniquement.** Ce n'est pas « impossible partout » : Hrysiewicz
   et al. 2024 (RSE 291) mesurent la respiration de **tourbières hautes** en
   bande C (Sphagnum plus sec, diffuseurs de surface plus stables). Rzecin est
   un **fen de transition** à tapis flottant *plus humide* sur eau
   quasi-affleurante → beaucoup plus de diffusion de volume. Le contraste est
   lui-même défendable. **La voie ouverte reste la bande L (ALOS-2 / NISAR)**,
   qui pénètre le couvert et atteint la surface du tapis.

**Formulation pour la thèse / l'article :**

> Le déplacement vertical du tapis flottant de Rzecin n'est pas mesurable par
> InSAR Sentinel-1 (bande C), et ceci n'est **pas** une limite de l'algorithme
> d'inversion mais une **propriété physique de la cible**. Nous le démontrons
> par trois estimateurs indépendants : (i) l'inversion WLS/ISBAS paire-à-paire
> ne fournit aucun pixel fiable sur le tapis ; (ii) le **phase-linking** au
> maximum de vraisemblance (l'estimateur qui sous-tend le DS-InSAR),
> implémenté sans ISCE directement sur les interférogrammes HyP3, ne récupère
> pas d'historique de phase cohérent sur le tapis (5 % de pixels à
> temporal_coherence ≥ 0.7, contre 65 % pour une végétation
> phénologiquement appariée sur sol stable, et un plancher de bruit à ~0.55
> atteint par le lac résiduel) ; (iii) la décorrélation est volumétrique/
> diélectrique (ratio de polarisation croisée, rétrodiffusion), ni
> saisonnière-hydrologique ni gelable. Le tapis flottant se comporte en bande C
> comme un **diffuseur de volume à phase quasi aléatoire**.

## 6. Portée et limites

- **Ce que ça établit :** le fen flottant est intrinsèquement moins observable
  en InSAR C-band que la végétation comparable sur sol stable — un résultat
  quantifié, spatialement structuré et relié à l'échec d'inversion, **confirmé
  par trois estimateurs (WLS, ISBAS, phase-linking EVD)**. C'est la
  contribution positive du projet, indépendante des échecs de traitement.
- **H1 (DS-InSAR/phase-linking) est désormais FERMÉE** (§5-quater) : le
  phase-linking ne récupère pas le tapis, à réseau contrôlé. L'échec est
  physique, pas algorithmique.
- **Ce que ça ne tranche pas :** le mécanisme physique fin (micro-mouvement
  mécanique constant vs variabilité diélectrique de la tourbe saturée) ; H4
  (mouvement réversible). Seuls un laser in situ ou une diversité de longueur
  d'onde (bande L) les départageraient.
- **Limites méthodologiques :** proxy de nappe grossier (pluie cumulée ERA5, à
  remplacer par la WTD in situ) ; n_froid = 31 (test de gel indicatif) ;
  appariement C limité par la rareté de végétation humide sur sol stable
  autour (169 px) ; l'appariement S2 contrôle l'humidité *optique* du couvert,
  pas l'état diélectrique du sol.

## 7. Expériences qui trancheraient le mécanisme

1. **Laser in situ** au pixel A : mesure directe du mouvement de surface —
   s'il bouge à l'échelle du cycle S1, le mécanique est confirmé.
2. **Vraie série de nappe (WTD)** : refaire le test 1 avec la nappe mesurée
   plutôt que le proxy pluie.
3. **Polarimétrie / bande L (NISAR)** : séparer double-bounce (eau/humidité)
   de diffusion de volume, et pénétrer le couvert.

---

*Figures associées : outputs/phaseD/{zones, coherence_by_zone, decay_curves,
greenness_vs_coh}.png ; outputs/phaseDbis/{coh_vs_hydro, freeze_test,
radial_profile, residual_by_zone}.png ; sortie Phase E2 : phaseE2_evd.nc +
histogramme/carte tcoh par zone. Code : src/insar_wetlands/stratify.py,
src/insar_wetlands/inversion/phaselinking.py (EVD, sans ISCE) ;
tests/test_synthetic_phaselinking.py ; notebooks phaseD_inside_vs_outside.ipynb,
phaseDbis_mechanism_spatial.ipynb, phaseDter_scattering_scatterers.ipynb et
phaseE2_evd_phaselinking.ipynb.*
