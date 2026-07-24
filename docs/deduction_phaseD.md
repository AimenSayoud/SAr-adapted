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

- [INT — robuste] **H1 est fermée pour le réseau HyP3 disponible actuel.** Le
  « 0 pixel fiable » de la Phase A n'était **pas** une faiblesse du WLS :
  l'estimateur du maximum de vraisemblance échoue aussi sur le tapis, sur le
  **même réseau de 356 paires**. L'échec d'inversion est une **propriété
  physique de la cible**, pas de l'algorithme. (Portée exacte : un travail futur
  avec **toutes les SLC** — au-delà des 356 paires HyP3 — pourrait modifier les
  *valeurs absolues* de tcoh ; mais l'**écart relatif A vs C**, mesuré à réseau
  égal, devrait persister — voir garde-fou 1.)
- [INT — robuste] **Ce n'est pas « la végétation en bande C » en général.** C
  est phénologiquement/couvert-apparié à A et se récupère à 65 % **sur le même
  réseau** → la sparsité du réseau est *contrôlée* (elle handicape A et C de la
  même façon) et n'explique pas l'écart. La cause est **propre au tapis**.
  Réponse Phase D (A ≪ C) confirmée par une **seconde méthode indépendante**.
- [INT] **Convergence des quatre angles** : Phase D (Δcoh −0.069, p=2×10⁻⁴⁹) +
  D-bis (ni hydro-saisonnier ni gelable) + D-ter (diffusion de volume VH/VV,
  σ0 +1.3 dB, pas de double-bounce) + **E2 (phase-linking échoue sur A, réussit
  sur C)**. Ce que ces quatre angles établissent **solidement**, c'est le rejet
  d'un **mouvement de tapis en corps rigide couplé à la nappe** (aucun indice
  positif : ni couplage hydrologique, ni stabilisation au gel, ni double-bounce).
  Ce qui **reste** est **cohérent avec une décorrélation volumétrique/
  diélectrique** d'un couvert humide sur substrat saturé — c'est l'hypothèse la
  mieux étayée. [HYP — non exclue] **Mais éliminer le corps rigide ne prouve pas
  le diélectrique** : une **troisième famille** demeure possible, un
  **micro-mouvement NON rigide** (flexion locale, tassement différentiel
  infra-pixel — ni déplacement en bloc, ni pur effet diélectrique), qui
  produirait aussi une décorrélation de phase sans signature mécanique de
  corps rigide. Seul un laser in situ la départagerait du diélectrique. La carte
  de tcoh confirme visuellement l'unité distincte : le polygone tourbière est
  uniformément bas, les taches à tcoh élevée sont les structures stables/champs
  — pas d'artefact visuel.

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
> atteint par le lac résiduel) ; (iii) la signature radar (ratio de
> polarisation croisée, rétrodiffusion) est cohérente avec une décorrélation
> **volumétrique/diélectrique** et exclut un mouvement de tapis en corps rigide
> (ni saisonnière-hydrologique, ni gelable, ni double-bounce). Le tapis flottant
> se comporte en bande C comme un **diffuseur de volume à phase quasi
> aléatoire**. La nature exacte de l'instabilité résiduelle — état diélectrique
> de la tourbe saturée vs micro-mouvement non rigide (flexion/tassement
> infra-pixel) — ne peut être tranchée par ces seules observations radar et
> demande une mesure in situ (laser) ou une diversité de longueur d'onde.

## 5-quinquies. Phases G et H — agrégation spatiale et prédiction de l'échec

### G — Changement d'OBSERVABLE (pas une 7e inversion)

Toutes les phases précédentes cherchaient une carte **pixel par pixel**. G
mesure **un seul nombre pour tout le tapis** : la moyenne complexe sur N pixels
divise le bruit de phase par √N_eff (~1.5 rad par pixel → ~1 mm sur 499 px),
donc *le signal n'est pas sous le plancher de bruit, il est sous le plancher de
bruit PAR PIXEL*. Validé en synthétique : sur des données identiques, l'inversion
par pixel donne −13.7 mm/an (36 % de pixels utilisables) là où l'agrégat donne
−19.8 mm/an pour une vérité de −20.

- **[FAIT] G1** — phase commune agrégée |R| : C=0.569, D=0.426, B=0.395,
  **A=0.234**. Le tapis est le plus bas de toutes les zones. [INT] Cohérent avec
  D→E2. *Réserve* : |R| contient l'atmosphère (commune à tous les pixels), donc
  seul le **classement relatif** est informatif.
- **[FAIT] G2** — double différence agrégée A−C vs **contrôle nul** (deux bandes
  adjacentes de sol stable) : réseau complet −1.53 vs nul **−1.50** mm/an ;
  baselines ≤60 j −3.87 vs nul **−4.88** mm/an. **Signal indiscernable du nul.**
  [INT] Aucun mouvement différentiel détectable ; le résultat publiable est une
  **borne supérieure**, pas une mesure. Que signal et nul coïncident à
  ~−1.5 mm/an indique un **résidu systématique commun** (rampe résiduelle), donc
  bien un plancher méthodologique.
- **[INT — correction de méthode] G2 testait la mauvaise observable.** La
  respiration d'une tourbière est **saisonnière** (10-40 mm, Hrysiewicz 2024),
  pas une dérive linéaire : une régression de vitesse sur un signal périodique
  donne ~0 *par construction* et n'a donc **aucune puissance** sur la physique
  recherchée. D'où `seasonal_amplitude` (G2-bis), qui ajuste un cycle annuel et
  compare son amplitude a celle du nul.
- **[FAIT] G2-bis — détection apparente, NON encore confirmée.** Amplitude
  saisonnière A−C = **3.29 mm** contre 0.57 mm pour le nul (×5.8) ;
  r²_saisonnier 0.30 vs 0.055 ; maximum vers le **jour 104** (mi-avril), ce qui
  est physiquement cohérent avec un pic de gonflement printanier a nappe haute.
- **[FAIT] G2-ter — LA DÉTECTION EST CONFIRMÉE, p = 0.021.** Contre 46 nuls
  **appariés en taille** (n_A=499, n_C=398 exactement) : nul médian 0.82 mm,
  p95 1.93 mm, **observé 3.29 mm** — aucun nul n'atteint l'observé.
  *Réserve technique* : p=0.0213 = 1/(1+46) est le **plancher** atteignable avec
  46 tirages ; passer a 100-200 tirages resserrerait (p ≥ 0.0099 a n=100). La
  détection est réelle, sa p-value est bornée par le nombre de tirages.
- **[LIMITE — LA PLUS IMPORTANTE] Signal saisonnier ≠ MOUVEMENT saisonnier.**
  3.3 mm ne valent que **0.75 rad** en bande C : un **cycle annuel d'humidité**
  (profondeur de pénétration variable) produit exactement ce signal de phase
  **sans aucun déplacement de surface**. Et le pic mi-avril est compatible avec
  les **deux** hypothèses — gonflement de la tourbe *et* humidité maximale
  suivent tous deux la nappe — donc **le timing ne discrimine rien**. Ce qui est
  établi est un **signal de phase différentiel saisonnier**, pas une respiration
  mécanique.
- **[FAIT — DÉCISIF] Le test du lac tranche : le signal est DIÉLECTRIQUE.** La
  zone B **ne peut pas respirer mécaniquement**. Or elle oscille : **2.63 mm,
  p = 0.036**, phase jour 94.7 — soit **80 % de l'amplitude du tapis (3.29 mm,
  p = 0.011, jour 104.2), a 10 jours près**. Tapis et lac partagent donc un
  **cycle saisonnier commun** qui n'est pas mécanique : un cycle annuel
  d'humidité / de profondeur de pénétration propre aux surfaces **saturées**,
  par contraste avec la prairie minérale sèche (C).
- **[FAIT — argument indépendant] L'ordre de grandeur exclut la mécanique.**

  | | amplitude attendue |
  |---|---|
  | tapis flottant libre suivant une nappe ±10 cm | ~100 mm |
  | respiration publiée (tourbières hautes, Hrysiewicz 2024) | 10-40 mm |
  | **mesuré ici** | **3.3 mm** |

  Le signal est **3 a 60× trop petit** pour une flottaison ou une respiration
  mécanique. Deux arguments indépendants (lac + magnitude) convergent.
- **[FAIT — la triangulation se referme] Tapis moins lac : le signal
  S'ANNULE.** En référençant le tapis **au lac** (au lieu de la prairie), le
  cycle commun disparaît :

  | comparaison | amplitude | phase (jour) | r²_sais | p |
  |---|---|---|---|---|
  | A − C (tapis vs prairie) | 3.29 mm | **104** | 0.30 | 0.011 |
  | B − C (lac vs prairie) | 2.63 mm | **95** | 0.11 | 0.036 |
  | **A − B (tapis vs lac)** | **0.90 mm** | 146 *(aléatoire)* | **0.05** | **0.448** |

  Les deux premières s'accordent en phase (mi-avril) ; la troisième tombe **au
  niveau du nul** (médiane 0.83 mm). [INT] Tapis et lac sont saisonnièrement
  **indiscernables** : si le tapis respirait mécaniquement et pas le lac, A−B
  révélerait cette respiration. Il ne révèle rien.
- **[FAIT — LE CHIFFRE À PUBLIER] Borne supérieure sur le mouvement propre du
  tapis : < 2 mm** (p95 du nul apparié pour A−B = 2.0 mm ; observé 0.90 mm).
  À comparer aux **10-40 mm** rapportés sur tourbières hautes (Hrysiewicz 2024).
  C'est une contrainte quantitative et falsifiable, bien plus forte qu'un
  « nous n'avons rien détecté ».
- **[INT] Ce que ça établit, et ce que ça n'établit pas.** Nous mesurons le
  **contraste diélectrique saisonnier** entre une tourbière saturée et une
  prairie sèche — un résultat en soi, et la première mesure positive de
  l'étude. Mais **ce n'est pas la respiration du tapis** et ne doit pas être
  présenté comme telle.
  ⚠️ *Distinction a maintenir* : ceci établit que le **signal saisonnier** est
  diélectrique ; cela ne dit **rien** sur la nature du mécanisme de
  **décorrélation**, qui reste (a) ou (b). Deux questions différentes.
- **[FAIT] G3 — biais de phase de fermeture : PRÉDICTION FALSIFIÉE.** Nous
  attendions qu'en passant de 304 a ~3000 triplets le biais de A devienne
  significatif (~5 σ). Or le réseau ne contient que **518 triplets fermés**, et
  a 518 le biais a *diminué* : −0.090 rad a **1.6 σ** (contre −0.136 a 1.77 σ
  sur 304). Une estimation qui régresse vers zéro quand n augmente est le
  comportement d'une **fluctuation**. **Aucun biais diélectrique systématique
  n'est détecté.**
- **[FAIT] G3 — la dispersion, elle, est robuste** : |fermeture| médiane
  A=0.683 et B=0.777 contre C=0.212 et D=0.210 (**×3.2**, stable entre les deux
  runs ; π/2≈1.57 = triplets purement aléatoires, donc A est *partiellement*
  cohérent, pas du bruit pur).
- **[INT — important] G3 ne tranche PAS (a) vs (b).** Une forte dispersion
  **sans biais de signe** = reconfiguration **aléatoire** du diffuseur, et non
  une dérive diélectrique monotone. Or des fluctuations d'humidité **comme** un
  micro-mouvement non rigide produisent exactement cette signature. G3 mesure
  donc le **degré** de non-stationnarité, pas sa **nature** : la question du
  mécanisme reste ouverte (laser ou bande L).

### H — Prédire OÙ Sentinel-1 échoue (analyse INTRA-tapis)

Passage de « le tapis décorrèle » a « la décorrélation survient sous telles
conditions ». Cible = temporal_coherence E2 en **continu** (le seuil 0.7 est une
convention arbitraire qui réduit 499 pixels a un seul nombre) ; analyse de la
variabilité **interne** de A, pas A vs C.

- **[FAIT] Courbe multi-seuils** : A et C ne diffèrent pas par un simple
  décalage — a 0.50 elles sont quasi identiques (0.974 vs 0.995), l'écart se
  creuse vers 0.65-0.75 (A 0.236→0.006, C 0.794→0.446). Le « 5.4 % vs 64.7 % »
  de E2 n'est qu'**un point** d'une courbe qui sépare surtout dans la **queue
  haute**.
- **[FAIT] Modèle intra-tapis** (n=499, 12 covariables dont la polarimétrie) :
  **R² validé croisé = 0.238 ± 0.029** (forêt aléatoire : 0.337), contre 0.127
  sans les covariables radar — le gain vient donc bien de la polarimétrie.
  Prédicteur dominant : **σ0_VV** (ρ = −0.379, 1ʳᵉ importance en forêt
  aléatoire).
- **[INT] Pouvoir prédictif RÉEL mais MODESTE** : ~24 % de la variance
  intra-tapis (33 % en non linéaire) ; **les trois-quarts restent inexpliqués**.
  Image physique (coefficients **nettoyés**) : la cohérence se dégrade la ou
  c'est plus **brillant** (σ0 −0.275), plus **diffusant en volume**
  (RVI −0.251), plus **humide** (−0.240) et plus **dynamique** phénologiquement
  (−0.107) ; `dist_edge_m` +0.101 indique une cohérence **meilleure près du
  bord** et se dégradant vers le **centre** du tapis — cohérent avec le profil
  radial de la Phase D-bis. Formulation défendable : *« S1 échoue davantage sur
  les parties du tapis les plus humides, les plus brillantes et les plus
  dynamiques »*.
- **[CORRECTION] L'altitude va dans le sens INVERSE de ce qui avait été écrit** :
  ρ = −0.168 et coefficient −0.131, donc plus haut = **moins** cohérent (et non
  l'inverse). Réserve : sur un tapis flottant quasi plat, le relief du MNT
  (Copernicus 30 m) est de l'ordre du bruit — cette variable est probablement un
  **proxy de position**, a ne pas interpréter physiquement.
- **[INT] La verdure moyenne est un PROXY, pas un moteur** : Spearman +0.320
  (2ᵉ rang) mais coefficient partiel **+0.029** (dernier rang) — elle ne prédit
  rien une fois σ0 et l'humidité pris en compte. Illustration directe du piège
  du Spearman marginal.
- **[INT] RVI : interprétable mais fragile.** Coefficient nettoyé **−0.251**
  (plus de volume → moins de cohérence), ce qui soutient l'hypothèse volumique ;
  *mais* son importance en forêt aléatoire est ~0.001, donc il n'apporte presque
  rien au-dela de σ0. À présenter comme **cohérent avec**, pas comme un
  prédicteur indépendant. Le seul prédicteur robuste dans les deux modèles est
  **σ0_VV**.
- **[FAIT — le plus frappant] Les moteurs S'INVERSENT entre tapis et prairie.**
  σ0_VV : ρ=**−0.379** dans A contre −0.008 dans C (écart −0.371) ; verdure
  moyenne : **+0.320** contre −0.009 (écart +0.329) ; altitude : −0.168 contre
  **+0.430** (écart **−0.598**). [INT] Les mêmes covariables agissent **en sens
  opposé** selon la zone : le tapis a une **physique propre**, il ne se comporte
  pas comme de la végétation ordinaire. Argument INDÉPENDANT en faveur de la
  Conclusion 1 (§4.1), par une voie que ni D, ni E2, ni G n'empruntaient.
- **[CORRIGÉ — piège de colinéarité]** Au 1er run, `rvi` et `vh_vv_db` avaient
  un Spearman **identique** (−0.1848) : ce sont deux transformations
  **monotones** du même rapport r = VH/VV (RVI = 4r/(1+r) ;
  vh_vv_db = 10·log₁₀ r), donc de rangs identiques et quasi parfaitement
  colinéaires. En régression, cela produisait deux coefficients **énormes et de
  signes opposés** (−1.21 et +0.95) qui n'ajustaient que le bruit entre deux
  variables identiques — **non publiables**. Ajout d'un diagnostic **VIF**
  (`collinearity_report`) et d'un retrait itératif (`drop_redundant`) ; seuls
  les coefficients du modèle **nettoyé** doivent être lus.
- **[MÉTHODE] RVI dual-pol** = 4·VH/(VV+VH) en puissance : préférable au ratio
  VH/VV brut de la Phase D-ter car **normalisé par la puissance totale**, donc
  insensible a un biais de calibration commun aux deux polarisations.

## 5-septies. Phase I — le signal diélectrique EST un capteur d'humidité

**Question :** le signal saisonnier de la Phase G suit-il l'hydrologie ?

### Validation préalable des zones (six vues, dont deux objectives)

- **[FAIT] Contrôle de surface** : A+B = **90.24 ha** contre **89.7 ha**
  documentés → **+0.6 %**. Vérification **numérique** du géoréférencement du
  masque, qu'aucune inspection visuelle ne peut fournir.
- **[FAIT]** Le polygone ressort dans **cinq capteurs indépendants** (cohérence,
  σ0, RVI, humidité S2, temporal_coherence), avec une **marche nette** au bord
  dans les trois profils radiaux, et des **distributions par zone disjointes**.
- **[FAIT]** Humidité S2 : **A = −0.513, C = −0.522** (jumeaux phénologiques
  confirmés) contre **B = +0.185** et σ0 = **−15.4 dB** (eau libre sans
  ambiguïté). Médianes par zone : σ0 A −10.09 / B −15.41 / C −11.22 / D −9.28 ;
  RVI A 0.914 / B 0.993 / C 0.881 / D 1.045.

### Le piège saisonnier, et sa résolution

Le balayage brut donnait `s2_wetness` r = +0.576 (lag 54 j) **et** `t2m_c`
r = −0.509 (lag 72 j), tous deux « significatifs ». [INT] Inexploitable en
l'état : température et humidité sont fortement anti-corrélées saisonnièrement,
`api_mm` (le proxy hydrologique le plus **direct**) échouait, et deux signaux a
cycle annuel corrèlent **toujours** a un certain décalage — le balayage ne fait
qu'aligner les phases (54 j = 53° de phase annuelle, pas un délai physique).

**Test décisif** : retirer l'harmonique annuelle des **deux** séries et ne
corréler que les **anomalies**. Validé en synthétique : sur deux cycles annuels
*indépendants*, |r| passe de 0.84 a 0.24 ; un couplage réel sur anomalies, lui,
survit (> 0.7).

### [FAIT] Résultat sur les anomalies

| forçage | r saisonnier | p | **r ANOMALIES** | **p** | lag |
|---|---|---|---|---|---|
| **s2_wetness** | 0.576 | 0.021 | **+0.450** | **0.011** | **12 j** |
| t2m_c | −0.509 | 0.021 | 0.224 | 0.581 | 78 j |
| api_mm | 0.230 | 0.426 | 0.293 | 0.172 | 6 j |
| precip_mm | 0.191 | 0.681 | −0.225 | 0.312 | 66 j |

### [INT] Trois faits convergents

1. **La température s'effondre** (0.509 → 0.224 ; p 0.021 → 0.581) : sa
   corrélation n'était **que** le cycle annuel partagé. Le confondant
   température/humidité est **levé** — on peut désormais attribuer le signal a
   **l'eau**.
2. **L'humidité survit** (r = +0.450, p ≤ 0.011 contre 92 nuls appariés en
   taille), et provient d'un **capteur optique totalement indépendant** du radar.
3. **Le décalage chute de 54 j a 12 j** — soit **un cycle de revisite**, donc
   *instantané* a notre résolution d'échantillonnage. C'était le critère posé
   d'avance : réponse **diélectrique** quasi instantanée vs tassement
   **mécanique** retardé de plusieurs semaines. **Le lag confirme le
   diélectrique par une voie indépendante du test du lac (Phase G).**

Le **signe** est cohérent : plus humide → pénétration moindre → centre de phase
plus haut → **soulèvement apparent** (r positif). (Le signe seul ne discrimine
pas, un gonflement mécanique donnerait le même ; c'est la magnitude et le lac
qui excluent le mécanique.)

### ⚠️ [LIMITE — DÉCISIVE] La détection est MARGINALE et le test de robustesse ÉCHOUE

Deux faits imposent la prudence, et interdisent pour l'instant de revendiquer un
capteur d'humidité :

1. **La détection est marginale.** Observé r = 0.450 contre un **nul p95 =
   0.404** : seulement **11 % au-dessus de la queue** de la distribution nulle.
   Le p = 0.0108 est de plus le **plancher** 1/(1+92) — avec davantage de
   tirages, des nuls pourraient dépasser l'observé.
2. **Le forçage différentiel échoue et le signe s'inverse.**
   `s2_wetness_diff` = NDWI(A) − NDWI(C) donne r = **−0.316**, p = **0.15**. Or
   la série InSAR est elle-même différentielle : si le lien était « le contraste
   d'humidité A−C pilote le contraste de phase A−C », ce forçage devrait
   **renforcer** la corrélation, pas la détruire.

**Modèle proposé pour l'expliquer — a posteriori, donc a tester et non a
croire.** Une humidité **régionale** commune M(t) piloterait φ(A) avec une
sensibilité k_A et φ(C) avec k_C, k_A > k_C (tourbe saturée vs prairie
minérale) ; alors φ(A)−φ(C) = (k_A−k_C)·M(t) : la phase **différentielle**
suivrait l'humidité **absolue**, via un contraste de *sensibilité* et non un
contraste d'humidité — et NDWI(A)−NDWI(C) ≈ 0 + bruit, les deux surfaces
répondant optiquement de façon similaire.

**Prédiction falsifiable** (Phase I, cellule 3.5) : si ce modèle est correct, le
NDWI de **C** et celui de **D** — proxys du même M(t) — doivent aussi corréler
**positivement**, avec un r du même ordre. Si **seul** le NDWI de A corrèle, le
modèle est **réfuté** et la corrélation est propre a la zone A, voire fortuite.

### [INT] Formulation prudente (en l'état actuel)

> Sur les anomalies (cycle annuel retiré), la phase agrégée co-varie avec
> l'humidité optique de surface (r = 0.45) a décalage quasi nul (12 j), alors
> que la température ne survit pas a la désaisonnalisation (0.51 → 0.22,
> p = 0.58). **Ce résultat est toutefois marginal** (nul p95 = 0.404) et **n'est
> pas confirmé par un forçage différentiel** ; il constitue un **indice**
> d'une sensibilité a l'humidité, cohérent avec le mécanisme diélectrique
> établi en Phase G, mais **pas une capacité de mesure démontrée**.

### [ACQUIS ROBUSTE de la Phase I] — indépendant du point précédent

L'élimination de la **température** est solide : sa corrélation apparente
(r = −0.509) disparaît a la désaisonnalisation (0.224, p = 0.58). Cela **écarte
un artefact thermique** et vaut en soi. De même, le **décalage quasi nul**
(12 j = un cycle de revisite) reste cohérent avec une réponse **diélectrique
instantanée** plutôt qu'un tassement mécanique retardé.

### [LIMITE] Autres réserves

- La **vraie nappe (WTD) in situ** reste supérieure au NDWI comme mesure de
  l'état hydrique, et sa structure temporelle diffère de celle de la
  température — c'est elle qui trancherait proprement.
- Augmenter le nombre de tirages nuls (200-500) pour sortir du plancher de
  p-value.

## 5-sexies. SYNTHÈSE GLOBALE — où nous en sommes après D→H

### La question a changé, et c'est ce qui a débloqué le travail

Question initiale : *« Sentinel-1 peut-il mesurer le déplacement vertical du
tapis flottant de Rzecin ? »* — question fermée, dont la réponse est non.

Question effective : ***« Que mesure réellement Sentinel-1 au-dessus d'une
tourbière flottante, et pourquoi ? »*** — question de recherche, dont nous avons
maintenant une réponse en quatre points.

### Les quatre acquis

**1. L'échec n'est pas algorithmique — c'est une propriété de la cible.**
Six estimateurs aux hypothèses mathématiques différentes (SBAS, ISBAS, paires
annuelles, réseau hybride, phase-linking EVD, WLS) échouent de la même façon.
Le phase-linking, estimateur du maximum de vraisemblance, ne récupère que 5.4 %
des pixels du tapis a tcoh ≥ 0.7 contre 64.7 % pour une végétation
**appariée** sur sol stable, **avec le même réseau**. `H1 fermée pour le réseau
HyP3 actuel.`

**2. Le tapis est une unité radar distincte, avec une physique PROPRE.**
Δcohérence = −0.069 a couvert apparié (Wilcoxon p = 2×10⁻⁴⁹, jackknife par date
stable) ; frontière nette au bord du polygone ; dispersion de fermeture ×3.2 ;
|R| agrégé le plus bas de toutes les zones. Et surtout, en Phase H, **les
prédicteurs s'inversent** entre tapis et prairie (σ0 : −0.379 vs −0.008 ;
altitude : −0.168 vs +0.430) — le tapis ne se comporte pas comme de la
végétation ordinaire.

**3. Ce que Sentinel-1 mesure ici est un signal DIÉLECTRIQUE, pas un
mouvement.** Amplitude 3.3 mm (p = 0.011), maximum mi-avril — mais le **lac
oscille identiquement** (2.63 mm, même phase) alors qu'il ne peut pas respirer,
et **A−B s'annule** (0.90 mm, p = 0.45). Trois arguments indépendants (contrôle
du lac, annulation A−B, ordre de grandeur 3-60× trop petit) convergent.

**4. Un INDICE de sensibilité a l'humidité (Phase I) — non démontré.** Sur les
**anomalies**, l'humidité optique S2 co-varie avec la phase agrégée
(r = +0.450, lag **12 j** = un cycle de revisite = instantané), tandis que la
**température s'effondre** (0.509 → 0.224, p = 0.58) — ce dernier point est
**robuste** et écarte un artefact thermique. **Mais la détection d'humidité est
marginale** (nul p95 = 0.404 pour 0.450 observé) et **n'est pas confirmée par un
forçage différentiel** (r = −0.316, p = 0.15, signe inversé). À présenter comme
un **indice cohérent** avec le mécanisme diélectrique de la Phase G, **pas**
comme une capacité de mesure établie.

**5. Borne supérieure quantitative : le mouvement propre du tapis est
< 2 mm**, contre 10-40 mm publiés sur tourbières hautes. Résultat falsifiable,
comparatif, et bien plus fort qu'une absence de détection.

### Les deux contributions méthodologiques

- **Changement d'observable** : un nombre pour tout le tapis plutôt qu'une carte.
  Le bruit de phase décroît en 1/√N_eff ; le signal n'était pas sous le plancher
  de bruit, il était sous le plancher **par pixel**. Validé en synthétique
  (par pixel −13.7 mm/an, agrégat −19.8 pour une vérité de −20).
- **Protocole de test de signal faible** : nul **apparié en taille** (le
  plancher dépend de N) + N tirages → **p-value empirique**. Sans lui, deux
  fausses conclusions auraient été publiées (un nul 4× trop grand, puis un test
  de vitesse sur un signal périodique).

### Les erreurs corrigées en cours de route (a documenter comme telles)

| erreur | conséquence évitée |
|---|---|
| `ratio_to_floor` comparé entre zones de N_eff différents | classement faux des zones |
| nul non apparié en taille (2200 px vs 499) | fausse détection saisonnière |
| test de **vitesse** sur un signal **périodique** | puissance nulle sur la physique cherchée |
| colinéarité `rvi`/`vh_vv_db` (VIF ≈ 240) | coefficients ininterprétables (−1.21 / +0.95) |
| prédiction « biais de fermeture a 5 σ » | **falsifiée** par les données (1.6 σ) |
| corrélation naïve sur deux cycles annuels (lag 54 j) | fausse attribution a l'eau — levée par désaisonnalisation |
| prédiction « aucun forçage ne survivra » | partiellement falsifiée : l'humidité survit, mais **marginalement** |
| forçage différentiel présenté comme « plus propre » | **échoue** (r=−0.316, p=0.15) → conclusion Phase I revue a la baisse |

### Ce qui reste ouvert

- **Le mécanisme de la DÉCORRÉLATION** — (a) diélectrique vs (b) micro-mouvement
  non rigide. Le biais de fermeture ne tranche pas (dispersion sans biais de
  signe = reconfiguration aléatoire, que les deux produisent). *Distinct de la
  nature du signal saisonnier, qui, elle, est établie comme diélectrique.*
- **Validation in situ** : le laser reste la seule mesure directe capable de
  confirmer la borne < 2 mm et de trancher (a) vs (b).
- **La bande L** (NISAR, Phase F prête) : teste si la longueur d'onde débloque
  la cible.
- **Transposabilité** du modèle prédictif (R²cv 0.24) a d'autres tourbières.

## 6. Portée et limites

- **Ce que ça établit :** le fen flottant est intrinsèquement moins observable
  en InSAR C-band que la végétation comparable sur sol stable — un résultat
  quantifié, spatialement structuré et relié à l'échec d'inversion, **confirmé
  par trois estimateurs (WLS, ISBAS, phase-linking EVD)**. C'est la
  contribution positive du projet, indépendante des échecs de traitement.
- **H1 (DS-InSAR/phase-linking) est fermée pour le réseau HyP3 actuel**
  (§5-quater) : le phase-linking ne récupère pas le tapis, à réseau contrôlé.
  L'échec est physique, pas algorithmique. (Un phase-linking full-SLC pourrait
  changer les valeurs absolues, pas l'écart relatif A vs C.)
- **Ce que ça ne tranche pas :** le mécanisme physique fin — **trois** familles
  restent en jeu, pas deux : (a) variabilité **diélectrique** de la tourbe
  saturée, (b) **micro-mouvement non rigide** (flexion locale, tassement
  différentiel infra-pixel), (c) — désormais **écartée** — mouvement de **corps
  rigide** couplé à la nappe. (a) et (b) restent toutes deux compatibles avec
  les observations ; ainsi que H4 (mouvement réversible). Le test de **biais de
  fermeture** (Phase G3), qui devait départager (a) de (b), **ne détecte aucun
  biais systématique** : la signature est une non-stationnarité **aléatoire**,
  compatible avec les deux. Seuls un laser in situ ou une diversité de longueur
  d'onde (bande L) départageraient (a) de (b).
- **Signal saisonnier agrégé (Phase G) — LE RÉSULTAT POSITIF DU PROJET, mais
  DIÉLECTRIQUE :** l'agrégation spatiale des 499 px du tapis, référencée au sol
  stable adjacent, détecte une **amplitude saisonnière de 3.3 mm** (max
  mi-avril, **p = 0.011** contre 92 nuls appariés en taille). C'est la
  **première mesure positive** de l'étude, obtenue la ou six inversions par
  pixel échouaient — parce que l'observable a changé (un nombre pour le tapis
  entier, pas une carte) et que la quantité testée est devenue la bonne
  (amplitude **saisonnière**, pas vitesse : régresser une vitesse sur un signal
  périodique donne ~0 par construction).
  **Mais ce n'est PAS une respiration mécanique** : le **lac oscille aussi**
  (2.63 mm, 80 % de l'amplitude du tapis, même phase) alors qu'il ne peut pas
  respirer, et 3.3 mm est **3 a 60× trop petit** pour une flottaison. Il s'agit
  du **contraste diélectrique saisonnier** entre tourbière saturée et prairie
  sèche. La **borne supérieure sur la respiration mécanique du tapis est donc
  plus contraignante que 3.3 mm**.
- **Pouvoir prédictif partiel (Phase H) :** ~24 % (R² validé croisé) a 34 %
  (non linéaire) de la variance intra-tapis s'explique par la rétrodiffusion
  σ0_VV, l'humidité, la verdure et la dynamique phénologique — réel mais
  modeste ; les trois-quarts restent inexpliqués.
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
3. **Bande L — NISAR (l'expérience décisive, désormais faisable) :** [FAIT,
   mise à jour 2026] les données L-band de **NISAR** (NASA-ISRO, λ = **24 cm**,
   couverture **globale** dont la Pologne, cycle 12 j) sont **publiques depuis
   le 20 juillet 2026** via ASF DAAC / Earthdata (record complet attendu fin
   2026 ; observations à partir de juin 2026 — donc test **prospectif**, non
   rétroactif sur 2022-2024). La bande L à 24 cm **pénètre la canopée** et voit
   la surface du tapis, là où le C-band (5.5 cm) décorrèle par diffusion de
   volume. **C'est le test qui départage (a) diélectrique de (b) micro-mouvement
   non rigide** : si le L-band récupère une cohérence exploitable sur la zone A
   (là où le C-band est au plancher), le facteur limitant est la **longueur
   d'onde/pénétration**, et un signal de déplacement du tapis devient mesurable.
   Alternative historique 2015-2024 : **ALOS-2 / PALSAR-2** (JAXA, L-band 24 cm)
   — même longueur d'onde mais **accès restreint/payant** (pas gratuit comme
   ASF). *NB : Sentinel-1C/1D sont **C-band** comme tout Sentinel-1 — ils
   rétablissent le cycle 6 j (moindre décorrélation temporelle) mais n'apportent
   PAS de bande L ; le 6 j n'adresse pas la décorrélation volumétrique/
   diélectrique, qui dépend de λ, pas de Δt.*

### 7-bis. Positionnement méthodologique (état de l'art opérationnel)

Notre phase-linking EVD (Phase E2) n'est pas une bricole isolée : le produit
opérationnel **OPERA DISP-S1** (NASA/JPL, via ASF) réalise exactement
« *phase linking on the sample coherence matrix* » dans une approche hybride
**PS + DS** — le même cœur algorithmique que notre EVD, mais en chaîne
opérationnelle. DISP-S1 valide donc la **direction** de la Phase E2. Il n'est
toutefois **pas utilisable sur Rzecin** : (i) couverture **Amérique du Nord
uniquement** (USA, territoires ≤ 200 km de la frontière, Canada→Panama) ;
(ii) **C-band** — il donnerait le même plancher de bruit sur le tapis. À citer
comme **référence méthodologique** confirmant que l'approche testée est celle de
l'état de l'art, pas comme une solution au site.

---

*Figures associées : outputs/phaseD/{zones, coherence_by_zone, decay_curves,
greenness_vs_coh}.png ; outputs/phaseDbis/{coh_vs_hydro, freeze_test,
radial_profile, residual_by_zone}.png ; sortie Phase E2 : phaseE2_evd.nc +
histogramme/carte tcoh par zone. Code : src/insar_wetlands/stratify.py,
src/insar_wetlands/inversion/phaselinking.py (EVD, sans ISCE) ;
tests/test_synthetic_phaselinking.py ; notebooks phaseD_inside_vs_outside.ipynb,
phaseDbis_mechanism_spatial.ipynb, phaseDter_scattering_scatterers.ipynb et
phaseE2_evd_phaselinking.ipynb.*
