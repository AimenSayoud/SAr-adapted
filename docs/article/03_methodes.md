# 3. Méthodes

**Statut :** rédigé · **Sources :** `src/insar_wetlands/{inversion,aggregate,
stratify,predict_failure,hydro_link}.py` · **Figures :** F4 (schéma du protocole)

---

## 3.1 Estimateurs d'inversion comparés (H1)

Six approches, reposant sur des hypothèses **mathématiquement distinctes** :

| # | Méthode | Hypothèse propre |
|---|---|---|
| 1 | SBAS (MintPy) | réseau à petites baselines, inversion globale |
| 2 | ISBAS | tolère les pixels intermittents (sous-réseau par pixel) |
| 3 | Paires annuelles | contourne la décorrélation saisonnière par appariement d'état |
| 4 | Réseau hybride | combine baselines courtes et longues |
| 5 | Moindres carrés pondérés | pondération par cohérence, paire à paire |
| 6 | **Phase-linking (EVD)** | **maximum de vraisemblance** sur la matrice de cohérence complète |

Le sixième est le plus puissant : il exploite **simultanément** toutes les
paires via la matrice de cohérence complexe N × N par pixel, dont l'historique
de phase est le vecteur propre dominant. Implémenté en numpy directement sur les
interférogrammes HyP3 (qui *sont* les entrées de cette matrice), sans ISCE ni
SLC. Qualité = **cohérence temporelle**, analogue de `MiaplPy.temporalCoherence`.

**Plancher de bruit.** Sur un réseau de redondance ~4 (356 paires / 89 dates),
un pixel *purement décorrélé* rend une cohérence temporelle ≈ **0.55** (établi
par simulation). Toute valeur voisine de 0.55 doit donc être lue comme « bruit »,
non comme « signal faible ».

## 3.2 Comparaison à couvert apparié (H2)

Test **apparié par interférogramme** : chaque paire est observée dans A *et*
dans C, ce qui neutralise ligne de base perpendiculaire et atmosphère du jour.

**Statistique.** Wilcoxon signé sur les différences coh_A − coh_C. Les 356
paires n'étant **pas indépendantes** (elles partagent ~90 dates), un bootstrap
par paire sous-estimerait l'incertitude ; nous rapportons donc un **jackknife
par date** (retrait de chaque acquisition et de toutes ses paires) : le résultat
n'est retenu que si le signe est stable pour *toute* date retirée.

## 3.3 Agrégation spatiale — le changement d'observable (H3)

**Motivation quantitative.** L'écart-type de phase d'un pixel à γ = 0.4 est
~1.5 rad. La moyenne complexe sur N pixels le divise par √N_eff : sur les 499
pixels du tapis, ~0.07 rad (~0.3 mm) ; même avec N_eff = 50, ~1 mm — bien
en-dessous de la respiration de 10-40 mm recherchée. **Le signal n'est pas sous
le plancher de bruit : il est sous le plancher de bruit *par pixel*.**

**Hypothèse physique qui l'autorise :** le tapis est une unité hydrologique ; il
respire en bloc. Moyenner ne détruit donc pas le signal (contrairement à un
champ de déformation hétérogène) mais élimine la composante aléatoire.

**Trois observables :**

1. **|R|** — module du vecteur résultant moyen, `R = Σ w·exp(iφ) / Σ w` avec
   w = cohérence. Phases aléatoires → |R| ≈ 1/√N_eff. Calculé sur la phase
   **enroulée**, donc insensible aux erreurs de déroulement.
2. **Double différence agrégée A − C** — mêmes paire et échelle spatiale
   (~1 km) : atmosphère et orbite s'annulent, il reste le mouvement
   différentiel.
3. **Amplitude saisonnière** — ajustement de `y = c + d·t + a·cos(2πt) +
   b·sin(2πt)`, amplitude = √(a²+b²). **C'est la bonne observable** : la
   respiration d'une tourbière est saisonnière, pas une dérive ; régresser une
   *vitesse* sur un signal périodique donne ~0 par construction.

## 3.4 Protocole de test de signal faible

C'est l'apport méthodologique central, et il a **invalidé deux conclusions
intermédiaires**.

**Contrôle nul apparié en taille.** Le bruit d'un agrégat décroît en 1/√N. Un
nul construit sur 2 200 pixels alors que la zone testée en compte 499 possède
~2× moins de bruit et **sous-estime le plancher**, fabriquant de fausses
détections. Les nuls sont donc des **taches compactes adjacentes de sol stable,
d'effectifs exactement égaux** à ceux des zones comparées.

**p-value empirique.** Une réalisation nulle unique n'est pas un test. On tire
N nuls indépendants et on calcule p = (1 + #{nul ≥ observé}) / (1 + N).
⚠️ Cette p-value a un **plancher** 1/(1+N) : avec 92 tirages, p ne peut
descendre sous 0.0108. Toute valeur au plancher doit s'écrire « p ≤ ».

**Traitement identique du nul.** Lorsque la statistique observée résulte d'une
sélection (meilleur |r| sur ~16 décalages), le nul subit **le même balayage**,
sans quoi la comparaison est biaisée.

## 3.5 Discrimination du mécanisme (H3)

- **Contrôle du lac.** Le lac résiduel **ne peut pas respirer mécaniquement**.
  S'il présente la même amplitude saisonnière que le tapis, le signal est
  d'origine diélectrique.
- **Différence tapis − lac.** Référencer A à B annule tout cycle commun aux
  surfaces saturées ; il ne subsiste que le mouvement propre du tapis.
- **Ordre de grandeur.** Un tapis flottant libre suivant une nappe de ±10 cm se
  déplacerait de ~100 mm ; la respiration publiée est de 10-40 mm.
- **Biais de phase de fermeture.** Un déplacement, même non rigide, ferme les
  triplets à zéro ; une dérive diélectrique monotone les biaiserait
  (De Zan 2015 ; Ansari 2021). Calculé sur phase **enroulée** — distinct de la
  détection d'erreurs de déroulement (multiples de 2π).

## 3.6 Modèle prédictif intra-tapis (H2/H4)

Cible = cohérence temporelle **en continu** (le seuil 0.7 est une convention
arbitraire qui réduirait 499 pixels à un seul nombre) ; analyse de la
variabilité **interne** de A.

- **Spearman** (monotone, robuste) pour le classement marginal.
- **Régression multiple standardisée** pour les effets *partiels*.
- **VIF** obligatoire avant toute lecture de coefficient : `RVI = 4r/(1+r)` et
  `VH/VV (dB) = 10·log₁₀ r` sont deux transformations **monotones** du même
  rapport, donc quasi parfaitement colinéaires (VIF ≈ 240) — d'où des
  coefficients aberrants et de signes opposés si l'on n'y prend pas garde.
- **R² validé croisé** (5 blocs), jamais le R² d'ajustement.
- Forêt aléatoire en complément non linéaire, jamais comme argument principal.

## 3.7 Lien hydrologique (H4)

**Le piège de l'autocorrélation.** Série InSAR et forçage hydrologique sont
fortement autocorrélés : une p-value de corrélation naïve est massivement trop
optimiste. On réutilise donc les **séries nulles appariées en taille**, qui
partagent la même structure temporelle — la distribution nulle absorbe
l'autocorrélation par construction.

**Le piège saisonnier.** Deux signaux à cycle annuel corrèlent **toujours**
fortement à *un* décalage : le balayage ne fait qu'aligner les phases (54 j =
53° de phase annuelle, pas un délai physique). La conclusion causale n'est donc
tirée que **sur les anomalies désaisonnalisées** (harmonique annuelle retirée
des deux séries).

**Interprétation du décalage résiduel.** Une réponse diélectrique à l'humidité
est quasi instantanée ; un tassement mécanique suivrait la nappe avec des
semaines de retard.
