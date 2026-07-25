# 8. Discussion

**Statut :** rédigé · **Sources :** synthèse toutes phases + littérature

---

## 8.1 Pourquoi la bande C réussit sur tourbières hautes et échoue ici

Hrysiewicz et al. (2024) mesurent la respiration de tourbières hautes
irlandaises en InSAR Sentinel-1 avec r = 0.8-0.9. Nos résultats ne contredisent
pas les leurs : ils délimitent le **domaine de validité** de la méthode.

| | tourbière haute (*raised bog*) | fen flottant (Rzecin) |
|---|---|---|
| Couvert | *Sphagnum* relativement **sec** | *Sphagnum* + cypéracées, **saturé** |
| Nappe | plus profonde, plus variable | **quasi affleurante**, stable |
| Diffuseurs | de **surface**, stables | **volumiques**, non stationnaires |
| Substrat | tourbe consolidée | tapis reposant sur l'eau |
| Cohérence C-band | exploitable | 5.4 % des pixels ≥ 0.7 |

[INT] Le facteur discriminant n'est pas « tourbière ou non » mais **l'humidité
du couvert et la stationnarité des diffuseurs**. Sur un tapis saturé, le centre
de phase est dominé par un volume végétal humide dont la configuration change
entre deux passages, ce qui décorrèle en 12 j à λ = 5.5 cm.

**Contribution de portée générale :** le succès de la bande C sur tourbières ne
se transpose **pas** automatiquement aux tourbières flottantes, qui sont
pourtant celles où le signal attendu est le plus fort.

## 8.1-bis Comparaison avec les tourbières drainées : deux grandeurs différentes

Patil et al. (2026, *RSASE* 41) mesurent en InSAR Sentinel-1 (2015-2025) la
subsidence du **Great Fen** (Cambridgeshire, RU), tourbière basse **drainée**
en cours de restauration : **0.48 à 1.40 cm/an**. Il est tentant de confronter
ces valeurs à notre « < 2 mm ». **Ce serait une erreur de grandeur** : les deux
chiffres ne mesurent pas la même chose.

| | Patil et al. 2026 | Cette étude |
|---|---|---|
| Grandeur | **vitesse** de subsidence | **amplitude** d'un cycle annuel |
| Unité | mm **par an** | mm (crête à moyenne) |
| Valeur | 4.8 à 14.0 mm/an | < 2 mm |
| Projection | verticale (déclarée) | **LOS** (voir ci-dessous) |
| Durée | 10 ans | 3 ans |

### Deux précisions indispensables

**(a) Nous rapportons la ligne de visée (LOS), pas la verticale.** À un angle
d'incidence Sentinel-1 IW typique (~39°), la conversion sous hypothèse de
mouvement purement vertical est d_vert = d_LOS / cos(θ) ≈ **1.29 × d_LOS**.
Notre borne devient donc **< 2.6 mm en vertical**. Nous conservons le LOS dans
le texte parce que l'hypothèse de verticalité pure n'est pas vérifiable ici.

**(b) La comparaison pertinente porte sur la VITESSE, que nous avons aussi
mesurée.** Notre vitesse différentielle A−C est de **−1.53 mm/an** — mais le
contrôle nul donne **−1.50 mm/an** : le signal est **indiscernable du bruit**.
Notre plancher de détection en vitesse est donc de l'ordre de **1.5 à 5 mm/an**
selon le réseau.

### Ce que la comparaison enseigne réellement

| Site | État hydrologique | Subsidence |
|---|---|---|
| Fermes restaurées tardivement (Great Fen) | drainé, en restauration | 14.0 mm/an |
| Fermes restaurées tôt (Great Fen) | drainé, en restauration | 11.7 mm/an |
| Réserves naturelles (Holme/Woodwalton) | **conservé, plus humide** | **4.8 mm/an** |
| **Rzecin (cette étude)** | **naturel, saturé, non drainé** | **non détecté** (< plancher) |

[INT] Le gradient est cohérent et **hydrologiquement interprétable** : plus le
site est drainé, plus il s'affaisse. Rzecin, tourbière naturelle à nappe
affleurante et jamais drainée, se situe **en-dessous du site le moins affaissé**
de leur série — ce qui est exactement l'attendu écologique, et non un échec de
mesure. Notre non-détection est donc **cohérente avec la littérature**, pas en
contradiction avec elle.

[INT] Cette comparaison délimite aussi le **domaine d'application** de la
méthode : l'InSAR bande C mesure la subsidence des tourbières **drainées**
(surfaces agricoles, cohérence élevée, signal de plusieurs cm/an) et échoue sur
les tourbières **saturées à tapis flottant** (couvert humide, faible cohérence,
signal millimétrique).

### Une mise en garde méthodologique que nos résultats suggèrent

Patil et al. rapportent que les fluctuations saisonnières s'alignent sur les
anomalies d'humidité du sol et de précipitations, et y voient un « fort contrôle
hydrologique du mouvement de surface de la tourbe ».

[INT] Nos résultats invitent à la prudence sur ce type d'inférence : **une
oscillation saisonnière corrélée à l'humidité n'est pas automatiquement un
mouvement**. Sur notre site, un signal de 3.3 mm parfaitement corrélé à
l'humidité s'est révélé **diélectrique** — le lac résiduel, qui ne peut pas
respirer, oscillait de la même amplitude et à la même phase. Un **contrôle sur
une surface en eau** (ou toute cible dont le mouvement est physiquement exclu)
est peu coûteux et permet de séparer le déplacement réel de l'effet de
profondeur de pénétration. Nous suggérons de l'intégrer systématiquement aux
études de mouvement de surface en tourbière, en particulier lorsque le signal
rapporté est de l'ordre du millimètre à la dizaine de millimètres.

## 8.2 Deux apports méthodologiques transposables

### (a) Changement d'observable

Le bruit de phase par pixel (~1.5 rad à γ = 0.4) décroît en 1/√N_eff sous
agrégation. Sur 499 pixels, il tombe à ~0.3 mm (~1 mm à N_eff = 50), très
en-dessous du signal cherché. **Le signal n'était pas sous le plancher de bruit,
il était sous le plancher de bruit *par pixel*.**

Ce raisonnement s'applique à toute cible **spatialement cohérente mais
temporellement décorrélée** : tourbières, glaciers rocheux, zones humides,
cultures. La condition est que la cible se déforme comme une unité — hypothèse
à justifier physiquement, et testable en découpant la zone en sous-régions.

### (b) Protocole de test de signal faible

Trois règles, qui ont chacune invalidé une conclusion intermédiaire :

1. **Nul apparié en taille.** Le bruit d'un agrégat décroît en 1/√N : un nul
   4× plus grand a ~2× moins de bruit et **fabrique de fausses détections**.
2. **Distribution nulle, pas réalisation unique.** Une seule réalisation n'est
   pas un test ; N tirages donnent une p-value empirique — dont il faut rappeler
   le **plancher** 1/(1+N).
3. **Traitement identique du nul.** Si l'observé résulte d'une sélection
   (meilleur |r| sur 16 décalages), le nul doit subir le même balayage.

Ces règles sont peu coûteuses et devraient accompagner toute revendication de
signal faible en terrain décorrélé.

## 8.3 Erreurs corrigées en cours d'analyse

Nous les documentons délibérément : elles montrent que le protocole a résisté
aux attentes des auteurs.

| Erreur | Conséquence évitée |
|---|---|
| Rapport au plancher comparé entre zones de N_eff différents | classement erroné des zones |
| Contrôle nul non apparié en taille (2 200 px contre 499) | fausse détection saisonnière |
| Test d'une **vitesse** sur un signal **périodique** | puissance nulle sur la physique cherchée |
| Colinéarité RVI / VH-VV (VIF ≈ 240) | coefficients ininterprétables (−1.21 / +0.95) |
| Corrélation naïve entre deux cycles annuels | fausse attribution à l'eau |
| **Prédiction** « biais de fermeture à ~5 σ » | **falsifiée** par les données (1.6 σ) |
| **Prédiction** « aucun forçage ne survivra » | **falsifiée** (l'humidité survit) |

[INT] Deux prédictions explicites des auteurs ont été **réfutées par leurs
propres données**, et une explication *a posteriori* a été soumise à un test
falsifiable avant d'être retenue. C'est le fonctionnement attendu d'un protocole
robuste.

## 8.4 Ce qui reste ouvert

**Le mécanisme de la décorrélation.** Deux familles restent compatibles :
(a) variabilité **diélectrique** de la tourbe saturée ; (b) **micro-mouvement
non rigide** (flexion locale, tassement infra-pixel). Une troisième — mouvement
de **corps rigide** couplé à la nappe — est écartée (absence de couplage
hydrologique de la cohérence, absence de stabilisation au gel, absence de
signature double-bounce). Le biais de fermeture, qui devait départager (a) de
(b), ne détecte aucun biais systématique : une dispersion élevée **sans biais de
signe** est compatible avec les deux.

⚠️ Cette question est **distincte** de celle du signal saisonnier, dont
l'origine diélectrique est établie (§6).

**Validation in situ.** Un **laser** au sol reste la seule mesure directe
capable de confirmer la borne < 2 mm et de trancher (a) vs (b).

**La nappe mesurée.** Rzecin est un site instrumenté ; une série WTD in situ
transformerait le test H4, sa structure temporelle différant de celle de la
température.

## 8.5 Perspectives instrumentales

**La bande L est la voie la plus directe.** À λ = 24 cm, le radar **pénètre la
canopée** et atteint la surface du tapis, là où la bande C (5.5 cm) décorrèle par
diffusion de volume. **NISAR** (NASA-ISRO) fournit désormais des données L-band
**globales et gratuites** (publiques depuis juillet 2026), dont le produit
**GUNW** est l'analogue direct de nos interférogrammes — l'ensemble de la chaîne
développée ici s'y rebranche sans modification. Test **prospectif** : l'archive
ne couvre pas rétroactivement 2022-2024.

**Le retour à 6 jours ne suffira pas.** Sentinel-1C et 1D rétablissent la cadence
à deux satellites, ce qui réduit la décorrélation **temporelle**. Mais tout
Sentinel-1 reste en **bande C** : la décorrélation **volumique** dépend de λ, pas
de Δt. Nous ne nous attendons donc pas à ce que le 6 jours débloque ce site.

**Alternative historique.** ALOS-2/PALSAR-2 (L-band, 2015-2024) permettrait un
test rétroactif sur notre fenêtre, mais l'accès est restreint.

## 8.6 Limites de l'étude

- **Un seul site, une seule trace, une seule polarisation** (VV) : la
  transposabilité du modèle prédictif (R²cv 0.24) à d'autres tourbières reste à
  démontrer.
- **Aucune validation in situ** (ni laser, ni WTD) : la borne < 2 mm est
  interne à l'InSAR.
- **Fenêtre S1A seul** (12 j) : cadence dégradée par rapport à ce qui sera
  disponible.
- **p-values au plancher** (1/(1+N)) pour plusieurs tests : augmenter le nombre
  de tirages nuls resserrerait les bornes.
- **Hypothèse d'unité** de l'agrégation : le tapis est supposé se déformer en
  bloc ; à vérifier par découpage en sous-régions si un signal mécanique
  apparaissait.
