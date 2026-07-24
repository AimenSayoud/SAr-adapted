# 6. Résultats — H3 : mouvement de surface ou effet diélectrique ?

**Statut :** rédigé · **Sources :** `phaseG_aggregation` · **Figures :**
F9 (série agrégée vs nul), F10 (test de significativité), F11 (fermeture)

> **H3.** *Le signal résiduel, invisible pixel par pixel, traduit un mouvement
> vertical du tapis.*

---

## 6.1 Le changement d'observable

Toutes les analyses précédentes cherchaient une **carte** pixel par pixel. Nous
mesurons ici **un seul nombre pour tout le tapis**.

**[FAIT — validation synthétique]** Sur des données simulées identiques,
l'inversion par pixel rend **−13.7 mm/an** (36 % de pixels utilisables) là où
l'agrégat rend **−19.8 mm/an** pour une vérité de **−20**.

**[INT]** Le signal n'était pas sous le plancher de bruit : il était sous le
plancher de bruit **par pixel**.

## 6.2 Existe-t-il une phase commune ? (|R|)

**[FAIT]** Module du vecteur résultant agrégé, par zone :

| Zone | \|R\| médian |
|---|---|
| C — prairie | **0.569** |
| D — autres | 0.426 |
| B — lac | 0.395 |
| **A — tapis** | **0.234** |

**[INT]** Le tapis a la phase commune **la plus faible de toutes les zones**,
cohérent avec H2.

⚠️ **Deux précautions.** (i) |R| contient l'**atmosphère**, commune à tous les
pixels : dépasser le plancher n'est donc pas une preuve de signal de sol — seul
le **classement relatif** est informatif. (ii) Le rapport au plancher
1/√N_eff **n'est pas comparable entre zones**, car N_eff varie de ~16 (B) à
~2 700 (D).

## 6.3 La vitesse était la mauvaise observable

**[FAIT]** Vitesse différentielle A − C contre contrôle nul apparié :

| réseau | signal A−C | nul (plancher) |
|---|---|---|
| complet | −1.53 mm/an | **−1.50 mm/an** |
| baselines ≤ 60 j | −3.87 mm/an | **−4.88 mm/an** |

**[INT]** Signal **indiscernable du nul**. Mais surtout : la respiration d'une
tourbière est **saisonnière**, pas une dérive. Régresser une **vitesse** sur un
signal périodique donne ~0 *par construction* — ce test n'avait **aucune
puissance** sur la physique recherchée. C'est une erreur de conception que le
contrôle nul a permis d'identifier.

## 6.4 Amplitude saisonnière : détection

**[FAIT]** Ajustement d'un cycle annuel sur la série agrégée A − C, testé contre
**92 nuls appariés en taille** :

| | amplitude | phase (jour) | r² saisonnier | p |
|---|---|---|---|---|
| **A − C** | **3.29 mm** | **104** (mi-avril) | 0.30 | **≤ 0.011** |

Nul : médiane 0.82 mm, p95 1.93 mm. Aucun des 92 nuls n'atteint l'observé.

**[INT]** Détection réelle, obtenue là où six inversions par pixel ne voyaient
rien. Le maximum de mi-avril est cohérent avec un gonflement printanier à nappe
haute.

⚠️ p = 0.0108 = 1/(1+92) est le **plancher** atteignable avec 92 tirages.

## 6.5 Trois arguments indépendants excluent le mouvement

### (a) Le lac oscille aussi

**[FAIT]** Le lac résiduel **ne peut pas respirer mécaniquement**. Or :

| | amplitude | phase (jour) | p |
|---|---|---|---|
| A − C (tapis) | 3.29 mm | 104 | ≤ 0.011 |
| **B − C (lac)** | **2.63 mm** | **95** | 0.036 |

Soit **80 % de l'amplitude du tapis, à 10 jours près de la même phase**.

### (b) La différence tapis − lac s'annule

**[FAIT]** En référençant A **au lac** plutôt qu'à la prairie :

| | amplitude | phase | r² sais. | p |
|---|---|---|---|---|
| **A − B** | **0.90 mm** | 146 *(aléatoire)* | 0.05 | **0.448** |

Nul médian : 0.83 mm. **Le signal tombe exactement au plancher.**

**[INT]** Tapis et lac sont saisonnièrement **indiscernables**. Si le tapis
respirait mécaniquement et pas le lac, A − B révélerait cette respiration. Il ne
révèle rien.

### (c) L'ordre de grandeur

| | amplitude attendue |
|---|---|
| tapis flottant libre suivant une nappe de ±10 cm | ~100 mm |
| respiration publiée, tourbières hautes (Hrysiewicz 2024) | 10-40 mm |
| **mesuré ici** | **3.3 mm** |

**[INT]** Le signal est **3 à 60× trop faible** pour une flottaison ou une
respiration mécanique.

## 6.6 Le biais de fermeture ne discrimine pas

**[FAIT]** Sur les 518 triplets fermés du réseau :

| Zone | biais moyen (rad) | σ | \|closure\| médiane |
|---|---|---|---|
| A | −0.090 | 1.6 | **0.683** |
| B | −0.088 | 1.5 | **0.777** |
| C | +0.027 | 1.1 | 0.212 |
| D | −0.021 | 1.1 | 0.210 |

**[INT]** **Aucun biais systématique n'est détecté.** (Une prédiction préalable
annonçait ~5 σ en augmentant le nombre de triplets ; le réseau n'en contient que
518, et à 518 l'estimation a *diminué* — comportement d'une fluctuation. La
prédiction est **falsifiée**.)

Ce qui reste robuste est la **dispersion** (×3.2, stable entre runs ; π/2 ≈ 1.57
correspondrait à des triplets purement aléatoires, donc A reste *partiellement*
cohérent). Mais une forte dispersion **sans biais de signe** traduit une
reconfiguration **aléatoire** du diffuseur, que produisent *aussi bien* des
fluctuations d'humidité qu'un micro-mouvement non rigide. **Ce test mesure le
degré de non-stationnarité, pas sa nature.**

## 6.7 Conclusion sur H3

> **[INT — robuste] H3 est rejetée.** Le signal saisonnier détecté (3.3 mm,
> p ≤ 0.011) est d'origine **diélectrique** : le lac, qui ne peut pas respirer,
> oscille identiquement ; la différence tapis − lac s'annule ; et l'amplitude est
> 3 à 60× trop faible pour une flottaison. Nous mesurons un **contraste
> d'humidité saisonnier** entre surfaces saturées et prairie sèche.

> **[FAIT] Borne supérieure : le mouvement propre du tapis est < 2 mm**
> (p95 du nul apparié pour A − B), contre 10-40 mm publiés sur tourbières
> hautes. Résultat falsifiable et comparatif, bien plus fort qu'une absence de
> détection.

⚠️ **Distinction à maintenir.** Ceci établit que le **signal saisonnier** est
diélectrique. Cela ne dit rien sur la nature du mécanisme de **décorrélation**,
qui reste indéterminé entre variabilité diélectrique et micro-mouvement non
rigide (voir §8).
