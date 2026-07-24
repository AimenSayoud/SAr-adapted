# 1. Introduction

**Statut :** rédigé · **Sources :** revue de littérature, `revue_critique_litterature.md`

---

## 1.1 Enjeu

Les tourbières couvrent ~3 % des terres émergées mais stockent environ un tiers
du carbone organique des sols. Leur bilan carbone dépend étroitement du niveau
de la nappe, lui-même lié à un mouvement vertical de la surface : la tourbe
gonfle quand elle se sature et s'affaisse quand elle s'assèche
(*Mooratmung* / « respiration de la tourbière »). Mesurer ce mouvement
renseigne donc indirectement sur l'état hydrologique et, à terme, sur le
fonctionnement carbone.

L'interférométrie radar satellitaire offre en principe la mesure idéale :
gratuite, continue, spatialisée, avec une sensibilité millimétrique. Sentinel-1,
en particulier, fournit une couverture systématique depuis 2014.

## 1.2 État de l'art — et la tension qu'il contient

**La bande C fonctionne sur tourbières hautes.** Hrysiewicz et al. (2024, *RSE*
291) mesurent la respiration de tourbières hautes irlandaises en InSAR
Sentinel-1, avec r = 0.8-0.9 contre données in situ, RMSE ~7 mm/an, et
récupèrent des amplitudes de 10-40 mm. Ce résultat **réfute** l'idée répandue
que la bande C serait inutilisable sur tourbière.

**Mais toutes les tourbières ne se ressemblent pas.** Les tourbières hautes
(*raised bogs*) portent un couvert de *Sphagnum* relativement sec, à diffuseurs
de surface stables. Les **tourbières flottantes** (*Schwingmoor*, fens de
transition) présentent un tapis végétal reposant sur une nappe quasi
affleurante : couvert bas dense, humidité permanente, substrat mobile. Rien ne
garantit que le succès obtenu sur les premières se transpose aux secondes.

**La question n'a pas été tranchée frontalement.** La littérature sur InSAR et
zones humides porte majoritairement sur la détection d'inondation par
double-bounce, la cohérence comme indicateur de couvert, ou les tourbières
drainées/exploitées. Le cas du tapis flottant — celui où le signal
géomorphologique attendu est le plus fort — est peu documenté.

## 1.3 Le piège méthodologique que nous voulons éviter

Un travail qui échoue à mesurer un déplacement peut toujours être attribué à un
mauvais choix de traitement. La littérature InSAR offre un large éventail
d'algorithmes (SBAS, PS, *phase-linking*/DS, réseaux à baselines variées) et il
est tentant, face à un échec, d'en essayer un de plus indéfiniment.

Nous adoptons la démarche inverse : **si des estimateurs reposant sur des
hypothèses mathématiques différentes échouent de la même manière, la charge de
la preuve se déplace** — l'explication n'est probablement pas algorithmique. Le
travail consiste alors à caractériser la cible, pas à chercher un septième
algorithme.

## 1.4 Reformulation de la question

Nous ne demandons donc pas :

> *L'InSAR Sentinel-1 peut-il mesurer le déplacement vertical du tapis flottant
> de Rzecin ?*

— question fermée dont la réponse serait un simple oui/non — mais :

> **Que mesure réellement l'InSAR Sentinel-1 en bande C au-dessus d'une
> tourbière flottante, et pourquoi ?**

## 1.5 Structure de l'étude

Nous organisons l'analyse autour de **quatre hypothèses concurrentes**,
testables et falsifiables :

| | Hypothèse | Test principal | Section |
|---|---|---|---|
| **H1** | L'échec vient de l'**algorithme** d'inversion | six estimateurs indépendants sur le même réseau | [§4](04_resultats_H1_algorithme.md) |
| **H2** | Le **tapis** est une cible radar distincte | comparaison à couvert apparié + validation multi-capteurs | [§5](05_resultats_H2_cible.md) |
| **H3** | Le signal résiduel est un **mouvement** de surface | agrégation spatiale + contrôles (lac, nul apparié, magnitude) | [§6](06_resultats_H3_mecanisme.md) |
| **H4** | Le signal mesuré traduit l'**état hydrique** | corrélation sur anomalies désaisonnalisées | [§7](07_resultats_H4_hydrologie.md) |

## 1.6 Contributions

1. Une **démonstration multi-méthodes** que la limitation est physique et non
   algorithmique, sur un site où le signal attendu est maximal.
2. Une **borne supérieure quantitative** (< 2 mm) sur le mouvement de surface,
   directement comparable aux 10-40 mm publiés sur tourbières hautes.
3. La démonstration que le signal saisonnier détectable est **diélectrique**,
   établie par trois arguments indépendants.
4. Deux apports méthodologiques transposables : le **changement d'observable**
   (agrégation spatiale) et un **protocole de test de signal faible**
   (nul apparié en taille + p-value empirique).

## Références à citer ici

- Hrysiewicz, A., et al. (2024). *Remote Sensing of Environment*, 291 — InSAR
  C-band sur tourbières hautes.
- Juszczak, R., et al. (2013) — hydrologie et flux du site de Rzecin.
- Fornaro, G., et al. — estimateurs de phase-linking (EMI/EVD).
- Ferretti, A., et al. (2011) — SqueeSAR / diffuseurs distribués.
- De Zan, F., et al. (2015) ; Ansari, H., et al. (2021) — phase de fermeture et
  biais lié à l'humidité du sol.
- Kellndorfer, J., et al. (2022). *Scientific Data* — cohérence globale S1.
