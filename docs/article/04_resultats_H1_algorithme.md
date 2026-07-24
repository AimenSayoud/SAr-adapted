# 4. Résultats — H1 : l'échec est-il algorithmique ?

**Statut :** rédigé · **Sources :** `phase08`-`phase15`, `phaseA_hybrid`,
`phaseE2_evd_phaselinking` · **Figures :** F5 (tcoh par zone), F6 (carte tcoh)

> **H1.** *L'incapacité à mesurer un déplacement provient du choix de
> l'algorithme d'inversion.*

---

## 4.1 Six estimateurs, un même échec

| Méthode | Résultat sur le tapis |
|---|---|
| SBAS (MintPy) | aucun pixel exploitable |
| ISBAS | idem (pixels intermittents non sauvés) |
| Paires annuelles | idem |
| Réseau hybride | 14 238 pixels « résolus » sur l'emprise, **0 fiable** sur l'AOI (résidu médian 2.5 rad) |
| Moindres carrés pondérés | résidu médian **2.46 rad** (A) contre **1.92 rad** (C), à nombre de paires valides équivalent |
| **Phase-linking (EVD)** | voir §4.2 |

[INT] Le cas du réseau hybride est instructif : le gain de couverture apparent
(×78) était **creux** — les pixels « résolus » ne l'étaient pas sur la zone
d'intérêt. Un critère de couverture sans critère de fiabilité induit en erreur.

## 4.2 Le test décisif : le phase-linking

Le *phase-linking* est l'estimateur du **maximum de vraisemblance** sous modèle
gaussien circulaire : il pondère simultanément toutes les paires via la matrice
de cohérence complexe. S'il échoue, aucun estimateur du même réseau ne réussira.

**[FAIT] Cohérence temporelle par zone (356 paires, ~90 dates)**

| Zone | médiane | p25-p75 | **% pixels ≥ 0.7** |
|---|---|---|---|
| **C** — prairie appariée | **0.734** | 0.671-0.803 | **64.7 %** |
| D — autres couverts | 0.639 | 0.597-0.693 | 23.1 % |
| **A** — tapis flottant | **0.604** | 0.566-0.647 | **5.4 %** |
| B — lac résiduel | 0.584 | 0.542-0.630 | 1.5 % |

**[INT] Lecture par le plancher de bruit.** À cette redondance, un pixel
purement décorrélé rend ≈ 0.55.

- **B (lac) = 0.584** tombe **au plancher** — validation interne : on *sait* que
  l'eau végétalisée est décorrélée, et la méthode la classe correctement. La
  chaîne se comporte donc correctement sur la cible-témoin.
- **A (tapis) = 0.604** n'est qu'à ~0.05 au-dessus du plancher, juste à côté du
  lac. Ce n'est pas « un peu moins bon » que C : **A ≈ bruit, C ≈ signal**.
- **C = 0.734** est franchement au-dessus, avec 65 % de pixels de qualité DS —
  **sur exactement le même réseau**.

## 4.3 Analyse multi-seuils

Le seuil 0.7 est une convention ; la courbe complète est plus informative.

| seuil | A | B | C | D |
|---|---|---|---|---|
| 0.50 | 0.974 | 0.938 | 0.995 | 0.983 |
| 0.55 | 0.840 | 0.692 | 0.976 | 0.931 |
| 0.60 | 0.531 | 0.431 | 0.882 | 0.734 |
| 0.65 | 0.236 | 0.185 | 0.794 | 0.429 |
| **0.70** | **0.054** | 0.015 | **0.647** | 0.232 |
| 0.75 | 0.006 | 0.000 | 0.446 | 0.125 |
| 0.80 | 0.000 | 0.000 | 0.267 | 0.062 |

[INT] A et C sont quasi indiscernables à 0.50 (0.974 vs 0.995) : l'écart se
creuse dans la **queue haute** (≥ 0.65). Le tapis n'est donc pas globalement
décalé — il est **privé de ses meilleurs pixels**, ce qui est précisément ce qui
interdit l'inversion.

## 4.4 Conclusion sur H1

> **[INT — robuste] H1 est rejetée pour le réseau disponible.** Six estimateurs
> aux hypothèses mathématiques distinctes échouent identiquement, dont
> l'estimateur du maximum de vraisemblance. L'échec d'inversion est une
> **propriété physique de la cible**, pas un défaut d'algorithme.

**Portée exacte.** Un travail futur exploitant **toutes les SLC** (au-delà des
356 paires HyP3) pourrait modifier les *valeurs absolues* de cohérence
temporelle. Mais l'**écart relatif A vs C**, mesuré à réseau strictement égal,
devrait persister : C réussit là où A échoue, avec les mêmes paires, la même
sparsité et le même estimateur.

**Positionnement.** L'approche testée est celle de l'état de l'art
opérationnel : le produit **OPERA DISP-S1** (NASA/JPL) réalise un
*phase-linking* hybride PS+DS sur la matrice de cohérence — même cœur
algorithmique. Il n'est pas applicable ici (couverture Amérique du Nord, et
bande C, donc même limite), mais il confirme que la voie explorée n'est pas
marginale.
