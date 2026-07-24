# 5. Résultats — H2 : le tapis est-il une cible radar distincte ?

**Statut :** rédigé · **Sources :** `phaseD`, `phaseDbis`, `phaseDter`,
`phaseH`, `phaseI` (§2) · **Figures :** F2 (zones), F3 (profils radiaux),
F7 (distributions), F8 (prédicteurs A vs C)

> **H2.** *Le tapis flottant se comporte comme une cible radar particulière,
> distincte d'une végétation équivalente sur sol stable.*

---

## 5.1 Déficit de cohérence à couvert apparié

**[FAIT]** Comparaison appariée par interférogramme, A vs C :

| Mesure | Valeur |
|---|---|
| Cohérence moyenne A / C | 0.408 / 0.492 |
| **Δ apparié (coh_A − coh_C)** | **−0.069** |
| Wilcoxon signé | **p = 2.2 × 10⁻⁴⁹** |
| Fraction de paires de signe négatif | 88 % |
| **Jackknife par date** | Δ ∈ [−0.0705, −0.0652], **signe stable pour toute date retirée** |
| Temps de décorrélation τ | 21 j (A) contre 32 j (C) |

**[INT]** Le déficit ne dépend d'aucune acquisition particulière, et subsiste
après contrôle de la ligne de base, de l'atmosphère (appariement), de la pente
(MNT) et de l'humidité optique du couvert (features S2).

## 5.2 L'appariement de couvert est effectif

**[FAIT]** A et C sont des **jumeaux phénologiques** : humidité optique S2
médiane **−0.513** (A) contre **−0.522** (C) ; même classe ESA WorldCover
(prairie) ; verdure et amplitude saisonnière dans les mêmes gammes.

**[INT]** Le Δcohérence n'est donc **pas un artefact de couvert végétal**. La
différence vient d'une propriété de surface **non optique**, à l'échelle radar.

## 5.3 Une unité spatialement délimitée

**[FAIT]**
- **Profil radial** : plateau bas et plat (~0.40) sur tout l'intérieur, **marche
  nette au passage du bord**, pic juste à l'extérieur (~0.47). Même
  discontinuité en σ0 et en RVI.
- **Cinq capteurs indépendants** font ressortir le polygone : cohérence, σ0 VV,
  RVI, humidité S2, cohérence temporelle.
- **Surface** : A + B = 90.24 ha contre 89.7 ha documentés (**+0.6 %**).
- **Distributions par zone disjointes** (voir F7).

**[INT]** Ce n'est pas un bruit diffus mais une **unité physique délimitée**,
dont les contours — tracés à partir de sources vectorielles et optiques —
coïncident avec une structure visible dans des champs **radar indépendants**.

## 5.4 Signature de diffusion

**[FAIT]**

| Mesure | A (tapis) | C (prairie) | Lecture |
|---|---|---|---|
| σ0 VV médian | **−10.09 dB** | −11.22 dB | A **plus brillant** (+1.1 dB) |
| RVI (dual-pol) | **0.914** | 0.881 | A **plus dépolarisant** |
| Dispersion d'amplitude D_A | 0.243 (59 % < 0.25) | 0.238 (69 %) | A **n'est pas radar-sombre** |
| Dispersion de fermeture (\|closure\| médiane) | **0.683 rad** | 0.212 rad | **×3.2** |

**[INT]** Le tapis se comporte en bande C comme un **volume diffusant plus dense
et plus humide** que la prairie sèche, malgré une phénologie optique identique.
Le RVI plus élevé exclut un double-bounce de type eau libre. La dispersion de
fermeture ×3.2 mesure directement une **non-stationnarité du diffuseur** :
les triplets du tapis ne ferment pas, ceux du sol stable oui.

Point important : **A n'est pas dépourvu de cible** (59 % des pixels ont
D_A < 0.25). Son problème n'est pas l'absence de rétrodiffusion mais une
**phase instable** — ce qui justifiait a priori l'essai du phase-linking (§4).

## 5.5 Les prédicteurs de cohérence s'inversent entre zones

Analyse **intra-zone** (variabilité interne, pas A vs C) :

**[FAIT]** Modèle sur les 499 pixels du tapis, 12 covariables :
R² validé croisé = **0.239 ± 0.022** (forêt aléatoire 0.326), contre 0.127 sans
les covariables radar.

| covariable | ρ dans **A** (tapis) | ρ dans **C** (prairie) | écart |
|---|---|---|---|
| σ0 VV | **−0.379** | −0.008 | −0.371 |
| verdure moyenne | **+0.320** | −0.009 | +0.329 |
| altitude | −0.168 | **+0.430** | **−0.598** |

**[INT]** Les mêmes variables agissent **en sens opposé** selon la zone. C'est un
argument **indépendant** — obtenu par une voie que ni la comparaison appariée ni
le phase-linking n'empruntaient — que le tapis possède une **physique propre**
et ne se comporte pas comme de la végétation ordinaire.

**Précautions de lecture.**
- La colinéarité `RVI` / `VH-VV` (VIF ≈ 240, transformations monotones du même
  rapport) produisait des coefficients aberrants (−1.21 et +0.95) ; seuls les
  coefficients du modèle **nettoyé** sont interprétables (σ0 −0.275, RVI −0.251,
  humidité −0.240).
- La **verdure est un proxy** : ρ = +0.320 en marginal mais coefficient partiel
  +0.029 — elle ne prédit rien une fois σ0 et l'humidité pris en compte.
- L'**altitude** est probablement un proxy de position : sur un tapis quasi
  plat, le relief du MNT (30 m) est de l'ordre du bruit. À ne pas interpréter
  physiquement.
- Le seul prédicteur robuste dans les deux modèles (linéaire et non linéaire)
  est **σ0 VV**.

## 5.6 Conclusion sur H2

> **[INT — robuste] H2 est confirmée.** À couvert apparié et après contrôle de
> la ligne de base, de l'atmosphère, de la pente et de l'humidité optique, le
> tapis flottant présente une cohérence significativement inférieure, une
> frontière nette, une signature de diffusion volumique et une
> non-stationnarité ×3.2 — et les facteurs qui gouvernent sa cohérence
> **s'inversent** par rapport à la prairie. C'est une **unité radar distincte**,
> et non « de la végétation en bande C ».

**Pouvoir prédictif partiel.** ~24 % (33 % en non linéaire) de la variance
intra-tapis s'explique : la cohérence se dégrade là où c'est plus brillant, plus
diffusant en volume, plus humide et plus dynamique phénologiquement. Les
trois-quarts restent inexpliqués — le modèle identifie une direction, pas une
loi.
