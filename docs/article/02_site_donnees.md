# 2. Site d'étude et données

**Statut :** rédigé · **Sources :** `phase01`, `phase02`, `phaseD`, `phaseI`
· **Figures :** F1 (localisation), F2 (zones A/B/C/D)

---

## 2.1 Site

**Tourbière de Rzecin**, Grande-Pologne (52.7632 °N, 16.3098 °E), **89.7 ha**.

Fen de transition (*poor fen* / tourbière de transition) à **tapis flottant de
*Sphagnum*** (*Schwingmoor*), comportant un **lac résiduel** en cours de
terrestrialisation. Caractéristiques pertinentes pour le radar :

- **Nappe quasi affleurante** : 0-30 cm sous la surface, hydrologiquement
  stable (Juszczak et al. 2013).
- **Couvert bas dense** : *Sphagnum*, cypéracées, éricacées — volume diffusant
  important en bande C.
- **Substrat mobile** : le tapis repose sur l'eau, d'où l'hypothèse de départ
  d'un mouvement vertical marqué.

C'est donc, a priori, le site où le **signal géomorphologique attendu est
maximal** — ce qui en fait un cas-test exigeant.

## 2.2 Données radar

| | |
|---|---|
| Capteur | Sentinel-1A, mode IW, SLC bursts |
| Polarisation | VV (interférométrie) ; VV+VH (rétrodiffusion) |
| Trace | 175, orbite ascendante |
| Burst | `175_374052_IW1` (couverture AOI vérifiée) |
| Période | 2022-01-01 → 2024-12-31 |
| Cadence | 12 j (période **S1A seul** : S1B hors service, S1C pas encore opérationnel) |
| Interférogrammes | **356 paires**, ~90 dates |
| Traitement | ASF **HyP3** `INSAR_ISCE_BURST` (corégistration, déroulement SNAPHU, géocodage) |
| Grille | UTM, ~40 m, emprise recadrée 129 × 138 px |

> **Note de portée.** La fenêtre 2022-2024 correspond à l'ère **S1A seul**
> (12 j). Le retour d'une constellation à deux satellites (S1C lancé fin 2024,
> S1D en 2025) rétablit la cadence 6 j : la décorrélation temporelle sera
> moindre pour de futures études. Cela ne change en revanche **pas la longueur
> d'onde**, dont dépend la décorrélation volumique (voir §8).

## 2.3 Données auxiliaires

| Source | Usage | Volume |
|---|---|---|
| **Sentinel-2 L2A** | NDWI/MNDWI → verdure, humidité de surface, appariement de couvert | 68-69 dates |
| **Sentinel-1 RTC** (Microsoft Planetary Computer) | σ0 VV/VH, ratio, **RVI**, dispersion d'amplitude | 85 dates |
| **ERA5** | précipitations, T2m → API, test de gel | journalier |
| **ESA WorldCover 10 m** | classe de couvert pour l'appariement | 1 tuile |
| **MNT Copernicus** (via HyP3) | pente (contrôle), altitude | statique |

## 2.4 Stratification en zones

Quatre zones sont définies sur la grille radar :

| Zone | Définition | n px | Surface |
|---|---|---|---|
| **A** | tapis végétalisé — intérieur du polygone, hors eau | 499 | 79.84 ha |
| **B** | lac résiduel — intérieur, fraction inondée > 0.30 | 65 | 10.40 ha |
| **C** | **prairie appariée** — extérieur, même classe WorldCover que A, features S2 dans [p10, p90] de A, pente < 5° | 398 | 63.68 ha |
| **D** | autres couverts extérieurs (contexte, réservoir de nuls) | 10 750 | 1 720 ha |

**La zone C est le cœur du dispositif** : c'est un témoin *apparié* en couvert
végétal et en phénologie, ce qui permet d'isoler ce qui est propre au tapis de
ce qui relève de « la végétation en bande C » en général.

### Validation objective du masque

- **A + B = 90.24 ha** contre **89.7 ha** documentés → **écart de +0.6 %**.
  Vérification *numérique* du géoréférencement, qu'aucune inspection visuelle ne
  peut fournir.
- **A et C sont des jumeaux phénologiques** : humidité optique S2 médiane
  **−0.513** (A) contre **−0.522** (C). L'appariement est donc effectif.
- **B est sans ambiguïté de l'eau** : σ0 VV = **−15.41 dB** (spéculaire) et
  humidité S2 = **+0.185**, très au-dessus de toutes les autres zones.

### Signature multi-capteurs des zones

| Zone | cohérence | σ0 VV (dB) | RVI | humidité S2 | coh. temporelle |
|---|---|---|---|---|---|
| A (tapis) | 0.408 | −10.09 | 0.914 | −0.513 | 0.604 |
| B (lac) | 0.396 | **−15.41** | 0.993 | **+0.185** | 0.584 |
| C (prairie) | **0.492** | −11.22 | 0.881 | −0.522 | **0.734** |
| D (autres) | 0.438 | −9.28 | 1.045 | −0.440 | 0.639 |

Le polygone ressort dans **cinq capteurs indépendants**, avec des distributions
par zone disjointes et une **marche nette** au passage du bord dans les profils
radiaux (voir [figures.md](figures.md), F2-F3).
