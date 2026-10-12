# Phase C — Tester le phase-linking / DS-InSAR (la méthode jamais essayée)

## Pourquoi

Tout notre traitement (HyP3 + MintPy + ISBAS maison) repose sur des
interférogrammes **déjà multi-vus**, traités indépendamment. Nous n'avons
JAMAIS testé le **phase-linking / DS-InSAR** (SqueeSAR-like), la classe de
méthode conçue précisément pour les **cibles distribuées à faible cohérence**
(végétation, tourbe) : elle estime la phase optimale de chaque pixel à partir
de la matrice de cohérence de TOUTES les paires du stack SLC, récupérant des
pixels que le SBAS classique abandonne. C'est le principal levier « méthode »
non exploré (hypothèse H1 du doc verification_et_challenge.md).

## Distinction importante (honnête)

| Outil | Nature | Entrée requise | Effort | Ce que ça teste |
|---|---|---|---|---|
| **LiCSBAS** | SBAS + fermeture de boucle, masquage cohérence | Produits **LiCSAR** (interférogrammes prêts) | **Faible-moyen** | Un SBAS mieux outillé (déroulement, masquage) récupère-t-il le signal ? **PAS du phase-linking.** |
| **MiaplPy** | **Phase-linking / DS-InSAR** vrai | Stack SLC **corégistré** (format ISCE) | **Élevé** | La vraie méthode DS. Nécessite de reprocesser les SLC. |
| **FRInGE** | Phase-linking (EVD/EMI) | Stack SLC corégistré (ISCE) | **Élevé** | Idem, alternatif. |

**Notre stack HyP3 = interférogrammes multi-vus, PAS des SLC.** Donc :
- **LiCSBAS est faisable rapidement** (produits LiCSAR tout prêts) → à faire
  d'abord ; il teste l'hypothèse « meilleur déroulement/masquage » sans
  reprocesser.
- **Le vrai phase-linking (MiaplPy/FRInGE) exige un stack SLC corégistré**
  (ISCE2 `stackSentinel.py` sur les SLC bruts) → lourd (compute + stockage),
  à planifier comme étape suivante si LiCSBAS est prometteur.

## Étape C1 — LiCSBAS (à faire en premier)

### C1.1 Trouver le frame LiCSAR couvrant Rzecin
- Portail COMET LiCSAR : **https://comet.nerc.ac.uk/comet-lics-portal/**
- Chercher par coordonnées (52.763, 16.312). Noter le(s) **frame ID(s)**
  ascendant(s) — idéalement la même trace que nous (orbite 175 asc.), ex.
  `NNN A_FFFF_...`. LiCSAR fournit interférogrammes + cohérence + E-W/U geo.

### C1.2 Installer et lancer LiCSBAS (Colab)
```bash
pip install LiCSBAS   # ou: git clone https://github.com/comet-licsar/LiCSBAS
# Workflow standard :
LiCSBAS01_get_geotiff.py -f <FRAME_ID> -s 20220101 -e 20241231   # télécharge le frame
LiCSBAS02_ml_prep.py -i GEOC -o GEOCml10 -n 10                    # multilook 10
LiCSBAS03op_GACOS.py  # (optionnel, si GACOS déposé) correction atmosphérique
LiCSBAS04op_mask.py   # masquage
LiCSBAS05op_clip.py -g <lon1/lon2/lat1/lat2>                      # clip sur Rzecin
LiCSBAS11_check_unw.py    # qualité de déroulement
LiCSBAS12_loop_closure.py # FERMETURE DE BOUCLE -> corrige les erreurs de déroulement
LiCSBAS13_sb_inv.py       # inversion SBAS -> série temporelle + vitesse
LiCSBAS16_filt_ts.py      # filtrage
```
**Points clés vs notre pipeline :** LiCSBAS applique la **fermeture de boucle**
(LiCSBAS12) pour corriger les erreurs de déroulement — exactement le problème
diagnostiqué chez nous (changement de signe des paires 2 ans). Il intègre aussi
la correction GACOS nativement (LiCSBAS03).

### C1.3 Comparer à notre résultat
- Extraire la vitesse LiCSBAS sur l'AOI (même geojson) et au pixel de référence.
- Comparer : couverture (pixels valides), vitesse médiane AOI, et surtout le
  **résidu / RMS de la série** vs nos 2–2.9 rad.
- Verdict : si LiCSBAS (avec fermeture de boucle) donne une série cohérente là
  où nous avons 0 pixel fiable → **c'était le déroulement/la méthode (H1)**,
  pas le site. Sinon → indice fort pour « site difficile » (H3), et on passe au
  vrai phase-linking (C2) comme test ultime.

## Étape C2 — Vrai phase-linking (MiaplPy), si C1 est prometteur ou non concluant

### C2.1 Prérequis (lourd)
- Télécharger les **SLC Sentinel-1** (SLC IW, VV) de la trace 175 asc.,
  2022–2024 (~90 acquisitions) depuis ASF/CDSE.
- Corégistrer en stack ISCE2 : `stackSentinel.py -s SLC -d DEM -a AUX -o ORBITS
  -w work -c 1 -p vv` → génère le stack corégistré + baselines.

### C2.2 MiaplPy
```bash
pip install miaplpy
miaplpyApp.py miaplpy.cfg   # phase-linking (EMI/EVD) sur le stack SLC
# -> estime la phase DS optimale, puis inversion type MintPy
```
`miaplpy.cfg` : pointer sur le stack ISCE, définir la fenêtre de phase-linking
(ex. 15×15), le seuil de similarité temporelle. Sortie : série temporelle DS
+ masque de qualité de phase-linking (temporal coherence).

### C2.3 Le test décisif de la méthode
Comparer la **couverture DS** et le **résidu** à notre ISBAS et à LiCSBAS. Le
phase-linking est censé récupérer les pixels distribués végétalisés : si LUI
non plus ne donne rien de cohérent sur le tapis, alors l'hypothèse « site
physiquement difficile » (H3) devient la plus robuste — un résultat fort et
publiable (« même le DS-InSAR échoue sur ce fen flottant »).

## Ordre recommandé et coût

1. **EGMS** (Phase B) — gratuit, immédiat, contexte externe.
2. **LiCSBAS** (C1) — faible-moyen, produits prêts, teste la fermeture de
   boucle + GACOS. **Le meilleur rapport info/effort.**
3. **MiaplPy** (C2) — élevé (reprocesser les SLC), seulement si nécessaire ;
   c'est le test ultime « méthode vs site ».

Chaque étape est conçue pour FALSIFIER « site difficile ». Si les trois
échouent de façon cohérente, la conclusion devient robuste ; si l'une réussit,
on a notre signal — et la réponse à la question depuis le début.
