# Guide — obtenir les corrections GACOS pour Rzecin

GACOS (*Generic Atmospheric Correction Online Service for InSAR*) fournit un
délai troposphérique zénithal total (ZTD) **spatialement variable** (0.125°,
composante stratifiée + turbulente) — la seule correction atmosphérique
susceptible d'améliorer notre petit site plat (l'ERA5 scalaire est annulé par
la calibration au point de référence).

## 1. Où

Formulaire en ligne : **http://www.gacos.net** (ou le miroir
http://ceg-research.ncl.ac.uk/v2/gacos/). Service gratuit, sur inscription
e-mail.

## 2. Paramètres à renseigner pour Rzecin

- **Region of interest (bounding box)** — englober le crop (buffer 2 km) :
  - North: 52.79, South: 52.74, West: 16.28, East: 16.34
  - (élargir légèrement ne coûte rien ; garder EPSG:4326 / degrés)
- **Acquisition UTC time : 16:36** — la trace S1 (orbite 175 ascendante) passe
  à ~16:36 UTC (cf. noms de granules `..._T1636xx_...`). NE PAS mettre 05:00.
- **Dates** : une par date S1 utilisée. Pour le réseau hybride, il faut TOUTES
  les dates qui apparaissent dans `hybrid_pairs.csv` (colonnes ref_date +
  sec_date, dédupliquées) — soit ~90 dates. Le formulaire accepte une liste
  de dates ou un fichier ; on peut aussi ne demander que les ~5 dates
  annuelles d'avril si on veut d'abord corriger les paires annuelles fragiles.

Pour extraire la liste des dates depuis le dépôt :
```python
import pandas as pd
p = pd.read_csv("hybrid_pairs.csv")
dates = sorted(set(p.ref_date) | set(p.sec_date))
for d in dates: print(pd.Timestamp(d).strftime("%Y-%m-%d"))
```

## 3. Ce que GACOS renvoie

Par e-mail, des liens de téléchargement : pour chaque date, un fichier
`YYYYMMDD.ztd` (binaire float32 little-endian, délai en **mètres**) + un
`YYYYMMDD.ztd.rsc` (géo-référencement : WIDTH, FILE_LENGTH, X_FIRST, Y_FIRST,
X_STEP, Y_STEP). GACOS indique aussi si la correction est **recommandée** pour
chaque date (qualité du modèle) — noter cette info.

## 4. Où les mettre + activation automatique

Déposer tous les `.ztd` + `.ztd.rsc` dans :
```
/content/drive/MyDrive/insar_rzecin/gacos/
```
Le notebook Phase A détecte automatiquement ce dossier (`USE_GACOS`) et
`atmosphere.apply_gacos_tropo` applique alors la correction spatiale par paire
(différence sec−ref, projetée en LOS). Aucune autre manip : relancer les
cellules §4→§5.

## 5. Attendu

Sur un site plat de <5 km, GACOS restera assez lisse spatialement (résolution
0.125° ≈ 14 km) — il corrige surtout la composante temporelle/stratifiée, pas
la turbulence sub-kilométrique. On DOCUMENTE donc son effet (comparaison
résidu avec/sans, sur pixels stables) plutôt que d'en attendre un miracle. Si
le gain est marginal, c'est un résultat en soi : l'atmosphère n'était pas le
facteur limitant ici (le réseau et la décorrélation le sont).
