# 9. Conclusion

**Statut :** rédigé · **Sources :** synthèse

---

Nous avons testé systématiquement l'applicabilité de l'InSAR Sentinel-1 en bande
C au suivi du déplacement vertical d'une **tourbière flottante** — le cas où le
signal géomorphologique attendu est le plus fort, et où la méthode n'avait pas
été évaluée frontalement.

## Quatre résultats

**1. La limitation est physique, pas algorithmique.** Six estimateurs aux
hypothèses mathématiques distinctes — jusqu'au *phase-linking* par maximum de
vraisemblance — échouent identiquement. Sur le même réseau de 356 paires, le
tapis rend **5.4 %** de pixels à cohérence temporelle ≥ 0.7 contre **64.7 %**
pour une végétation **appariée** sur sol stable. La charge de la preuve se
déplace : il n'y a pas lieu de chercher un septième algorithme.

**2. Le tapis est une cible radar distincte.** À couvert et phénologie appariés,
sa cohérence est significativement inférieure (Δ = −0.069, p = 2 × 10⁻⁴⁹,
jackknife par date stable), sa frontière est nette, sa diffusion est volumique,
sa non-stationnarité de fermeture est **3.2×** supérieure — et les facteurs qui
gouvernent sa cohérence **s'inversent** par rapport à la prairie voisine.

**3. Le signal saisonnier mesurable est diélectrique, et le mouvement est
borné.** Une agrégation spatiale extrait une amplitude saisonnière de **3.3 mm**
(p ≤ 0.011) invisible pixel par pixel. Mais le lac résiduel — qui ne peut pas
respirer — oscille identiquement (2.63 mm, même phase), la différence
tapis − lac s'annule (0.90 mm, p = 0.45), et l'amplitude est **3 à 60× trop
faible** pour une flottaison. **Le mouvement propre du tapis est < 2 mm**,
contre 10-40 mm publiés sur tourbières hautes.

**4. Ce que l'InSAR mesure ici est une sensibilité à l'humidité de surface.**
Sur les anomalies désaisonnalisées, la phase agrégée co-varie avec l'humidité
optique Sentinel-2 (r = 0.42-0.45, p ≤ 0.022) à **décalage quasi nul** (12 j),
alors que la température ne survit pas à la désaisonnalisation. Le couplage
opère via un **contraste de sensibilité** entre tourbe saturée et prairie
minérale — prédiction confirmée par l'échec attendu d'un forçage différentiel.

## Réponse à la question posée

> **Que mesure réellement l'InSAR Sentinel-1 en bande C au-dessus d'une
> tourbière flottante ?**
>
> Pas le mouvement de la surface, qui reste sous 2 mm — mais l'**état hydrique**
> de cette surface, au travers d'un effet de **profondeur de pénétration**. La
> tourbière flottante se comporte en bande C comme un **volume diffusant humide
> à phase quasi aléatoire**, dont le centre de phase se déplace avec l'humidité
> et non avec le substrat.

## Apports méthodologiques

- Le **changement d'observable** : agréger spatialement une cible qui se déforme
  comme une unité fait tomber le bruit en 1/√N_eff et révèle un signal
  inaccessible pixel par pixel.
- Un **protocole de test de signal faible** — nul apparié en taille,
  distribution nulle plutôt que réalisation unique, traitement identique du nul
  — qui a invalidé deux conclusions intermédiaires et réfuté deux prédictions
  des auteurs.
- Une implémentation **légère du phase-linking**, applicable directement à des
  produits interférométriques standards sans chaîne SLC.

## Perspectives

La voie la plus directe est la **bande L**. À λ = 24 cm, le radar pénètre la
canopée et atteint la surface du tapis. Les données **NISAR** sont désormais
globales et gratuites, et leur produit GUNW est l'analogue direct de nos
interférogrammes : la chaîne développée ici s'y applique sans modification. Le
retour de la constellation Sentinel-1 à une cadence de 6 jours réduira la
décorrélation temporelle mais **ne changera pas la longueur d'onde**, dont
dépend la décorrélation volumique.

Une validation **in situ** — laser de surface et nappe mesurée — reste
nécessaire pour confirmer la borne de 2 mm et pour départager, dans le mécanisme
de décorrélation, la variabilité diélectrique du micro-mouvement non rigide.
