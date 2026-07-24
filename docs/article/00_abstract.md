# Résumé

**Statut :** rédigé · **Sources :** toutes phases · **Cible :** ~250 mots

---

## Résumé (version longue de travail)

Les tourbières stockent un carbone considérable et leur surface se déforme avec
le niveau de nappe (« respiration » de 10-40 mm/an documentée sur tourbières
hautes). L'InSAR Sentinel-1 permettrait un suivi gratuit et continu de ce
signal, mais son applicabilité aux **tourbières flottantes** (*Schwingmoor*),
plus humides et à couvert bas dense, n'est pas établie.

Nous testons systématiquement l'InSAR bande C sur le fen flottant de Rzecin
(Pologne, 89.7 ha) sur 356 interférogrammes Sentinel-1 (2022-2024), en
organisant l'analyse autour de quatre hypothèses concurrentes.

**(H1) L'échec n'est pas algorithmique.** Six estimateurs aux hypothèses
distinctes — SBAS, ISBAS, paires annuelles, réseau hybride, moindres carrés
pondérés, et *phase-linking* par décomposition en valeurs propres — échouent
identiquement. Le phase-linking, estimateur du maximum de vraisemblance, ne
récupère que **5.4 %** des pixels du tapis à cohérence temporelle ≥ 0.7 contre
**64.7 %** pour une végétation appariée sur sol stable, **avec le même réseau**.

**(H2) Le tapis est une cible radar distincte.** À couvert apparié
(ESA WorldCover + phénologie Sentinel-2), sa cohérence est significativement
inférieure (Δ = **−0.069**, Wilcoxon p = 2×10⁻⁴⁹, jackknife par date stable), la
frontière est nette, et les prédicteurs de cohérence **s'inversent** entre
tapis et prairie.

**(H3) Le signal saisonnier est diélectrique, pas mécanique.** Une agrégation
spatiale (moyenne complexe sur 499 pixels) extrait une amplitude saisonnière de
**3.3 mm** (p ≤ 0.011 contre nuls appariés en taille) que six inversions par
pixel ne voyaient pas. Mais le **lac résiduel oscille identiquement** (2.63 mm,
même phase) alors qu'il ne peut pas respirer, la différence tapis−lac
**s'annule** (0.90 mm, p = 0.45), et l'amplitude est **3 à 60× trop faible**
pour une flottaison. Borne supérieure sur le mouvement propre du tapis :
**< 2 mm**.

**(H4) Ce que l'InSAR mesure est une sensibilité à l'humidité.** Sur les
anomalies désaisonnalisées, la phase agrégée co-varie avec l'humidité optique
Sentinel-2 (r = 0.42-0.45, p ≤ 0.022) à **décalage quasi nul** (12 j), tandis
que la température ne survit pas à la désaisonnalisation. Le couplage opère via
un **contraste de sensibilité** entre tourbe saturée et prairie minérale —
prédiction confirmée par l'échec attendu d'un forçage différentiel.

Nous concluons que la limitation est **physique et non méthodologique**, et
proposons deux apports transposables : le **changement d'observable** (agrégat
plutôt que carte, le signal étant sous le plancher de bruit *par pixel*) et un
**protocole de test de signal faible** (nul apparié en taille + p-value
empirique) qui a invalidé deux conclusions intermédiaires.

---

## Points clés (bullet highlights)

- Six méthodes d'inversion InSAR échouent identiquement sur une tourbière
  flottante : la limitation est physique, pas algorithmique.
- L'agrégation spatiale extrait un signal saisonnier de 3.3 mm invisible pixel
  par pixel, mais trois tests indépendants montrent qu'il est **diélectrique**.
- Borne supérieure de **2 mm** sur le mouvement de surface du tapis, contre
  10-40 mm publiés sur tourbières hautes.
- La phase agrégée répond aux **anomalies d'humidité** à décalage quasi nul, via
  un contraste de sensibilité entre surfaces.
- Protocole reproductible de test de signal faible en terrain décorrélé.

---

## Mots-clés

InSAR ; Sentinel-1 ; bande C ; tourbière flottante ; *Schwingmoor* ;
décorrélation ; phase-linking ; diffuseurs distribués ; humidité de surface ;
agrégation spatiale.
