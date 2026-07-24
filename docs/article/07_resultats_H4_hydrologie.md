# 7. Résultats — H4 : que mesure-t-on alors ?

**Statut :** rédigé · **Sources :** `phaseI_hydro_sensor` · **Figures :**
F12 (série InSAR vs humidité), F13 (significativité), F14 (saisonnier vs anomalies)

> **H4.** *Le signal diélectrique détecté porte de l'information sur l'état
> hydrique de surface.*

---

## 7.1 Le piège saisonnier

**[FAIT]** Balayage brut (série non désaisonnalisée) :

| forçage | r | lag | p |
|---|---|---|---|
| s2_wetness | +0.576 | 54 j | ≤ 0.021 |
| t2m_c | −0.509 | 72 j | ≤ 0.021 |
| api_mm | +0.230 | 6 j | 0.426 |
| precip_mm | +0.191 | 72 j | 0.681 |

**[INT] Inexploitable en l'état**, pour trois raisons :

1. La **température corrèle presque autant** que l'humidité, or les deux sont
   fortement anti-corrélées saisonnièrement → attribution impossible.
2. L'**API échoue** — c'est pourtant le proxy hydrologique le plus *direct*
   (mémoire des précipitations). Qu'il échoue pendant que les variables purement
   saisonnières passent est un signal d'alarme.
3. Les **décalages de 54-72 j ne sont pas des délais physiques** : deux signaux à
   cycle annuel corrèlent toujours à *un* décalage ; le balayage aligne les
   phases (54 j = 53° de phase annuelle).

## 7.2 Le test décisif : les anomalies

On retire l'**harmonique annuelle des deux séries** ; il ne reste que les
anomalies inter-annuelles et événementielles.

*(Validé en simulation : sur deux cycles annuels **indépendants**, |r| passe de
0.84 à 0.24 ; un couplage réel sur anomalies, lui, survit au-delà de 0.7.)*

**[FAIT]** Corrélations sur anomalies, contre 92 nuls appariés en taille :

| forçage | r saisonnier | **r ANOMALIES** | lag | **p** |
|---|---|---|---|---|
| **NDWI zone A** | 0.576 | **+0.450** | **12 j** | **≤ 0.011** |
| **NDWI zone C** | 0.490 | **+0.427** | **12 j** | 0.022 |
| **NDWI zone D** | 0.519 | **+0.424** | 42 j | ≤ 0.011 |
| NDWI(A) − NDWI(C) | 0.395 | −0.316 | 6 j | 0.150 |
| api_mm | 0.230 | 0.293 | 6 j | 0.172 |
| precip_mm | 0.191 | −0.225 | 66 j | 0.312 |
| **t2m_c** | **−0.509** | **0.224** | 78 j | **0.581** |

## 7.3 Trois faits convergents

**(a) La température s'effondre.** −0.509 → 0.224, p de 0.021 à 0.581. Sa
corrélation n'était **que** le cycle annuel partagé. Le confondant
température/humidité est **levé** : on peut désormais attribuer le signal à
l'eau. Ce point écarte aussi un **artefact thermique**.

**(b) L'humidité survit**, et provient d'un **capteur optique totalement
indépendant du radar** (Sentinel-2 : autre plateforme, autre physique de mesure).

**(c) Le décalage chute de 54 j à 12 j** — soit **un cycle de revisite**, donc
*instantané* à notre résolution d'échantillonnage. C'était le critère posé
**a priori** : réponse diélectrique quasi instantanée vs tassement mécanique
retardé de plusieurs semaines. **Le lag confirme le mécanisme diélectrique par
une voie indépendante du test du lac (§6).**

**Cohérence du signe.** r positif : plus humide → pénétration moindre → centre
de phase plus haut → soulèvement apparent. (Le signe seul ne discrimine pas — un
gonflement mécanique donnerait le même — mais il est cohérent.)

## 7.4 Le mécanisme : un contraste de sensibilité

L'échec du forçage **différentiel** (NDWI_A − NDWI_C : r = −0.316, p = 0.15) est
d'abord surprenant, la série InSAR étant elle-même différentielle.

**Modèle proposé.** Une humidité **régionale** commune M(t) pilote φ(A) avec une
sensibilité k_A et φ(C) avec k_C, avec **k_A > k_C** (tourbe saturée plus
réactive que prairie minérale). Alors :

> φ(A) − φ(C) = (k_A − k_C) · M(t)

La phase **différentielle** suit donc l'humidité **absolue**, via un contraste de
*sensibilité* et non un contraste d'humidité. Et NDWI(A) − NDWI(C) ≈ 0 + bruit,
les deux surfaces répondant optiquement de façon similaire.

**Prédiction falsifiable** (énoncée **avant** le test) : le NDWI de **C** et
celui de **D**, proxys du même M(t), doivent aussi corréler **positivement**,
avec un r du même ordre.

**[FAIT] Prédiction confirmée** : NDWI(A) +0.450, NDWI(C) +0.427, NDWI(D) +0.424
— amplitudes quasi identiques — tandis que le différentiel échoue.

**[INT]** Le modèle est donc validé par une **prédiction préalable**, non par une
rationalisation *a posteriori*. **L'échec du forçage différentiel n'infirme pas
le lien hydrologique : il en est une conséquence attendue.**

## 7.5 Conclusion sur H4

> **[INT] H4 est confirmée, avec un effet modéré.** Sur les anomalies, la phase
> agrégée co-varie avec l'humidité de surface (r = 0.42-0.45 selon la zone de
> référence, p ≤ 0.022) à **décalage quasi nul**, tandis que la température ne
> survit pas à la désaisonnalisation. Le couplage opère par un **contraste de
> sensibilité** entre tourbe saturée et prairie minérale.

## 7.6 Limites

- **Effet modéré** : 0.450 contre un nul p95 de 0.404 — à présenter comme une
  **sensibilité mesurable**, non comme un produit hydrologique opérationnel.
- Les trois NDWI zonaux **ne sont pas indépendants** (ils mesurent le même M(t)
  régional) : ce n'est pas une triple confirmation, mais **une** confirmation
  d'un motif prédit — dont l'échec du différentiel est la partie réfutable.
- Les p-values sont au **plancher** 1/(1+92) : écrire « p ≤ ».
- La **vraie nappe (WTD) in situ** reste supérieure au NDWI, et sa structure
  temporelle diffère de celle de la température — c'est elle qui trancherait
  proprement. Rzecin étant un site instrumenté, cette donnée existe
  probablement.
