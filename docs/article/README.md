# Article — Rzecin / InSAR bande C sur tourbière flottante

**Titre de travail**

> *What does Sentinel-1 C-band InSAR measure over a floating peatland?
> Multi-method evidence for a dielectric-dominated signal and an upper bound on
> mat motion*

Version française : *Que mesure l'InSAR Sentinel-1 en bande C au-dessus d'une
tourbière flottante ? Faisceau multi-méthodes pour un signal d'origine
diélectrique et borne supérieure sur le mouvement du tapis.*

---

## Principe d'organisation

L'article est structuré **par hypothèses concurrentes**, pas chronologiquement.
Un lecteur ne doit pas revivre l'ordre des phases : il doit voir quelles
explications ont été testées et lesquelles ont survécu.

| Fichier | Contenu | Statut |
|---|---|---|
| [`00_abstract.md`](00_abstract.md) | Résumé + points clés | rédigé |
| [`01_introduction.md`](01_introduction.md) | Contexte, littérature, question | rédigé |
| [`02_site_donnees.md`](02_site_donnees.md) | Site, jeux de données, zones | rédigé |
| [`03_methodes.md`](03_methodes.md) | Traitement, agrégation, protocole statistique | rédigé |
| [`04_resultats_H1_algorithme.md`](04_resultats_H1_algorithme.md) | H1 — l'échec est-il algorithmique ? | rédigé |
| [`05_resultats_H2_cible.md`](05_resultats_H2_cible.md) | H2 — le tapis est-il une cible distincte ? | rédigé |
| [`06_resultats_H3_mecanisme.md`](06_resultats_H3_mecanisme.md) | H3 — mouvement ou diélectrique ? | rédigé |
| [`07_resultats_H4_hydrologie.md`](07_resultats_H4_hydrologie.md) | H4 — que mesure-t-on alors ? | rédigé |
| [`08_discussion.md`](08_discussion.md) | Comparaison littérature, portée, limites | rédigé |
| [`09_conclusion.md`](09_conclusion.md) | Conclusion + perspectives | rédigé |
| [`figures.md`](figures.md) | Inventaire des figures (source, légende, statut) | à produire |
| [`traceabilite.md`](traceabilite.md) | Phase → résultat → section → figure | rédigé |

## Convention de rédaction

- **Chaque chiffre cité doit être traçable** à un notebook via
  [`traceabilite.md`](traceabilite.md). Aucun nombre « de mémoire ».
- Étiquetage hérité de `../deduction_phaseD.md` : **[FAIT]** mesuré ·
  **[INT]** interprétation étayée · **[HYP]** hypothèse · **[SPEC]** spéculation.
  Ces étiquettes disparaissent à la mise au propre finale mais structurent la
  rédaction.
- **Les résultats négatifs et les prédictions falsifiées sont conservés**
  (§ « erreurs corrigées » de la discussion). C'est un atout en évaluation, pas
  une faiblesse.
- **Toujours distinguer AMPLITUDE (mm) et VITESSE (mm/an)**, et préciser
  **LOS vs vertical**. Notre borne de 2 mm est une *amplitude saisonnière en
  LOS* ; la littérature sur tourbières drainées rapporte des *vitesses en
  vertical* (cm/an). Confondre les deux est l'erreur la plus facile a commettre.
- Les p-values empiriques ont un **plancher** 1/(1+n_tirages) : toujours écrire
  « p ≤ x » quand la valeur est au plancher.

## Mise à jour après un nouveau run

1. Reporter les chiffres dans le fichier de section concerné.
2. Mettre à jour [`traceabilite.md`](traceabilite.md) (phase → chiffre).
3. Si une figure change, mettre à jour [`figures.md`](figures.md).
4. Répercuter dans `../deduction_phaseD.md` (document de travail détaillé).

## Message principal en une phrase

> Sur une tourbière flottante, l'échec de l'InSAR Sentinel-1 n'est **pas**
> algorithmique mais **physique** ; le signal mesurable n'est pas un mouvement
> (< 2 mm) mais une **sensibilité diélectrique à l'humidité de surface**.
