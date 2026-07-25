# 10. Falsification — explications alternatives et tests

**Statut :** rédigé (tests 1-3 et 6 à exécuter) · **Sources :** transverses

> *« Pour chaque conclusion majeure : quelle observation prouverait que j'ai
> tort ? »* Chaque hypothèse alternative éliminée renforce l'explication
> retenue ; celles qui ne peuvent pas l'être doivent être déclarées.

---

## 10.1 Le mécanisme proposé, en une chaîne causale

```
        Nappe quasi affleurante  (état hydrologique)
                  │
                  ▼
   Teneur en eau du couvert et de la tourbe de surface
                  │
                  ▼
   Constante diélectrique élevée et VARIABLE (ε)
                  │
        ┌─────────┴─────────┐
        ▼                   ▼
Profondeur de        Diffusion de VOLUME
pénétration           dominante (RVI ↑, σ0 ↑)
variable                     │
        │                    ▼
        │        Centre de phase instable entre passages
        │                    │
        │                    ▼
        │            COHÉRENCE BASSE  (Δ = −0.069)
        │                    │
        │                    ▼
        │      Échec de l'inversion PAR PIXEL (6 méthodes)
        │                    │
        │                    ▼
        │      Agrégation spatiale → bruit ÷ √N_eff
        ▼                    │
   Décalage du centre        ▼
   de phase avec ε  ──►  Signal saisonnier de 3.3 mm
                             │
                             ▼
              Corrélé aux ANOMALIES d'humidité (r ≈ 0.45, lag 12 j)
                             │
                             ▼
                  ⇒ signal DIÉLECTRIQUE, non mécanique
                    (lac oscille pareil ; A−B s'annule)
```

Toutes les observations pointent vers **un seul mécanisme**. La section suivante
teste si un autre mécanisme pourrait produire les mêmes observations.

---

## 10.2 Explications alternatives

| # | Hypothèse alternative | Statut | Test |
|---|---|---|---|
| 1 | **Neige / gel** | ⏳ à exécuter | exclure les paires hivernales |
| 2 | **Atmosphère** | ✅ largement écartée | double différence + nul apparié |
| 3 | **Géométrie / incidence** | ⏳ à quantifier | Δincidence entre A et C |
| 4 | **Phénologie seule** | ✅ écartée | A et C jumeaux phénologiques |
| 5 | **Erreurs de déroulement** | ✅ écartée | |R| sur phase enroulée ; filtre baseline |
| 6 | **Corrélation spatiale (N_eff)** | ⏳ à mesurer | longueur de corrélation empirique |
| 7 | **Couvert / classe erronée** | ✅ écartée | WorldCover + appariement S2 + surface +0.6 % |
| 8 | **Le tapis ET le lac bougent ensemble** | ❌ **non écartée** | nécessite le laser |

### 1. Neige et gel — *test prêt, à exécuter*

**Pourquoi c'est sérieux.** Un manteau neigeux modifie fortement la
rétrodiffusion et la cohérence, affecte différemment une tourbière saturée et
une prairie drainée, et possède **un cycle annuel** — donc une explication
alternative complète du signal saisonnier de 3.3 mm.

**Test.** Recalculer l'amplitude saisonnière en **excluant les paires dont une
date tombe en décembre-février** (`filter_pairs(..., exclude_months=(12,1,2))`).

- Signal **survit** → neige et gel écartés ; le signal est porté par la saison
  végétative.
- Signal **disparaît** → il était hivernal ; l'interprétation diélectrique
  estivale tombe.

*Indice existant* : le test de gel (Phase D-bis) montrait que le tapis gagne
**moins** de cohérence au gel (+0.028) que la prairie (+0.077) — le tapis ne
gèle pas comme un sol stable. Cela n'exclut pas un effet neigeux sur la phase.

### 2. Atmosphère — *écartée par construction*

L'observable est une **double différence** entre deux zones distantes de ~1 km
vues dans **la même paire** : l'écran atmosphérique à cette échelle est commun et
s'annule au premier ordre. Surtout, le **contrôle nul** est construit sur des
zones de sol stable **soumises au même écran** : toute contribution
atmosphérique résiduelle apparaît donc dans la distribution nulle et est
absorbée par la p-value empirique.

⚠️ Réserve : une composante atmosphérique **corrélée à la topographie** ne
s'annulerait pas parfaitement. Le relief étant très faible ici (< quelques m),
l'effet attendu est négligeable — mais non mesuré.

### 3. Géométrie et angle d'incidence — *à quantifier*

A et C appartiennent au **même burst**, distants de ~1 km en portée : la
variation d'incidence sur cette distance est de l'ordre de **0.05-0.1°** sur une
fauchée de 250 km. Un tel écart ne peut pas produire un signal saisonnier.

**Test.** Extraire l'incidence médiane de A et de C depuis la couche `lv_theta`
et rapporter l'écart. Une ligne suffit à clore la question.

### 4. Phénologie seule — *écartée*

A et C sont des **jumeaux phénologiques** : humidité optique médiane −0.513 vs
−0.522, même classe WorldCover, verdure et amplitude saisonnière appariées. Si
la phénologie seule pilotait le signal, la double différence A − C l'annulerait.

⚠️ Réserve : l'appariement porte sur l'humidité **optique** du couvert, pas sur
l'état diélectrique du **sol** — c'est précisément la variable que nous
invoquons.

### 5. Erreurs de déroulement — *écartée*

|R| et le biais de fermeture sont calculés sur la phase **enroulée**, donc
insensibles aux sauts de 2π. Le filtrage des baselines > 60 j (paires annuelles,
±25 mm de dispersion) teste la sensibilité de l'inversion agrégée : le résultat
saisonnier n'en dépend pas qualitativement.

### 6. Corrélation spatiale et N_eff — *à mesurer*

**La critique est fondée.** L'argument « le bruit décroît en 1/√N » suppose des
pixels indépendants. Ils ne le sont pas.

**Réponse en deux temps :**

**(a) Les résultats principaux n'en dépendent pas.** La significativité de
l'amplitude saisonnière et des corrélations repose sur des **nuls empiriques
appariés en taille** — construits sur des zones réelles, avec la corrélation
spatiale réelle. Aucune valeur de N_eff n'entre dans ces p-values. Le facteur
1/√N n'est qu'une **motivation** de la démarche, pas une étape du calcul.

**(b) Là où N_eff intervient (plancher indicatif de |R|), il doit être mesuré.**
`correlation_length()` estime la portée par autocorrélation empirique (seuil
1/e) et `effective_looks(..., field=...)` en déduit N_eff sur les données au
lieu de le supposer. À exécuter et à reporter.

### 7. Couvert mal attribué — *écartée*

Trois vérifications convergentes : classe ESA WorldCover, appariement sur
features Sentinel-2, et surtout **contrôle de surface** (A+B = 90.24 ha contre
89.7 ha documentés, **+0.6 %**) qui valide numériquement le géoréférencement.

### 8. Tapis et lac bougeant ensemble — ❌ **NON écartée**

C'est la **faiblesse principale** de la borne affinée (§6.8, niveau 2). Le lac et
le tapis flottent sur la même nappe : un mouvement **commun** produirait la même
annulation A − B qu'une **absence** de mouvement.

**Ce qui limite la portée du problème** : le plafond robuste du niveau 1
(≤ 4.2 mm, A − C contre sol stable) ne dépend pas du lac et exclut déjà la
flottaison libre.

**Ce qui la lèverait** : le **laser** in situ — mesure directe et absolue du
mouvement du tapis.

---

## 10.3 Ce que le laser et le drone doivent tester

Le relecteur a raison : leur rôle n'est **pas** de valider un déplacement que
nous ne prétendons pas mesurer, mais de **tester le mécanisme**.

| Instrument | Question posée | Issue et lecture |
|---|---|---|
| **Laser** | Le tapis bouge-t-il, et de combien ? | > 5 mm → notre borne est fausse, chercher l'erreur ; < 4 mm → borne confirmée, et l'ambiguïté n°8 est levée |
| **Laser** | Le mouvement est-il en phase avec notre signal de 3.3 mm ? | en phase → composante mécanique réelle ; déphasé ou absent → confirme le diélectrique |
| **Laser + WTD** | Le mouvement suit-il la nappe ? | oui → flottaison partielle (tapis ancré) ; non → tapis contraint |
| **Drone** | Microtopographie buttes/dépressions | teste directement le modèle de la Phase H (σ0 élevé = dépressions humides = échec) |
| **Drone** | Hétérogénéité interne du tapis | teste l'hypothèse d'**unité** de l'agrégation |
| **WTD** | Forçage à structure temporelle ≠ température | lève la limite principale de H4 |

**Le cas le plus intéressant n'est pas la confirmation.** Si le laser montre un
mouvement réel de 20 mm alors que l'InSAR n'en voit que 3.3 mm dont l'essentiel
est diélectrique, cela **quantifie directement l'insensibilité de la bande C** au
mouvement de cette surface — un résultat plus fort que n'importe quelle
validation croisée.

---

## 10.4 Deux questions distinctes

Nos résultats répondent en réalité à **deux questions séparées**, qui méritent
des figures et des discussions distinctes :

| | Question | Réponse | Sections |
|---|---|---|---|
| **Q1** | Sentinel-1 peut-il mesurer le déplacement vertical ? | **Non**, pas de façon fiable ; borne ≤ 4 mm | §4, §5, §6 |
| **Q2** | Sentinel-1 peut-il renseigner l'état hydrologique saisonnier ? | **Possiblement oui**, effet modéré mais mesurable | §6.5, §7 |

Q1 est un résultat de **limite instrumentale** ; Q2 est un résultat de
**capacité**. Les confondre affaiblirait les deux.
