
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rEAG <a href='https://github.com/JHuguenin/rEAG'><img src="https://raw.githubusercontent.com/JHuguenin/rEAG/master/inst/img/logo_rEAG_1.png" align="right" height="138"/></a>

<!-- badges: start -->
<!-- badges: end -->

rEAG is a R package for import and analyse EAG data.

## Installation

You can install the development version of rEAD from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JHuguenin/rEAG")
library(rEAG)
```

## Utilisation

Le package rEAG est conçu pour importer et traiter les donnees EAG,
issues du setup du CEFE (UMR 5175) de Montpellier composé d’un appareil
Syntech,

### Plan d’experience

Le plan d’experience doit etre concu sous cette forme :

<img src="https://raw.githubusercontent.com/JHuguenin/rEAG/master/inst/img/expdes_table.png" align="center"/>

Le fichier peut etre genere vide grace a la commande suivante :

``` r
# create experimental design
create.empty.exp.design(VOC = c("L","B","N"), control = "T", nS = 10)
```

L’argument *VOC* est la sequence de molecules testes. Cette sequence est
encadree par un *control*. *nS* represente le nombre d’echantillons.
Chaque molecule sequence *VOC* doit etre separee par un espace mais peut
etre composee d’une ou plusieurs lettres. A partir de ce fichier vide,
des colonnes peuvent etre ajoutee librement a droite, selon vos besoins.

Une fois cree, le plan d’experience doit etre importe :

``` r
# import experimental design
expdes <- import.exp.design("Experimental_design")
```

### Import des EAG

Un fichier EAG doit exister pour chaque ligne du plan d’experience. Ce
fichier doit avoir le meme nom que celui de la colonne *File*
correspondante. Les donnees d’un fichier EAG sont disposees en colonnes.
Les deux premieres sont les moyennes du reste du fichier. Ensuite,
chaque groupe de trois colonne correspond a une impulsion de l’EAG. Et
les impulsions de VOC suivent evidemment la sequence precisee dans le
plan d’experience.

Pour chaque implusion, il y a trois colonnes (EAG, FID, DIG suivi d’un
numero) :  
- EAG : l’electroantennograme, en mV.  
- FID : inutile dans cette analyse.  
- DIG : le marqueur d’impulsion (+/- 0.5, sans unite).

Importez chaque fichier :

``` r
# import data
S_A1 <- eag.import(Sname = "Sample_A1")
S_A2 <- eag.import(Sname = "Sample_A2")
S_A3 <- eag.import(Sname = "Sample_A3")
S_A4 <- eag.import(Sname = "Sample_A4")
S_B1 <- eag.import(Sname = "Sample_B1")
S_B2 <- eag.import(Sname = "Sample_B2")
S_B3 <- eag.import(Sname = "Sample_B3")
S_B4 <- eag.import(Sname = "Sample_B4")
```

A noter, quatre variables peuvent etre precisee :  
- wd (= NULL) : le repertoire de travail. Si *NULL*, wd prend la valeur
de *getwd()*.  
- control (= “T”) : le nom du control. C’est le VOC qui permet de
calculer la ligne de base.  
- tmP (= 50) : la duree de l’impulsion. Depend du setup.  
- tmD (= 2\*tmP) : la duree de la depolarisation. De base, 2xPulse mais
peut augmenter si les depolarisation sont trop fortes. Dans l’ideal, il
doit etre equivalent pour tout les echantillons. - ws (= 25) : la
largeur de la fenetre pour le lissage.

Vous obtenez un objet S4 “eag” avec les electrogrammes dans la matrice
@eag et les donnees calculees dans la matrice @depol :  
- Idp : Intensite de la depolarisation - Tdp : Temps de depolarisation -
con : Concentration du VOC - seq : Nom ou abreviation du VOC - Idp\_norm
: Intensite normee - variable : nom racourci visible sur le graphe -
Idp\_adj : Intensite absolue - Idp\_norm\_adj : Intensite normee absolut

### Regroupement des echantillons

Une fois vos echantillons correctement importes, rassemblez les dans un
objet meag, pour “multiple EAG”.

``` r
## rassembler les echantillons #
Sname <- c("S_A1", "S_A2", "S_A3", "S_A4",
           "S_B1", "S_B2", "S_B3", "S_B4")
Meag <- eag.merge(S_A1, S_A2, S_A3, S_A4,
                  S_B1, S_B2, S_B3, S_B4,
                  eag_names = Sname)
View(Meag@depol) # to view results
```

Cette fonction ne fait que rassembler les resultats. Le tableau
*<Meag@depol>* peut etre utilise directement pour vos analyses.

### Graphiques

Vous pouvez generer des graphes supplementaires grace a la fonction
*eag.print()*. Vous pouvez grouper vos echantillons grace a l’argument
*moda* qui peut être soit “con”, soit “seq” soit un autre titre de votre
tableau du plan d’experience. Si vous souhaitez supprimer des VOC ou des
concentrations, utilisez *delet\_seq* ou *delet\_con*.

``` r
eag.print(Meag, moda = "seq", delet_seq = "T")
eag.print(Meag, moda = "seq", delet_con = c("0","0.1","1")
```

Les bonnes idees, les remonter de bugs, les demandes de developpement,
d’aide ou de colaboration peuvent etre envoyees a l’adresse renseignee
dans mon profil.
