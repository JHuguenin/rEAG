
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
```

## Utilisation

Les deux premieres etapes, “Import” et “Regroupement” d’echantillons,
doivent etre effectue selon la metodologie suivante.

### Import

Le package rEAD est conçu pour importer et traiter les donnees EAD du
setup du CEFE (UMR 5175) de Montpellier. Ce setup est compose d’un EAG
de la marque Syntech.  
Les donnees de l’EAG sont disposees en colonnes structures de groupe de
vecteurs subcompose le spectre FID et l’autre l’electroantennogramme
EAD. Ces vecteurs ne sont pas indexes sur le temps. Les deux spectres
FID respectifs permettent d’effectuer l’alignement et l’etalonnage.
Trois difficultes supplementaire s’ajoutent a l’analyse :

<img src="https://raw.githubusercontent.com/JHuguenin/rEAG/master/inst/img/expdes_table.PNG" align="center"/>

et des trucs

``` r
library(rEAD)

## gestion de repertoir de travail # 
wd = "C:/Users/huguenin/Documents/R/rEAD_moustiques" # working directory
setwd(wd)  # changement du repertoire dans R

## 1er import #
a01 <- import.GC.EAD(file_csv = "data/sample_01_EAD.csv")
```

### Regroupement des echantillons

Une fois vos echantillons correctement importes, rassemblez les dans un
objet mead, pour “multiple EAD”.

``` r
## rassembler les echantillons #
list_test <- gcead.merge(a01, a02, a03, a04, a05, # tous les echantillons importes
                         gcead_names = c("a01","a02","a03","a04","a05"), # leurs noms reduits
                         RT_limits = c(6,18), # le temps tronques des parties inutiles
                         gap_FID = 1) # un facteur multiplicatif pour le signal FID
```

Essentiellement, la fonction *gcead.merge* utilise les signaux FID de
chaque echantillon pour aligner toutes les donnes. L’argument
*gcead\_names* permet de simplifier le nom des echantillons. Utilisez
des noms courts (et sans regex, ou alors ne comptez pas sur moi pour
venir vous aider). Le *RT\_limits* permet de zoomer directement sur les
zones d’analyses. Le *gap\_FID* applique un facteur multiplicatif au
signal FID afin d’augmenter la lisibilite des figures. Plusieus
parametres detailles dans l’aide permettre d’affinier les pretraitements
effectues sur les donnees FID. Les donnees peuvent etres visualisees
grace a la fonction *mead.graph()*.

Enfin, la fonction *m.EAD.norm* permet de mettre en forme les signaux
EAD pour permettre l’analyse.

``` r
list_test <- m.EAD.norm(shift = -0.5, amplitude = 5, overlay = FALSE, pk_mat = list_test) # une mise en forme des signaux
```

### Analyses

``` r
## measure #
depol_param <- pk.measure(list_test)  # calcul des intensites et de la duree de depolarisation des pics EAD
# View(depol_param[[1]]) pour la moyenne
# write.csv2(depol_param[[1]], "depol_param.csv") # pour sauvegarder le fichier csv
```

Les bonnes idees, les remonter de bugs, les demandes de developpement,
d’aide ou de colaboration peuvent etre envoyer a l’adresse renseignee
dans mon profil.
