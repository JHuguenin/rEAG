% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rEAG.R
\name{create.empty.exp.design}
\alias{create.empty.exp.design}
\title{Create a empty experimental design}
\usage{
create.empty.exp.design(
  VOC = c("L", "B", "N"),
  control = "T",
  nS = 10,
  wd = NULL
)
}
\arguments{
\item{VOC}{VOC names}

\item{control}{control name}

\item{nS}{number of pre-drawn sequences}

\item{wd}{working directory}
}
\value{
a csv file
}
\description{
Genere un fichier .csv avec les sequences de VOC pretires. De colonnes et des
lignes peuvent etre ajoutees a ce fichier si necessaire. Les nouvelles colonnes
doivent toutes etre ajoutees sur la droite.
}
\examples{
create.empty.exp.design(VOC = c("L","B","N"), control = "T", nS = 10)
}
