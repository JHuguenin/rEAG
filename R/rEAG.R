#' rEAG : read EAG data
#'
#' import, pretraitement and tools for analyze data of EAG
#'
#' @docType package
#' @name rEAG
#'
#' @import stringr
#' @import utils
#' @import stats
#' @import methods
#' @import magrittr
#' @importFrom htmlwidgets saveWidget
#' @importFrom data.table melt
#' @importFrom data.table setDT
#' @importFrom plotly add_markers
#' @importFrom plotly plot_ly
#' @importFrom pracma savgol
#' @importFrom RColorBrewer brewer.pal
NULL

# Library ####

# library(plotly)
# library(stringr)
# library(pracma)
# library(RColorBrewer)
# library(data.table)

# Import Experimental Design ####

#' Import experimental design
#'
#' Import un plan d'experience, les concentrations des VOC, les sequences,
#' modalites et commentaires
#'
#' @param ExpDes file name (without extension)
#' @param wd working directory
#'
#' @return a data frame
#' @export
#'
#' @examples
#' # expdes <- import.exp.design("Experimental_design", wd = wod)
import.exp.design <- function(ExpDes = "Experimental_Design_Name", wd = NULL){

  # check ####
  if (is.null(wd) == TRUE) wd <- getwd()
  if (!is.character(wd)) stop("'wd' must be character")
  if (!is.character(ExpDes)) stop("'ExpDes' must be character")
  if (length(ExpDes) != 1) stop("Length of 'ExpDes' must be 1.")

  # import
  ExpDes <- str_remove(ExpDes,".csv") # suppresion du l'extension
  expdes <- read.csv2(paste0(wd,"/",ExpDes,".csv")) # plan d'experience

  # verification des colonnes
  if (length(colnames(expdes)) < 4) stop("Experimental Design csv file must respect the format.
                                         Use 'create.empty.exp.design()' for a good practise")
  fmr <- colnames(expdes)[1:4] %>% str_flatten()
  if (fmr != "DateVOCconcentrationVOCseqFile" ) stop("Experimental Design csv file must respect the format.
  First columns must be 'Date' 'VOCconcentration' 'VOCseq' and 'File'. Use 'create.empty.exp.design()' for a good practise")

  # verification des formats
  expdes$VOCconcentration <- as.numeric(expdes$VOCconcentration) # numeric
  expdes$File <- str_remove_all(expdes$File,".csv") # suppresion du l'extension

  return(expdes)
}

#' Create a empty experimental design
#'
#' Genere un fichier .csv avec les sequences de VOC pretires. De colonnes et des
#' lignes peuvent etre ajoutees a ce fichier si necessaire. Les nouvelles colonnes
#' doivent toutes etre ajoutees sur la droite.
#'
#' @param VOC VOC names
#' @param control control name
#' @param nS number of pre-drawn sequences
#' @param wd working directory
#'
#' @return a csv file
#' @export
#'
#' @examples
#' create.empty.exp.design(VOC = c("L","B","N"), control = "T", nS = 10)
create.empty.exp.design <- function(VOC = c("L","B","N"), control = "T", nS = 10, wd = NULL){

  Vs <- NULL
  for(i in 1:nS) Vs <- c(Vs,str_flatten(c(control,sample(VOC),control)," "))

  fmr <- data.frame(File = "",
             VOCconcentration = "",
             VOCseq =  Vs)
  write.csv2(fmr,"Exp_Design_Empty.csv",quote = FALSE,row.names = FALSE)
}

# Import EAG ####

#' Import EAG
#'
#' Importe un fichier csv, detecte les pulses de COV puis centre chaque pulse sur un T0.
#' On calcule l'intensite et le temps de depolarisation pour chaque pulse.
#'
#' @param Sname name of file.csv
#' @param control name of control (concentration zero)
#' @param tmP duration of pulse (in centi second)
#' @param wd working directory
#' @param ws width selection
#'
#' @return a S4 eag object
#' @export
#'
#' @examples
#' # B_137 <- eag.import("Bombus137_C_0ppb_1h", control = "T", tmP = 50, wd = NULL)
eag.import <- function(Sname = filenames[1], control = "T", tmP = 50, wd = NULL, ws =25){

  # check ####
  if (is.null(wd) == TRUE) wd <- getwd()
  if (!is.character(wd)) stop("'wd' must be character")
  if (("figures" %in% dir(wd))==FALSE) dir.create(paste0(wd,"/figures"))
  if (!is.character(Sname)) stop("'Sname' must be character")
  if (!is.numeric(ws)) stop("'ws' must be numeric")
  if (!is.numeric(tmP)) stop("'tmP' must be numeric")
  if (!is.character(control)) stop("'control' must be a character")
  if (length(ws) != 1) stop("Length of 'ws' must be 1.")
  if (length(tmP) != 1) stop("Length of 'tmP' must be 1.")
  if (length(control) != 1) stop("Length of 'control' must be 1.")

  # Importation du fichier EAG ####
  eag <- read.csv(paste0(wd,"/",Sname,".csv"))[,-(1:2)] # import du fichier EAG
  neag <- ncol(eag)/3  # nombre de signaux EAG

  fmr <- which(Sname == expdes$File) # index dans le plan d'experience
  Cvoc <- expdes$VOCconcentration[fmr] # Concentrations associees
  Svoc <- expdes$VOCseq[fmr] # sequences VOC associees

  # Traitement de chaque signal ####
  resM <- NULL # matrice des resultats
  for(k in 1:neag){ # k=1

    # index effectif du fichier signal
    i <- seq(length.out = neag,by = 3)[k]

    # La sequence VOC
    VOCseq <- str_squish(Svoc[k]) %>% str_split(" ") %>% unlist()

    # recherche des pulses
    pulse <- eag[,i+2] - c(eag[-1,i+2],-0.5)
    dP <- data.frame(str = which(pulse == -1),
                     end = which(pulse == -1)+tmP) # position des pulses start et end
    if(length(VOCseq) != nrow(dP) ) stop("The 'VOC sequence' of experimental design does not correspond to the number of pulses")

    # name pulse
    rownames(dP) <- paste0(VOCseq, ave(VOCseq, VOCseq, FUN = seq), "_c",Cvoc[k])

    # calcules preliminaires
    dP$X1 <- dP$str - 2*tmP # un blanc anterieur au pulse d'une duree deux fois plus longue que le pulse
    dP$X2 <- dP$end + 2*tmP # idem pour le blanc ulterieur
    dP$Idp <- rep(NA,nrow(dP)) # NA pour l'intensite de la depolarisation
    dP$Tdp <- rep(NA,nrow(dP)) # NA pour la duree de depolaristation

    matP <- matrix(NA,nrow(dP),length(dP$X1[1]:dP$X2[1])) # une matrix pour l'EAG brut du pulse elargie
    matS <- matP # la meme pour le signal EAG lissee

    # chaque peuf
    for(j in 1:nrow(dP)) { #j=1 # pour chaque pulse de la sequence

      # preparation des matrix signals
      IndP <- dP$X1[j]:dP$X2[j] # les index du pulse elargie
      matS[j,] <- savgol(eag[IndP,i],ws) # lissage
      fmr <- matS[j,2*tmP+1] # le shift a T0
      matS[j,] <- matS[j,] - fmr # suppression du shift
      matP[j,] <- eag[IndP,i] - fmr # idem pour l'EAG non lissee

      # caracterisation de l'impuslion
      Ips <- (2*tmP):(2*tmP + 2*tmP) # de T0 ?? T0+2x la duree du pulse du VOC
      dP$Idp[j] <- min(matS[j,Ips]) # l'intensite max de la depol
      dP$Tdp[j] <- (which.min(matS[j,Ips])-2)/100 # sa duree
    }

    # Nettoyage des donnees
    matS[,1:((ws-1)/2)] <- NA # suppression des vagues induites par le lissage
    matS[,(ncol(matS)-((ws-1)/2)):ncol(matS)] <- NA # idem pour la fin du signal

    Tp <- dP$X1[j]:dP$X2[j] - dP$str[j]# le temps recentre

    # # Graphe control
    # matplot(Tp/100,t(matS), type = "l", lwd = 2, lty = 1, col = rainbow(nrow(dP)),
    #         xlab = "Time (s)", ylab = "EAG (mA)", main = paste(Sname,"\nC =",Cvoc[k]))
    # matplot(Tp/100,t(matP), type = "l", lwd = 1, lty = 1, col = rainbow(nrow(dP)), add = TRUE)
    # abline(v = c(0,tmP/100), col = "orange", lwd = 2, lty = 2)
    # points(dP$Tdp, dP$Idp, col = "black", pch = 16)
    # legend("bottomleft", legend = rownames(dP), lty = 1, lwd = 2,col = rainbow(nrow(dP)))

    # Union des matrix
    dP$con <- as.numeric(Cvoc[k])
    dP$seq <- VOCseq

    # mise en forme
    colnames(matS) <- paste0("T_",Tp) # le temps
    fmr <- !is.na(matS[1,]) # pour le calcul du temps lors de la partie eag ci-dessus ... ouais c'est pas ouf
    matS <- matS[,!is.na(matS[1,])] # suppression des vagues du au lissage
    resM <- rbind(resM,cbind(dP[,5:8], matS)) # maj de la matrice
  }

  # Mise en forme finale ####

  depol <- resM[,1:4] # on garde l'intensit?? et le temps de depol, ainsi que la concentration et le VOC
  depol$con[which(depol$seq == control)] <- 0 # mise a zero des temoins
  depol$con <- as.character(depol$con) %>% as.factor() # les numeriques posent preobleme pour le graphe
  depol$seq <- as.factor(depol$seq) # facteur
  depol$variable <- rownames(depol) # pour le graphe

  eag <- data.frame(time = Tp[fmr]/100, t(resM[,-(1:4)])) # data
  setDT(eag) # transformation en data.table (???)
  Teag <- melt(eag, id.vars = "time") # mise en forme pour plotly
  Teag$seq <- Teag$variable # ajout de la sequence des VOC
  levels(Teag$seq) <- depol$seq
  Teag$con <- Teag$variable # ajout des concentrations
  levels(Teag$con) <- depol$con

  # Graphes ####

  ### by Concentration
  fig <- plot_ly(Teag, type = "scatter", mode = "lines", x = ~time, y = ~value,
                 color = ~con, name = ~variable)
  fig <- plotly::layout(fig, title = Sname,
                shapes = list(list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                                   x0 = 0, x1 = 0, line = list(color = "orange", dash="dot")),
                              list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                                   x0 = tmP/100, x1 = tmP/100, line = list(color = "orange", dash="dot"))),
                xaxis = list(title = 'Time (sec)'),
                yaxis = list(title = 'EAG (mA)'),
                legend = list(title=list(text='<b> VOC_concentration </b>'), x = 0.02, y = 0.9))
  fig <- add_markers(fig, x = ~Tdp, y = ~Idp, data = depol, showlegend= FALSE)

  # print(fig)
  htmlwidgets::saveWidget(fig, paste0(wd,"/figures/options_graph_dy.html"), selfcontained = TRUE)
  file.rename(from = paste0(wd,"/figures/options_graph_dy.html"),
              to = paste0(wd,"/figures/",Sname,"_by_concentration.html"))

  ### by VOC
  fig <- plot_ly(Teag, type = "scatter", mode = "lines", x = ~time, y = ~value,
                 color = ~seq, name = ~variable)
  fig <- plotly::layout(fig, title = Sname,
           shapes = list(list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                              x0 = 0, x1 = 0, line = list(color = "orange", dash="dot")),
                         list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                              x0 = tmP/100, x1 = tmP/100, line = list(color = "orange", dash="dot"))),
           xaxis = list(title = 'Time (sec)'),
           yaxis = list(title = 'EAG (mA)'),
           legend = list(title=list(text='<b> VOC_concentration </b>'),
                         x = 0.02, y = 0.9))
  fig <- add_markers(fig, x = ~Tdp, y = ~Idp, data = depol, showlegend= FALSE)

  # print(fig)
  htmlwidgets::saveWidget(fig, paste0(wd,"/figures/options_graph_dy.html"), selfcontained = TRUE)
  file.rename(from = paste0(wd,"/figures/options_graph_dy.html"),
              to = paste0(wd,"/figures/",Sname,"_by_VOC.html"))

  # Export ####

  return(new(Class = "eag", depol = depol, eag = eag, tmP = tmP, wd = wd, names = Sname))
}

# Merge EAG ####

#' Merge EAG
#'
#' @param ... list of eag
#' @param eag_names reduce name of samples
#' @param print.graph logical. Print or not a graph for view depolarisation, group by seq and by concentration
#' @param tmP duration of pulse (in cs)
#' @param wd working directory
#'
#' @return a S4 eag object
#' @export
#'
#' @examples
#' # B_137 <- eag.import("Bombus137_C_0ppb_1h")
#' # B_170 <- eag.import("Bombus170_C_0ppb_1h")
#' # B_171 <- eag.import("Bombus171_C_0ppb_1h")
#' # B_172 <- eag.import("Bombus172_C_0ppb_1h")
#' #
#' # Sname <- c("B_137","B_170","B_171","B_172")
#' # Meag <- eag.merge(B_137, B_170, B_171, B_172, eag_names = Sname)
eag.merge <- function(..., eag_names = NULL, tmP = NULL, wd = NULL, print.graph = FALSE){

  ls_eag <- list(...) # ls_eag <- list(B_137, B_170, B_171, B_172)

  # check
  sapply(ls_eag, function(X) if (class(X) != "eag") stop("variables must be a eag S4 object"))
  if (is.null(eag_names)) eag_names <- sapply(ls_eag, function(X) return(X@names))
  if (!is.character(eag_names))  stop("eag_names must be a character")
  if (length(eag_names) != length(ls_eag)) stop("Length of eag_names must be egal to length of eag objects")
  if (is.null(wd) == TRUE){
    wd <- lapply(ls_eag, function(X) return(X@wd)) %>% unlist() %>% unique()
    if (length(wd) != 1) print(paste("Several 'wd' have been detected. The current 'wd' is :",wd[1]))
    wd <- wd[1]
  }
  if (!is.character(wd)) stop("'wd' must be character")
  if (("figures" %in% dir(wd))==FALSE) dir.create(paste0(wd,"/figures"))
  if (is.null(tmP) == TRUE){
    tmP <- lapply(ls_eag, function(X) return(X@tmP)) %>% unlist() %>% unique()
    if (length(tmP) != 1) print(paste("Several 'tmP' have been detected. The current 'tmP' is :",tmP[1]))
    tmP <- tmP[1]
  }
  if (!is.numeric(tmP)) stop("'tmP' must be numeric")
  if (length(tmP) != 1) stop("Length of 'tmP' must be 1.")
  if (!is.logical(print.graph)) stop("'print.graph must be a logical")

  # mise en forme
  leag <- length(ls_eag)
  for(i in 1:leag){ # i=1
    rownames(ls_eag[[i]]@depol) <- paste0(eag_names[i],"_",rownames(ls_eag[[i]]@depol))
    ls_eag[[i]]@depol$variable <- rownames(ls_eag[[i]]@depol)
    colnames(ls_eag[[i]]@eag) <- c("time",rownames(ls_eag[[i]]@depol))
  }

  # Merge ####
  ls_depol <- lapply(ls_eag, function(X) return(X@depol))
  depol <- do.call(rbind,ls_depol)

  list_eag <- lapply(ls_eag, function(X) return(X@eag[,-1]))
  eag <- do.call(cbind,list_eag) %>% cbind("time" = ls_eag[[1]]@eag$time,.)

  # Mise en forme ####
  Teag <- melt(eag, id.vars = "time") # mise en forme pour plotly
  Teag$seq <- Teag$variable # ajout de la sequence des VOC
  levels(Teag$seq) <- depol$seq
  Teag$con <- Teag$variable # ajout des concentrations
  levels(Teag$con) <- depol$con

  # Graphes ####
  if(print.graph == TRUE){
    ### by Concentration
    fig <- plot_ly(Teag, type = "scatter", mode = "lines", x = ~time, y = ~value,
                   color = ~con, name = ~variable)
    fig <- plotly::layout(fig, title = str_flatten(string = eag_names,collapse = " "),
                  shapes = list(list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                                     x0 = 0, x1 = 0, line = list(color = "orange", dash="dot")),
                                list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                                     x0 = tmP/100, x1 = tmP/100, line = list(color = "orange", dash="dot"))),
                  xaxis = list(title = 'Time (sec)'),
                  yaxis = list(title = 'EAG (mA)'),
                  legend = list(title=list(text='<b> VOC_concentration </b>'),
                                x = 0.02, y = 0.9))
    fig <- add_markers(fig, x = ~Tdp, y = ~Idp, data = depol, showlegend= FALSE)

    # print(fig)
    htmlwidgets::saveWidget(fig, paste0(wd,"/figures/options_graph_dy.html"), selfcontained = TRUE)
    file.rename(from = paste0(wd,"/figures/options_graph_dy.html"),
                to = paste0(wd,"/figures/all_EAG_by_concentration.html"))

    ### by VOC
    fig <- plot_ly(Teag, type = "scatter", mode = "lines", x = ~time, y = ~value,
                   color = ~seq, name = ~variable)
    fig <- plotly::layout(fig, title = str_flatten(eag_names," "),
             shapes = list(list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                                x0 = 0, x1 = 0, line = list(color = "orange", dash="dot")),
                           list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                                x0 = tmP/100, x1 = tmP/100, line = list(color = "orange", dash="dot"))),
             xaxis = list(title = 'Time (sec)'),
             yaxis = list(title = 'EAG (mA)'),
             legend = list(title=list(text='<b> VOC_concentration </b>'),
                           x = 0.02, y = 0.9))
    fig <- add_markers(fig, x = ~Tdp, y = ~Idp, data = depol, showlegend= FALSE)

    # print(fig)
    htmlwidgets::saveWidget(fig, paste0(wd,"/figures/options_graph_dy.html"), selfcontained = TRUE)
    file.rename(from = paste0(wd,"/figures/options_graph_dy.html"),
                to = paste0(wd,"/figures/all_EAG_by_VOC.html"))
  }

  # export
  return(new(Class = "eag", depol = depol, eag = eag, tmP = tmP, wd = wd, names = eag_names))
}

# Print EAG ####

#' Print EAG
#'
#' a simple function for print results of EAG
#'
#' @param eag a eag objet
#' @param moda a modality for group samples
#'
#' @return a figure
#' @export
#'
#' @examples
#' # Meagdepol$con_seq <- paste0(Meag@depol$seq,"_",Meag@depol$con)
#' # eag.print(Meag, "con_seq")
eag.print <- function(eag, moda = "seq"){

  # check
  if (class(eag) != "eag") stop("eag must be a eag S4 object")
  if (!is.character(moda)) stop("'moda' must be character")
  if (("figures" %in% dir(wd))==FALSE) dir.create(paste0(wd,"/figures"))
  if (length(moda) != 1) stop("Length of 'moda' must be 1.")
  im <- match(moda,colnames(eag@depol))
  if (is.na(im)) stop("moda must be a vector of eag@depol")

  # Mise en forme ####
  Teag <- melt(eag@eag, id.vars = "time") # mise en forme pour plotly
  Teag$moda <- Teag$variable # ajout de la sequence des VOC
  levels(Teag$moda) <- eag@depol[,im]
  eag@depol$moda <- eag@depol[,im]

  # Graphe ####
  nM <- length(levels(Teag$moda))
  dcol <- c(brewer.pal(8,"Accent"),brewer.pal(8,"Dark2"),brewer.pal(8,"Set2"))[1:nM]

  ### by Modality
  fig <- plot_ly(Teag, type = "scatter", mode = "lines", x = ~time, y = ~value,
                 color = ~moda, name = ~variable, colors = dcol)
  fig <- plotly::layout(fig, title = paste(str_flatten(eag@names, " "),"by",moda),
           shapes = list(list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                              x0 = 0, x1 = 0, line = list(color = "orange", dash="dot")),
                         list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                              x0 = tmP/100, x1 = tmP/100, line = list(color = "orange", dash="dot"))),
           xaxis = list(title = 'Time (sec)'),
           yaxis = list(title = 'EAG (mA)'),
           legend = list(title=list(text='<b> VOC_concentration </b>'), x = 0.02, y = 0.9))
  fig <- add_markers(fig, x = ~Tdp, y = ~Idp, data = eag@depol, showlegend = FALSE)

  print(fig)

  htmlwidgets::saveWidget(fig, paste0(wd,"/figures/options_graph_dy.html"), selfcontained = TRUE)
  file.rename(from = paste0(wd,"/figures/options_graph_dy.html"),
              to = paste0(wd,"/figures/eag_by_",moda,".html"))
}

# Gestion des 'S4' ####

#' object used for rEAG package
#'
#' @slot depol data.frame.
#' @slot eag data.table.
#' @slot tmP numeric.
#' @slot wd character.
#' @slot names character.
#'
#' @return a eag object
#' @export
#'
#' @examples
#' #class(eag)
setClass("eag",representation(depol = "data.frame",
                              eag = "data.table",
                              tmP = "numeric",
                              names = "character",
                              wd = "character"))
