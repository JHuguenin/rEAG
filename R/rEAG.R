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

# Experimental Design ####

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
  fmr <- colnames(expdes)[1:3] %>% str_flatten()
  if (fmr != "FileNumEAGVOCconcentration" ) stop("Experimental Design csv file must respect the format.
  First columns must be 'File','NumEAG','VOCconcentration' then 'VOCseq'. Use 'create.empty.exp.design()' for a good practise")

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

  fmr <- data.frame(File = "", NumEAG = rep(1,nS), VOCconcentration = "", VOCseq =  Vs)
  write.csv2(fmr,"Exp_Design_Empty.csv",quote = FALSE,row.names = FALSE)
}

# Import EAG ####

#' Import EAG
#'
#' Importe un fichier csv, detecte les pulses de COV puis centre chaque pulse sur un T0.
#' On calcule l'intensite et le temps de depolarisation pour chaque pulse.
#'
#' @param Sname name of file.csv
#' @param expdes name od data.fram Experimental Design
#' @param control name of control (concentration zero)
#' @param tmP duration of pulse (in centi second)
#' @param tmD duration of depolarisation (in centi second)
#' @param wd working directory
#' @param ws width selection
#'
#' @return a S4 eag object
#' @export
#'
#' @examples
#' # B_137 <- eag.import("Bombus137_C_0ppb_1h", control = "T", tmP = 50, wd = NULL)
eag.import <- function(Sname, expdes = ExDs, control = "T", tmP = 50, tmD = NULL, wd = NULL, ws =25){

  # check ####
  if (is.null(wd) == TRUE) wd <- getwd()
  if (is.null(tmD) == TRUE) tmD <- 2*tmP
  if (!is.character(wd)) stop("'wd' must be character")
  if (("figures" %in% dir(wd))==FALSE) dir.create(paste0(wd,"/figures"))
  if (!is.character(Sname)) stop("'Sname' must be character")
  if (!is.numeric(ws)) stop("'ws' must be numeric")
  if (!is.numeric(tmP)) stop("'tmP' must be numeric")
  if (!is.numeric(tmD)) stop("'tmD' must be numeric")
  if (!is.character(control)) stop("'control' must be a character")
  if (length(ws) != 1) stop("Length of 'ws' must be 1.")
  if (length(tmP) != 1) stop("Length of 'tmP' must be 1.")
  if (length(tmD) != 1) stop("Length of 'tmD' must be 1.")
  if (length(control) != 1) stop("Length of 'control' must be 1.")
  if (!is.data.frame(expdes)) stop("Experimental Design must be a data frame")

  # Importation du fichier EAG ####
  eag <- read.csv(paste0(wd,"/",Sname,".csv")) # import du fichier EAG
  fmr <- grep("AVE",colnames(eag)) # detection des colonnes AVE
  if(length(fmr)>0)  eag <- eag[, -fmr] # suppresion des colonnes AVE
  neag <- ncol(eag)/3  # nombre de signaux EAG

  fmr <- str_split(Sname,"/",simplify = TRUE)
  Sname <- as.character(fmr[length(fmr)])

  fmr <- which(Sname == expdes$File) # index dans le plan d'experience
  Cvoc <- expdes$VOCconcentration[fmr] # Concentrations associees
  Svoc <- expdes$VOCseq[fmr] # sequences VOC associees
  NumEAG <- expdes$NumEAG[fmr] # ID of associated sequence
  if(neag != length(NumEAG)) stop(paste("the number of lines named",
    Sname,"in the experimental design does not correspond to the number of EAGs
    in the",Sname,"file."))

  # Figure des EAG brutes ####
  dt_eag <- eag  # data
  dt_eag$time <- (1:nrow(eag))/100
  dt_eag <- dt_eag[,-grep("FID",colnames(dt_eag))]
  setDT(dt_eag) # transformation en data.table (???)
  Teag <- melt(dt_eag, id.vars = "time") # mise en forme pour plotly
  Teag$num <- Teag$variable # ajout des concentrations
  levels(Teag$num) <- str_replace(levels(Teag$num),"DIG","EAD")
  fig <- plot_ly(Teag, type = "scatter", mode = "lines", x = ~time, y = ~value,
                 color = ~num, name = ~variable)
  fig <- plotly::layout(fig, title = Sname,
                        xaxis = list(title = 'Time (sec)'),
                        yaxis = list(title = 'EAG (mV)'))
  htmlwidgets::saveWidget(fig, paste0(wd,"/figures/options_graph_dy.html"), selfcontained = TRUE)
  file.rename(from = paste0(wd,"/figures/options_graph_dy.html"),
              to = paste0(wd,"/figures/",Sname,"_raw.html"))

  # Traitement de chaque signal ####
  resM <- NULL # matrice des resultats
  for(k in NumEAG){ # k=1

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
    dP$X2 <- dP$end + tmD # idem pour le blanc ulterieur
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
      Ips <- (2*tmP):(2*tmP + tmD) # de T0 à T0+TmD la duree de depolarisation du VOC
      dP$Idp[j] <- min(matS[j,Ips]) # l'intensite max de la depol
      dP$Tdp[j] <- (which.min(matS[j,Ips])-2)/100 # sa duree
    }

    # Nettoyage des donnees
    matS[,1:((ws-1)/2)] <- NA # suppression des vagues induites par le lissage
    matS[,(ncol(matS)-((ws-1)/2)):ncol(matS)] <- NA # idem pour la fin du signal

    Tp <- dP$X1[j]:dP$X2[j] - dP$str[j]# le temps recentre

    # # Graphe control
    # matplot(Tp/100,t(matS), type = "l", lwd = 2, lty = 1, col = rainbow(nrow(dP)),
    #         xlab = "Time (s)", ylab = "EAG (mV)", main = paste(Sname,"\nC =",Cvoc[k]))
    # matplot(Tp/100,t(matP), type = "l", lwd = 1, lty = 1, col = rainbow(nrow(dP)), add = TRUE)
    # abline(v = c(0,tmP/100), col = "orange", lwd = 2, lty = 2)
    # points(dP$Tdp, dP$Idp, col = "black", pch = 16)
    # legend("bottomleft", legend = rownames(dP), lty = 1, lwd = 2,col = rainbow(nrow(dP)))

    # Union des matrix
    dP$con <- as.numeric(Cvoc[k])
    dP$seq <- VOCseq
    ict <- which(dP$seq == control)
    dP$Idp_norm <- dP$Idp - mean(dP$Idp[ict])

    # mise en forme
    colnames(matS) <- paste0("T_",Tp) # le temps
    fmr <- !is.na(matS[1,]) # pour le calcul du temps lors de la partie eag ci-dessus ... ouais c'est pas ouf
    matS <- matS[,!is.na(matS[1,])] # suppression des vagues du au lissage
    resM <- rbind(resM,cbind(dP[,5:9], matS)) # maj de la matrice
  }

  # Mise en forme finale ####
  depol <- resM[,1:5] # on garde l'intensité et le temps de depol, ainsi que la concentration et le VOC
  i_ctrl <- which(depol$seq == control) # l'index des controls

  depol$con[i_ctrl] <- 0 # mise a zero des temoins
  depol$con <- as.character(depol$con) %>% as.factor() # les numeriques posent probleme pour le graphe
  depol$seq <- as.factor(depol$seq) # facteur
  depol$variable <- rownames(depol) # pour le graphe

  depol$Idp_adj <- abs(depol$Idp)
  depol$Idp_adj[i_ctrl] <- 0
  depol$Idp_norm_adj <- abs(depol$Idp_norm)
  depol$Idp_norm_adj[i_ctrl] <- 0

  eag <- data.frame(time = Tp[fmr]/100, t(resM[,-(1:5)])) # data
  setDT(eag) # transformation en data.table
  Teag <- melt(eag, id.vars = "time") # mise en forme pour plotly
  Teag$seq <- Teag$variable # ajout de la sequence des VOC
  levels(Teag$seq) <- depol$seq
  Teag$con <- Teag$variable # ajout des concentrations
  levels(Teag$con) <- depol$con

  # Graphes ####

  ### by Concentration ...
  # plotly for interactive plot
  fig <- plot_ly(Teag, type = "scatter", mode = "lines", x = ~time, y = ~value,
                 color = ~con, name = ~variable)
  fig <- plotly::layout(fig, title = Sname,
                shapes = list(list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                                   x0 = 0, x1 = 0, line = list(color = "orange", dash="dot")),
                              list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                                   x0 = tmP/100, x1 = tmP/100, line = list(color = "orange", dash="dot"))),
                xaxis = list(title = 'Time (sec)'),
                yaxis = list(title = 'EAG (mV)'),
                legend = list(title=list(text='<b> VOC_concentration </b>'), x = 0.02, y = 0.9))
  fig <- add_markers(fig, x = ~Tdp, y = ~Idp, data = depol, showlegend= FALSE)
  htmlwidgets::saveWidget(fig, paste0(wd,"/figures/options_graph_dy.html"), selfcontained = TRUE)
  file.rename(from = paste0(wd,"/figures/options_graph_dy.html"),
              to = paste0(wd,"/figures/",Sname,"_by_concentration.html"))

  # matplot for fixed plot
  depol$col_C <- depol$con
  levels(depol$col_C) <- RColorBrewer::brewer.pal(length(unique(depol$con)),"Dark2")

  tiff(filename = paste0(wd,"/figures/",Sname,"_by_concentration.tiff"), width = 1000, height = 580)
    par(mar = c(5,5,3,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))

      matplot(eag$time,eag[,-1], type = "l", lty = 1, col = as.character(depol$col_C), lwd = 2,
              main = Sname, xlab = "Time (sec)", ylab = "EAG (mV)")
      abline(h = 0)
      abline(v = c(0,tmP/100), lty = 2, col = "orange")
      points(depol$Tdp, depol$Idp, pch = 16, col = as.character(depol$col_C), cex = 1.5)
      legend("bottomleft", legend = rownames(depol),col = as.character(depol$col_C), lty = 1)
  dev.off()

  ### by VOC
  fig <- plot_ly(Teag, type = "scatter", mode = "lines", x = ~time, y = ~value,
                 color = ~seq, name = ~variable)
  fig <- plotly::layout(fig, title = Sname,
           shapes = list(list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                              x0 = 0, x1 = 0, line = list(color = "orange", dash="dot")),
                         list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                              x0 = tmP/100, x1 = tmP/100, line = list(color = "orange", dash="dot"))),
           xaxis = list(title = 'Time (sec)'),
           yaxis = list(title = 'EAG (mV)'),
           legend = list(title=list(text='<b> VOC_concentration </b>'),
                         x = 0.02, y = 0.9))
  fig <- add_markers(fig, x = ~Tdp, y = ~Idp, data = depol, showlegend= FALSE)
  htmlwidgets::saveWidget(fig, paste0(wd,"/figures/options_graph_dy.html"), selfcontained = TRUE)
  file.rename(from = paste0(wd,"/figures/options_graph_dy.html"),
              to = paste0(wd,"/figures/",Sname,"_by_VOC.html"))

  # matplot for fixed plot
  depol$col_S <- depol$seq
  levels(depol$col_S) <- RColorBrewer::brewer.pal(length(unique(depol$seq)),"Accent")

  tiff(filename = paste0(wd,"/figures/",Sname,"_by_sequence.tiff"), width = 1000, height = 580)
    par(mar = c(5,5,3,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))

      matplot(eag$time,eag[,-1], type = "l", lty = 1, col = as.character(depol$col_S), lwd = 2,
              main = Sname, xlab = "Time (sec)", ylab = "EAG (mV)")
      abline(h = 0)
      abline(v = c(0,tmP/100), lty = 2, col = "orange")
      points(depol$Tdp, depol$Idp, pch = 16, col = as.character(depol$col_S), cex = 1.5)
      legend("bottomleft", legend = rownames(depol),col = as.character(depol$col_S), lty = 1)
  dev.off()

  # Export ####
  depol$col_C <- NULL
  depol$col_S <- NULL

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
                  yaxis = list(title = 'EAG (mV)'),
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
             yaxis = list(title = 'EAG (mV)'),
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
#' @param delet_con character. Concentration to be deleted
#' @param delet_seq character. Name of the sequences to be deleted
#'
#' @return a figure
#' @export
#'
#' @examples
#' # Meagdepol$con_seq <- paste0(Meag@depol$seq,"_",Meag@depol$con)
#' # eag.print(Meag, "con_seq")
eag.print <- function(eag, moda = "seq", delet_seq = NULL, delet_con = NULL){

  # check
  if (class(eag) != "eag") stop("eag must be a eag S4 object")
  if (!is.character(moda)) stop("'moda' must be character")
  if (("figures" %in% dir(wd))==FALSE) dir.create(paste0(wd,"/figures"))
  if (length(moda) != 1) stop("Length of 'moda' must be 1.")
  im <- match(moda,colnames(eag@depol))
  if (is.na(im)) stop("moda must be a vector of eag@depol")
  if (!is.null(delet_con)) if (!is.character(delet_con)) stop("'delet_con' must be null or a character")
  if (!is.null(delet_seq)) if (!is.character(delet_seq)) stop("'delet_seq' must be null or a character")
  if (is.character(delet_con)) if(mean(delet_con %in% levels(eag@depol$con)) != 1) stop("'delet_con' isn't a level of experimental design")
  if (is.character(delet_seq)) if(mean(delet_seq %in% levels(eag@depol$seq)) != 1) stop("'delet_seq' isn't a level of experimental design")

  # suppresion des echantillon non desires ###
  i_del <- c(which(match(eag@depol$seq,delet_seq) > 0),
             which(match(eag@depol$con,delet_con) > 0)) %>% unique() %>% sort() # recherche des indices

  eag@depol <- eag@depol[-i_del,] # suppression dans depol
  i_del <- i_del + 1             # le +1 compense la colonne time
  eag@eag <- eag@eag[,-..i_del] # suppression dans eag

  # Mise en forme ####
  Teag <- melt(eag@eag, id.vars = "time") # mise en forme pour plotly
  Teag$moda <- Teag$variable # ajout de la sequence des VOC
  levels(Teag$moda) <- eag@depol[,im]
  eag@depol$moda <- eag@depol[,im]

  ### Graph by modality
  fig <- plot_ly(Teag, type = "scatter", mode = "lines", x = ~time, y = ~value,
                 color = ~moda, name = ~variable)
  fig <- plotly::layout(fig, title = paste(str_flatten(eag@names, " "),"by",moda),
           shapes = list(list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                              x0 = 0, x1 = 0, line = list(color = "orange", dash="dot")),
                         list(type = "line", y0 = 0, y1 = 1, yref = "paper",
                              x0 = eag@tmP/100, x1 = eag@tmP/100, line = list(color = "orange", dash="dot"))),
           xaxis = list(title = 'Time (sec)'),
           yaxis = list(title = 'EAG (mV)'),
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
