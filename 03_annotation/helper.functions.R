readSAMGenomeIndex <- function(index = NULL,
                               filterSize = NULL){
  # read ingroup chromosome information
  chroms <- read.table(index, 
                       as.is = T,
                       header = F)
  chroms <- chroms[,c(1,2)]
  colnames(chroms) <- c("chromosome","length")
  chroms <- chroms[order(chroms$length,decreasing = T),]
  if(is.null(filterSize) == F){
    chroms <- chroms[chroms$length > filterSize,]    
  }
  return(chroms)
}

# make a function to plot the chromosomes of the outgroup
plotFocalChroms <- function(genomeIndex = NULL,
                            minSize = NULL,
                            axisLab = NULL,
                            sp = NULL,
                            alpha = 80,
                            chromNum = NULL){
  
  # read genome index of the outgroup
  outgroup.chroms <- readSAMGenomeIndex(index = genomeIndex)
  # when multiple scaffolds are present lets keep everything greater than 1 Mb
  # these mostly represent chromosomes
  mb <- 1000000
  if(is.null(chromNum)){
    nChroms <- sum(outgroup.chroms$length > minSize)
    if(nChroms >= 13){
      reduceChroms <- T
      cat("Adjustung minumum chromosome size to keep a maximum of 12 chromosomes \n")
      while (reduceChroms) {
        minSize <- minSize + 100000
        nChroms <- sum(outgroup.chroms$length > minSize)
        if(nChroms <= 12){
          reduceChroms <- F
        }
      }
    }
  }else{
    nChroms <- chromNum
  }
  if(outgroup.chroms$length[nChroms] < mb){
    cat("size of the smallest chromosome is less than 1 MB.
        Which might not be a true chromosome. Re adjustung the chromosome numbers
        to filter smaller scaffolds/contigs less than Megabase size")
    reduceChroms <- T
    while (reduceChroms) {
      nChroms <- nChroms -1
      if(outgroup.chroms$length[nChroms] >= mb){
        reduceChroms <- F
      }
    }
  }else{
    cat("Minimum chromosome size set to", minSize, "bp \n",sep = " ")
  }
  
  # this will assign a gray colour to all scaffolds
  outgroup.chroms$col <- "gray40"
  cols <- RColorBrewer::brewer.pal(n = nChroms, name = "Paired")
  # set transparency
  cols <- paste(cols,alpha, sep = "")
  outgroup.chroms$col[1:nChroms] <- cols
  #plot outgroup
  plot(x = NULL,
       y = NULL,
       xlim = c(0,nChroms) ,
       ylim = c(0,max(outgroup.chroms$length[1:nChroms])/ 1000000),axes = F,
       xlab = "Scaffold",
       ylab = "Size (Mb)",
       main = sp)
  for(j in 1:nChroms){
    rect(xleft = (j-0.25),xright = (j+0.25),
         ybottom = 0,
         ytop = outgroup.chroms$length[j]/ mb,
         col = outgroup.chroms$col[j],
         border = "gray80",
         lwd = 0.5)
  }
  # plot axis
  if(is.null(axisLab)){
    axisLab <- 1:(nChroms)
  }
  axis(side = 2)
  axis(side = 1,
       at = 1:nChroms,
       labels = axisLab,
       tick = F)
  return(outgroup.chroms)
}


plotSynteny <- function(ingroup.chroms = NULL,
                        dat = NULL,
                        outgroup.chroms = NULL,
                        labels = NULL,
                        minSize = NULL,
                        sp = NULL){
  mb <- minSize
  # keep only chromosomes that are greater than 1 mb
  #plot
  plot(x = NULL,
       y = NULL,
       xlim = c(0,nrow(ingroup.chroms)) ,
       ylim = c(0,max(ingroup.chroms$length)/ mb),axes = F,
       xlab = "Scaffold",
       ylab = "Size (Mb)",
       main = sp)
  for(j in 1:nrow(ingroup.chroms)){
    chrom1 <- dat[dat$V4 %in% ingroup.chroms$chromosome[j],]
    # paint using majority syntenic chromosome
    syn_lengths <- chrom1$V6 - chrom1$V5
    #names(syn_lengths) <- chrom1$V6
    # minimum alignment length to 100
    chrom1 <- chrom1[syn_lengths >= 1,]
    if(nrow(chrom1) == 0){
      rect(xleft = (j-0.3),xright = (j+0.3),
           ybottom = 0,
           ytop = ingroup.chroms$length[j]/ mb,
           col = NA,border = "gray80",
           lwd = 0.5)
      next
    }
    for(i in 1:nrow(chrom1)){
      rect(xleft = (j-0.3),xright = (j+0.3),
           ybottom = chrom1$V5[i]/ mb,
           ytop = chrom1$V6[i]/ mb,
           col = outgroup.chroms$col[outgroup.chroms$chromosome %in% chrom1$V1[i]],
           border = outgroup.chroms$col[outgroup.chroms$chromosome %in% chrom1$V1[i]],
           lwd = 0)
    }
    rect(xleft = (j-0.3),xright = (j+0.3),
         ybottom = 0,
         ytop = ingroup.chroms$length[j]/ mb,
         col = NA,border = "gray80",
         lwd = 0.5)
    #box()
    chrom.list <- unique(chrom1$V1)
    synteny.length <- sum((chrom1$V6 - chrom1$V5)+1)
    hit <- c()
    for(k in 1:length(chrom.list)){
      hit[k] <- sum((chrom1$V6[chrom1$V1 %in% chrom.list[k]] - chrom1$V5[chrom1$V1 %in% chrom.list[k]])+1)
    }
    names(hit) <- (hit / synteny.length) * 100
  }
  
  if(is.null(labels)){
    labels <- 1:nrow(ingroup.chroms) 
  }
  
  axis(side = 2)
  axis(side = 1,
       at = 1:nrow(ingroup.chroms),
       labels = labels,
       tick = F,las = 2)
}

plotKimuraDistance <- function(DivsumTab = NULL,
                               genomesSize = NULL,
                               plotUnkowns = NULL,
                               repClassPlot = F,
                               savePlot = F,
                               zoom = F,
                               zoomPar = list(xmin = NULL,
                                              xmax = NULL,
                                              ymin = NULL,
                                              ymax = NULL,
                                              plotX1 = NULL,
                                              plotX2 = NULL,
                                              plotY1 = NULL,
                                              plotY2 = NULL)){
  
  # read kimura distance data
  KimuraDistance <- read.csv(DivsumTab,sep=" ")
  #add here the genome size in bp
  genomes_size <- genomesSize
  
  # get repeat class
  for(i in 1:ncol(KimuraDistance)){
    colnames(KimuraDistance)[i] <- unlist(strsplit(colnames(KimuraDistance)[i],split = ".",fixed = T))[1]
  }
  
  # subset data based on the repeat class
  DNA <- KimuraDistance[,colnames(KimuraDistance) == "DNA"]
  if(is.data.frame(DNA)){
    DNA_sum <- rowSums(DNA)
  }else{
    DNA_sum <- DNA  
  }
  
  LINE <- KimuraDistance[,colnames(KimuraDistance) == "LINE"]
  if(is.data.frame(LINE)){
    LINE_sum <- rowSums(LINE)
  }else{
    LINE_sum <- LINE  
  }
  
  LTR <- KimuraDistance[,colnames(KimuraDistance) == "LTR"]
  if(is.data.frame(LTR)){
    LTR_sum <- rowSums(LTR)
  }else{
    LTR_sum <- LTR  
  }
  
  SINE <- KimuraDistance[,colnames(KimuraDistance) == "SINE"]
  if(is.data.frame(SINE)){
    SINE_sum <- rowSums(SINE)
  }else{
    SINE_sum <- SINE  
  }
  
  Satellite <- KimuraDistance[,colnames(KimuraDistance) == "Satellite"]
  if(is.data.frame(Satellite)){
    Satellite_sum <- rowSums(Satellite)
  }else{
    Satellite_sum <- Satellite  
  }
  
  Simple_repeat <- KimuraDistance[,colnames(KimuraDistance) == "Simple_repeat"]
  if(is.data.frame(Simple_repeat)){
    Simple_repeat_sum <- rowSums(Simple_repeat)
  }else{
    Simple_repeat_sum <- Simple_repeat  
  }
  
  Unkown <- KimuraDistance[,colnames(KimuraDistance) == "Unknown"]
  if(is.data.frame(Unkown)){
    Unkown_sum <- rowSums(Unkown)
  }else{
    Unkown_sum <- Unkown  
  }
  
  if(plotUnkowns == T){
    
    # assign data to a new data fram for plotting
    reps <- as.data.frame(matrix(nrow = nrow(KimuraDistance),ncol = 8))
    colnames(reps) <- c("Div","DNA","LINE","LTR","SINE","Satellite","SSR", "Unkown")
    reps$Div <- KimuraDistance$Div
    reps$DNA <- DNA_sum
    reps$LINE <- LINE_sum
    reps$LTR <- LTR_sum
    reps$SINE <- SINE_sum
    reps$Satellite <- Satellite_sum
    reps$SSR <- Simple_repeat_sum
    reps$Unkown <- Unkown_sum
    
    # set colours
    col <-c("#d73027",
            "#fc8d59",
            "#fee090",
            "#e0f3f8",
            "#91bfdb",
            "#4575b4",
            "gray")
    
    # calculate the pecentage of the repeat sequences in the genome
    kd_melt = melt(reps,id="Div")
    kd_melt$norm = kd_melt$value/genomes_size * 100
    
  }
  
  if(plotUnkowns == F){
    # plot without unknowns
    # assign data to a new data fram for plotting
    reps <- as.data.frame(matrix(nrow = nrow(KimuraDistance),ncol = 7))
    colnames(reps) <- c("Div","DNA","LINE","LTR","SINE","Satellite","SSR")
    reps$Div <- KimuraDistance$Div
    reps$DNA <- DNA_sum
    reps$LINE <- LINE_sum
    reps$LTR <- LTR_sum
    reps$SINE <- SINE_sum
    reps$Satellite <- Satellite_sum
    reps$SSR <- Simple_repeat_sum
    
    # set colours
    col <-c("#d73027",
            "#fc8d59",
            "#fee090",
            "#e0f3f8",
            "#91bfdb",
            "#4575b4",
            "gray")
    
    # calculate the pecentage of the repeat sequences in the genome
    kd_melt = melt(reps,id="Div")
    kd_melt$norm = kd_melt$value/genomes_size * 100
  }
  
  # plot
  rPlot <- ggplot(kd_melt, aes(fill=variable, y=norm, x=Div)) + 
    geom_bar(position="stack", stat="identity",color="black", width= 0.75, linewidth = 0.05) +
    scale_fill_viridis(discrete = T) +
    scale_fill_manual("Repeat class",values=col) +
    theme_classic() +
    xlab("Kimura substitution level (CpG adjusted)") +
    ylab("Percent of the genome") + 
    labs(fill = "") +
    coord_cartesian(xlim = c(0, 50)) +
    theme(axis.text=element_text(size=12),
          axis.title =element_text(size=14),
          panel.grid.major.y = element_line(colour = "gray30",linewidth = 0.05,linetype = 5),
          panel.grid.minor.y = element_line(colour = "gray70",linewidth = 0.05,linetype = 5),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14,face = "bold"))
  
  finalPlot <- rPlot
  
  if(repClassPlot == T){
    repClass <- ggplot() + 
      geom_line(data = kd_melt,mapping = aes(colour=variable, y=norm, x=Div),linewidth = 1) +
      xlab("Kimura substitution level (CpG adjusted)") +
      ylab("Percent of the genome") + 
      labs(fill = "") +
      coord_cartesian(xlim = c(0, 50)) +
      scale_color_discrete("Repeat class",type = col) +
      theme_classic() +
      theme(axis.text=element_text(size=12),
            axis.title =element_text(size=14),
            panel.grid.major.y = element_line(colour = "gray30",linewidth = 0.05,linetype = 5),
            panel.grid.minor.y = element_line(colour = "gray70",linewidth = 0.05,linetype = 5),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14,face = "bold"))
    
    finalPlot <- ggarrange(rPlot, repClass, nrow = 1, labels = c("A","B"))
  }
  
  if(zoom == T){
    p1 <- ggplot() +
      geom_line(data = kd_melt,mapping = aes(colour=variable, y=norm, x=Div),linewidth = 1) +
      xlab("Kimura substitution level (CpG adjusted)") +
      ylab("Percent of the genome") + 
      labs(fill = "") +
      coord_cartesian(xlim = c(0, 50)) +
      geom_rect(aes(xmin = zoomPar$xmin,
                    xmax = zoomPar$xmax, 
                    ymin = zoomPar$ymin, 
                    ymax = zoomPar$ymax), color = "black", alpha = 0)+
      scale_color_discrete("Repeat class",type = col) +
      theme_classic() +
      theme(axis.text=element_text(size=12),
            axis.title =element_text(size=14),
            panel.grid.major.y = element_line(colour = "gray30",linewidth = 0.05,linetype = 5),
            panel.grid.minor.y = element_line(colour = "gray70",linewidth = 0.05,linetype = 5),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14,face = "bold"))
    
    p2 <- ggplot() +
      geom_line(data = kd_melt,mapping = aes(colour=variable, y=norm, x=Div),linewidth = 1) +
      scale_color_discrete("Repeat class",type = col) +
      theme_classic() +
      xlab("") +
      ylab("") + 
      labs(fill = "") +
      coord_cartesian(xlim = c(zoomPar$xmin,
                               zoomPar$xmax),
                      ylim = c(zoomPar$ymin,
                               zoomPar$ymax)) +
      theme(axis.text=element_text(size=9),
            axis.title =element_text(size=9),
            panel.grid.major.y = element_line(colour = "gray30",linewidth = 0.05,linetype = 0),
            panel.grid.minor.y = element_line(colour = "gray70",linewidth = 0.05,linetype = 5),
            legend.position = "none",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14))
    
    data <-data.frame(x= c(zoomPar$xmin,
                           zoomPar$plotX1,
                           zoomPar$xmax,
                           zoomPar$plotX2),
                      y=c(zoomPar$ymax,
                          zoomPar$plotY1,
                          zoomPar$ymax,
                          zoomPar$plotY1),
                      grp=c(1,1,2,2))
    
    p3 <- p1 + 
      annotation_custom(ggplotGrob(p2), 
                        xmin = zoomPar$plotX1,
                        xmax = zoomPar$plotX2,
                        ymin = zoomPar$plotY1,
                        ymax = zoomPar$plotY2) +
      geom_rect(aes(xmin = zoomPar$plotX1,
                    xmax = zoomPar$plotX2,
                    ymin = zoomPar$plotY1,
                    ymax = zoomPar$plotY2), 
                color='black', linetype='dashed', alpha=0) +
      geom_path(data = data,aes(x = x, y = y, group = grp),linetype='dashed')
    
    finalPlot <- ggarrange(rPlot, p3, nrow = 1, labels = c("A","B"))
  }
  
  return(finalPlot)
  
  if(savePlot == T){
    # save plot
    ggsave(
      "repeatLandscapePlot.pdf",
      plot = finalPlot,
      device = "pdf",
      path = getwd(),
      scale = 1,
      width = 12,
      height = 5,
      units = "in",
      dpi = 300,
      limitsize = TRUE,
      bg = NULL
    )  
  }
}