#################################
###                           ###
### Changing PCOA plot colors ###
###      By Jessika Neves     ###
###         March 2021        ###
###                           ###
#################################
#
#This script was developed under the dartR version 1.1.6 available at https://cran.r-project.org/src/contrib/Archive/dartR/
#Run this function before run the main script
#
#Load the packages
library(dartR)
library(ggplot2)
library(adegenet)
library(directlabels)
#
gl.pcoa.plot <- function(glPca, data, scale=FALSE, ellipse=FALSE, p=0.95, labels="pop", hadjust=1.5, 
                         vadjust=1, xaxis=1, yaxis=2) {
  
  if(class(glPca)!="glPca" | class(data)!="genlight") {
    cat("Fatal Error: glPca and genlight objects required for glPca and data parameters respectively!\n"); stop()
  }
  
  # Tidy up the parameters
  #  if (labels=="smart") { hadjust <- 0; vadjust <- 0 }
  
  # Create a dataframe to hold the required scores
  m <- cbind(glPca$scores[,xaxis],glPca$scores[,yaxis])
  df <- data.frame(m)
  
  # Convert the eigenvalues to percentages
  s <- sum(glPca$eig)
  e <- round(glPca$eig*100/s,1)
  
  # Labels for the axes
  xlab <- paste("PCoA Axis", xaxis, "(",e[xaxis],"%)")
  ylab <- paste("PCoA Axis", yaxis, "(",e[yaxis],"%)")
  
  
  # If labels = none
  
  if (labels == "none" | labels==FALSE) {
    cat("Plotting points with no labels\n")
    pop <- factor(pop(data))
    df <- cbind(df,pop)
    colnames(df) <- c("PCoAx","PCoAy","pop")
    
    # Plot
    p <- ggplot(df, aes(x=df$PCoAx, y=df$PCoAy,colour=pop)) +
      geom_point(shape=1, size=5,colour=c("darkorange2", "deepskyblue2", "darkviolet","darkorchid1", "deeppink2", "chartreuse3", "blue2")[pop],aes(colour=pop)) +
      #geom_dl(aes(label=ind),method="first.points") +
      #ggtitle(paste("PCoA Plot")) +
      theme(axis.title=element_text(face="bold.italic",size="20", color="black"),
            axis.text.x  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            axis.text.y  = element_text(face="bold",angle=0, vjust=0.5, size=10),
            legend.title = element_text(colour="black", size=18, face="bold"),
            legend.text = element_text(colour="black", size = 16, face="bold")
      ) +
      labs(x=xlab, y=ylab) +
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0)+
      theme(legend.position="left")
    # Scale the axes in proportion to % explained, if requested
    if(scale==TRUE) { p <- p + coord_fixed(ratio=e[yaxis]/e[xaxis]) }
    # Add ellipses if requested
    if(ellipse==TRUE) {p <- p + stat_ellipse(aes(colour=pop), type="norm", level=0.95)}
  }
  
  return (p)
}
###############################################################################################################################
#
#jessika.neves@icbs.ufal.br

