#written by noah friedman
#This script contains multiple functionalities:
#1. the function find_and_save_all_distributions will take a file with 
#cancer types and tmbs and for each one performs a clustering using Ckmeans.1d.dp and parameters you set
#then it writes the results to a directory you specify if desired
#2. the function load_and_plot_all_distributions which plots all distributions together using the function plot_distribution

#note please refer to the attached ipython script for functionality for how I change classifications to highMutation burden for certain cancer types
#please refer to the attached flow chart for a visual explanation of how this method works

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(stringr)

library("mclust")
library(mclust, quietly=TRUE)
library(gridBase)
library(Ckmeans.1d.dp)

plot_distribution <- function(df, binW = 1, title='', hideLegend=TRUE){
  
  emptyTheme <- theme(axis.line = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank())
  
  
  plt<- ggplot(df, aes(x=Nmut_Mb))+
    
    geom_histogram(binwidth = binW, aes(fill=factor(hypermutantClassification, levels=c('highMutationBurden', 'Hypermutated', 'Indeterminate', 'Normal'))))+
    #geom_density(aes(y=binW * ..count..), alpha=0.5)+
    
    theme(axis.text.x=element_text(size=3, face="bold"))+
    scale_fill_manual(values = c('orange', 'maroon', '#2a2a2a', 'gray'), drop=FALSE)+
    coord_cartesian(ylim=c(0,100), xlim=c(0,500))+
    ggtitle(title)+
    scale_y_continuous(breaks=c(0,25,50,100), labels=c('0', '25', '50', '>100'))+
    theme(axis.text.x = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=20),
          plot.title = element_text(size=30)
          )+
    ylab('N cases')+
    xlab('Mutations/Megabase')+
    emptyTheme+
    guides(fill=guide_legend(title="Hypermutant Classification"))
  if(hideLegend){
    plt <- plt + theme(legend.position = 'none') 
  }
  return(plt)
}

#New method written by Noah Friedman on 11/19
fit_and_analyze_dist <- function(dataFrame,
                                 quantileBufferSize=.1, #the quantile determines how many cases we throw out from the hypermutant cluster
                                 minimumHypermutantThresh = 20, #this is the number below which we well never consider a case hypermutated
                                 maxNClusters = 10
){
  
  dataFrame = dataFrame[!is.na(dataFrame$tmb),]
  
  clus <- Ckmeans.1d.dp(dataFrame$tmb, k=c(1,maxNClusters)) #WE cluster using ckmeans 1d on log(mutations per megabase).
  
  #create a dataframe with the clustering results
  df <- data.frame(c(clus$cluster), dataFrame$Tumor_Sample_Barcode, dataFrame$tmb)
  colnames(df) <- c("cluster", "Tumor_Sample_Barcode", "Nmut_Mb")

  nClusters <- length(unique(clus$cluster))
  minimumHypermutatedCluster = nClusters + 1
  for(i in seq(1, nClusters)){
    if(clus$centers[[i]] > minimumHypermutantThresh){
      minimumHypermutatedCluster = i
      break
    }
  }
  #IF SOMEHOW WE DIDNT FIND ANY CLUSTERS MARK THAT DOWN
  if(nClusters == 1){
    df$HypermutantThresh = 'Only one cluster found'
    df$NormalThresh = 'Only one cluster found'
    df$hypermutantClassification = 'Normal-No Clusters Found'
    return(df)
  }
  else{
    hypermutantCluster = df[df$cluster >= minimumHypermutatedCluster,] #the numerically largest cluster always has the biggest mean
    notHypermutantCluster = df[df$cluster < minimumHypermutatedCluster,]
    
    #Create the bounds for the indeterimante clusters using this information and the parameters passed in
    indeterminateUpperBound <- max(quantile(hypermutantCluster$Nmut_Mb, c(quantileBufferSize)), minimumHypermutantThresh)
    indeterminateLowerBound <- min(minimumHypermutantThresh, quantile(notHypermutantCluster$Nmut_Mb, c(1)))

    #make a dataframe with this info
    df <- mutate(df,
                 hypermutantClassification = ifelse(Nmut_Mb>=indeterminateUpperBound, "Hypermutated",
                                                    ifelse(Nmut_Mb>=indeterminateLowerBound, "Indeterminate",
                                                           ifelse(Nmut_Mb<indeterminateLowerBound, "Normal", "No_Nmut_Mb_Info"))))
    #FIX classifications that end up as NA
    df$hypermutantClassification <- sapply(df$hypermutantClassification, function(x) if(is.na(x)) 'Normal' else x)
    
    df$HypermutantThresh = indeterminateUpperBound
    df$NormalThresh = indeterminateLowerBound
    return(df)
  }
}

#a big omnibus function that does the following:
#finds hypermutants based on clustering
#simulataneously if you specify the save option with your own save path it will save them
#relies on you to have previousy adjust the $cancerType column beforehand
find_and_save_all_distributions <- function(tmbData, minNSamples=250,
  quantileBuffer=.1, #the quantile determines how many cases we throw out from the hypermutant cluster
  minHypermutantThresh = 20, #this is the number below which we well never consider a case hypermutated
  maxClusters = 10, #NOTE that the choice of between 1 and 3 clusters is arbitrary its just what work
  writeData=TRUE, writePath ='~/Desktop/mnt/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds/'){
  
  cancerTypes <- unique(tmbData$cancerType)
  for (cancerType in cancerTypes){
    print(cancerType)
    analyzeD = tmbData[tmbData$cancerType == cancerType,]
    if(dim(analyzeD)[[1]] >= minNSamples){
      data <- fit_and_analyze_dist(analyzeD, quantileBufferSize=quantileBuffer,
                                       minimumHypermutantThresh = minHypermutantThresh,
                                       maxNClusters = maxClusters)
      
      if(writeData == TRUE){
        if(nchar(cancerType) > 1){
          print(paste('writingData for ', cancerType))
          cancerTypeAdj <- gsub(" ", "_", cancerType)
          fullPath <- paste(writePath, cancerTypeAdj, '.tsv', sep='')
          write.table(data, file=fullPath, sep='\t')
        }
      }
    }
  }
}

#load the distributions (they may have been changed by my python script/any other steps the user deems necessary)
load_and_plot_all_distributions <-  function(fileDir = '~/Desktop/mnt/juno/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds/'){
  l <- list()
  i <- 1
  for(file in list.files(fileDir)){
    df <- read.table(paste(fileDir, file, sep = ""), sep='\t', header=TRUE)
    
    cancerType = str_replace(file, ".tsv", "")
    p <- plot_distribution(df, title=cancerType)
    l[[i]] <- p
    i <- i + 1
  }
  #dummy plot to get a legend
  legend <- get_legend(plot_distribution(df, title=cancerType, hideLegend=FALSE))
  print(legend)
  l[[i]] <- legend
  return(l)
}

tmbData <- read.table('~/Desktop/mnt/juno/work/taylorlab/friedman/myAdjustedDataFiles/mutations_TMB_and_MSI_stats.txt', sep='\t', header=TRUE)

#FIND THE DISTRIBUTIONS AND SAVE THEM IF NEEDED
find_and_save_all_distributions(tmbData, minNSamples=100, 
    writeData=TRUE, writePath ='~/Desktop/mnt/juno/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds/')

#MAKE THE BIG COMINED PLOT
listOfPlots <- load_and_plot_all_distributions()
p <- plot_grid(plotlist=listOfPlots)
plotWithTitleAndCaption <- plot_grid(ggplot()+ggtitle('Hypermutation Classification')+theme(plot.title=element_text(size=50)),
                                    p, ggplot()+labs(caption='plotAndDefineHypermutationThreshold.R\ngenerate_mut_classification_figure#1a#.ipynb'),
                                    nrow=3, rel_heights = c(.1,1,.1))
                                    
ggsave('~/Desktop/plot.pdf', plot=plotWithTitleAndCaption,  width = 60, height = 20, units = c("in"), limitsize = FALSE)

#make the smaller plot as figure 1a.

endometrialData <- read.table('~/Desktop/mnt/juno/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds/Endometrial_Cancer.tsv', sep='\t', header=TRUE)
colorectalData <- read.table('~/Desktop/mnt/juno/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds/Colorectal_Cancer.tsv', sep='\t', header=TRUE)
bladderData <- read.table('~/Desktop/mnt/juno/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds/Bladder_Cancer.tsv', sep='\t', header=TRUE)

pEndo <- plot_distribution(endometrialData, title='Endometrial Cancer')
pColo <- plot_distribution(colorectalData, title='Colorectal Cancer')
pBladder <- plot_distribution(bladderData, title='Bladder Cancer')
legend <- get_legend(plot_distribution(endometrialData, title='endo', hideLegend=FALSE))
alignedPlot <- plot_grid(pEndo, pColo, pBladder, nrow=3)
alignedPlotWithLegend <- plot_grid(alignedPlot, legend, ncol=2, rel_widths = c(1,.5))

finalPlot <- plot_grid(alignedPlotWithLegend, ggplot()+labs(caption='plotAndDefineHypermutationThreshold.R\ngenerate_mut_classification_figure#1a#.ipynb'), nrow=2, rel_heights=c(1,.1))
ggsave('~/Desktop/plot.pdf', plot=finalPlot,  width = 16, height = 12, units = c("in"), limitsize = FALSE)
