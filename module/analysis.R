draw_clust<-function(cancertype,feature){
  #-read CNV
  GBM_CNV <- read.csv("./data/GBM/GBM_CNV_core.txt", stringsAsFactors = FALSE,sep="\t")
  KIRC_CNV <- read.csv("./data/KIRC/KIRC_CNV_core.txt", stringsAsFactors = FALSE,sep="\t")
  OV_CNV <- read.csv("./data/OV/OV_CNV_core.txt", stringsAsFactors = FALSE,sep="\t")
  #chrom_posGBM=""
  #chrom_posKIRC=""
  #chrom_posOV=""
  
  LUSC_CNV <- read.csv("./data/LUSC/LUSC_CNV_core.txt", stringsAsFactors = FALSE,sep="\t")
  LUSC_CNV<-LUSC_CNV[,1:(ncol(LUSC_CNV)-1)]
  chrom_posLUSC<-colnames(LUSC_CNV[,-1])
  #-read methylation
  GBM_methyl <- read.csv("./data/GBM/GBM_methylation_core.txt", stringsAsFactors = FALSE,sep="\t")
  KIRC_methyl <- read.csv("./data/KIRC/KIRC_methylation_core.txt", stringsAsFactors = FALSE,sep="\t")
  OV_methyl <- read.csv("./data/OV/OV_methylation_core.txt", stringsAsFactors = FALSE,sep="\t")
  
  methyl_array<-read.csv("./data/methyl_array.txt", stringsAsFactors = FALSE,sep="\t")
  methyl_array_ref<-subset(methyl_array,methyl_array[,2]!="")
  
  library(cluster)
  library(fpc)
  ccnt=0
  fcnt=0
  for(i in 1:length(cancertype)){
    #cat(cancertype[i])
    switch(cancertype[i],
          "1"={
            cat("1")
            ccnt=ccnt+1
            for(j in 1:length(feature)){
              #cat(feature[i])
              switch(feature[i],
                     "1"={
                       
                       GBM_CNV<-GBM_CNV[,1:(ncol(GBM_CNV)-1)]
                       #chrom_posGBM<-colnames(GBM_CNV[,-1])
                     }
                )
              }
            
          },
          "2"={
                cat("2")
                ccnt=ccnt+2
                for(j in 1:length(feature)){ 
                  switch(feature[i],  
                      "1"={ 
                        
                        KIRC_CNV<-KIRC_CNV[,1:(ncol(KIRC_CNV)-1)]
                        #chrom_posKIRC<-colnames(KIRC_CNV[,-1])
                      }
                  )
                }
          },
          "3"={
            cat("3")
            ccnt=ccnt+3
            for(j in 1:length(feature)){ 
              switch(feature[i],  
                     "1"={ 
                       
                       OV_CNV<-OV_CNV[,1:(ncol(OV_CNV)-1)]
                       #chrom_posOV<-colnames(OV_CNV[,-1])
                     }
              )
            }            
               
          }
    )
    
  }
  
  for(j in 1:length(feature)){
    fcnt=fcnt+as.numeric(feature[j])
  }
  
  if(ccnt==1){ #brain cancer
    if(fcnt==1){ #CNV
        chrom_posGBM<-colnames(GBM_CNV[,-1])
        common_index<-chrom_posGBM
        com_GBM_CNV<-GBM_CNV[,common_index]
        rownames(com_GBM_CNV)<-paste("GBM_CNV",seq(1,nrow(com_GBM_CNV)),sep="_")
        View(com_GBM_CNV)
        k.means.fit<-kmeans(com_GBM_CNV,1)
        clusplot(cnv_com_total[,1:41],k.means.fit$cluster,color=TRUE,labels=4,xlab = "", ylab = "",main="Clustering Plot(Copy Number Variants)",sub ="Group 1=GBM")
        
    }
    if(fcnt==2){ #DNV meth
        GBM_fe<-GBM_methyl[,methyl_index]
        k.means_met.fit<-kmeans(GBM_fe,1)
        clusplot(GBM_fe[,1:300],k.means_met.fit$cluster,color=TRUE,labels=4,xlab = "", ylab = "",main="Clustering Plot(DNA methylation)",sub ="Group 1=GBM")    }
    if(fcnt==3){ # both
        chrom_posGBM<-colnames(GBM_CNV[,-1])
        com_GBM_CNV<-GBM_CNV[,common_index]
        rownames(com_GBM_CNV)<-paste("GBM_CNV",seq(1,nrow(com_GBM_CNV)),sep="_")
        GBM_fe<-GBM_methyl[,methyl_index]
        GBM_fe[is.na(GBM_fe)]<-0
        GBM_total_fe<-cbind(com_GBM_CNV,GBM_fe)
        cnv_methyl_total<-rbind(GBM_total_fe[,1:1709])
        #cnv_methyl_total[is.na(cnv_methyl_total)]<-0
        k.means_mix.fit<-kmeans(cnv_methyl_total,3)
        
        clusplot(cnv_methyl_total[,1:341],k.means_mix.fit$cluster,color=TRUE,labels=4,xlab = "", ylab = "",main="Clustering Plot(CNV+DNA methylation)",sub ="Group 1=GBM")
        
    }
    
  }
  
  if(ccnt==6){
    if(fcnt==1){
      chrom_posGBM<-colnames(GBM_CNV[,-1])  #chromosome position
      #chrom_posLUSC<-colnames(LUSC_CNV[,-1]) 
      chrom_posKIRC<-colnames(KIRC_CNV[,-1])
      chrom_posOV<-colnames(OV_CNV[,-1])
      common_index<-Reduce(intersect, list(chrom_posGBM,chrom_posKIRC,chrom_posOV)) #tupling in same DNA
      
      com_GBM_CNV<-GBM_CNV[,common_index]
      rownames(com_GBM_CNV)<-paste("GBM_CNV",seq(1,nrow(com_GBM_CNV)),sep="_")
      #com_LUSC_CNV<-LUSC_CNV[,common_index]
      #rownames(com_LUSC_CNV)<-paste("LUSC_CNV",seq(1,nrow(com_LUSC_CNV)),sep="_")
      com_KIRC_CNV<-KIRC_CNV[,common_index]
      rownames(com_KIRC_CNV)<-paste("KIRC_CNV",seq(1,nrow(com_KIRC_CNV)),sep="_")
      com_OV_CNV<-OV_CNV[,common_index]
      rownames(com_OV_CNV)<-paste("OV_CNV",seq(1,nrow(com_OV_CNV)),sep="_")
      
      cnv_com_total<-rbind(com_GBM_CNV,com_KIRC_CNV,com_OV_CNV)
      View(cnv_com_total)
      k.means.fit<-kmeans(cnv_com_total,3)
      
      clusplot(cnv_com_total[,1:41],k.means.fit$cluster,color=TRUE,labels=4,xlab = "", ylab = "",main="Clustering Plot(Copy Number Variants)",sub ="Group 1=GBM, Group 2=KIRC, Group 3=OV")
    }
    if(fcnt==2){
      methyl_index<-Reduce(intersect, list(colnames(GBM_methyl[,2:(ncol(GBM_methyl)-1)]),colnames(KIRC_methyl[,2:(ncol(KIRC_methyl)-1)]),colnames(OV_methyl[,2:(ncol(OV_methyl)-1)])))
      GBM_fe<-GBM_methyl[,methyl_index]
      KIRC_fe<-KIRC_methyl[,methyl_index]
      OV_fe<-OV_methyl[,methyl_index]
      
      methyl_com_total<-rbind(GBM_fe,KIRC_fe,OV_fe)
      methyl_com_total[is.na(methyl_com_total)]<-0
      k.means_met.fit<-kmeans(methyl_com_total,3)
      clusplot(methyl_com_total[,1:300],k.means_met.fit$cluster,color=TRUE,labels=4,xlab = "", ylab = "",main="Clustering Plot(DNA methylation)",sub ="Group 1=GBM, Group 2=KIRC, Group 3=OV")
    }
    if(fcnt==3){
      #CNV
      chrom_posGBM<-colnames(GBM_CNV[,-1])
      chrom_posLUSC<-colnames(LUSC_CNV[,-1])
      chrom_posKIRC<-colnames(KIRC_CNV[,-1])
      chrom_posOV<-colnames(OV_CNV[,-1])
      common_index<-Reduce(intersect, list(chrom_posGBM,chrom_posKIRC,chrom_posLUSC,chrom_posOV))
      
      com_GBM_CNV<-GBM_CNV[,common_index]
      rownames(com_GBM_CNV)<-paste("GBM_CNV",seq(1,nrow(com_GBM_CNV)),sep="_")
      com_KIRC_CNV<-KIRC_CNV[,common_index]
      rownames(com_KIRC_CNV)<-paste("KIRC_CNV",seq(1,nrow(com_KIRC_CNV)),sep="_")
      com_OV_CNV<-OV_CNV[,common_index]
      rownames(com_OV_CNV)<-paste("OV_CNV",seq(1,nrow(com_OV_CNV)),sep="_")
      
      #DNA methylation
      methyl_index<-Reduce(intersect, list(colnames(GBM_methyl[,2:(ncol(GBM_methyl)-1)]),colnames(KIRC_methyl[,2:(ncol(KIRC_methyl)-1)]),colnames(OV_methyl[,2:(ncol(OV_methyl)-1)])))
      GBM_fe<-GBM_methyl[,methyl_index]
      GBM_fe[is.na(GBM_fe)]<-0
      KIRC_fe<-KIRC_methyl[,methyl_index]
      KIRC_fe[is.na(KIRC_fe)]<-0
      OV_fe<-OV_methyl[,methyl_index]
      OV_fe[is.na(OV_fe)]<-0
      #--combine
      GBM_total_fe<-cbind(com_GBM_CNV,GBM_fe)
      KIRC_total_fe<-cbind(com_KIRC_CNV,KIRC_fe)
      OV_total_fe<-cbind(com_OV_CNV,OV_fe)
      
      cnv_methyl_total<-rbind(GBM_total_fe[,1:1709],KIRC_total_fe[,1:1709],OV_total_fe[,1:1709])
      #cnv_methyl_total[is.na(cnv_methyl_total)]<-0
      k.means_mix.fit<-kmeans(cnv_methyl_total,3)
      
      clusplot(cnv_methyl_total[,1:341],k.means_mix.fit$cluster,color=TRUE,labels=4,xlab = "", ylab = "",main="Clustering Plot(CNV+DNA methylation)",sub ="Group 1=GBM, Group 2=KIRC, Group 3=OV")
      
    }
    
  }
  
}