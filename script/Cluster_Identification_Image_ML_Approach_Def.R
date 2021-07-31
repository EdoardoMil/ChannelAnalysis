require(raster)
require(reshape2)
require(plyr)
require(png)
require(Matrix)
require(cluster)
require(REdaS)
require(OpenImageR)
require(opencv)
require(ggplot2)
require(tiff)

# path definition
#dir_script <- "/Users/edo/Documents/Progetti/ImageGPCR/script/"
dir_im <- "../TIFF/"
dir_output_ML <- "../results_temp_singleIm/"
dir_output_Zern <- "../Zernike_temp_singleIm/"
dir_BestFigC4 <- "../Figure_Best_C4/"
dir_BestZernC4 <- "../Zernike_Best_C4/"

#dir_basic <- "/Users/edo/Documents/Progetti/ImageGPCR/"

system(paste("rm -r"," ", dir_BestFigC4, sep=""))
system(paste("mkdir"," ", dir_BestFigC4, sep=""))

system(paste("rm -r"," ", dir_BestZernC4, sep=""))
system(paste("mkdir"," ", dir_BestZernC4, sep=""))


FL <- list.files(path = dir_im, pattern = ".tif")

#Parameter definition
SizeMax <- 2000
conta <- 0
SizeSide <- 250
SizeSide_IM <- 80
NoiseAdd <- 2
SizeMin <- 150
SizeMax <- 2000
LatoProbe <- 3
ing <- 3
edge <- 25

#Starting the loop on the images 
for(f in 1:length(FL)){
  
  str_name<-paste(dir_im,FL[f], sep="")
  name_file <- unlist(strsplit(FL[f],".tiff"))
  
  a <- readTIFF(str_name, as.is = TRUE)
  a <- round((a*255)/max(a))
  
  rainbowcols <- topo.colors(20)
  rainbowcols <- c(rainbowcols, "white")
  rainbowcols <- colorRampPalette(c("black","gray","white"))
  rainbowcols <- rainbowcols(120)
  
  #If you want to see the image or part of it.
  #plot(raster(a), col=rainbowcols)
  #plot(raster(a[1:1500,1:1500]), col=rainbowcols)
  
  #smoothing 1
  c <- a
  RhoMatIm <- c
  LatoProbe <- 10
  for(i in (1+LatoProbe):(nrow(c)-LatoProbe)){
    print(i)
    for(j in (1+LatoProbe):(ncol(c)-LatoProbe)){
      matAus <-c[(i-LatoProbe):(i+LatoProbe),(j-LatoProbe):(j+LatoProbe)]
      rho <- mean(matAus)
      RhoMatIm[i,j] <- rho
    }
  }
  plot(raster(RhoMatIm), col=rainbowcols)
  plot(raster(RhoMatIm[1:1500,1:1500]), col=rainbowcols)
  
  #smoothing 2
  b <- a
  b[b>mean(b)] <- 255
  
  RhoMat <- b
  LatoProbe <- 5
  for(i in (1+LatoProbe):(nrow(b)-LatoProbe)){
    print(i)
    for(j in (1+LatoProbe):(ncol(b)-LatoProbe)){
      matAus <-b[(i-LatoProbe):(i+LatoProbe),(j-LatoProbe):(j+LatoProbe)]
      rho <- mean(matAus)
      RhoMat[i,j] <- rho
    }
  }
  
  #If you want to see the image or part of it.
  #plot(raster(RhoMat[1:1500,1:1500]), col=rainbowcols)
  #plot(raster(RhoMat[700:1300,20:700]), col=rainbowcols)
  
  
  #binarization on the filter, eliminating values below a certain threshold
  RhoMatBin <- RhoMatIm
  #hist(RhoMatIm, 100)
  cutoffBin <- mean(RhoMatIm)-2*sd(RhoMatIm)
  cutoffBin <- mean(as.vector(RhoMatIm)[as.vector(RhoMatIm)!=255])-1.8*sd(as.vector(RhoMatIm)[as.vector(RhoMatIm)!=255])
  #abline(v=cutoffBin, col="red", lwd=2)
  
  RhoMatBin[RhoMatBin <= cutoffBin] <- 0.1
  RhoMatBin[RhoMatBin > cutoffBin] <- 0
  RhoMatBin <- RhoMatBin*10
  
  #If you want to see the image or part of it.
  #plot(raster(RhoMatBin), col=rainbowcols)
  #plot(raster(RhoMatBin[1:1500,1:1500]), col=rainbowcols)
  
  RhoMatBin_Melt <- melt(RhoMatBin)
  sm = RhoMatBin_Melt[RhoMatBin_Melt[,3]==1,1:2]
  rownames(sm) <- NULL
  colnames(sm) <- c("i","j")
  dim(sm)
  
  SizeMin <- 150
  MaxClusterCalcul <- 30000
  MaxClCalAus_inf <- 1
  MaxClCalAus_sup <- 30000
  Ncl_aus <- 0
  SizeCluster <- 50 
  LatoProbe <- 15
  
  NumStep <- round(dim(sm)[1]/MaxClusterCalcul)
  
  RhoMat_Cluster_Melt <- data.frame()
  for(p in 1:round(nrow(sm)/MaxClusterCalcul)){
    print(paste("File Number: ", f, " -- Step Sub-Region: ", p, " -- Total Step: ", NumStep, sep=""))
    
    if(MaxClCalAus_sup <= nrow(sm)){
      sm_aus <- sm[MaxClCalAus_inf:MaxClCalAus_sup,]
      d = dist(sm_aus, "manhattan")
      gr = cutree(hclust(d, "single"), h = 1)
      RhoMat_Cluster_aus <- sparseMatrix(i = sm_aus[, "i"], j = sm_aus[, "j"], x = gr)
      RhoMat_Cluster_Melt_aus <- summary(RhoMat_Cluster_aus)
      
      names_largeSize <- names(table(summary(RhoMat_Cluster_aus)[,3])[table(summary(RhoMat_Cluster_aus)[,3])>=SizeCluster])
      RhoMat_Cluster_Melt_aus <- RhoMat_Cluster_Melt_aus[RhoMat_Cluster_Melt_aus$x %in% as.numeric(names_largeSize),]
      
      RhoMat_Cluster_Melt_aus$x <- RhoMat_Cluster_Melt_aus$x + Ncl_aus
      RhoMat_Cluster_Melt <- as.matrix(rbind(RhoMat_Cluster_Melt, RhoMat_Cluster_Melt_aus))
      
      MaxClCalAus_inf <- MaxClCalAus_inf + MaxClusterCalcul
      MaxClCalAus_sup <- MaxClCalAus_sup + MaxClusterCalcul
      
      if(length( max(as.numeric(names_largeSize))) > 0){
        Ncl_aus <- Ncl_aus + max(as.numeric(names_largeSize))
      }else{
        Ncl_aus <- Ncl_aus + 0
      }
      
    }else{
      
      sm_aus <- sm[MaxClCalAus_inf:nrow(sm),]
      d = dist(sm_aus, "manhattan")
      gr = cutree(hclust(d, "single"), h = 1)
      RhoMat_Cluster_aus <- sparseMatrix(i = sm_aus[, "i"], j = sm_aus[, "j"], x = gr)
      RhoMat_Cluster_Melt_aus <- summary(RhoMat_Cluster_aus)
      
      names_largeSize <- names(table(summary(RhoMat_Cluster_aus)[,3])[table(summary(RhoMat_Cluster_aus)[,3])>=SizeCluster])
      RhoMat_Cluster_Melt_aus <- RhoMat_Cluster_Melt_aus[RhoMat_Cluster_Melt_aus$x %in% as.numeric(names_largeSize),]
      
      RhoMat_Cluster_Melt_aus$x <- RhoMat_Cluster_Melt_aus$x + Ncl_aus
      RhoMat_Cluster_Melt <- as.matrix(rbind(RhoMat_Cluster_Melt, RhoMat_Cluster_Melt_aus))
      
    }
  }
  
  RhoMat_Cluster_Melt <- as.data.frame(RhoMat_Cluster_Melt)
  colnames(RhoMat_Cluster_Melt) <- c("i","j","x")
  
  #Calculate the size of each cluster 
  Ngroups <- length(unique(RhoMat_Cluster_Melt$x))
  SizeCluster <- c()
  for(i in 1:Ngroups){
    print(i)
    aus <- unique(RhoMat_Cluster_Melt$x)[i]
    size_aus <- nrow(RhoMat_Cluster_Melt[RhoMat_Cluster_Melt$x==aus,])
    SizeCluster <- c(SizeCluster, size_aus)
  }
  
  
  #Cluster distribution
  RhoMat_Cluster_Melt$SizeCluster <- 0
  for(i in 1:Ngroups){
    aus <- unique(RhoMat_Cluster_Melt$x)[i]
    RhoMat_Cluster_Melt[RhoMat_Cluster_Melt$x==aus,"SizeCluster"] <- SizeCluster[i]
  }
  
  #hist(SizeCluster,100)
  #abline(v=SizeMin, col="red", lwd=2)
  
  RhoMat_Cluster_Melt_Red <- RhoMat_Cluster_Melt[RhoMat_Cluster_Melt$SizeCluster >= SizeMin & RhoMat_Cluster_Melt$SizeCluster <= SizeMax,]
  Ngroups_red <- length(unique(RhoMat_Cluster_Melt_Red$x))
  groups_red <- unique(RhoMat_Cluster_Melt_Red$x)
  
  RhoMatBinRed <- matrix(0, nrow = nrow(RhoMatBin), ncol = ncol(RhoMatBin))
  for(i in 1:nrow(RhoMat_Cluster_Melt_Red)){
    print(i)
    i_aus <- RhoMat_Cluster_Melt_Red[i,"i"]
    j_aus <- RhoMat_Cluster_Melt_Red[i,"j"]
    val_aus <- RhoMat_Cluster_Melt_Red[i,"x"]
    val_aus <- 1
    RhoMatBinRed[i_aus,j_aus] <- val_aus
  }
  
  #plot(raster(RhoMatBin), col=rainbowcols)
  #plot(raster(RhoMatBinRed), col=rainbowcols)
  #plot(raster(RhoMatBinRed[1:800,1:800]), col=rainbowcols)
  
  RhoMatBinRed_fig <- RhoMatBinRed
  RhoMatBinRed_fig[RhoMatBinRed_fig==0] <- 2
  RhoMatBinRed_fig[RhoMatBinRed_fig==1] <- 0
  RhoMatBinRed_fig[RhoMatBinRed_fig==2] <- 1
  #plot(raster(RhoMatBinRed_fig[1:1500,1:1500]), col=rainbowcols)
  
  
  #Mapping each cluster identified in the original matrix to obtain the image of the single cluster
  VetPercCutoff <- seq(from = 0, to = 1.5, by = 0.1)
  Ncut <- length(VetPercCutoff)

  par(mfrow=c(3,3))
  for(n in 1:Ngroups_red){
    #setwd(dir_basic)
    
    system(paste("rm -r"," ", dir_output_ML, sep=""))
    system(paste("mkdir"," ", dir_output_ML, sep=""))
    
    system(paste("rm -r"," ", dir_output_Zern, sep=""))
    system(paste("mkdir"," ", dir_output_Zern, sep=""))
    
    conta_ML <- 0
    for(c in 1:Ncut){

      cutoffBin <- mean(RhoMat) - VetPercCutoff[c]*sd(RhoMat)
      
      df_prova <- RhoMat_Cluster_Melt_Red[RhoMat_Cluster_Melt_Red$x==groups_red[n],c("i","j")]
      hInf <- min(df_prova[,"i"]) - edge
      hSup <- max(df_prova[,"i"]) + edge
      lInf <- min(df_prova[,"j"]) - edge
      lSup <- max(df_prova[,"j"]) + edge
    
      if(hInf < 0){
        hInf=0
      }
      if(hSup < 0){
        hSup=0
      }
      if(lInf < 0){
        lInf=0
      }
      if(lSup < 0){
        lSup=0
      }
      
      if(lSup <= ncol(RhoMat) & hSup <= nrow(RhoMat)){
        
        #Attention: here you can insert RhoMatIm (with a lesser smoothing) or RhoMat with a greater smoothing.
        RhoMatCiclo <- RhoMat #RhoMat2
        RhoMat_spec <- RhoMatCiclo[(hInf):(hSup),(lInf):(lSup)]
        #plot(raster(RhoMat_spec), axes=FALSE, box=FALSE, legend=F,col=rainbowcols, main=groups_red[n])
        
        centralPoin_y <- round(hInf + (hSup - hInf)/2)
        centralPoin_x <- round(lInf + (lSup - lInf)/2)
        centralPoin_y <- nrow(RhoMat) - centralPoin_y
        
        Cluster_Rho_NoNoise <- RhoMatCiclo[hInf:hSup,lInf:lSup]
        Cluster_Rho_NoNoise[Cluster_Rho_NoNoise > cutoffBin] <- max(Cluster_Rho_NoNoise)
        #Cluster_Rho_NoNoise[Cluster_Rho_NoNoise < cutoffBin] <- min(Cluster_Rho_NoNoise)
        
        if(lSup <= ncol(RhoMat) & hSup <= nrow(RhoMat)){
          
          
          #Clustering to remove small islands
          Cluster_Rho_NoNoise_Bin <- round(Cluster_Rho_NoNoise)
          bg <- round(as.numeric(names(table(as.vector(Cluster_Rho_NoNoise_Bin))[table(as.vector(Cluster_Rho_NoNoise_Bin))==max(table(as.vector(Cluster_Rho_NoNoise_Bin)))])))
          Cluster_Rho_NoNoise_Bin[Cluster_Rho_NoNoise_Bin==bg] <- 255
          Cluster_Rho_NoNoise_Bin[Cluster_Rho_NoNoise_Bin!=255] <- 1
          Cluster_Rho_NoNoise_Bin[Cluster_Rho_NoNoise_Bin==255] <- 0
          
          RhoMat_Bin_Melt <- melt(Cluster_Rho_NoNoise_Bin)
          
          sm = RhoMat_Bin_Melt[RhoMat_Bin_Melt[,3]==1,1:2]
          rownames(sm) <- NULL
          colnames(sm) <- c("i","j")
          
          if(dim(sm)[1]>0){
            d = dist(sm, "manhattan")
            gr = cutree(hclust(d, "single"), h = 1)
            RhoMat_Cluster_aus <- sparseMatrix(i = sm[, "i"], j = sm[, "j"], x = gr)
            RhoMat_Cluster_Melt_aus <- summary(RhoMat_Cluster_aus)
            clMax <- as.numeric(names(table(RhoMat_Cluster_Melt_aus$x)[table(RhoMat_Cluster_Melt_aus$x)==max(table(RhoMat_Cluster_Melt_aus$x))]))
            RhoMat_Cluster_Melt_aus_ClMax <- RhoMat_Cluster_Melt_aus[RhoMat_Cluster_Melt_aus$x==clMax[1],] # da rivedere perch?? possono esserci valori multipli
            RhoMat_Cluster_ClMax <- matrix(255, nrow = nrow(Cluster_Rho_NoNoise_Bin), ncol = ncol(Cluster_Rho_NoNoise_Bin))
            for(i in 1:nrow(RhoMat_Cluster_Melt_aus_ClMax)){
              i_aus <- RhoMat_Cluster_Melt_aus_ClMax[i,"i"]
              j_aus <- RhoMat_Cluster_Melt_aus_ClMax[i,"j"]
              val_aus <- Cluster_Rho_NoNoise[i_aus,j_aus]
              RhoMat_Cluster_ClMax[i_aus,j_aus] <- val_aus
            }
          }
          
          #Check that the edges are not cut
          control_1 <- length(unique(RhoMat_Cluster_ClMax[,ncol(RhoMat_Cluster_ClMax)])) == 1
          control_2 <- length(unique(RhoMat_Cluster_ClMax[,1])) == 1
          control_3 <- length(unique(RhoMat_Cluster_ClMax[nrow(RhoMat_Cluster_ClMax),])) == 1
          control_4 <- length(unique(RhoMat_Cluster_ClMax[1,])) == 1
          
            StatCol <- min(which(apply(RhoMat_Cluster_ClMax,2,mean)!=255))
            EndCol <- max(which(apply(RhoMat_Cluster_ClMax,2,mean)!=255))
            StatRow <- min(which(apply(RhoMat_Cluster_ClMax,1,mean)!=255))
            EndRow <- max(which(apply(RhoMat_Cluster_ClMax,1,mean)!=255))
            ImRow <- (EndRow - StatRow)+1
            ImCol <- (EndCol - StatCol)+1
            
            SizeAddMat <- round((SizeSide_IM-ImCol)/2)
            SizeAddMat <- SizeAddMat + NoiseAdd
            
            if(SizeAddMat>0){
              RhoMat_Cluster_ClMax_red <- RhoMat_Cluster_ClMax[,apply(RhoMat_Cluster_ClMax,2,mean)!=255]
              if(SizeAddMat %% 2 == 0){
                AddCol <- matrix(255, ncol = SizeAddMat, nrow = nrow(RhoMat_Cluster_ClMax_red))
                RhoMat_Cluster_ClMax_ausCol <- cbind(AddCol,RhoMat_Cluster_ClMax_red,AddCol)
              }else{
                SizeAddMat <- round(SizeAddMat)
                AddCol1 <- matrix(255, ncol = SizeAddMat, nrow = nrow(RhoMat_Cluster_ClMax_red))
                AddCol2 <- matrix(255, ncol = (SizeAddMat-1), nrow = nrow(RhoMat_Cluster_ClMax_red))
                RhoMat_Cluster_ClMax_ausCol <- cbind(AddCol1,RhoMat_Cluster_ClMax_red,AddCol2)
              }
              RhoMat_Cluster_ClMax_red_red <- RhoMat_Cluster_ClMax_ausCol[apply(RhoMat_Cluster_ClMax_ausCol,1,mean)!=255,]
              
              SizeAddMat <- round((SizeSide_IM-ImRow)/2)
              SizeAddMat <- SizeAddMat + NoiseAdd
              
              if(SizeAddMat>0){
                
                if(SizeAddMat %% 2 == 0){
                  AddRow <- matrix(255, ncol = ncol(RhoMat_Cluster_ClMax_red_red), nrow = SizeAddMat)
                  RhoMat_Cluster_ClMax_ausRow <- rbind(AddRow,RhoMat_Cluster_ClMax_red_red,AddRow)
                }else{
                  SizeAddMat <- round(SizeAddMat)
                  AddRow1 <- matrix(255, nrow = SizeAddMat, ncol = ncol(RhoMat_Cluster_ClMax_ausCol))
                  AddRow2 <- matrix(255, nrow = (SizeAddMat-1), ncol = ncol(RhoMat_Cluster_ClMax_ausCol))
                  RhoMat_Cluster_ClMax_ausRow <- rbind(AddRow1,RhoMat_Cluster_ClMax_red_red,AddRow2)
                }
                
              }
              
              RhoMat_Cluster_ClMax_ausRow <- RhoMat_Cluster_ClMax_ausRow[1:SizeSide_IM, 1:SizeSide_IM]
            }else{
              RhoMat_Cluster_ClMax_ausRow <- matrix(0, ncol = SizeSide_IM, nrow = SizeSide_IM)
              
            } 
        
            #Double the number of pixels in the image to increase the resolution.
            RhoMat_Cluster_ClMax_ausRow_new <- RhoMat_Cluster_ClMax_ausRow
            RhoMat_Cluster_ClMax_ausRow_HR <- matrix(255, ncol = ncol(RhoMat_Cluster_ClMax_ausRow_new)*ing, nrow = nrow(RhoMat_Cluster_ClMax_ausRow_new)*ing)
            RhoMat_Cluster_ClMax_ausRow_control <- matrix(0, ncol = ncol(RhoMat_Cluster_ClMax_ausRow_new)*ing, nrow = nrow(RhoMat_Cluster_ClMax_ausRow_new)*ing)
            for(righe in 1:(nrow(RhoMat_Cluster_ClMax_ausRow_HR)-ing)){
              for(colonne in 1:(ncol(RhoMat_Cluster_ClMax_ausRow_HR)-ing)){
                RhoMat_Cluster_ClMax_ausRow_HR[righe,colonne] <- RhoMat_Cluster_ClMax_ausRow_new[(round(righe/ing)+1),(round(colonne/ing)+1)]
              }
            }
            
            RhoMat_Cluster_ClMax_ausRow <- RhoMat_Cluster_ClMax_ausRow_HR
            #plot(raster(RhoMat_Cluster_ClMax_ausRow_HR), axes=FALSE, box=FALSE, legend=F,col=rainbowcols)
            
         
            #Smoothing.
            RhoMatIm_1 <- RhoMat_Cluster_ClMax_ausRow
            
            for(i in (1+LatoProbe):(nrow(RhoMatIm_1)-LatoProbe)){
              print(i)
              for(j in (1+LatoProbe):(ncol(RhoMatIm_1)-LatoProbe)){
                matAus_1 <-RhoMat_Cluster_ClMax_ausRow[(i-LatoProbe):(i+LatoProbe),(j-LatoProbe):(j+LatoProbe)]
                rho_1 <- mean(matAus_1)
                RhoMatIm_1[i,j] <- rho_1
              }
            }
            
            RhoMat_Cluster_ClMax_ausRow <- RhoMatIm_1
            
            #plot(raster(RhoMat_Cluster_ClMax_ausRow), axes=FALSE, box=FALSE, legend=F,col=rainbowcols, main=paste("file number: ", f, " --- ","cutoff:", VetPercCutoff[c],"--", sep=""))
            #plot(raster(RhoMat_spec), axes=FALSE, box=FALSE, legend=F,col=rainbowcols, main=groups_red[n])
            #plot(raster(Cluster_Rho_NoNoise), axes=FALSE, box=FALSE, legend=F,col=rainbowcols, main=groups_red[n])
            #plot(raster(RhoMat_Cluster_ClMax_ausRow), axes=FALSE, box=FALSE, legend=F,col=rainbowcols, main=paste("file number: ", f, sep=""))
            
            RhoMat_Cluster_ClMax_ausRow <- abs(RhoMat_Cluster_ClMax_ausRow-255)
            #plot(raster(RhoMat_Cluster_ClMax_ausRow), axes=FALSE, box=FALSE, legend=F,col=rainbowcols, main=groups_red[n])
            #plot(raster(Cluster_Rho_NoNoise_Circle), yaxt="n", xaxt="n",col=rainbowcols, main=groups_red[n])
            
            #Delete the image considered twice
            Control_all_black <- sum(RhoMat_Cluster_ClMax_ausRow==0)/(nrow(RhoMat_Cluster_ClMax_ausRow)*ncol(RhoMat_Cluster_ClMax_ausRow))
            Control_all_black <- round(Control_all_black,2)
            if(Control_all_black > 0.4 & Control_all_black < 0.9){
              plot(raster(RhoMat_Cluster_ClMax_ausRow), axes=FALSE, box=FALSE, legend=F,col=rainbowcols, main=paste("file number: ", f, " --- ","cutoff:", VetPercCutoff[c],"--", sep=""))

              conta_ML <- conta_ML + 1
              writePNG(RhoMat_Cluster_ClMax_ausRow/255, target = paste(dir_output_ML,"GPCR_Number_",conta_ML,".png", sep=""))
            }
        }
      }  
    }
    
    #########################################################################
    #########################################################################
    #Minimization of the coefficients relating to C4 symmetry
    system(paste("python3", "ZernikeMoment2D_SingleImage_All_Def.py", sep = " "))
    #system(paste("/Users/edo/anaconda3/bin/python3", "ZernikeMoment2D_SingleImage_All_Def.py", sep = " "))
    
    coord_name <- paste("X_", centralPoin_x, "_Y_", centralPoin_y, sep="")

    #system(paste("mv"," ", "/Users/edo/Documents/Progetti/ImageGPCR/Zernike_Best_C4/Zernike_GPCR.txt", " /Users/edo/Documents/Progetti/ImageGPCR/Zernike_Best_C4/File_", name_file,"_Cluster_", coord_name, ".txt", sep=""))
    system(paste("mv"," ", dir_BestFigC4, "FIG_GPCR.png", " ", dir_BestFigC4, "Fig_", name_file,"_Cluster_", coord_name, ".png", sep=""))
    system(paste("mv"," ", dir_BestZernC4, "Zernike_GPCR.txt", " ", dir_BestZernC4, "File_", name_file,"_Cluster_", coord_name, ".txt", sep=""))
    
    #########################################################################
    #########################################################################
    
  }
  
}#End of loop on files ...
