# For ellipse plot generation
library(GOplot)
# For processing directory names
library(stringr)
library(pracma)
# For storing data in excel file 
library(readxl)
# For determining ellipse formation
library(car)
# For matchpt
library(Biobase)
# For mutate
library(dplyr)
# For write_xlsx
library(writexl)
# For generating quantiles
library(binr)

## K-means clustering analysis
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
#library(HH)

## Density-based Cluster analysis
library(factoextra)
library(fpc)

## SPlit directories for storage
library(DescTools)

# https://uc-r.github.io/kmeans_clustering


# Organize all data quant files with:
# Inside parent directories with separate excel files and following naming
# 1) Stain - Marker Expression to evaluate for vector colocalization and localization (e.g. in organoid)
# 2) Attempt - Replicates
# 3) Condition - Experimental parameter
# Will generate data based on distance input variable (e.g. default = 100)

## File nameing convention is stain_attempt#_cond_#
# ex. CD133_attempt1_A_1, CD133_attempt1_A_2, CD133_attempt2_A_1, 
# CD133_attempt1_B_1, CD133_attempt2_B_1, etc.

## Please have study file have individual conditions with subdirectories 
# for each image for that condition (e.g. replicates)
# these subdirectories will store the imaris datafiles of the same name

organoid_analysis_Partitioned_v3 <- function(stain="CD133", cond=c("A", "B"), eps_v=0.02, prox=2, mpoints=2, quant=5)  {
  
  ### Variables Description ###
  
  # stain - name of stain (e.g. ch3)
  # cond - conditions (e.g. pGFP and pCMV-LonP1-IRES-GFP)
  # dist - distance to define quantiles for binning as measure from periphery of ellipse
  # eps - fraction distance (e.g. total organoid size) for defining the limits of each potential cluster
  # prox - multiplier distance (e.g. total cluster size) from cluster center for correlation and measurements taken for each cluster
  
  # Collect all file directories for each sample 
  
  # Selecting directory where study files are located
  rep_nm <- list()
  cond_nm <- list()
  
  # Store directory name where study files are located
  dir_nm <- choose.dir(default = "", caption = paste("Study Directory for ", cond[x]))
  
  # List subdirectories containing individual statistics representing replicates 
  for (x in 1:length(cond)) {
    cond_nm[x] <- paste(dir_nm, "\\", cond[x], sep="")
    rep_nm[[cond[x]]] <- paste(cond_nm[x], "\\", list.files(cond_nm[[x]], pattern=NULL, all.files=TRUE,
                                    full.names=FALSE, no.. = TRUE), sep="")
  }
  
  
  
  ### Process Statistic files to list_object
  # Store correlations and covariances for those in 'prox' around "Class B" cells 
  rep_cor <- list()
  rep_cov <- list()
  rep_cor_stat <- list()
  rep_cor_dist_stat <- list()
  rep_data_A_avg <- list()
  rep_data_B_avg <- list()
  # Max size of organoid
  prox_cluster_max_tl<-0
  
  for (x in 1:length(cond)) {
    
    # Initialize statistic storage 
    temp_data <- data.frame()
    temp_data_A_avg <- data.frame()
    temp_data_B_avg <- data.frame()
    tmp_data_stat <- list()
    tmp_data_dist_stat <- list()
    d<-1
    
    for (y in 1:length(rep_nm[[x]]))  {
      print(paste("Starting Processing:", x, " ; Replicate:", y, sep=""))
      
      ### Determining files to access for relevant data
      # Gather different statistic file dir names for condition x and replicate y
      rep_file_list<- paste(rep_nm[[cond[x]]][[y]], "\\", list.files(rep_nm[[cond[x]]][[y]], pattern=NULL, all.files=TRUE,
                 full.names=FALSE, no.. = TRUE), sep="")
      
      # Upload position data
      postion_file <- rep_file_list[grepl("Position.csv", rep_file_list)]
      tmp_pos <- read.csv(postion_file, skip=3)
      
      # Upload intensity data for channel 3 (e.g. channel used to determine marker expression)
      int_ch3_file <- rep_file_list[grepl("Intensity_Mean_Ch=3_Img=1.csv", rep_file_list)]
      tmp_int_ch3 <- read.csv(int_ch3_file, skip=3)
      
      # Determining conf_levels
      conf_levels <- dataEllipse(as.vector(unlist(tmp_pos[,1][tmp_pos["Set.1"]=="Class A"])), as.vector(unlist(tmp_pos[,2][tmp_pos["Set.1"]=="Class A"])),
                                 levels=0.95, ellipse.label=0.95, lty=2, fill=TRUE, fill.alpha=0.1, robust=TRUE)
      
      # Determine the diameter of the ellipse for ascertaining cluster size based on eps_v
      org_diameter <- sqrt((max(conf_levels[,1]) - min(conf_levels[,1]))^2 + (max(conf_levels[,2]) - min(conf_levels[,2]))^2)
      
      # Store maximum size organoid
      if (prox_cluster_max_tl<org_diameter) {
        prox_cluster_max_tl<-org_diameter
      }
      
      ### Processing position data relative to ellipse perimeter
      # Determine ellipse center position
      center_x <- (max(conf_levels[,1]) - min(conf_levels[,1]))/2 + min(conf_levels[,1])
      center_y <- (max(conf_levels[,2]) - min(conf_levels[,2]))/2 + min(conf_levels[,2])
      
      tmp_x <- as.vector(unlist(tmp_pos[,1]))
      tmp_y <- as.vector(unlist(tmp_pos[,2]))
      
      # For each id determine distance from nearest point on ellipse perimeter (e.g. distance from organoid edge)
      coordinates <- data.frame(x=tmp_x, y=tmp_y, id=as.vector(unlist(tmp_pos[,9])))
      
      index_perim<-matchpt(as.matrix(coordinates[,1:2]), conf_levels)
      
      coordinates <- cbind(coordinates, index_perim[,2])
      
      ### Incorporating intensity values for each cell
      # Add corresponding intensities values and class-based cutoff values (e.g. threshold for GFP expression)
      coordinates <- cbind(coordinates, tmp_int_ch3["Intensity.Mean"], tmp_int_ch3["Set.1"])
      
      # Create empty column in coordinates dataframe for input of quantile in terms of distance from edge of ellipse
      coordinates["dist_quantile"]<-NA
      
      ### Determine distance from perimeter for each cell
      # y determines number of quantiles based on max distance from periphery divided by desired dist (e.g for each quantile)
      for (z in 1:ceiling(max(coordinates[,4])/dist))  {
        # Determine values that fall into specific quantiles and assign specific 'dist_quantile' values
        coordinates[coordinates[,4]<(z*dist) & coordinates[,4]>((z-1)*dist),7]<-z
      }
      
      ### Store values for replicates into list object
      list_object[[cond[[x]]]][[y]] <- coordinates
      
      ### Determine proximity correlation of intensity around "Class B"
      # Assumes "Class B" represents the GFP labeled cells and maps expression relative to proxmity of 'prox'
      # Create new masked dataframe based on "Class B" denomination
      B_coordinates <- coordinates[coordinates[,6]=="Class B",]
      A_coordinates <- coordinates[coordinates[,6]=="Class A",]
      
      ### Loop through all cells in B_coordinates
      # Store in dataframe with rbind
      stor_dataframe <- data.frame()
      stor_dataframe_b <- data.frame()
      
      ######################## TEST CODE #################################
      
      ### Generate DBSCAN cluster ellipse perimeter and save plot
      print("Generating Plot...")
      db <- fpc::dbscan(B_coordinates[,1:2], eps = org_diameter*eps_v, MinPts = mpoints)
      
      heatcols <- heat.colors(max(db[["cluster"]]))
      
      jpeg(file=paste(dir_nm, "\\cond_", cond[x], "_rep_", y, "_clust", ".jpeg", sep=""))
      
      if (max(db[["cluster"]])!=0) {
        plot(fviz_cluster(db, data = B_coordinates[,1:2], stand = FALSE,
                          ellipse = TRUE, show.clust.cent = TRUE,
                          geom = "point",palette = heatcols, ggtheme = theme_classic()))
      }
      
      dev.off()
      print("DBSCAN cluster Plot Generated. Processing Data...")
      
      ### Generating Proximity Plot based on clusters
      # Determine center positions of proximity plot
      
      if (max(db[["cluster"]])!=0)  {
        centers <- list()
        prox_cluster_max<-0
        for (z in 1:max(db[["cluster"]]))  {
          tmp <- B_coordinates[,1:2][db[["cluster"]]==z,]
          
          prox_cluster <- sqrt((max(tmp[,1]) - min(tmp[,1]))^2 + (max(tmp[,2]) - min(tmp[,2]))^2)
          
          center_x <- mean(tmp[,1])
          center_y <- mean(tmp[,2])
          
          centers[[z]]<-c(center_x, center_y, prox_cluster)
          
          # Store max cluster size
          if (prox_cluster>prox_cluster_max) {
            prox_cluster_max <- prox_cluster
          }
          
        }
        
        # Determine intensities around proximity plot 
        stor_dataframe <- data.frame()
        stor_list <- list()
        quant_prox_int <- list()
        
        # Initialize quantile list
        for (t in 1:quant)  {
          quant_prox_int[[t]]<-1
        }
        
        
        for (z in 1:length(centers))  {
          
          # Create mask and create stor value with proximity distance matched with intensity value
          prox_val <- sqrt((A_coordinates["x"] - centers[[z]][1])^2 + (A_coordinates["y"] - centers[[z]][2])^2)
          prox_val_b <- sqrt((B_coordinates["x"] - centers[[z]][1])^2 + (B_coordinates["y"] - centers[[z]][2])^2)
          
          prox_mask <- (prox*centers[[z]][3])>prox_val & prox_val>0
          prox_mask_b <- (centers[[z]][3])>prox_val_b
            
          stor <- as.data.frame(cbind(prox_val[prox_mask], A_coordinates[prox_mask,5]))
          stor_b <- as.data.frame(cbind(prox_val_b[prox_mask_b], B_coordinates[prox_mask_b,5]))
          
          # Storage as list
          stor_list[[z]] <- stor
          
          # Appending storage format
          
          # Store values in temporary dataframe
          stor_dataframe <- rbind(stor_dataframe, stor)
          stor_dataframe_b <- rbind(stor_dataframe_b, stor_b)
          
          #org_diameter*eps_v
          # Assign quantiles
          tmp <- as.data.frame(data.matrix(stor_dataframe)) %>% mutate(quantile = ntile(stor_dataframe[,1], quant))
          
          
          # Determine the intensity statistics with respect to quantile distances of 'quant' number from cluster centers
          for (l in 1:quant)  {
            quant_prox_int[[l]][[z]] <- mean(tmp[tmp["quantile"]==l,2])
          }
        }
        
        # Generating Proximity Plot 
        print("Generating Proximity Plot...")

        if (!isempty(unlist(stor_dataframe))) {
 
          jpeg(file=paste(dir_nm, "\\cond_", cond[x], "_rep_", y, "_proximity_", prox, ".jpeg", sep=""))
          plot(stor_dataframe, ylim = c(1,200), xlim = c(1,prox_cluster_max*prox))
          points(axes=FALSE, stor_dataframe_b[,1], stor_dataframe_b[,2], xlab = "", ylab = "", col = "green", pch = 2)
          
          dev.off()
          
          # Generating Proximity-Quantile Plot for specific organoid
          # Determine statistics of quantiles 
          tmp_stat_vec <- data.frame(Mean=NA, Std=NA)
          
          for (t in 1:quant)  {
            tmp_stat_vec[t,"Mean"] <- mean(na.omit(quant_prox_int[[t]]))
            tmp_stat_vec[t,"Std"] <- std(na.omit(quant_prox_int[[t]]))
          }
          
          print("Generating Proximity Quantile Plot...")
          'if (x==2 & y==2)  {
            print("Here")
          }'
          jpeg(file=paste(dir_nm, "\\cond_", cond[x], "_rep_", y, "_proximity-dist-quantile__", prox, ".jpeg", sep=""))
          barplot(unlist(tmp_stat_vec["Mean"]), main="Cluster Quantile Intensity Distance Plot",
                  xlab="Quantiles 1-5", ylim=c(min(tmp_stat_vec["Mean"])-5, max(tmp_stat_vec["Mean"])+5))
          
          barCenters <- barplot(height = unlist(tmp_stat_vec["Mean"]),
                                beside = true, las = 2,
                                ylim = c(min(tmp_stat_vec["Mean"])*0.8, max(tmp_stat_vec["Mean"])*1.2),
                                cex.names = 0.75, xaxt = "n",
                                main = "test",
                                ylab = "test2",
                                border = "black", axes = TRUE)
          
          segments(barCenters, unlist(tmp_stat_vec["Mean"]) - unlist(tmp_stat_vec["Std"]) * 2, barCenters,
                   unlist(tmp_stat_vec["Mean"]) + unlist(tmp_stat_vec["Std"]) * 2, lwd = 1.5)
          
          dev.off()
          
        }
        
        ### Stor correlation and covariance
        
        temp_data[y ,"cor"] <- cor(stor_dataframe)[1,2]
        temp_data[y, "cov"] <- cov(stor_dataframe)[1,2]
        temp_data_A_avg[y ,"mean"] <- mean(stor_dataframe[,2])
        temp_data_B_avg[y ,"mean"] <- mean(stor_dataframe_b[,2])
        
        # For all images process correlation
        for (z in 1:length(stor_list))  {
          tmp_data_stat[[d]] <- cor(stor_list[[z]])[1,2]
          tmp_data_dist_stat[[d]] <- matchpt(t(as.matrix(centers[[z]][1:2])), conf_levels)[2]
          d<-d+1
        }
        
        # Reset stor_dataframe for next replicate
        stor_dataframe <- data.frame()
      }
      
      
      ######################## TEST CODE #################################
      
      ### Generate ellipse perimeter and save plot
      print("Generating Plot...")
      jpeg(file=paste(dir_nm, "\\cond_", cond[x], "_rep_", y, ".jpeg", sep=""))
      conf_levels <- dataEllipse(as.vector(unlist(tmp_pos[,1][tmp_pos["Set.1"]=="Class A"])), as.vector(unlist(tmp_pos[,2][tmp_pos["Set.1"]=="Class A"])),
                                 levels=0.95, ellipse.label=0.95, lty=2, fill=TRUE, fill.alpha=0.1, robust=TRUE)
      points(axes=FALSE, tmp_pos[,1][tmp_pos["Set.1"]=="Class B"], tmp_pos[,2][tmp_pos["Set.1"]=="Class B"], xlab = "", ylab = "", col = "green", pch = 2)
      #points(axes=FALSE, test_clust[["centers"]], xlab = "", ylab = "", col = "red", pch = 2)
      
      dev.off()
      print("Plot Generated. Processing Data...")
      
    }
    
    # Store in respective condition
    rep_cor[[cond[x]]] <- as.data.frame(temp_data)
    rep_cor_stat[[cond[x]]] <- tmp_data_stat
    rep_cor_dist_stat[[cond[x]]] <- tmp_data_dist_stat
    rep_data_A_avg[[cond[x]]] <- temp_data_A_avg
    rep_data_B_avg[[cond[x]]] <- temp_data_B_avg
    # Reinitialize as empty for next condition
    temp_data <- data.frame()
    
  }
  
  ### Plots for correlation for each cluster by distance from organoid periphery 
  # For condition A
  jpeg(file=paste(dir_nm, "\\cond_A_correlation_by_distance.jpeg", sep=""))
  plot(main="Cond: A; Correlations of Clusters per Distance from Organoid Border", 
       xlab="Distance (uM)", ylab="Correlation", ylim=c(-1, 1), xlim=c(1, prox_cluster_max_tl/2), 
       na.omit(t(rbind(unlist(rep_cor_dist_stat[["A"]]), unlist(rep_cor_stat[["A"]])))))
  
  # Data matrix from plot
  tmp_matrix <- as.data.frame(na.omit(t(rbind(unlist(rep_cor_dist_stat[["A"]]), unlist(rep_cor_stat[["A"]])))))
  
  #Save summary of the linear model
  sum_lm<-lm(tmp_matrix[,2]~tmp_matrix[,1])
  
  #Get coefficients
  coef_lm<-sum_lm$coefficients
  
  #Set the abline as the coefficients obtained by your linear model
  abline(a = coef_lm[1], 
         b = coef_lm[2], 
         col = "lightblue",
         lwd = 2)
  
  dev.off()
  
  # For condition B
  jpeg(file=paste(dir_nm, "\\cond_B_correlation_by_distance.jpeg", sep=""))
  plot(main="Cond: B; Correlations of Clusters per Distance from Organoid Border", 
       xlab="Distance (uM)", ylab="Correlation", ylim=c(-1, 1), xlim=c(1, prox_cluster_max_tl/2), 
       na.omit(t(rbind(unlist(rep_cor_dist_stat[["B"]]), unlist(rep_cor_stat[["B"]])))))
  
  # Data matrix from plot
  tmp_matrix <- as.data.frame(na.omit(t(rbind(unlist(rep_cor_dist_stat[["B"]]), unlist(rep_cor_stat[["B"]])))))
  
  #Save summary of the linear model
  sum_lm<-lm(tmp_matrix[,2]~tmp_matrix[,1])
  
  #Get coefficients
  coef_lm<-sum_lm$coefficients
  
  #Set the abline as the coefficients obtained by your linear model
  abline(a = coef_lm[1], 
         b = coef_lm[2], 
         col = "lightblue",
         lwd = 2)
  dev.off()
  
  ### Key statistics that are saved in new datafile in directory
  data_summary<-list()
  attr(data_summary, "exp")<-paste("Experiment:", SplitPath(dir_nm)[["filename"]], sep="")
  # For condition A (e.g. LonP1 1:100): mean and std of all correlations for each cluster
  data_summary[["Cond A: Mean Correlations"]]<-mean(rep_cor[[cond[1]]][,1])
  data_summary[["Cond A: Std Correlations"]]<-std(rep_cor[[cond[1]]][,1])
  data_summary[["Cond B: Mean Correlations"]]<-mean(rep_cor[[cond[2]]][,1])
  data_summary[["Cond B: Std Correlations"]]<-std(rep_cor[[cond[2]]][,1])
  
  # Correlations for correlations based on distance
  data_summary[["Cond A: Distance Correlation"]]<-cor(na.omit(t(rbind(unlist(rep_cor_dist_stat[["A"]]), unlist(rep_cor_stat[["A"]])))))[1,2]
  
  data_summary[["Cond B: Distance Correlation"]]<-cor(na.omit(t(rbind(unlist(rep_cor_dist_stat[["B"]]), unlist(rep_cor_stat[["B"]])))))[1,2]
  
  # For condition A: Class A (e.g. non-GFP) - ch3 average and sd intensity in clusters
  data_summary[["Cond A: Class A: Mean Intensity"]]<-mean(unlist(rep_data_A_avg[["A"]]))
  data_summary[["Cond A: Class A: Std Intensity"]]<-std(unlist(rep_data_A_avg[["A"]]))
  
  # For condition A: Class B (e.g. GFP) - ch3 average and sd intensity in clusters
  data_summary[["Cond A: Class B: Mean Intensity"]]<-mean(unlist(rep_data_A_avg[["B"]]))
  data_summary[["Cond A: Class B: Std Intensity"]]<-std(unlist(rep_data_A_avg[["B"]]))
  
  # For condition B: Class A (e.g. non-GFP) - ch3 average and sd intensity in clusters
  data_summary[["Cond B: Class A: Mean Intensity"]]<-mean(unlist(rep_data_B_avg[["A"]]))
  data_summary[["Cond B: Class A: Std Intensity"]]<-std(unlist(rep_data_B_avg[["A"]]))
  
  # For condition B: Class B (e.g. GFP) - ch3 average and sd intensity in clusters
  data_summary[["Cond B: Class B: Mean Intensity"]]<-mean(unlist(rep_data_B_avg[["B"]]))
  data_summary[["Cond B: Class B: Std Intensity"]]<-std(unlist(rep_data_B_avg[["B"]]))
  
  ### Store processed data for data_summary table  file in study directory
  write_xlsx(as.data.frame(data_summary),
             paste(dir_nm, "\\", "Data_Summary_Clusters.xlsx", sep="" ))
  
  ### Sample code for loading data
  #   test<-read_xlsx(paste(dir_nm, "\\", "Data_Summary_Clusters.xlsx", sep="" ))

  ## Generate ellipse for every object to determine borders and reassign displacement
  ## as distance from 1 (e.g. edge) to 0 (e.g. center) and then quantile

}

