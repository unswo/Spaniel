#library(imager)
#library(sp)
#library(magick)
#library(tidyverse)

#library(doParallel)

######## TO DO!
### convert all of the images to black and white
### create a for loop to add annotations
### add annotations to sce object and and integrated seurat
### find markers for:
# annotate heart sections as “base” and “ventricles” ​
# gene-by-gene comparison pooled spots of XR v CTRL “bases” ​
# gene-by-gene comparison pooled spots of XR and CTRL “ventricles”​
# identify up/down-regulated genes in “base” only
# create find domain function for Spaniel

findDomain <-  function(imgFile, cln = 3, fll = 12){
  ## domain coordinates
  im <- load.image(imgFile)
  im <- grayscale(im)
  #plot(im)
  
  im <- clean(im,cln) %>% imager::fill(fll)
  #plot(im)
  plot.new()
  hl <- highlight(im)
  
  
  coords <- lapply(hl, function(coord){data.frame(coord$x, coord$y)})
  p1s <- lapply(coords, Polygon)  %>% Polygons(3)
  sps <- SpatialPolygons(list(p1s))
  return(sps)}


findDomainsAnno <- function(anno_sample, domain_list, sce, annotationPath){
  domain_name <- anno_sample$domain
  domain_list[[domain_name]] <- findDomain(file.path(annotationPath,anno_sample$jpg), 
                                           fll = 7, cln = 3)
  return(domain_list)
}


sceDomainAnno <- function(sce, domain_name, domain_list, points){
  colData(sce)[,paste0("domain_", domain_name)] <- over(points, domain_list[[domain_name]]) %>% 
    is.na() %>% 
    ifelse(yes = "", no = domain_name)
  return(sce)
}


annotateFig <- function(sce, sampleInformation = sampleInfo, scaleFactorsPath = scaleFacPath, annotationPath = annoPath){
  sampleInformation <- sampleInformation[sampleInformation$id == sce$id[1],]
  scaleFactorsPath <- file.path(scaleFactorsPath,
                                sampleInformation$Fname,
                                "spatial/scalefactors_json.json")
  annotations <- list.files(annotationPath)
  annotations_df <- data.frame(id_short = gsub("_[^_]{3}.*.jpg", "", annotations) ,
                               domain = annotations %>% gsub("\\.jpg", "", .), 
                               jpg = annotations)
  annotations_df <- annotations_df %>% left_join(sampleInfo)
  annotations_df <- annotations_df[annotations_df$id_short == sce$id[1],]
  scaleFactors <- jsonlite::fromJSON(txt = file.path(scaleFactorsPath))
  spots <- data.frame(Image_Y = sce$Image_Y, Image_X = sce$Image_X)
  spots$pixel_x <- spots$Image_X * scaleFactors$tissue_hires_scalef
  spots$pixel_y <- spots$Image_Y * scaleFactors$tissue_hires_scalef
  points <- spots %>% select(pixel_x, pixel_y) %>% as.matrix() %>%  SpatialPoints
  
  domain_list <- list()
  #id_short <- sampleInformation$id_short
  annotations_sample <- annotations_df#[annotations_df$id_short == id_short,]

  for ( i in 1:length(annotations_sample$id_short)){
    cat(i, "\n")
    domain_list <- findDomainsAnno(annotations_sample[i,], domain_list,  sce, annotationPath)
    sce <- sceDomainAnno(sce, annotations_sample$domain[i], domain_list, points)
  }
  # collapse domain lists
  domains <- colData(sce) %>% 
    data.frame() %>% 
    select(grep(paste0("domain_", sampleInformation$id_short, "_"),colnames(colData(sce)))) %>% 
    apply(1 , paste, collapse = "")
  ## remove overlapping domains
  domains[!domains %in% annotations_sample$domain] <- NA
  domains <- gsub(paste0("^", sampleInformation$id_short, "_"), "", domains)
  
  sce$domain <- domains
  names(sce) <- sce$id[1]
  return(sce)
}


spanielAnnotation <- function(sceList, sampleInformation = sampleInfo, scaleFactorsPath = scaleFacPath, annotationPath = annoPath, par = F){
  annotated_sce <- list()
  if (par !=F){
    registerDoParallel(par)
    annotated_sce <- foreach(i = 1:length(sceList), .export = c(ls(.GlobalEnv), "dplyr", "Spaniel")) %dopar% annotateFig(sceList[[i]])
    stopImplicitCluster()
      } else{
    for (j in 1:length(sceList)){
      annotated_sce[[j]] <- annotateFig(sceList[[j]])
    }
  }
  names(annotated_sce) <- names(sceList)
  system("mkdir -p images/")
  
  for (i in 1:length(annotated_sce)){
    imgOut <- paste0("images/", annotated_sce[[i]]$id[1], "_", annotated_sce[[i]]$tissue[1], "_annotations.tiff")
    p1 <- spanielPlot(object = annotated_sce[[i]],
                      plotType = "Cluster",
                      clusterRes = "domain",
                      showFilter = NULL,
                      techType = "Visium",
                      ptSizeMax = 1,
                      customTitle = "Domains")
    p1 %>% plot()
    ggsave(imgOut, width = 10, height = 8)
  }
  return(annotated_sce)
}

# ## read data
# seurat_list <- readRDS("/data/rachel/Spartan/Spartan_008/rObjects/seurat_list.rds")
# sObj_irradiated <- readRDS("/data/rachel/Spartan/Spartan_008/rObjects/irradiated_integrated.rds")
# 
# ### spot coordinates
# 
# sampleInfo <- readRDS("/data/rachel/Spartan/Spartan_008/rObjects/sampleInfo.rds")
# sampleInfo$id_short <- gsub("_.*", "", sampleInfo$id)
# sce_list <- readRDS("/data/rachel/Spartan/Spartan_008/rObjects/sce_list.rds")
# 
# sampleInformation <- sampleInfo
# scaleFacPath <- "/data/rachel/Spartan/Spartan_008/Results_nostromo/"
# annoPath <- "/data/rachel/Spartan/Spartan_008/annotations"
# 
# #system.time(anno_sce <- spanielAnnotation(sce_list))
# # user  system elapsed
# # 43.476   6.865  50.214
# 
# system.time(anno_sce <- spanielAnnotation(sce_list, par = 4))
# user  system elapsed
# 21.539  10.667  17.642
# 
# 
# 
# All_integrated <- readRDS("/data/rachel/Spartan/Spartan_008/rObjects/All_integrated.rds")
# 
# All_integrated$domain <- NA
# 
# for (id in names(annotated_sce)){
#   All_integrated$domain[All_integrated$id == id] <- annotated_sce[[id]]$domain
#   
# }
# 
# ## correct typos
# All_integrated$domain[All_integrated$domain == "left-venrticle"] <- "left-ventricle"
# sce_list[[3]]$domain[sce_list[[3]]$domain == "left-venrticle"] <- "left-ventricle"
# # saveRDS(All_integrated, "rObjects/All_integrated.rds")
# # saveRDS(sce_list, "rObjects/sce_list_annotated.rds")
# 
# DimPlot(All_integrated, group.by = "domain", label = TRUE) #+ ggsave("/data/rachel/Spartan/Spartan_008/images/umap_by_domain.tiff")
# 
# 
# All_integrated[[]] %>% select(c(id, condition, domain, seurat_clusters)) %>% 
#   group_by(domain, condition, seurat_clusters) %>% 
#   summarise(no_spots = n()) %>% ggplot(aes(x=seurat_clusters,no_spots, fill = condition )) + 
#   geom_bar(stat = "identity") + 
#   facet_wrap(~domain, scales = "free_y") #+ ggsave("/data/rachel/Spartan/Spartan_008/images/numbers_clusters_domains.png")




