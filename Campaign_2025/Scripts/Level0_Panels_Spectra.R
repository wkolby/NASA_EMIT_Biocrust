setwd("/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust/Campaign_2025")
github_dir <- "/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust/Campaign_2025/"
library(asdreader)
library(reshape2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)

#####################################Horse Thief Panels#######################################
master <- read.csv2(paste(github_dir,'Data/Level0_Field_Data/Data_Sheets/Panels/Datasheet_USGS2_Panels_H_05062025.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'Data/Level0_Field_Data/HorseThief/Panels',sep=''), full.names = T, pattern = paste("Panel",master$File[i],".asd",sep=''))
  print(file)
  md<-get_metadata(file)
  temp_data<-get_spectra(file)
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1001)]-data[which(wvl==1000)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)]+c1,data[which(wvl==1001):which(wvl==1800)],data[which(wvl==1801):which(wvl==2500)]+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  site=rep("HT",length(wvl))
  type=rep(master$Panel[i],length(wvl))
  rep=rep(master$Rep[i],length(wvl))
  sample=rep(master$Sample[i],length(wvl))
  
  ds=cbind(site,wvl,data_c,type,sample,rep)
  dataset <- rbind(dataset,ds)
}

#Name the columns
colnames(dataset) <- c("site","wavelength","reflectance","type","sample","rep")
dataset=transform(dataset,wavelength = as.numeric(wavelength))
dataset=transform(dataset,reflectance = as.numeric(reflectance))
dataset=transform(dataset,rep = as.numeric(rep))
dataset=transform(dataset,sample = as.numeric(sample))

#####ASD Reflectance Plots#################
###Species level plots (all samples)###
panels<-c('White','Grey','Black')
for(i in 1:length(panels)){
  ggplot(subset(dataset, type %in% panels[i]), aes(x=wavelength,y=reflectance, group=rep, color=rep)) +
    geom_line(show.legend = T,linewidth=.5,linetype="solid") +
    facet_wrap(~sample)+
    scale_y_continuous("Reflectance") +
    scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,400)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(text = element_text(size = 18))
  ggsave(paste(github_dir,'/Figures/Level0_Panel_Figures/USGS_',panels[i],'_Panel.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% panels[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/Level1_Field_Data/Panels_',panels[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
}

#####################################Horse Thief Tarp#######################################
master <- read.csv2(paste(github_dir,'Data/Level0_Field_Data/Data_Sheets/Panels/Datasheet_NUSO_Tarp_H_05062025.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'Data/Level0_Field_Data/HorseThief/Panels',sep=''), full.names = T, pattern = paste("Panel",master$File[i],".asd",sep=''))
  print(file)
  md<-get_metadata(file)
  temp_data<-get_spectra(file)
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1001)]-data[which(wvl==1000)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)]+c1,data[which(wvl==1001):which(wvl==1800)],data[which(wvl==1801):which(wvl==2500)]+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  site=rep("HT",length(wvl))
  type=rep(master$Panel[i],length(wvl))
  rep=rep(master$Rep[i],length(wvl))
  sample=rep(master$Sample[i],length(wvl))
  
  ds=cbind(site,wvl,data_c,type,sample,rep)
  dataset <- rbind(dataset,ds)
}

#Name the columns
colnames(dataset) <- c("site","wavelength","reflectance","type","sample","rep")
dataset=transform(dataset,wavelength = as.numeric(wavelength))
dataset=transform(dataset,reflectance = as.numeric(reflectance))
dataset=transform(dataset,rep = as.numeric(rep))
dataset=transform(dataset,sample = as.numeric(sample))

#####ASD Reflectance Plots#################
###Species level plots (all samples)###
panels<-c('White','Grey','Black')
for(i in 1:length(panels)){
  ggplot(subset(dataset, type %in% panels[i]), aes(x=wavelength,y=reflectance, group=rep, color=rep)) +
    geom_line(show.legend = T,linewidth=.5,linetype="solid") +
    facet_wrap(~sample)+
    scale_y_continuous("Reflectance") +
    scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,400)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(text = element_text(size = 18))
  ggsave(paste(github_dir,'/Figures/Level0_Panel_Figures/NUSO_',panels[i],'_Tarp.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% panels[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/Level1_Field_Data/Tarp_',panels[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
}

