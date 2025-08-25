setwd("/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust/Campaign_2025")
github_dir <- "/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust/Campaign_2025/"
library(asdreader)
library(reshape2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)

#####################################TEST#######################################
dataset <- c()
file_list <- list.files(paste(github_dir,'/Data/Level0_Calibration/Wavelength_Calibration_2025/UA_ASD_FieldSpec3/',sep=''), full.names = T, pattern = ".asd")
for(i in 1:length(file_list)){
  print(file_list[i])
  md<-get_metadata(file_list[1])
  temp_data<-get_spectra(file_list[1])
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1000)]-data[which(wvl==1001)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)],data[which(wvl==1001):which(wvl==1800)]+c1,data[which(wvl==1801):which(wvl==2500)]+c1+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  type=rep(i,length(wvl))
  id=rep("UA_ASD3",length(wvl))
  
  ds=cbind(wvl,data,type,id)
  dataset <- rbind(dataset,ds)
}
file_list <- list.files(paste(github_dir,'/Data/Level0_Calibration/Wavelength_Calibration_2025/JPL_ASD_FieldSpec4/',sep=''), full.names = T, pattern = ".asd")
for(i in 1:length(file_list)){
  print(file_list[i])
  md<-get_metadata(file_list[1])
  temp_data<-get_spectra(file_list[1])
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1000)]-data[which(wvl==1001)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)],data[which(wvl==1001):which(wvl==1800)]+c1,data[which(wvl==1801):which(wvl==2500)]+c1+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  type=rep(i,length(wvl))
  id=rep("JPL_ASD4",length(wvl))
  
  
  ds=cbind(wvl,data,type,id)
  dataset <- rbind(dataset,ds)
}
file_list <- list.files(paste(github_dir,'/Data/Level0_Calibration/Wavelength_Calibration_2025/USGS_ASD_FieldSpec4/',sep=''), full.names = T, pattern = ".asd")
for(i in 1:length(file_list)){
  print(file_list[i])
  md<-get_metadata(file_list[1])
  temp_data<-get_spectra(file_list[1])
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1000)]-data[which(wvl==1001)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)],data[which(wvl==1001):which(wvl==1800)]+c1,data[which(wvl==1801):which(wvl==2500)]+c1+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  type=rep(i,length(wvl))
  id=rep("USGS_ASD4",length(wvl))
  
  
  ds=cbind(wvl,data,type,id)
  dataset <- rbind(dataset,ds)
}

#Name the columns
colnames(dataset) <- c("wavelength","reflectance","rep","id")
dataset=transform(dataset,wavelength = as.numeric(wavelength))
dataset=transform(dataset,reflectance = as.numeric(reflectance))

#####ASD Reflectance Plots#################
ggplot(dataset, aes(x=wavelength,y=reflectance,group=id,color=id)) +
  geom_line(show.legend = T,linewidth=1,linetype="dashed") +
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(300,2500), breaks = seq(300,2500,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/Figures/Level0_Cal_Figures/Wavelength_Calibration_ALL_ASDs.png',sep=''),dpi=300,width=180,height=120,units='mm')
