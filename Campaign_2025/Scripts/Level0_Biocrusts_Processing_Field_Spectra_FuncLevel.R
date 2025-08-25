setwd("/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust/Campaign_2025")
github_dir <- "/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust/Campaign_2025/"
library(asdreader)
library(reshape2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)

###FUNCTIONS##################################################################################
#Calculates vegetation simple ratio (a:b) / (c:d)
calc_SR <- function(df, a, b, c, d){
  b1_1 <- which.min(abs(a - as.numeric(names(df), options(warn=-1))))   #Identify first band and start range
  b1_2 <- which.min(abs(b - as.numeric(names(df), options(warn=-1))))   #Identify first band and end range
  b2_1 <- which.min(abs(c - as.numeric(names(df), options(warn=-1))))   #Identify second band and start range
  b2_2 <- which.min(abs(d - as.numeric(names(df), options(warn=-1))))   #Identify second band and end range
  (rowMeans(df_wide[b1_1:b1_2])/(rowMeans(df_wide[b2_1:b2_2])))
}
#Calculates normalized index using (a:b - c:d) / (a:b + c:d)
calc_VI <- function(df, a, b, c, d){
  b1_1 <- which.min(abs(a - as.numeric(names(df), options(warn=-1))))   #Identify first band and start range
  b1_2 <- which.min(abs(b - as.numeric(names(df), options(warn=-1))))   #Identify first band and end range
  b2_1 <- which.min(abs(c - as.numeric(names(df), options(warn=-1))))   #Identify second band and start range
  b2_2 <- which.min(abs(d - as.numeric(names(df), options(warn=-1))))   #Identify second band and end range
  ((rowMeans(df_wide[b1_1:b1_2]) - rowMeans(df_wide[b2_1:b2_2])) /      #vegetation index equation
      (rowMeans(df_wide[b1_1:b1_2]) + rowMeans(df_wide[b2_1:b2_2])))
}

#####################################Horse Thief#######################################
master <- read.csv2(paste(github_dir,'Data/Level0_Field_Data/Data_Sheets/Biocrust/Datasheet_Biocrust_H_05062025.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'/Data/Level0_Field_Data/HorseThief/Biocrust/',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
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
  type=rep(master$Biocrust[i],length(wvl))
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
species<-c('BGR','LCY','DCY','LCNy','COTE','MOSS')
for(i in 1:length(species)){
  ggplot(subset(dataset, type %in% species[i]), aes(x=wavelength,y=reflectance, group=rep, color=rep)) +
    geom_line(show.legend = T,linewidth=.5,linetype="solid") +
    facet_wrap(~sample)+
    scale_y_continuous("Reflectance") +
    scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,400)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(text = element_text(size = 18))
  ggsave(paste(github_dir,'/Figures/Level0_Field_Figures/HorseThief_',species[i],'_Biocrust.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% species[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/Level1_Field_Data/HorseThief_',species[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
}

###Species level plots (mean +/- sd)###
#soil biocrust
biocrust<-subset(dataset, type %in% c('BGR','LCY','DCY','LCNy','COTE','MOSS'))
biocrust_mean_std <- biocrust %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

biocrust_mean_std$type <- factor(biocrust_mean_std$type,levels=c('BGR','LCY','DCY','LCNy','COTE','MOSS')) #Set the order for manual fill
ggplot(biocrust_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  scale_color_manual(values=c('sandybrown','red2','red4','yellow3','blue2','green4'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  scale_fill_manual(values=c('sandybrown','red2','red4','yellow3','blue2','green4'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(350,2500), breaks = seq(300,2500,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/Figures/Level1_Field_Figures/HorseThief_Full_Soil_Biocrusts.png',sep=''),dpi=300,width=180,height=120,units='mm')

#####################################Salt Valley#######################################
master <- read.csv2(paste(github_dir,'Data/Level0_Field_Data/Data_Sheets/Biocrust/Datasheet_Biocrust_S_05062025.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'/Data/Level0_Field_Data/SaltValley/Biocrust/',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
  print(file)
  md<-get_metadata(file)
  temp_data<-get_spectra(file)
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1001)]-data[which(wvl==1000)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)]+c1,data[which(wvl==1001):which(wvl==1800)],data[which(wvl==1801):which(wvl==2500)]+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  site=rep("SV",length(wvl))
  type=rep(master$Biocrust[i],length(wvl))
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
species<-c('BGR','LCY','MOSS','LCNy','LCNp')
for(i in 1:length(species)){
  ggplot(subset(dataset, type %in% species[i]), aes(x=wavelength,y=reflectance, group=rep, color=rep)) +
    geom_line(show.legend = T,linewidth=.5,linetype="solid") +
    facet_wrap(~sample)+
    scale_y_continuous("Reflectance") +
    scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,400)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(text = element_text(size = 18))
  ggsave(paste(github_dir,'/Figures/Level0_Field_Figures/SaltValley_',species[i],'_Biocrust.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% species[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/Level1_Field_Data/SaltValley_',species[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
}

###Species level plots (mean +/- sd)###
#soil biocrust
biocrust<-subset(dataset, type %in% c('BGR','LCY','MOSS','LCNy','LCNp'))
biocrust_mean_std <- biocrust %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

biocrust_mean_std$type <- factor(biocrust_mean_std$type,levels=c('BGR','LCY','MOSS','LCNy','LCNp')) #Set the order for manual fill
ggplot(biocrust_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  scale_color_manual(values=c('sandybrown','red2','blue2','yellow4','pink2'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  scale_fill_manual(values=c('sandybrown','red2','blue2','yellow4','pink2'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(350,2500), breaks = seq(300,2500,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/Figures/Level1_Field_Figures/SaltValley_Full_Soil_Biocrusts.png',sep=''),dpi=300,width=180,height=120,units='mm')

#####################################Panorama#######################################
master <- read.csv2(paste(github_dir,'Data/Level0_Field_Data/Data_Sheets/Biocrust/Datasheet_Biocrust_P_05062025.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'/Data/Level0_Field_Data/Panorama/Biocrust/',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
  print(file)
  md<-get_metadata(file)
  temp_data<-get_spectra(file)
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1001)]-data[which(wvl==1000)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)]+c1,data[which(wvl==1001):which(wvl==1800)],data[which(wvl==1801):which(wvl==2500)]+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  site=rep("P",length(wvl))
  type=rep(master$Biocrust[i],length(wvl))
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
species<-c('BGR','LCY','MOSS')
for(i in 1:length(species)){
  ggplot(subset(dataset, type %in% species[i]), aes(x=wavelength,y=reflectance, group=rep, color=rep)) +
    geom_line(show.legend = T,linewidth=.5,linetype="solid") +
    facet_wrap(~sample)+
    scale_y_continuous("Reflectance") +
    scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,400)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(text = element_text(size = 18))
  ggsave(paste(github_dir,'/Figures/Level0_Field_Figures/Panorama_',species[i],'_Biocrust.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% species[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/Level1_Field_Data/Panorama_',species[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
}

###Species level plots (mean +/- sd)###
#soil biocrust
biocrust<-subset(dataset, type %in% c('BGR','LCY','MOSS','LCNy','LCNp'))
biocrust_mean_std <- biocrust %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

biocrust_mean_std$type <- factor(biocrust_mean_std$type,levels=c('BGR','LCY','MOSS')) #Set the order for manual fill
ggplot(biocrust_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  scale_color_manual(values=c('sandybrown','red2','blue2'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  scale_fill_manual(values=c('sandybrown','red2','blue2'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(350,2500), breaks = seq(300,2500,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/Figures/Level1_Field_Figures/Panorama_Full_Soil_Biocrusts.png',sep=''),dpi=300,width=180,height=120,units='mm')

