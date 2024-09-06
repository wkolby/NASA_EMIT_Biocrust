setwd("/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust")
github_dir <- "/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust/"
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

#####################################TEST#######################################
master <- read.csv2(paste(github_dir,'Data/NASA_EMIT_Campaign_062024/Level0/Data_Sheets/TEST/Datasheet_Test2_06052024.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level0/TEST2',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
  print(file)
  md<-get_metadata(file)
  temp_data<-get_spectra(file)
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1000)]-data[which(wvl==1001)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)],data[which(wvl==1001):which(wvl==1800)]+c1,data[which(wvl==1801):which(wvl==2500)]+c1+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  type=rep(master$Biocrust[i],length(wvl))
  
  ds=cbind(wvl,data_c,type)
  dataset <- rbind(dataset,ds)
}

#Name the columns
colnames(dataset) <- c("wavelength","reflectance","type")
dataset=transform(dataset,wavelength = as.numeric(wavelength))
dataset=transform(dataset,reflectance = as.numeric(reflectance))

#####ASD Reflectance Plots#################
###Isolate soil signature
soil<-subset(dataset, type %in% c('LCY','SOIL'))
soil_mean_std <- soil %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

ggplot(soil_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  scale_color_manual(values=c('red3','green3','cyan2','purple'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  scale_fill_manual(values=c('red3','green3','cyan2','purple'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/TEST2_Full_Soil_LCY.png',sep=''),dpi=300,width=180,height=120,units='mm')

###Isolate calibration panel rep1 signature
cal<-subset(dataset, type %in% c('CAL1','CAL2','CAL3'))
cal_mean_std <- cal %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

ggplot(cal_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  scale_color_manual(values=c('lightgrey','darkgrey','black'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .5) +
  scale_fill_manual(values=c('lightgrey','darkgrey','black'))+
  scale_y_continuous("Reflectance",limits=c(0,1)) +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/TEST2_Full_Calibration.png',sep=''),dpi=300,width=180,height=120,units='mm')

###Isolate calibration panel rep2 signature
cal<-subset(dataset, type %in% c('CAL1_REP2','CAL2_REP2','CAL3_REP2'))
cal_mean_std <- cal %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

ggplot(cal_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  scale_color_manual(values=c('lightgrey','darkgrey','black'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .5) +
  scale_fill_manual(values=c('lightgrey','darkgrey','black'))+
  scale_y_continuous("Reflectance",limits=c(0,1)) +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/TEST2_Full_Calibration_REP2.png',sep=''),dpi=300,width=180,height=120,units='mm')

#####################################SandFlats#######################################
master <- read.csv2(paste(github_dir,'Data/NASA_EMIT_Campaign_062024/Level0/Data_Sheets/Biocrust/Datasheet_Biocrust_SF_06102024.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level0/SF_soil',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
  print(file)
  md<-get_metadata(file)
  temp_data<-get_spectra(file)
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1000)]-data[which(wvl==1001)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)],data[which(wvl==1001):which(wvl==1800)]+c1,data[which(wvl==1801):which(wvl==2500)]+c1+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  type=rep(master$Biocrust[i],length(wvl))
  rep=rep(master$Rep[i],length(wvl))
  sample=rep(master$Sample[i],length(wvl))
  
  ds=cbind(wvl,data_c,type,sample,rep)
  dataset <- rbind(dataset,ds)
}

#Name the columns
colnames(dataset) <- c("wavelength","reflectance","type","sample","rep")
dataset=transform(dataset,wavelength = as.numeric(wavelength))
dataset=transform(dataset,reflectance = as.numeric(reflectance))
dataset=transform(dataset,rep = as.numeric(rep))
dataset=transform(dataset,sample = as.numeric(sample))

#####ASD Reflectance Plots#################
###Species level plots (all samples)###
species<-c('BGR','LCY','DCY','COTE','SQUA','Rock_MOSS','Rock_LCHk','Rock_LCHg','Rock_LCHy','Rock_LCHb')
for(i in 1:length(species)){
  ggplot(subset(dataset, type %in% species[i]), aes(x=wavelength,y=reflectance, group=rep, color=rep)) +
    geom_line(show.legend = T,linewidth=.5,linetype="solid") +
    facet_wrap(~sample)+
    scale_y_continuous("Reflectance") +
    scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,400)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(text = element_text(size = 18))
  ggsave(paste(github_dir,'/figures/Species/SF_',species[i],'_Biocrust.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% species[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level1/SF_',species[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
}

###Species level plots (mean +/- sd)###
#soil vs lcy
soil<-subset(dataset, type %in% c('LCY','BGR'))
soil_mean_std <- soil %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

ggplot(soil_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  #scale_color_manual(values=c('red3','green3','cyan2','purple'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  #scale_fill_manual(values=c('red3','green3','cyan2','purple'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/SF_Full_Soil_LCY.png',sep=''),dpi=300,width=180,height=120,units='mm')

#rock biocrusts
rock<-subset(dataset, type %in% c('Rock_MOSS','Rock_LCHk','Rock_LCHg','Rock_LCHy','Rock_LCHb'))
rock_mean_std <- rock %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

rock_mean_std$type <- factor(rock_mean_std$type,levels=c('Rock_MOSS','Rock_LCHk','Rock_LCHg','Rock_LCHy','Rock_LCHb')) #Set the order for manual fill
ggplot(rock_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  scale_color_manual(values=c('purple4','darkgrey','green3','yellow2','cyan3'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  scale_fill_manual(values=c('purple4','darkgrey','green3','yellow2','cyan3'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/SF_Full_Rock_Biocrusts.png',sep=''),dpi=300,width=180,height=120,units='mm')

#soil biocrusts
biocrust<-subset(dataset, type %in% c('BGR','LCY','DCY','COTE','SQUA'))
biocrust_mean_std <- biocrust %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

biocrust_mean_std$type <- factor(biocrust_mean_std$type,levels=c('BGR','LCY','DCY','COTE','SQUA')) #Set the order for manual fill
ggplot(biocrust_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  scale_color_manual(values=c('sandybrown','red2','red4','green4','cyan2'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  scale_fill_manual(values=c('sandybrown','red2','red4','green4','cyan2'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/SF_Full_Soil_Biocrusts.png',sep=''),dpi=300,width=180,height=120,units='mm')

#####################################IC2#######################################
master <- read.csv2(paste(github_dir,'Data/NASA_EMIT_Campaign_062024/Level0/Data_Sheets/Biocrust/Datasheet_Biocrust_IC2_06062024.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level0/IC2_soil',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
  print(file)
  md<-get_metadata(file)
  temp_data<-get_spectra(file)
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1000)]-data[which(wvl==1001)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)],data[which(wvl==1001):which(wvl==1800)]+c1,data[which(wvl==1801):which(wvl==2500)]+c1+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  type=rep(master$Biocrust[i],length(wvl))
  rep=rep(master$Rep[i],length(wvl))
  sample=rep(master$Sample[i],length(wvl))
  
  ds=cbind(wvl,data_c,type,sample,rep)
  dataset <- rbind(dataset,ds)
}

#Name the columns
colnames(dataset) <- c("wavelength","reflectance","type","sample","rep")
dataset=transform(dataset,wavelength = as.numeric(wavelength))
dataset=transform(dataset,reflectance = as.numeric(reflectance))
dataset=transform(dataset,rep = as.numeric(rep))
dataset=transform(dataset,sample = as.numeric(sample))

#####ASD Reflectance Plots#################
###Species level plots (all samples)###
species<-c('BGR','LCY','DCY','MOSS')
for(i in 1:length(species)){
  ggplot(subset(dataset, type %in% species[i]), aes(x=wavelength,y=reflectance, group=rep, color=rep)) +
    geom_line(show.legend = T,linewidth=.5,linetype="solid") +
    facet_wrap(~sample)+
    scale_y_continuous("Reflectance") +
    scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,400)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(text = element_text(size = 18))
  ggsave(paste(github_dir,'/figures/Species/IC2_',species[i],'_Biocrust.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% species[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level1/IC2_',species[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
}

###Species level plots (mean +/- sd)###
#lcy vs soil
soil<-subset(dataset, type %in% c('LCY','BGR'))
soil_mean_std <- soil %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

ggplot(soil_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  #scale_color_manual(values=c('red3','green3','cyan2','purple'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  #scale_fill_manual(values=c('red3','green3','cyan2','purple'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/IC2_Full_Soil_LCY.png',sep=''),dpi=300,width=180,height=120,units='mm')

#soil biocrust
biocrust<-subset(dataset, type %in% c('BGR','LCY','DCY','MOSS'))
biocrust_mean_std <- biocrust %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

biocrust_mean_std$type <- factor(biocrust_mean_std$type,levels=c('BGR','LCY','DCY','MOSS')) #Set the order for manual fill
ggplot(biocrust_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  scale_color_manual(values=c('sandybrown','red2','red4','purple4'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  scale_fill_manual(values=c('sandybrown','red2','red4','purple4'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/IC2_Full_Soil_Biocrusts.png',sep=''),dpi=300,width=180,height=120,units='mm')

#####################################IC1#######################################
master <- read.csv2(paste(github_dir,'Data/NASA_EMIT_Campaign_062024/Level0/Data_Sheets/Biocrust/Datasheet_Biocrust_IC1_06072024.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level0/IC1_soil',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
  print(file)
  md<-get_metadata(file)
  temp_data<-get_spectra(file)
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1000)]-data[which(wvl==1001)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)],data[which(wvl==1001):which(wvl==1800)]+c1,data[which(wvl==1801):which(wvl==2500)]+c1+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  type=rep(master$Biocrust[i],length(wvl))
  rep=rep(master$Rep[i],length(wvl))
  sample=rep(master$Sample[i],length(wvl))
  
  ds=cbind(wvl,data_c,type,sample,rep)
  dataset <- rbind(dataset,ds)
}

#Name the columns
colnames(dataset) <- c("wavelength","reflectance","type","sample","rep")
dataset=transform(dataset,wavelength = as.numeric(wavelength))
dataset=transform(dataset,reflectance = as.numeric(reflectance))
dataset=transform(dataset,rep = as.numeric(rep))
dataset=transform(dataset,sample = as.numeric(sample))

#####ASD Reflectance Plots#################
###Species level plots (all samples)###
species<-c('BGR','LCY','DCY','MOSS')
for(i in 1:length(species)){
  ggplot(subset(dataset, type %in% species[i]), aes(x=wavelength,y=reflectance, group=rep, color=rep)) +
    geom_line(show.legend = T,linewidth=.5,linetype="solid") +
    facet_wrap(~sample)+
    scale_y_continuous("Reflectance") +
    scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,400)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(text = element_text(size = 18))
  ggsave(paste(github_dir,'/figures/Species/IC1_',species[i],'_Biocrust.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% species[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level1/IC1_',species[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
}

###Species level plots (mean +/- sd)###
#lcy vs soil
soil<-subset(dataset, type %in% c('LCY','BGR'))
soil_mean_std <- soil %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

ggplot(soil_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  #scale_color_manual(values=c('red3','green3','cyan2','purple'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  #scale_fill_manual(values=c('red3','green3','cyan2','purple'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/IC1_Full_Soil_LCY.png',sep=''),dpi=300,width=180,height=120,units='mm')

#soil biocrust
biocrust<-subset(dataset, type %in% c('BGR','LCY','DCY','MOSS'))
biocrust_mean_std <- biocrust %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

biocrust_mean_std$type <- factor(biocrust_mean_std$type,levels=c('BGR','LCY','DCY','MOSS')) #Set the order for manual fill
ggplot(biocrust_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  scale_color_manual(values=c('sandybrown','red2','red4','purple4'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  scale_fill_manual(values=c('sandybrown','red2','red4','purple4'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/IC1_Full_Soil_Biocrusts.png',sep=''),dpi=300,width=180,height=120,units='mm')

#####################################IC3#######################################
master <- read.csv2(paste(github_dir,'Data/NASA_EMIT_Campaign_062024/Level0/Data_Sheets/Biocrust/Datasheet_Biocrust_IC3_06082024.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level0/IC3_soil',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
  print(file)
  md<-get_metadata(file)
  temp_data<-get_spectra(file)
  data<-as.numeric(temp_data)
  wvl=as.integer(colnames(temp_data))
  c1=data[which(wvl==1000)]-data[which(wvl==1001)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)],data[which(wvl==1001):which(wvl==1800)]+c1,data[which(wvl==1801):which(wvl==2500)]+c1+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  type=rep(master$Biocrust[i],length(wvl))
  rep=rep(master$Rep[i],length(wvl))
  sample=rep(master$Sample[i],length(wvl))
  
  ds=cbind(wvl,data_c,type,sample,rep)
  dataset <- rbind(dataset,ds)
}

#Name the columns
colnames(dataset) <- c("wavelength","reflectance","type","sample","rep")
dataset=transform(dataset,wavelength = as.numeric(wavelength))
dataset=transform(dataset,reflectance = as.numeric(reflectance))
dataset=transform(dataset,rep = as.numeric(rep))
dataset=transform(dataset,sample = as.numeric(sample))

#####ASD Reflectance Plots#################
###Species level plots (all samples)###
species<-c('BGR','LCY','DCY')
for(i in 1:length(species)){
  ggplot(subset(dataset, type %in% species[i]), aes(x=wavelength,y=reflectance, group=rep, color=rep)) +
    geom_line(show.legend = T,linewidth=.5,linetype="solid") +
    facet_wrap(~sample)+
    scale_y_continuous("Reflectance") +
    scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,400)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(text = element_text(size = 18))
  ggsave(paste(github_dir,'/figures/Species/IC3_',species[i],'_Biocrust.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% species[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level1/IC3_',species[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
}

###Species level plots (mean +/- sd)###
#lcy vs soil
soil<-subset(dataset, type %in% c('LCY','BGR'))
soil_mean_std <- soil %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

ggplot(soil_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  #scale_color_manual(values=c('red3','green3','cyan2','purple'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  #scale_fill_manual(values=c('red3','green3','cyan2','purple'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/IC3_Full_Soil_LCY.png',sep=''),dpi=300,width=180,height=120,units='mm')

#soil biocrust
biocrust<-subset(dataset, type %in% c('BGR','LCY','DCY'))
biocrust_mean_std <- biocrust %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

biocrust_mean_std$type <- factor(biocrust_mean_std$type,levels=c('BGR','LCY','DCY')) #Set the order for manual fill
ggplot(biocrust_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  scale_color_manual(values=c('sandybrown','red2','red4'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  scale_fill_manual(values=c('sandybrown','red2','red4'))+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/IC3_Full_Soil_Biocrusts.png',sep=''),dpi=300,width=180,height=120,units='mm')



