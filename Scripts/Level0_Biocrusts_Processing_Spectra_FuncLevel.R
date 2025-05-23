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
  c1=data[which(wvl==1001)]-data[which(wvl==1000)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)]+c1,data[which(wvl==1001):which(wvl==1800)],data[which(wvl==1801):which(wvl==2500)]+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  type=rep(master$Biocrust[i],length(wvl))
  site=rep('Test',length(wvl))
  
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

#write csv file
cal_mean <- cal %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean)) %>%
  as.data.frame()
cal_wide<- cal_mean %>% pivot_wider(names_from = type, values_from = reflectance) #samples as columns
write.csv(cal_wide,paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level1/Calibration_Rep1_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)

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

#write csv file
cal_wide<- cal %>% pivot_wider(names_from = type, values_from = reflectance) #samples as columns
write.csv(cal_wide,paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level1/Calibration_Rep2_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)

#####################################Sand Flats#######################################
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
  c1=data[which(wvl==1001)]-data[which(wvl==1000)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)]+c1,data[which(wvl==1001):which(wvl==1800)],data[which(wvl==1801):which(wvl==2500)]+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  site=rep('SF',length(wvl))
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

#Across site
biocrust_sf<-biocrust

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
  c1=data[which(wvl==1001)]-data[which(wvl==1000)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)]+c1,data[which(wvl==1001):which(wvl==1800)],data[which(wvl==1801):which(wvl==2500)]+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  site=rep('IC2',length(wvl))
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

#Across site
biocrust_ic2<-biocrust

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
  c1=data[which(wvl==1001)]-data[which(wvl==1000)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)]+c1,data[which(wvl==1001):which(wvl==1800)],data[which(wvl==1801):which(wvl==2500)]+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  site=rep('IC1',length(wvl))
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

#Across site
biocrust_ic1<-biocrust

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
  c1=data[which(wvl==1001)]-data[which(wvl==1000)]
  c2=data[which(wvl==1800)]-data[which(wvl==1801)]
  data_c=c(data[1:which(wvl==1000)]+c1,data[which(wvl==1001):which(wvl==1800)],data[which(wvl==1801):which(wvl==2500)]+c2)
  data_s=data_c/sqrt(sum(data_c^2))
  site=rep("IC3",length(wvl))
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

#Across site
biocrust_ic3<-biocrust

#####################################Horse Thief#######################################
master <- read.csv2(paste(github_dir,'Data/NASA_EMIT_Campaign_052025/Level0_Data/Data_Sheets/Biocrust/Datasheet_Biocrust_H_05062025.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'/Data/NASA_EMIT_Campaign_052025/Level0_Data/HorseThief/Biocrust/',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
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
  ggsave(paste(github_dir,'/Figures/Species/HorseThief_',species[i],'_Biocrust.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% species[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/NASA_EMIT_Campaign_052025/Level1/HorseThief_',species[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
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
ggsave(paste(github_dir,'/figures/HorseThief_Full_Soil_LCY.png',sep=''),dpi=300,width=180,height=120,units='mm')

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
ggsave(paste(github_dir,'/figures/HorseThief_Full_Soil_Biocrusts.png',sep=''),dpi=300,width=180,height=120,units='mm')

#Across site
biocrust_ht<-biocrust

#####################################Salt Valley#######################################
master <- read.csv2(paste(github_dir,'Data/NASA_EMIT_Campaign_052025/Level0_Data/Data_Sheets/Biocrust/Datasheet_Biocrust_S_05062025.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'/Data/NASA_EMIT_Campaign_052025/Level0_Data/SaltValley/Biocrust/',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
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
  ggsave(paste(github_dir,'/Figures/Species/SaltValley_',species[i],'_Biocrust.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% species[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/NASA_EMIT_Campaign_052025/Level1/SaltValley_',species[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
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
ggsave(paste(github_dir,'/figures/SaltValley_Full_Soil_LCY.png',sep=''),dpi=300,width=180,height=120,units='mm')

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
ggsave(paste(github_dir,'/figures/SaltValley_Full_Soil_Biocrusts.png',sep=''),dpi=300,width=180,height=120,units='mm')

#Across site
biocrust_sv<-biocrust

#####################################ACROSS SITE#######################################
all_sites<-rbind(biocrust_sf,biocrust_ic1,biocrust_ic2,biocrust_ic3,biocrust_ht,biocrust_sv)

#soil biocrust
all_sites<-subset(all_sites, type %in% c('BGR','LCY','DCY','LCNy','COTE','MOSS'))
all_mean_std <- all_sites %>%
  group_by(site,wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean=mean, sd=sd)) %>%
  as.data.frame()

all_mean_std$type <- factor(all_mean_std$type,levels=c('BGR','LCY','DCY','LCNy','COTE','MOSS')) #Set the order for manual fill
ggplot(all_mean_std, aes(x=wavelength,y=mean,group=site,color=site)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  #scale_color_manual(values=c())+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = site), alpha = .2) +
  #scale_fill_manual(values=c())+
  facet_wrap(~type)+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(300,2500), breaks = seq(350,2450,500)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/ALL_SITES_Full_Soil_Biocrusts_byType.png',sep=''),dpi=300,width=220,height=120,units='mm')

all_mean_std$site <- factor(all_mean_std$site,levels=c('HT','SV','IC1','IC2','IC3','SF')) #Set the order for manual fill
ggplot(all_mean_std, aes(x=wavelength,y=mean,group=type,color=type)) +
  geom_line(show.legend = T,linewidth=.5,linetype="solid") +
  #scale_color_manual(values=c('sandybrown','red2','red4','blue2','green4'))+
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type), alpha = .2) +
  #scale_fill_manual(values=c('sandybrown','red2','red4','blue2','green4'))+
  facet_wrap(~site)+
  scale_y_continuous("Reflectance") +
  scale_x_continuous("Wavelength (nm)",limits = c(300,2500), breaks = seq(350,2450,500)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text = element_text(size = 18))
ggsave(paste(github_dir,'/figures/ALL_SITES_Full_Soil_Biocrusts_bySite.png',sep=''),dpi=300,width=220,height=120,units='mm')

