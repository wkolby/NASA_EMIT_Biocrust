setwd("/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust")
github_dir <- "/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust/"
library(asdreader)
library(reshape2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)

#####################################TEST#######################################
master <- read.csv2(paste(github_dir,'Data/NASA_EMIT_Campaign_062024/Level0/Data_Sheets/TEST/Datasheet_Test2_06052024.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'Data/NASA_EMIT_Campaign_062024/Level0/TEST2',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
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
ggsave(paste(github_dir,'figures/TEST2_Full_Calibration.png',sep=''),dpi=300,width=180,height=120,units='mm')

#write csv file
cal_mean <- cal %>%
  group_by(wavelength,type) %>%
  summarise_at(vars(reflectance), list(mean)) %>%
  as.data.frame()
cal_wide<- cal_mean %>% pivot_wider(names_from = type, values_from = reflectance) #samples as columns
write.csv(cal_wide,paste(github_dir,'Data/NASA_EMIT_Campaign_062024/Level1/Calibration_Rep1_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)

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
ggsave(paste(github_dir,'figures/TEST2_Full_Calibration_REP2.png',sep=''),dpi=300,width=180,height=120,units='mm')

#write csv file
cal_wide<- cal %>% pivot_wider(names_from = type, values_from = reflectance) #samples as columns
write.csv(cal_wide,paste(github_dir,'Data/NASA_EMIT_Campaign_062024/Level1/Calibration_Rep2_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)

#####################################Horse Thief Panels#######################################
master <- read.csv2(paste(github_dir,'Data/NASA_EMIT_Campaign_052025/Level0_Data/Data_Sheets/Panels/Datasheet_USGS2_Panels_H_05062025.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'Data/NASA_EMIT_Campaign_052025/Level0_Data/HorseThief/Panels',sep=''), full.names = T, pattern = paste("Panel",master$File[i],".asd",sep=''))
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
  ggsave(paste(github_dir,'/Figures/Panels/USGS_',panels[i],'_Panel.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% panels[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/NASA_EMIT_Campaign_052025/Level1/Panels_',panels[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
}

#####################################Horse Thief Tarp#######################################
master <- read.csv2(paste(github_dir,'Data/NASA_EMIT_Campaign_052025/Level0_Data/Data_Sheets/Panels/Datasheet_NUSO_Tarp_H_05062025.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  file <- list.files(paste(github_dir,'Data/NASA_EMIT_Campaign_052025/Level0_Data/HorseThief/Panels',sep=''), full.names = T, pattern = paste("Panel",master$File[i],".asd",sep=''))
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
  ggsave(paste(github_dir,'/Figures/Panels/NUSO_',panels[i],'_Tarp.png',sep=''),dpi=300,width=180,height=120,units='mm')
  
  #write csv file
  dataset_rep_mean <- subset(dataset, type %in% panels[i]) %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
  dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
  write.csv(dataset_wide,paste(github_dir,'/Data/NASA_EMIT_Campaign_052025/Level1/Tarp_',panels[i],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
}

