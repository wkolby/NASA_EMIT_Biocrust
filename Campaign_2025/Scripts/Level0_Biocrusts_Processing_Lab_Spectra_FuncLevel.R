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

#####################################Lab Spectra#######################################
master <- read.csv2(paste(github_dir,'Data/Level0_Lab_Data/Data_Sheets/EMIT_metadata_hyperspec_campaign2025.csv',sep=''),sep=',',header=T)
master$spectrum..<-str_pad(master$spectrum.., 5, pad = "0")


###Horse Thief###
#sites=c('HorseThief')
#covers=c('M','LC','DC','B','L','Lnf')

###Salt VAlley###
#sites=c('SaltValley')
#covers=c('M','LC','B','L')

###Panorama###
sites=c('Panorama')
covers=c('M','LC','B')

###ALL###
treats=c('dry','wet')

for(m in 1:length(sites)){
  for(k in 1:length(covers)){
    for(j in 1:length(treats)){
      dataset <- data.frame()
      msub <- subset(master,site == sites[m] & treatment == treats[j] & cover.type == covers[k])
      for(i in 1:length(msub$spectrum..)){
        file <- list.files(paste(github_dir,'/Data/Level0_Lab_Data/Spectra/',sep=''), full.names = T, pattern = paste(msub$spectrum..[i],".asd",sep=''))
        print(file)
        md<-get_metadata(file)
        temp_data<-get_spectra(file)
        data<-as.numeric(temp_data)
        wvl=as.integer(colnames(temp_data))
        c1=data[which(wvl==1001)]-data[which(wvl==1000)]
        c2=data[which(wvl==1800)]-data[which(wvl==1801)]
        data_c=c(data[1:which(wvl==1000)]+c1,data[which(wvl==1001):which(wvl==1800)],data[which(wvl==1801):which(wvl==2500)]+c2)
        data_s=data_c/sqrt(sum(data_c^2))
        site=rep(msub$site,length(wvl))
        fgroup=rep(msub$functional.group[i],length(wvl))
        ctype=rep(msub$cover.type[i],length(wvl))
        treat=rep(msub$treatment[i],length(wvl))
        rep=rep(msub$rep[i],length(wvl))
        sample=rep(msub$sample[i],length(wvl))
        
        ds=cbind(site,wvl,data_c,fgroup,ctype,treat,sample,rep)
        dataset <- rbind(dataset,ds)
      }
    
    
      #Name the columns
      colnames(dataset) <- c("site","wavelength","reflectance","fgroup","ctype","treat","sample","rep")
      dataset=transform(dataset,wavelength = as.numeric(wavelength))
      dataset=transform(dataset,reflectance = as.numeric(reflectance))
      dataset=transform(dataset,rep = as.numeric(rep))
      dataset=transform(dataset,sample = as.numeric(sample))
    
      #####ASD Reflectance Plots#################
      ###Species level plots (all samples)###
      ggplot(dataset, aes(x=wavelength,y=reflectance, group=rep, color=rep)) +
        geom_line(show.legend = T,linewidth=.5,linetype="solid") +
        facet_wrap(~sample)+
        scale_y_continuous("Reflectance") +
        scale_x_continuous("Wavelength (nm)",limits = c(400,2400), breaks = seq(400,2400,400)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        theme(text = element_text(size = 18))
      ggsave(paste(github_dir,'/Figures/Level0_Lab_Figures/LabSpectra_',sites[m],'_',covers[k],'_',treats[j],'_Biocrust.png',sep=''),dpi=300,width=180,height=120,units='mm')
      
      #write csv file
      dataset_rep_mean <- dataset %>% group_by(wavelength,sample) %>% summarise_at(vars(reflectance), list(reflectance=mean)) %>% as.data.frame()
      dataset_wide<- dataset_rep_mean %>% pivot_wider(names_from = sample, values_from = reflectance) #samples as columns
      write.csv(dataset_wide,paste(github_dir,'/Data/Level1_Lab_Data/',sites[m],'_',covers[k],'_',treats[j],'_Spectra.csv',sep=''),row.names=FALSE,col.names=TRUE)
    
    }
  }
}