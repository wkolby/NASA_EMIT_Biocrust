setwd("/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust")
github_dir <- "/Users/wksmith/Documents/GitHub/NASA_EMIT_Biocrust/"

install.packages("readxl")
install.packages("gridExtra")
install.packages("asdreader")
install.packages("tidyverse")
install.packages("patchwork")

library(asdreader)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(dplyr)


master <- read.csv2(paste(github_dir,'Data/NASA_EMIT_Campaign_062024/Level0/Data_Sheets/TEST/Datasheet_Test2_06052024.csv',sep=''),sep=',',header=T)
master$File<-str_pad(master$File, 5, pad = "0")
dataset <- data.frame()

for(i in 1:length(master$File)){
  spectral_files <- list.files(paste(github_dir,'/Data/NASA_EMIT_Campaign_062024/Level0/TEST2',sep=''), full.names = T, pattern = paste(master$File[i],".asd",sep=''))
}

for(i in 1:length(spectral_files)){
  file = spectral_files[i]
  spectrum = get_spectra(file)
  reflectance = as.numeric(spectrum)
  wavelength = as.integer(colnames(spectrum))
  c1 = reflectance[which(wavelength==1000)]-reflectance[which(wavelength==1001)]
  c2 = reflectance[which(wavelength==1800)]-reflectance[which(wavelength==1801)]
  reflectance_c = c(reflectance[1:which(wavelength==1000)]-c1,reflectance[which(wavelength==1001):which(wavelength==1800)],reflectance[which(wavelength==1801):which(wavelength==2500)]+c2)
  spectrum_id = basename(file)
  
  ds=cbind(wavelength, reflectance, reflectance_c, spectrum_id)
  dataset <- rbind(dataset,ds)
}

# creating a data frame or data sheet for the spectral data that can then be joined with the existing excel data sheet
spectra_df <- purrr::map_dfr(files, function(file) {
  spectrum <- get_spectra(file)  # Get the spectrum data from the .asd file
  data.frame(
    wavelength = as.integer(colnames(spectrum)),   # Assuming wavelengths are sequential
    reflectance = as.numeric(spectrum), # Convert the spectrum to numeric values
    spectrum_id = basename(file)        # Extract the file name as spectrum_id
  )
})

# making sure that each data frame has the same header (spectrum_id)
head(spectra_df)
head(metadata)

#Joining the data frames
full_data <- left_join(spectra_df, metadata, by = "spectrum_id")
head(full_data)


# Plotting literally evertything against each other, not super helpful
ggplot(full_data, aes(x = wavelength, y = reflectance, color = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Spectra by Treatment",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Treatment") +
  theme_minimal()


# Plotting one sample based on treatment type. Dry vs. Wet
sample_731 <- full_data %>% 
  filter(`Sample ID` == "731")

ggplot(sample_731, aes(x = wavelength, y = reflectance, color = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Lichen (1) Spectra Dry vs. Wet",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Treatment") +
  theme_minimal()


sample_732 <- full_data %>% 
  filter(`Sample ID` == "732")

ggplot(sample_732, aes(x = wavelength, y = reflectance, color = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Lichen (2) Spectra Dry vs. Wet",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Treatment") +
  theme_minimal()


sample_728 <- full_data %>% 
  filter(`Sample ID` == "728")

ggplot(sample_728, aes(x = wavelength, y = reflectance, color = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "bare Spectra Dry vs. Wet",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Treatment") +
  theme_minimal()


## Plotting 732 while omitting the ones where the lights were on

sample_732_filtered <- full_data %>%
  filter(`Sample ID` == "732") %>%
  filter (`light` == "off")
print(sample_732_filtered)


ggplot(sample_732_filtered, aes(x = wavelength, y = reflectance, color = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Sample 732 Dry vs. Wet Omitting lights",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Treatment") +
  theme_minimal()


#Plotting lichen against bare

lichen_lc_comparison_dry <- full_data %>%
  filter(`treatment` == "dry" & `Sample ID` %in% c("731", "728"))
print(lichen_lc_comparison_dry)

plot1 <- ggplot(lichen_lc_comparison_dry, aes(x = wavelength, y = reflectance, color = `cover type`)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Comparison of Lichen and Bare in Dry Treatment",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Sample ID") +
  theme_minimal()



lichen_lc_comparison_wet <- full_data %>%
  filter(`treatment` == "wet" & `Sample ID` %in% c("731", "728"))

plot2 <- ggplot(lichen_lc_comparison_wet, aes(x = wavelength, y = reflectance, color = `cover type`)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Comparison of Lichen and Light Cyano in wet Treatment",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Sample ID") + 
  scale_y_continuous(limits = c(0,0.5)) +
  theme_minimal()

## Putting the Dry and Wet comparisons ^^ next to each other

grid.arrange(plot1, plot2, ncol = 2)

## Plotting lichen (731) vs Bare, ignoring treatment 

lichen_lc_comparison <- full_data %>%
  filter(`Sample ID` %in% c("731", "728"))

ggplot(lichen_lc_comparison, aes(x = wavelength, y = reflectance, color = `cover type`, linetype = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Comparison of Lichen and Bare",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Sample ID") +
  scale_linetype_manual(values = c("dashed", "solid"))+
  theme_minimal()


#Plotting sample 732 to compare lights on vs lights off spectra

sample_732_light <- full_data %>%
  filter(`Sample ID` == "732") %>%
  filter (`treatment` == "dry")
print(sample_732_light)


ggplot(sample_732_light, aes(x = wavelength, y = reflectance, color = light)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Lights on vs. Lights off",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Treatment") +
  theme_minimal()


#Plotting light cyano vs. bare
bare_lc_comparison_dry <- full_data %>%
  filter(`treatment` == "dry" & `cover type` %in% c("bare", "light cyano"))
print(bare_lc_comparison_dry)

lc_bare_plot_dry <- ggplot(bare_lc_comparison_dry, aes(x = wavelength, y = reflectance, color = `cover type`)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Comparison of Light Cyanobacteria and Bare in Dry Treatment",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Sample ID") +
  theme_minimal() 


lc_bare_plot_dry  


###### Wet Treatment 
bare_lc_comparison_wet <- full_data %>%
  filter(`treatment` == "wet" & `cover type` %in% c("bare", "light cyano"))


lc_bare_plot_wet <- ggplot(bare_lc_comparison_dry, aes(x = wavelength, y = reflectance, color = `cover type`)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Comparison of Light Cyanobacteria and Bare in Wet Treatment",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Sample ID") +
  theme_minimal() 


lc_bare_plot_wet


###both treatments

bare_lc_comparison <- full_data %>%
  filter(`cover type` %in% c("light cyano", "bare"))

ggplot(bare_lc_comparison, aes(x = wavelength, y = reflectance, color = `cover type`, linetype = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Comparison of Light Cyano and Bare",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Sample ID") +
  scale_linetype_manual(values = c("dashed", "solid"))+
  theme_minimal()

### Light vs Dark Cyano
dc_lc_comparison <- full_data %>%
  filter(`cover type` %in% c("light cyano", "dark cyano"))

ggplot(dc_lc_comparison, aes(x = wavelength, y = reflectance, color = `cover type`, linetype = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Comparison of Light Cyano and Dark Cyano",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Sample ID") +
  scale_linetype_manual(values = c("solid", "dashed"))+
  theme_minimal()

###Comparing the two species of lichen
licen_species_comparison <- full_data %>%
  filter(species %in% c("Squamarina", "Collema cocco"))

ggplot(licen_species_comparison, aes(x = wavelength, y = reflectance, color = species, linetype = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Collema cocco vs. Squamarina",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Sample ID") +
  scale_linetype_manual(values = c("solid", "dashed"))+
  theme_minimal()


##Dry treatment
licen_species_comparison_dry <- full_data %>%
  filter(species %in% c("Squamarina", "Collema cocco") & treatment %in% c("dry"))

ggplot(licen_species_comparison_dry, aes(x = wavelength, y = reflectance, color = species)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Comparison of Lichen Species",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Sample ID") +
  scale_linetype_manual(values = c("solid", "dashed"))+
  theme_minimal()


## Collema cocco dry vs. wet treatment 

Collema_cocco_plot <- full_data %>% 
  filter(species %in% c("Collema cocco"))

ggplot(Collema_cocco_plot, aes(x = wavelength, y = reflectance, color = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Collema cocco dry vs. wet treatment",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Treatment") +
  theme_minimal()

## Squamarina
Squamarina_plot <- full_data %>% 
  filter(species %in% c("Squamarina"))

ggplot(Squamarina_plot, aes(x = wavelength, y = reflectance, color = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Squamarina dry vs. wet treatment",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Treatment") +
  theme_minimal()

## Squamarina DRY
Squamarina_dry_plot <- full_data %>% 
  filter(species %in% c("Squamarina") & treatment %in% c("dry"))

ggplot(Squamarina_dry_plot, aes(x = wavelength, y = reflectance, color = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "Squamarina dry treatment",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Treatment") +
  theme_minimal()

## Bare dry vs. wet

bare_comparison <- full_data %>% 
  filter(`cover type` %in% c("bare"))

bare_comparison_plot <- ggplot(bare_comparison, aes(x = wavelength, y = reflectance, color = treatment)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs (title = "Bare Ground Dry vs. Wet Treatment",
        x = "Wavelength",
        y = "Reflectance",
        color = "Treatment") +
  theme_minimal()

# Light Cyano Dry vs. Wet
lightcyano_comparison <- full_data %>%
  filter(`cover type` %in% c("light cyano"))

lightcyano_plot <- ggplot(lightcyano_comparison, aes(x =wavelength, y = reflectance, color = treatment)) +
  geom_line(aes(group =spectrum_id), alpha = 0.7) +
  labs (title = "Light Cyano Dry vs. Wet Treatment",
        x = "Wavelength",
        y = "Reflectance", 
        color = "Treatment") +
  theme_minimal() +
  theme(legend.position = "none")
  
       
lightcyano_plot + bare_comparison_plot
    

head(full_data)

##All cover types

all_cover_comparison <- full_data %>%
  filter(species %in% c("Squamarina", "Collema cocco", "light cyano", "dark cyano", "bare") & treatment %in% c("dry"))

ggplot(all_cover_comparison, aes(x = wavelength, y = reflectance, color = species)) +
  geom_line(aes(group = spectrum_id), alpha = 0.7) +
  labs(title = "All Cover Type Comparison",
       x = "Wavelength (nm)",
       y = "Reflectance",
       color = "Cover Type") +
  theme_minimal()


###### Summarize spectra by cover type and treatment
summarized_spectra <- full_data %>%
  group_by(species, treatment, wavelength) %>%
  summarise(mean_reflectance = mean(reflectance, na.rm = TRUE), 
            n = n(),  # Sample size
            se_reflectance = sd(reflectance, na.rm = TRUE) / sqrt(n)) %>%  
  ungroup()

print(summarized_spectra)
head(summarized_spectra)

summarized_comparison_dry <- summarized_spectra %>% 
  filter(species %in% c("Squamarina", "Collema cocco", "light cyano", "dark cyano", "bare") & treatment %in% c("dry"))



# Plot summarized spectra
ggplot(summarized_comparison_dry, aes(x = wavelength, y = mean_reflectance, color = `species`, linetype = treatment)) +
  geom_line(size = 1) +  # Line for the mean reflectance
  labs(title = "Summarized Spectra by Cover Type and Treatment",
       x = "Wavelength (nm)",
       y = "Mean Reflectance",
       color = "Cover Type",
       linetype = "Treatment",
       fill = "Cover Type") +
  
  theme_minimal()

summarized_comparison <- summarized_spectra %>% 
  filter(species %in% c("Squamarina", "Collema cocco", "light cyano", "dark cyano", "bare"))

ggplot(summarized_comparison, aes(x = wavelength, y = mean_reflectance, color = `species`, linetype = treatment)) +
  geom_line(size = 1) +  
  labs(title = "Summarized Spectra by Cover Type and Treatment",
       x = "Wavelength (nm)",
       y = "Mean Reflectance",
       color = "Cover Type",
       linetype = "Treatment",
       fill = "Cover Type") +
  facet_wrap(~treatment) +
  theme_minimal()      
       
## with standard error
ggplot(summarized_comparison, aes(x = wavelength, y = mean_reflectance, color = `species`, linetype = treatment)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_reflectance - sd_reflectance, 
                  ymax = mean_reflectance + sd_reflectance, 
                  fill = species), alpha = 0.2, color = NA) +
  labs(title = "Summarized Spectra by Cover Type and Treatment",
       x = "Wavelength (nm)",
       y = "Mean Reflectance",
       color = "Cover Type",
       linetype = "Treatment",
       fill = "Cover Type") +
  facet_wrap(~treatment) +
  theme_minimal()  




