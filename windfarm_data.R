library(dplyr)
source("functions.R")
# ----------------------------------------------------------------
# Reduce windfarm data to only abundant species of fish
# i.e., species observed more than 4 times
# 29 July 2021
# Only keep stations that were observed before and after
# ----------------------------------------------------------------
dir <- "C:/Users/som/Desktop/Maeve/Projects/WindFarm/"
data_fold <- "Data/"
plots_folder <- "Plots/"
load(paste0(dir, data_fold, "wind_long.RData"))

data_tmp <- data_long %>%
  group_by(species) %>%
  mutate(totabund = sum(abundance))%>%
  ungroup() %>%
  rename(year = Year) %>%
  arrange(-totabund)

spp_order <- unique(data_tmp$species)
data_tmp <- data_tmp %>%
  mutate(species = factor(species, levels = spp_order),
         Station = as.factor(Station))

stat03 <- unique(data_tmp[data_tmp$year=="2003",]$Station)
stat10 <- unique(data_tmp[data_tmp$year=="2010",]$Station)

data_tmp <- data_tmp %>%
  filter(Station %in% stat03 )
data_tmp <- data_tmp %>%
  group_by(species) %>%
  mutate(obs = abundance > 0,
         n_ID = sum(obs))

spp_tot <- data_tmp %>%
  group_by(species) %>%
  summarise(totabund = totabund,
            n_ID = n_ID) %>%
  distinct()

# Remove species observed less than 5 times
datafit <- data_tmp %>%
  filter(n_ID > 4) %>%
  droplevels() %>%
  arrange(id) %>%
  mutate(Station = as.factor(Station),
         logAbund = log(abundance + 1)) %>%
  dplyr::select(- totabund, -obs, -n_ID)

datafit$Zone <- relevel(datafit$Zone, ref =c("WF") )
spp_order <- unique(datafit$species)
windfarm <- setNames(datafit, c("Year", "Zone", "Station", "Impact", "ID", "Area", "Species",
                                "abundance", "logAbund"))

save.image(file = "Data/windfarm_abund_data.RData")

