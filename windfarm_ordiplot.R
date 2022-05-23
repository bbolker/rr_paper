library(TMB)
library(dplyr)
library(lme4)
library(glmmTMB)
library(ggplot2)
library(RColorBrewer)
source("functions.R")
set.seed(666)
# ----------------------------------------------------------------
# Read in Data
# ----------------------------------------------------------------
dir <- "C:/Users/som/Desktop/Maeve/Projects/WindFarm/"
data_fold <- "Data/"
plots_folder <- "Plots/"
load("Results/wf_LRsimulation_results_20211011.RData")

wf.glmm.unc.fl.b <- extract.rr.fl.b(wf.glmm.unc)
unc.wf.b <- wf.glmm.unc.fl.b$b
unc.wf.fl <-  as.data.frame(wf.glmm.unc.fl.b$fl)

zone_col <- windfarm %>% ungroup() %>%
  select(ID, Zone, Year) %>% distinct() %>% mutate(ID = as.character(ID))
alpha <- 1
alpha_plot=sqrt(max(apply(unc.wf.b^2,1,max)))/sqrt(max(apply(unc.wf.fl^2,1,max)))*alpha
sp.col="blue"

## Adjust scaling so that both the site and species ordinations are on the same scale.
wf.unc.b2 <- as.data.frame(unc.wf.b*matrix(sqrt(colSums(unc.wf.fl^2)),nrow(unc.wf.b),2,byrow=T))
names(wf.unc.b2) <- paste0("LV", 1:ncol(unc.wf.fl))
wf.unc.b2$ID <- rownames(wf.unc.b2)
wf.unc.b2 <- left_join(wf.unc.b2, zone_col, by = "ID")
wf.unc.fact_load2 <- as.data.frame(unc.wf.fl*matrix(sqrt(colSums(unc.wf.fl^2)),nrow(unc.wf.fl),2,byrow=T))
cols <- brewer.pal(3, "Dark2")


listname <- "cond"
cnms <- wf.glmm.unc$modelInfo$reTrms[[listname]]$cnms
labels = cnms$ID
labels <- gsub("Species", "", labels)

unc.plot <- ggplot(wf.unc.b2) +
  geom_point(aes(x=LV1,y=LV2, color = Zone, shape = Year ) , size = 2) +
  geom_text(aes(x = V1, y = V2, label = labels), data = wf.unc.fact_load2, colour="dark grey") +
  theme_mine() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16)) +
  labs(title="(a) Unconstrained ordination", x = "Latent variable 1", y = "Latent variable 2") +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(-11, 10), ylim = c(-8, 8)) +
  scale_shape_manual(values=c(15,17,19)) # Change shapes
unc.plot

#### ---- Partial ordination plot

wf.glmm.fl.b <- extract.rr.fl.b(wf.glmm)
wf.b <- wf.glmm.fl.b$b
wf.fact_load <- wf.glmm.fl.b$fl

## Adjust scaling so that both the site and species ordinations are on the same scale.
wf.b2 <- as.data.frame(wf.b*matrix(sqrt(colSums(wf.fact_load^2)),nrow(wf.b),2,byrow=T))
names(wf.b2) <- c("LV1", "LV2")
wf.b2$ID <- rownames(wf.b2)
wf.b2 <- left_join(wf.b2, zone_col, by = "ID")
wf.fact_load2 <- as.data.frame(wf.fact_load*matrix(sqrt(colSums(wf.fact_load^2)),nrow(wf.fact_load),2,byrow=T))

residual.plot <- ggplot(wf.b2) +
  geom_point(aes(x=LV1,y=LV2, color = Zone, shape = Year ) , size = 2) +
  geom_text(aes(x = V1, y = V2, label = labels), data = wf.fact_load2, colour="dark grey") +
  theme_mine() +
  theme(legend.position =  c(0.9,0.8),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4),
        plot.title = element_text(size = 16)) +
  labs(title="(b) Residual ordination", x = "Latent variable 1", y = "Latent variable 2") +
  scale_color_brewer(palette = "Dark2")+
  coord_cartesian(xlim = c(-11, 10), ylim = c(-8, 8)) +
  scale_shape_manual(values=c(15,17,19)) # Change shapes

residual.plot <- residual.plot +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
unc.plot <- unc.plot +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

  require(gridExtra)
  grid.arrange(unc.plot, residual.plot, ncol=2)
#


  ###
  unc.wf.ll <- (unc.wf.fl%*%t(unc.wf.fl))
  wf.ll <- (wf.fact_load%*%t(wf.fact_load))
  sum(diag(unc.wf.ll))
  sum(diag(wf.ll))
  var.explained <- (sum(diag(unc.wf.ll))-sum(diag(wf.ll)))/sum(diag(unc.wf.ll))

