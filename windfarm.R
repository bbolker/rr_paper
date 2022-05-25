library(dplyr)
library(lme4)
library(glmmTMB)
library(ggplot2)
library(gridExtra)
source("functions.R")
set.seed(101)
# ----------------------------------------------------------------
# Read in Data
# ----------------------------------------------------------------
#current results
load(file = "Data/windfarm_abund_data.RData")
## Boxplot of abundance of species by year and zone
labs=c(0,1,10,100)
wf_boxplot <- ggplot(windfarm, aes(x = Year, y = log(abundance + 1), fill = Zone)) +
  geom_boxplot(outlier.shape = 21)+
  facet_wrap(~ Species, ncol = 3, scales = "free") +
  labs(x = "Year", color  = "Zone", y = "Abundance [log(y+1)-scale]") +
  scale_y_continuous( breaks = log(labs + 1), labels = labs)+
  facet_wrap(~ Species, ncol = 3, scales="free") +
  theme_mine() +
  theme(legend.justification = c("top"),
        # legend.position = c(.95, .95),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4)) +
  scale_fill_brewer(palette = "Set2")

plot.name <- "wf_boxplot"
ggsave(filename = paste0("Plots/", plot.name, ".png"), device = "png",
       plot  = wf_boxplot, height=8, width = 16 )

# ----------------------------------------------------------------
#### Fit models
# # ----------------------------------------------------------------
wf.glmm.null <- glmmTMB(abundance ~  Zone + Year + diag(Zone + Year|Species) + (1|Station) +
                          rr(Species + 0 | ID, d = 2),
                        family = "poisson", data = windfarm,
                        control = glmmTMBControl(start_method = list(method = "res")  ))

wf.glmm <- glmmTMB(abundance ~ Zone*Year + diag(Zone*Year|Species) + (1|Station) +
                     rr(Species + 0 | ID, d = 2),
                   family = "poisson", data = windfarm,
                   control = glmmTMBControl(start_method = list(method = "res")  ) )

wf.glmm.comp.sym <- glmmTMB(abundance ~ Zone*Year + diag(Zone*Year|Species) + (1|Station) +
                              (1 | ID/Species),
                   family = "poisson", data = windfarm)

AIC(wf.glmm); AIC(wf.glmm.comp.sym)

### Windfarm model without environmental covariates (for ordination plot)
wf.glmm.unc <- glmmTMB(abundance ~  Species + rr(Species + 0 | ID, d = 2),
                       family = "poisson", data = windfarm, control = glmmTMBControl(start_method = list(method = "res")  ))

LRobs <- 2*logLik(wf.glmm) - 2*logLik(wf.glmm.null) #observed test stat
### --------------------------------------
# Parametric bootstrap test
### --------------------------------------
nSims <- 1000
LR <- rep(NA,nSims)
for(i in 1:nSims){
  simAbund <- simulate(wf.glmm.null)$sim_1 #simulate data under the null
    null <- try(glmmTMB(simAbund ~  Zone + Year   +
                              diag(Zone + Year|Species) + (1|Station) +
                              rr(Species + 0 | ID, d = 2),
                            family = "poisson",
                            data = windfarm))
    alt <- try(glmmTMB(simAbund ~ Zone*Year   +
                         diag(Zone*Year|Species) + (1|Station) +
                         rr(Species + 0 | ID, d = 2),
                       family = "poisson",
                       data = windfarm))
  LR[i] <- tryCatch({2*logLik(alt) - 2*logLik(null)}, error = function(e) {NA})  #if model doesn't converge, LR = NA
}

p <- mean(LR > LRobs, na.rm=TRUE) #P-value

### --------------------------------------
### Plot random effects - diag
### --------------------------------------
wf.fl <- wf.glmm$obj$env$report(wf.glmm$fit$parfull)$fact_load[[3]]
wf_bs = ranef(wf.glmm)
block.plot <- 1
b <- as.matrix(wf_bs$cond[[block.plot]])
## Get the names of levels
cnms <- wf.glmm$modelInfo$reTrms[["cond"]]$cnms   ## list of (named) terms and X columns
flist <- wf.glmm$modelInfo$reTrms[["cond"]]$flist ## list of grouping variables
levs <- lapply(flist, levels)
diag.plot.long <- as.data.frame(t(b)) %>%
  rownames_to_column("Coefficient") %>%
  pivot_longer(cols = -Coefficient, names_to = c(names(levs[block.plot])), values_to = "b" )

s1 <- TMB::sdreport(wf.glmm$obj,  getJointPrecision = TRUE)
# subset first
s1.b <- s1$jointPrecision[rownames(s1$jointPrecision)=="b", colnames(s1$jointPrecision)=="b"]
nbseq = split.bseq(wf.glmm)
b.cols <- split(1:ncol(s1.b), nbseq)
b.cols.plot <-  b.cols[[block.plot]]
s1.b.plot <- s1.b[b.cols.plot, b.cols.plot]
s1.b.inv <- solve(s1.b.plot)

reStruc <- wf.glmm$modelInfo$reStruc[[paste0("cond", "ReStruc")]] ## random-effects structure
nc <- vapply(reStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block
nb <- vapply(reStruc, function(x) x$blockReps, numeric(1)) ## number of blocks per RE (may != nlevs in some cases)

sd.b <- data.frame(sd.b = sqrt( diag( s1.b.inv) ))
coef.nam <- rep(cnms[[block.plot]], nb[[block.plot]])
nam.grp <- rep(levs[[block.plot]], each = nc[[block.plot]])
sd.diag.re <- cbind(coef.nam, nam.grp, sd.b)

names(sd.diag.re) <- c("Coefficient", names(levs[block.plot]), "sd.b")

diag.plot.data <- left_join(diag.plot.long, sd.diag.re, by = intersect(names(diag.plot.long), names(sd.diag.re)))
diag.plot.data <- diag.plot.data %>%
  mutate(conf.low = b - 1.96*sd.b,
         conf.high =  b + 1.96*sd.b,
         not.zero = ifelse(conf.high < 0 | conf.low > 0, 1, 0),
         highlight = ifelse(Coefficient %in% c("ZoneN:Year2010", "ZoneS:Year2010") &  not.zero, 1, 0 ) )
cols <- c("1" = "red", "0" = "black")

coef_labs <- gsub("Zone", "", unique(diag.plot.data$Coefficient ))
coef_labs <- gsub("Year", "", coef_labs)
#
plot.wf.diag <- diag.plot.data %>%
  mutate(Coefficient = factor(Coefficient, levels = unique(Coefficient), labels = coef_labs)) %>%
  filter(Coefficient %in% c("N:2010", "S:2010")) %>%
  ggplot(aes(x = Coefficient, y = b, colour = as.factor(highlight), fill = as.factor(highlight))) +
  geom_point()+
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0, linetype  = "dashed") +
  facet_grid(~ factor(Species, levels = spp_order) , scales = "free") +
  xlab("Coefficient") + ylab("Estimate") +
  theme_mine() +
  theme(legend.position="none",
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14 , vjust = -2),
        axis.title.y = element_text(size=14 , vjust = 1 ),
        strip.text.x = element_text(size=14, angle=0, hjust = 0)) +
  scale_color_manual(values = c("black", "red"))
plot.wf.diag

plot.name <- "wf_diag_int"
ggsave(filename = paste0("Plots/", plot.name, ".png"), device = "png",
       plot  = plot.wf.diag, height=8, width = 16 )
### --------------------------------------
### Ordination plots
### --------------------------------------
wf.glmm.unc.fl.b <- extract.rr.fl.b(object = wf.glmm.unc)
unc.wf.b <- wf.glmm.unc.fl.b$b
unc.wf.fl <-  as.data.frame(wf.glmm.unc.fl.b$fl)

zone_col <- windfarm %>%
  ungroup() %>%
  select(ID, Zone, Year) %>%
  distinct() %>%
  mutate(ID = as.character(ID))

## Adjust scaling so that both the site and species ordinations are on the same scale.
wf.unc.b2 <- as.data.frame( unc.wf.b * matrix(sqrt(colSums(unc.wf.fl^2)),
                                              nrow(unc.wf.b), 2, byrow=T))
names(wf.unc.b2) <- paste0("LV", 1:ncol(unc.wf.fl))
wf.unc.b2$ID <- rownames(wf.unc.b2)
wf.unc.b2 <- left_join(wf.unc.b2, zone_col, by = "ID")
wf.unc.fact_load2 <- as.data.frame(unc.wf.fl * matrix( sqrt(colSums(unc.wf.fl^2)), nrow(unc.wf.fl), 2, byrow=T))

cnms <- wf.glmm.unc$modelInfo$reTrms[["cond"]]$cnms
labels <- gsub("Species", "", cnms$ID)

unc.plot <- ggplot(wf.unc.b2) +
  geom_point(aes(x = LV1, y = LV2, color = Zone, shape = Year ) , size = 2) +
  geom_text(data = wf.unc.fact_load2, aes(x = V1, y = V2, label = labels),  colour="dark grey") +
  theme_mine() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14)) +
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
wf.b2 <- as.data.frame(wf.b * matrix( sqrt( colSums(wf.fact_load^2) ), nrow(wf.b), 2, byrow=T))
names(wf.b2) <- c("LV1", "LV2")
wf.b2$ID <- rownames(wf.b2)
wf.b2 <- left_join(wf.b2, zone_col, by = "ID")
wf.fact_load2 <- as.data.frame( wf.fact_load * matrix( sqrt( colSums(wf.fact_load^2) ), nrow( wf.fact_load ), 2, byrow=T))

residual.plot <- ggplot(wf.b2) +
  geom_point(aes(x=LV1,y=LV2, color = Zone, shape = Year ) , size = 2) +
  geom_text(aes(x = V1, y = V2, label = labels), data = wf.fact_load2, colour="dark grey") +
  theme_mine() +
  theme(legend.position =  c(0.9,0.8),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4),
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14)) +
  labs(title="(b) Residual ordination", x = "Latent variable 1", y = "Latent variable 2") +
  scale_color_brewer(palette = "Dark2")+
  coord_cartesian(xlim = c(-11, 10), ylim = c(-8, 8)) +
  scale_shape_manual(values=c(15,17,19)) # Change shapes

ordiplots <- grid.arrange(unc.plot, residual.plot, ncol=2)

plot.name <- "wf_ordiplot_2"
ggsave(filename = paste0("Plots/", plot.name, ".png"), device = "png",
       plot  = ordiplots, height=8, width = 16)

# Variance explained
unc.wf.ll <- (as.matrix(unc.wf.fl) %*% as.matrix(t(unc.wf.fl)))
wf.ll <- (wf.fact_load%*%t(wf.fact_load))
sum(diag(unc.wf.ll))
sum(diag(wf.ll))
var.explained <- (sum(diag(unc.wf.ll))-sum(diag(wf.ll)))/sum(diag(unc.wf.ll))

# save.image("Results/wf_full_results_May.RData")
load("Results/wf_full_results_May.RData")
