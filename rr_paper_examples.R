library(dplyr)
library(glmmTMB)
library(ggplot2)
library(gridExtra)
library(foreign)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)

source("functions.R")
set.seed(1001)
#### ----------------------------------------------------------------
#### ----------------------------------------------------------------
#### Wind farm example
#### ----------------------------------------------------------------
#### ----------------------------------------------------------------
data_fold <- "Data/"
plots_folder <- "Plots/"
load(paste0(data_fold, "wind_long.RData"))

wf_tmp <- data_long %>%
  group_by(species) %>%
  mutate(totabund = sum(abundance))%>%
  ungroup() %>%
  mutate(Station = as.factor(Station))

stat03 <- unique(wf_tmp[wf_tmp$Year=="2003",]$Station)
stat10 <- unique(wf_tmp[wf_tmp$Year=="2010",]$Station)

wf_tmp <- wf_tmp %>%
  filter(Station %in% stat03 ) %>%
  group_by(species) %>%
  mutate(obs = abundance > 0,
         n_ID = sum(obs) )

spp_tot <- wf_tmp %>%
  group_by(species) %>%
  summarise(totabund = totabund,
            n_ID = n_ID) %>%
  distinct() %>%
  ungroup() %>%
  arrange(-totabund)

# Remove species observed less than 5 times
datafit <- wf_tmp %>%
  filter(n_ID > 4) %>%
  droplevels() %>%
  arrange(-totabund)

spp_order <- unique(datafit$species)
windfarm <- datafit %>%
  arrange(id) %>%
  mutate(species = as.factor(species),
         Zone = relevel(Zone, ref = c("WF"))) %>%
  select(Year, Zone, Station, id, species, abundance) %>%
  rename(ID = id, Species = species)

## Boxplot of abundance of species by year and zone
labs=c(0,1,10,100)
wf_boxplot <- ggplot(windfarm, aes(x = Year, y = log(abundance + 1), fill = Zone)) +
  geom_boxplot(outlier.shape = 21)+
  facet_wrap(~ Species, ncol = 3, scales = "free") +
  labs(x = "Year", color  = "Zone", y = "Abundance [log(y+1)-scale]") +
  scale_y_continuous( breaks = log(labs + 1), labels = labs)+
  facet_wrap(~ Species, ncol = 3, scales="free") +
  theme_mine() +
  theme(legend.position="none",
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18 , vjust = -2),
        axis.title.y = element_text(size=18 , vjust = 1 ),
        strip.text.x = element_text(size=16, angle=0, hjust = 0)) +
  scale_fill_brewer(palette = "Set2")
wf_boxplot
#### ----------------------------------------------------------------
#### Fit models
#### ----------------------------------------------------------------
# predictor of interest is Zone:Year
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

#### ----------------------------------------------------------------
#### Parametric bootstrap test
#### ----------------------------------------------------------------
LRobs <- 2*logLik(wf.glmm) - 2*logLik(wf.glmm.null) #observed test stat
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

#### ----------------------------------------------------------------
#### ----------------------------------------------------------------
### Plot of diagonal random effect
#### ----------------------------------------------------------------
#### ----------------------------------------------------------------
### names of random effect terms, e.g. (x | group)
cnms <- wf.glmm$modelInfo$reTrms[["cond"]]$cnms   ## list of terms, i.e. names of x
flist <- wf.glmm$modelInfo$reTrms[["cond"]]$flist ## list of grouping variables
levs <- lapply(flist, levels) # levels of group
reStruc <- wf.glmm$modelInfo$reStruc$condReStruc ## random-effects structure
nc <- vapply(reStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block
nb <- vapply(reStruc, function(x) x$blockReps, numeric(1)) ## (nb is the number of levels in group)

block.plot <- 1
wf_bs = ranef(wf.glmm)
b.plot <- as.matrix(wf_bs$cond[[block.plot]])
diag.plot.long <- as.data.frame(t(b.plot)) %>%
  rownames_to_column("Coefficient") %>%
  pivot_longer(cols = -Coefficient, names_to = c(names(levs[block.plot])), values_to = "b" )

## Get conditional standard error of b
s1 <- TMB::sdreport(wf.glmm$obj,  getJointPrecision = TRUE)
# subset first
s1.b <- s1$jointPrecision[rownames(s1$jointPrecision)=="b", colnames(s1$jointPrecision)=="b"]
nbseq = split.bseq(wf.glmm)
b.cols <- split(1:ncol(s1.b), nbseq)
b.cols.plot <-  b.cols[[block.plot]]
s1.b.plot <- s1.b[b.cols.plot, b.cols.plot]
s1.b.inv <- solve(s1.b.plot)
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

coef_labs <- gsub("Zone", "", unique(diag.plot.data$Coefficient ))
coef_labs <- gsub("Year", "", coef_labs)

plot.wf.diag <- diag.plot.data %>%
  mutate(Coefficient = factor(Coefficient, levels = unique(Coefficient), labels = coef_labs)) %>%
  filter(Coefficient %in% c("N:2010", "S:2010")) %>%
  ggplot(aes(x = Coefficient, y = b, colour = as.factor(highlight), fill = as.factor(highlight))) +
  geom_point()+
  geom_linerange(aes(ymin = conf.low, ymax = conf.high, alpha = 0.5), size = 2) +
  geom_hline(yintercept = 0, linetype  = "dashed") +
  facet_grid(~ factor(Species, levels = spp_order) , scales = "free") +
  xlab("Coefficient") + ylab("Estimate") +
  theme_mine() +
  theme(legend.position="none",
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=14 , vjust = -2),
        axis.title.y = element_text(size=14 , vjust = 1 ),
        strip.text.x = element_text(size=14, angle=0, hjust = 0)) +
  scale_color_manual(values = c("black", "red"))
plot.wf.diag

#### ----------------------------------------------------------------
### Ordination plots
#### ----------------------------------------------------------------
#### unconstrained ordination plot
#### --------------------------------------
wf.glmm.unc.fl.b <- extract.rr.fl.b(object = wf.glmm.unc)
unc.wf.b <- wf.glmm.unc.fl.b$b
unc.wf.fl <-  as.data.frame(wf.glmm.unc.fl.b$fl)

zone_col <- windfarm %>%
  ungroup() %>%
  select(ID, Zone, Year) %>%
  distinct() %>%
  mutate(ID = as.character(ID))

## Adjust scaling so that both the site (u) and species (fl) ordinations are on the same scale.
wf.unc.b2 <- as.data.frame( unc.wf.b * matrix(sqrt(colSums(unc.wf.fl^2)), nrow(unc.wf.b), 2, byrow=T))
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

#### ----------------------------------------------------------------
#### Partial ordination plot - after controlling for environmental covariates
#### ----------------------------------------------------------------
wf.glmm.fl.b <- extract.rr.fl.b(wf.glmm)
wf.b <- wf.glmm.fl.b$b
wf.fl <- wf.glmm.fl.b$fl
## Adjust scaling so that both the site and species ordinations are on the same scale.
wf.b2 <- as.data.frame(wf.b * matrix( sqrt( colSums(wf.fl^2) ), nrow(wf.b), 2, byrow=T))
names(wf.b2) <- c("LV1", "LV2")
wf.b2$ID <- rownames(wf.b2)
wf.b2 <- left_join(wf.b2, zone_col, by = "ID")
wf.fl2 <- as.data.frame( wf.fl * matrix( sqrt( colSums(wf.fl^2) ), nrow( wf.fl ), 2, byrow=T))

residual.plot <- ggplot(wf.b2) +
  geom_point(aes(x=LV1,y=LV2, color = Zone, shape = Year ) , size = 2) +
  geom_text(aes(x = V1, y = V2, label = labels), data = wf.fl2, colour="dark grey") +
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

#### ----------------------------------------------------------------
#### ----------------------------------------------------------------
#### PIRLS example
#### ----------------------------------------------------------------
#### ----------------------------------------------------------------
rm(list = ls())
source("functions.R")
load("Data/pirls_cat.RData")

pirls <- pirls %>%
  droplevels() %>%
  mutate(Eco_disad = relevel(Eco_disad, ref = "0 to 10%"),
         Size_lib = relevel(Size_lib, ref = "500 Book Titles or Fewer")) %>%
  select(Overall1, Student_ID, SchoolsInCountry, Country, Eco_disad, Size_lib)

pirls$Overall <- pirls$Overall1
pirls$School <- pirls$SchoolsInCountry

###----- Interaction: School location and eco disadvantage
nlv <- 1:4
fit_list <- lapply(nlv,
                   function(d) {
                     fit.rr <- glmmTMB(Overall ~  Size_lib*Eco_disad +
                                         (1 |School) +
                                         rr(Size_lib*Eco_disad | Country, d),
                                       data = pirls,
                                       family = gaussian(),
                                       control = glmmTMBControl(start_method = list(method = "res")))
                   })

aic_vec <- sapply(fit_list, AIC)
rank <- nlv[ which(aic_vec - min(aic_vec, na.rm = TRUE)==0) ]

fit.rr <- glmmTMB(Overall ~  Size_lib*Eco_disad +
                    (1 |School) +
                    rr(Size_lib*Eco_disad | Country, rank),
                  data = pirls,
                  family = gaussian(),
                  control = glmmTMBControl(optCtrl=list(iter.max=1e5,eval.max=1e3),
                                           start_method = list(method = "res")))

#### ----------------------------------------------------------------
#### ----------------------------------------------------------------
## Plot of rr estimates
#### ----------------------------------------------------------------
#### ----------------------------------------------------------------
object <- fit.rr
### names for random effect terms, e.g. (x | group)
cnms <- object$modelInfo$reTrms[["cond"]]$cnms   ## list of terms, i.e. names of x
flist <- object$modelInfo$reTrms[["cond"]]$flist ## list of grouping variables
levs <- lapply(flist, levels) # levels of group
# (nb is the number of levels in group)
nb <- vapply(object$modelInfo$reStruc$condReStruc, function(x) x$blockReps, numeric(1)) ## number of blocks per RE (may != nlevs in some cases)

## Note: b is the name for all latent variables in the model, i.e., all random effects in the model
## u.rr will refer to the latent variables of rr block
## Then the estimates of the rr random effect are b.rr = lambda * u.rr
## lambda is a (nc x d) matrix of factor loadings (nc is the number of levels in x, i.e. school vars)
## where U ~ N(0, 1) and u.rr is (nb x d) matrix
block.rr <- 2
fl <- object$obj$env$report(object$fit$parfull)$fact_load[[block.rr]]
rownames(fl) = cnms[[block.rr]]
u.rr <- ranef(fit.rr)$cond[[block.rr]][,1:ncol(fl)]
b.rr <- fl %*% t(u.rr)
b.rr.long <- as.data.frame(b.rr) %>%
  rownames_to_column("Coefficient") %>%
  pivot_longer(cols = -Coefficient, names_to = c("Country"), values_to = "lu" )

## Get standard deviations of params - Square root of the diagonal of the hessian matrix.
s1 <- TMB::sdreport(object$obj,  getJointPrecision = TRUE)
## s1 has all params (beta, theta etc), but we only want b
s1.b <- s1$jointPrecision[rownames(s1$jointPrecision)=="b", colnames(s1$jointPrecision)=="b"]
h.b.inv <- solve(s1.b) # inverse matrix of all bs
bseq <- split.bseq(object) ### provides split of bs by random effect
b.cols <- split( 1:ncol(h.b.inv), bseq) # column numbers split by random effects
h.inv.block.rr <- h.b.inv[b.cols[[block.rr]], b.cols[[block.rr]]] ## us for rr block
### u.rr has (nb x d) estimates so only need (nb x d) X (nb x d) matrix (the remaining are set to zero)
u.rr.dim <-  nrow(u.rr) * ncol(fl)
h.inv.u.rr <- h.inv.block.rr[1:u.rr.dim, 1:u.rr.dim]

## Get standard errors of the estimates of lambda * u by
## Var[Lu] = L u L' where L is matrix of lambdas i.e., factor loadings (fixed)
## Create L matrix, where L is (nc x nb) x (nb x d)
I <- diag( rep(1,  nb[[block.rr]] ) ) ## Identity matrix
L <- kronecker(I, (fl)) ## factor loadings on diagonal
rownames(L) <- rep( rownames(fl), nb[[block.rr]] )
H.l.u <- L %*% h.inv.u.rr %*% t(L)
sd.lu <- sqrt( diag(H.l.u) )
### There's probably a better way to do this
sd.lu <- data.frame( rownames(H.l.u), sd.lu)
nam.grp <- rep( levs$Country, each = nrow(fl) )
sd.re <- cbind( sd.lu, nam.grp )
names(sd.re) <- c("Coefficient", "sd.lu", "Country")
ranef.rr <- left_join(b.rr.long, sd.re, by = intersect(names(b.rr.long), names(sd.re))) ## include std errors to estimates

pirls.plot.rr <- ranef.rr ### data frame for plot of rr
## Fix labels, order and levels of Coefficient
### Change coefficient labels to be succinct
coef_labels <- c("(Intercept)", "More than 5,000 Books", "501-5,000 Books",  "No Library",
                 "11 to 25%", "26 to 50%",  "More than 50%", "More than 5,000 Books: 11 to 25%",
                 "501-5,000 Books: 11 to 25%", "No Library: 11 to 25%", "More than 5,000 Books: 26 to 50%",
                 "501-5,000 Books: 26 to 50%", "No Library: 26 to 50%", "More than 5,000 Books: More than 50%",
                 "501-5,000 Books: More than 50%", "No Library: More than 50%")
coef_order <- c("(Intercept)", "No Library",  "501-5,000 Books", "More than 5,000 Books",
                "11 to 25%", "26 to 50%",  "More than 50%", "No Library: 11 to 25%", "No Library: 26 to 50%", "No Library: More than 50%",
                "501-5,000 Books: 11 to 25%", "501-5,000 Books: 26 to 50%", "501-5,000 Books: More than 50%",
                "More than 5,000 Books: 11 to 25%", "More than 5,000 Books: 26 to 50%", "More than 5,000 Books: More than 50%")
coef_plot_labels <- c("(Intercept)", "No Library",  "501-5,000 Books", "More than 5,000 Books", "11 to 25%", "26 to 50%",  "More than 50%",
                      "No Library: \n 11 to 25%", "No Library: \n 26 to 50%", "No Library: \n More than 50%",
                      "501-5,000 Books: \n 11 to 25%", "501-5,000 Books: \n 26 to 50%", "501-5,000 Books: \n More than 50%",
                      "More than 5,000 Books: \n 11 to 25%", "More than 5,000 Books: \n 26 to 50%", "More than 5,000 Books: \n More than 50%")

### Highlight two countries
country1 <- "Bulgaria"
country2 <- "Georgia"
pirls.plot.rr$Coefficient_plot <- plyr::mapvalues(pirls.plot.rr$Coefficient, from = unique(pirls.plot.rr$Coefficient), to = coef_labels)
pirls.plot.rr <- pirls.plot.rr %>%
  mutate(Coefficient_plot = factor(Coefficient_plot, levels = coef_order, ordered = T),
         Country = factor(Country, levels = sort(unique(ranef.rr$Country)), ordered = T),
         highlight_c = case_when(Country == country1 ~ "colour 1",
                                 Country == country2 ~ "colour 2",
                                 TRUE ~ "colour"),
         highlight_c = factor(highlight_c, ordered = T))
pirls.plot.rr$Coefficient_plot_labels = plyr::mapvalues(pirls.plot.rr$Coefficient_plot, from = coef_order, to = coef_plot_labels)

plot.rr <- ggplot(pirls.plot.rr,
                  aes(x = reorder(Country, desc(Country))  ,
                      y = lu,
                      colour = (highlight_c))) +
  geom_point() +
  geom_linerange(aes(ymin = lu - 1.96*sd.lu, ymax = lu + 1.96*sd.lu)) +
  geom_hline(yintercept = 0, linetype  = "dashed") +
  facet_grid(~ factor(Coefficient_plot_labels, ordered = T), scales = "free") +
  coord_flip() + theme_mine() +
  xlab("Country")+
  ylab("Estimate") +
  theme( axis.text.y = element_text( size = 10),
         axis.text.x = element_text( size = 9),
         axis.title = element_text( size = 12),
         strip.text = element_text(size = 12),
         legend.position="none",
         strip.text.x = element_text(size=10, angle=90, hjust = 0, vjust =0.5))+
  scale_color_manual(values  = c("#000000", "#0033FF", "#0A7527" ))

#### ---------------------------------------------------
#### Prediction plot of interaction
#### ----------------------------------------------------
pirls_pred <- pirls %>%
  select(Country, Eco_disad, Size_lib) %>%
  mutate(Country = factor(Country, levels = sort(unique(ranef.rr$Country)), ordered = T)) %>%
  distinct() %>%
  na.omit() %>%
  mutate(School = NA)

pred1 <- predict(fit.rr, newdata = pirls_pred, allow.new.levels = T)
pirls_pred$pred_y <- pred1

beta <- fixef(object)$cond
## Include se of beta and rr ('ignoring' school random effect)
# Var[X beta] = X Var[beta] X'
s1.beta <- s1$jointPrecision[rownames(s1$jointPrecision)=="beta", colnames(s1$jointPrecision)=="beta"]
h.beta.inv <- solve(s1.beta)
### Now get X in the same order as fixed effects order
Sl_lev <- factor( na.omit( unique(pirls$Size_lib) ) )
Eco_lev <- factor( na.omit( unique(pirls$Eco_disad) ) )
x.data <- expand.grid(Size_lib = Sl_lev, Eco_disad = Eco_lev ) %>%
  mutate(Eco_disad = relevel(Eco_disad, ref = "0 to 10%"),
         Size_lib = relevel(Size_lib, ref = "500 Book Titles or Fewer"))
X.mat <- model.matrix( ~ 1 + Size_lib*Eco_disad, data = x.data)
if(!all(colnames(X.mat) == names(beta))) print("Warning! Check X model matrix")

### Get standard errors
covar.beta <- X.mat %*% h.beta.inv %*% t(X.mat)
sd_beta <- sqrt(diag(covar.beta))
names(sd_beta) <- names(beta)

fixed.effects = data.frame(beta = beta, sd_beta = sd_beta) %>%
  rownames_to_column() %>%
  dplyr::rename(Coefficient = rowname)

sd.error.pred <- left_join(ranef.rr, fixed.effects, by = "Coefficient") %>%
  select(Coefficient, Country, beta, sd_beta, u, sd.u) %>%
  mutate(sd_bu = sd_beta + sd.u)

### Now I need to match coefficient names with levels of categorical variables
#### May have done this way too complicated!!!
beta_names<- unique(sd.error.pred$Coefficient)
Ed_lev <- levels(pirls$Eco_disad)
sl_lev <- levels(pirls$Size_lib)
eco_disad_beta <- c(rep(Ed_lev[1], 4), Ed_lev[2:4], rep(Ed_lev[2:4], each = 3)) #Eco_lev
size_lib_beta <- c( sl_lev, rep(sl_lev[1], each = 3 ), rep( sl_lev[2:4], times = 3) ) #Sl_lev
##Check if they match
beta_check <- data.frame(cbind(beta_names, eco_disad_beta, size_lib_beta))

sd.error.pred$Eco_disad <- plyr::mapvalues(sd.error.pred$Coefficient, from = beta_names, to = eco_disad_beta)
sd.error.pred$Size_lib <- plyr::mapvalues(sd.error.pred$Coefficient, from = beta_names, to = size_lib_beta)
pirls_pred_data <- left_join(pirls_pred, sd.error.pred, by = c("Eco_disad", "Size_lib", "Country")) %>%
  arrange(Country, Eco_disad, Size_lib)
## ----------------------------------
## Plot of 2 countries only
## ----------------------------------
### Change labels of size library to look neater
size_lib_labels <- c("Less than 500 Books",  "More than 5,000 Books", "501-5,000 Books", "No School Library" )
size_lib_order <- c("No School Library", "Less than 500 Books", "501-5,000 Books",  "More than 5,000 Books"  )
pirls_pred_data$Size_lib_lab <- plyr::mapvalues(pirls_pred_data$Size_lib, from = sl_lev, to = size_lib_labels)

pirls_pred_data <- pirls_pred_data %>%
  distinct() %>%
  mutate(Size_lib_lab = factor(Size_lib_lab, levels = size_lib_order, ordered = T),
         LCL = pred_y - 1.96*sd_bu,
         UCL = pred_y + 1.96*sd_bu,
         label = interaction(Country, Size_lib_lab, sep="-", lex.order=TRUE))

pirls_colours <- pirls_pred_data %>%
  filter( (Country == country1 | Country == country2) & !is.na(pred_y)) %>%
  droplevels()

## Red and Green
colours1  <- scales::seq_gradient_pal(low = "#33CCFF", high = "#000099", space = "Lab")(1:4/4)
colours2  <- scales::seq_gradient_pal(low = "#14E44D", high = "#0A7527", space = "Lab")(1:4/4)
# create values for scale colour manual - so order remains, i.e., darker colour as books increase
values <- setNames(c(colours1, colours2), levels(pirls_colours$label))

plot.int <- ggplot(filter(pirls_pred_data,  !is.na(pred_y)),
                   aes(x=Eco_disad, y=pred_y, group = interaction(Size_lib_lab, Country), linetype = Size_lib_lab)) +
  geom_line(alpha=0.05) +
  geom_line(data=pirls_colours, aes(x = Eco_disad, y = pred_y, colour = label, linetype = Size_lib_lab),  size = 1) +
  geom_ribbon(data= pirls_colours, aes(ymin = LCL, ymax = UCL, x = Eco_disad, fill = label), alpha = 0.1) +
  theme_mine() +
  ylab("Literacy score")  +
  xlab("Economic disadvantage") +
  scale_x_discrete( ) +
  coord_cartesian(clip = 'off') +
  scale_colour_manual(values = values)+
  scale_fill_manual(values=values, name="fill") +
  scale_linetype_manual(values=c("dotted", "dotdash",  "dashed", "solid"), name = "Library Size") + # Change linetypes
  theme(axis.text.x = element_text(size= 12, angle=0, hjust = 0.4, vjust =0),
        axis.text = element_text( size = 12),
        axis.title.x = element_text(size = 12, margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(size = 12, margin = unit(c(0, 5,0 , 0), "mm")),
        plot.margin=unit(c(1,5,1.5,1.2),"cm"), legend.position= "none")  +
  # scale_y_continuous(limits=c(350,680), breaks=seq(350,680,40))+
  guides(colour="none") +
  geom_label_repel(data = filter(pirls_colours, Eco_disad == "More than 50%"),
                   aes(label = label , colour = label),
                   nudge_x = 3.5, nudge_y = 2,
                   force = 0.5, label.size = NA, #remove box
                   direction = "y", force_pull = 0,
                   xlim = c(4,5.2), na.rm = TRUE)

## ----------------------------------
### Save plot of individual countries
## ----------------------------------
country_levels <- unique(pirls_pred_data$Country)
for(plot.country in country_levels){
  # The palette with grey:
  data_1country <- filter(pirls_pred_data,  !is.na(pred_y) & Country == plot.country) %>% droplevels()
  cColours  <- scales::seq_gradient_pal(low = "#FF6600", high = "#990000", space = "Lab")(1:4/4)
  c1_values <- setNames(cColours, levels(data_1country$label))

  plot_country <-  ggplot(filter(pirls_pred_data,  !is.na(pred_y)), aes(x=Eco_disad, y=pred_y, group = interaction(Size_lib, Country), linetype=Size_lib)) +
    geom_line(alpha=0.05) +
    geom_line(data = data_1country, aes(x = Eco_disad, y = pred_y, colour = label, group = Size_lib, linetype = Size_lib)) +
    geom_ribbon(data =data_1country, aes(ymin = LCL, ymax = UCL, x = Eco_disad, fill = label), alpha = 0.1) +
    theme_classic() +
    ylab("Literacy score")  +
    xlab("Economic disadvantage") +
    scale_x_discrete( ) +
    coord_cartesian(clip = 'off') +
    scale_linetype_manual(values=c("dotted", "dotdash",  "dashed", "solid"), name = "Library Size") + # Change linetypes
    theme(axis.text.x = element_text(size= 12, angle=0, hjust = 0.4, vjust =0),
          axis.text.y = element_text(size= 12),
          axis.text = element_text( size = 12),
          axis.title.x = element_text(size = 12, margin = unit(c(5, 0, 0, 0), "mm")),
          axis.title.y = element_text(size = 12, margin = unit(c(0, 5,0 , 0), "mm")),
          plot.margin=unit(c(1,5,1.5,1.2),"cm"), legend.position= "none")  +
    guides(colour="none") +
    geom_label_repel(data = filter(data_1country, Eco_disad == "More than 50%"),
                     aes(label = label , colour = label),
                     nudge_x = 3.5, nudge_y = 2,
                     force = 0.5, label.size = NA, #remove box
                     direction = "y", force_pull = 0,
                     xlim = c(4,5.2), na.rm = TRUE) +
    scale_fill_manual(values = c1_values)+
    scale_colour_manual(values = c1_values)

  ggsave(filename = paste0("Plots/pirls_by_country/plot_",plot.country,".png"), device = "png", plot  = plot_country, height=8, width = 10 )
}
