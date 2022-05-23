library(glmmTMB)
library(foreign)
library(tidyverse)
library(dplyr)
library(TMB)
library(RColorBrewer)
library(grid)
library(directlabels)
library(ggrepel)
source("functions.R")

set.seed(101)
load("Data/pirls_cat.RData")

pirls <- pirls %>%
  droplevels()
pirls <- within(pirls, Eco_disad <- relevel(Eco_disad, ref = "0 to 10%"))
pirls <- within(pirls, Size_lib <- relevel(Size_lib, ref = "500 Book Titles or Fewer"))
pirls <- pirls %>%
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

save.image("Results/pirls_220521.RData")
load(("Results/pirls_220521.RData"))

## -------------------------------
## -------------------------------
## Plot of rr estimates
## -------------------------------
## -------------------------------
## Get estimates of the factor loadings
object = fit.rr
block.rr = 2
cnms <- object$modelInfo$reTrms[["cond"]]$cnms   ## list of (named) terms and X columns
flist <- object$modelInfo$reTrms[["cond"]]$flist ## list of grouping variables
levs <- lapply(flist, levels)

fl = object$obj$env$report(object$fit$parfull)$fact_load[[block.rr]]
rownames(fl) = cnms$Country ### CHANGE

## The variables for the reduced-rank block, b.rr = lambda * u where U ~ N(0, 1)
bs = ranef(fit.rr)
u = as.matrix(bs$cond[[block.rr]][,1:ncol(fl)])
b.rr = fl %*% t(u)
b.rr.long = as.data.frame(b.rr) %>%
  rownames_to_column("Coefficient") %>%
  pivot_longer(cols = -Coefficient, names_to = c("Country"), values_to = "lu" )

## Get standard deviations of params - Square root of the diagonal in the variance-covariance matrix.
s1 = TMB::sdreport(object$obj,  getJointPrecision = TRUE)
## subset b's first
s1.b = s1$jointPrecision[rownames(s1$jointPrecision)=="b", colnames(s1$jointPrecision)=="b"]
h.b.inv = solve(s1.b)
bseq = split.bseq(object)
b.cols = split(1:ncol(h.b.inv), bseq) ### provides split of bs by random effect
b.cols.rr =  b.cols[[block.rr]]
h.b.inv.tmp = h.b.inv[b.cols.rr, b.cols.rr]
## CHANGE!!! Will need to change how to get these
b0 = which( bs$cond[[block.rr]] == 0 )
h.b.rr.inv = h.b.inv.tmp[-b0, -b0]

## Get standard errors of the estimates of lambda * u by
## Var[Lu] = L u L' where L is matrix of lambdas i.e., factor loadings (fixed)
## Create L matrix
I = diag(rep(1, ncol(b.rr))) ## Identity matrix
L = kronecker(I, (fl)) ## Get factor loadings on diagonal
rownames(L) = rep( cnms$Country, ncol(b.rr) )
H.l.u = L %*% h.b.rr.inv %*% t(L)
sd.lu = sqrt(diag( H.l.u ))
### There's probably a better way to do this
sd.lu = data.frame( rownames(H.l.u), sd.lu)
nam.grp = rep( levs$Country, each = nrow(fl) )
sd.re = cbind( sd.lu, nam.grp )
names(sd.re) = c("Coefficient", "sd.u", "Country")

## include std errors to estimates
ranef.rr = left_join(b.rr.long, sd.re, by = intersect(names(b.rr.long), names(sd.re)))
ranef.rr$u = ranef.rr$lu

### Create data frame for plot of rr
pirls.plot.rr <- ranef.rr
## Fix labels, order and levels of Coefficient
### Change coefficient labels to be succinct
coef_labels = c("(Intercept)", "More than 5,000 Books", "501-5,000 Books",  "No Library",
                "11 to 25%", "26 to 50%",  "More than 50%", "More than 5,000 Books: 11 to 25%",
                "501-5,000 Books: 11 to 25%", "No Library: 11 to 25%", "More than 5,000 Books: 26 to 50%",
                "501-5,000 Books: 26 to 50%", "No Library: 26 to 50%", "More than 5,000 Books: More than 50%",
                "501-5,000 Books: More than 50%", "No Library: More than 50%")
coef_order = c("(Intercept)", "No Library",  "501-5,000 Books", "More than 5,000 Books",
               "11 to 25%", "26 to 50%",  "More than 50%", "No Library: 11 to 25%", "No Library: 26 to 50%", "No Library: More than 50%",
               "501-5,000 Books: 11 to 25%", "501-5,000 Books: 26 to 50%", "501-5,000 Books: More than 50%",
               "More than 5,000 Books: 11 to 25%", "More than 5,000 Books: 26 to 50%", "More than 5,000 Books: More than 50%")
coef_plot_labels = c("(Intercept)", "No Library",  "501-5,000 Books", "More than 5,000 Books", "11 to 25%", "26 to 50%",  "More than 50%",
                     "No Library: \n 11 to 25%", "No Library: \n 26 to 50%", "No Library: \n More than 50%",
                     "501-5,000 Books: \n 11 to 25%", "501-5,000 Books: \n 26 to 50%", "501-5,000 Books: \n More than 50%",
                     "More than 5,000 Books: \n 11 to 25%", "More than 5,000 Books: \n 26 to 50%", "More than 5,000 Books: \n More than 50%")

### Highlight two countries
country1 <- "Bulgaria"
country2 <- "Georgia"
pirls.plot.rr$Coefficient_plot = plysr::mapvalues(pirls.plot.rr$Coefficient, from = unique(pirls.plot.rr$Coefficient), to = coef_labels)
pirls.plot.rr = pirls.plot.rr %>%
  mutate(Coefficient_plot = factor(Coefficient_plot, levels = coef_order, ordered = T),
         Country = factor(Country, levels = sort(unique(ranef.rr$Country)), ordered = T),
         highlight_c = case_when(Country == country1 ~ "colour 1",
                                 Country == country2 ~ "colour 2",
                                 TRUE ~ "colour"),
         highlight_c = factor(highlight_c, ordered = T))
pirls.plot.rr$Coefficient_plot_labels = plyr::mapvalues(pirls.plot.rr$Coefficient_plot, from = coef_order, to = coef_plot_labels)

ggplot(pirls.plot.rr,
       aes(x = reorder(Country, desc(Country))  ,
           y = u,
           colour = (highlight_c))) +
  geom_point()+
  geom_linerange(aes(ymin = u - 1.96*sd.u, ymax = u + 1.96*sd.u)) +
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

beta = fixef(object)$cond
## Include se of beta and rr ('ignoring' school random effect)
# Var[X beta] = X Var[beta] X'
s1.beta = s1$jointPrecision[rownames(s1$jointPrecision)=="beta", colnames(s1$jointPrecision)=="beta"]
h.beta.inv = solve(s1.beta)
### Now get X in the same order as fixed effects order
Sl_lev <- factor( na.omit( unique(pirls$Size_lib) ) )
Eco_lev <- factor( na.omit( unique(pirls$Eco_disad) ) )
x.data <- expand.grid(Size_lib = Sl_lev ,
                      Eco_disad = Eco_lev ) %>%
  mutate(Eco_disad = relevel(Eco_disad, ref = "0 to 10%"),
         Size_lib = relevel(Size_lib, ref = "500 Book Titles or Fewer"))
X.mat <- model.matrix( ~ 1 + Size_lib*Eco_disad, data = x.data)
if(!all(colnames(X.mat) == names(beta))) print("Warning! Check X model matrix")

### Get standard errors
covar.beta = X.mat %*% h.beta.inv %*% t(X.mat)
sd_beta = sqrt(diag(covar.beta))
names(sd_beta) = names(beta)

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
# beta_check <- data.frame(cbind(beta_names, eco_disad_beta, size_lib_beta, levels(pirls_pred$Eco_disad), levels(pirls_pred$Size_lib)))
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
        axis.text.y = element_text(size= 12),
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

plot_name <- "pirls_interaction_with_sd2"
ggsave(filename = paste0("Plots/", plot_name,".png"), device = "png", plot  = plot.int, height=8, width = 14 )

level_sum <- pirls %>%
  group_by(Country, Eco_disad, Size_lib) %>%
  summarise(n = n()) %>%
  group_by(Country,Size_lib) %>%
  mutate(n_sl = sum(n)) %>%
  group_by(Country, Eco_disad) %>%
  mutate(n_eco = sum(n))

geor_sum <- level_sum[level_sum$Country == "Georgia",]
pirls_pred[pirls_pred$Country == "Georgia",]
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
