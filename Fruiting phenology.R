



rm(list=ls(all=TRUE)) 


# PART 1. TESTING PHYLOGENETIC SIGNAL IN FRUIT PHENOLOGY 


file.choose()
data_fruitdrop<-read.csv("/Users/kowiyou/Desktop/ACDB/students/PhD students/2026/phd/Kandolo/paper 1/data_f/fruit_drop_final.csv", header=T)
names(data_fruitdrop)


# 1. assembling the phylogeny

# Prepare species list

sp.list <- data_fruitdrop %>%
  select(species, genus, family) %>%
  distinct()

# Replace spaces with underscores in species names
sp.list$species <- gsub(" ", "_", sp.list$species)


# Load V.PhyloMaker data


data(GBOTB.extended)
data(nodes.info.1)

# update family names in my data to match the megatree

sp.list$family[sp.list$genus == "Aptandra"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Cordia"] <- "Boraginaceae"
sp.list$family[sp.list$genus == "Diogoa"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Heisteria"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Octoknema"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Okoubaka"] <- "Santalaceae"
sp.list$family[sp.list$genus == "Ongokea"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Strombosia"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Strombosiopsis"] <- "Olacaceae"


# replace NA with the correct genus name in the genus columns for 2 species

sp.list$genus[sp.list$species == "Chrysophyllum_claussensii"] <- "Chrysophyllum"

sp.list$genus[sp.list$species == "Randia_acuminata"] <- "Randia"


# Build phylogeny

phylo_result <- phylo.maker(
  sp.list = sp.list,
  tree = GBOTB.extended,
  nodes = nodes.info.1,
  scenarios = "S3"
)

# Extract the phylogenetic tree
phy <- phylo_result$scenario.3


# Plot tree

plot(
  phy,
  cex = 0.5,
  no.margin = TRUE
)

# Match dataset to tree

data_fruitdrop$species_tree <- gsub(" ", "_", data_fruitdrop$species)

df_phy <- data_fruitdrop %>%
  filter(species_tree %in% phy$tip.label)


# Lets Create phylogenetic covariance matrix for the modelling 


phy_cov <- vcv(phy, corr = TRUE)


# Save outputs


write.tree(phy, file = "tree_congo_basin")



## Phylogenetic signal test


# 1. for fruit drop onset

library(ape)
library(phytools)

# Species-level mean trait
trait_df <- aggregate(onset_doy ~ species,
                      data = data_fruitdrop,
                      FUN = mean,
                      na.rm = TRUE)

# Create named vector
trait <- trait_df$onset_doy
names(trait) <- trait_df$species

# Keep only species present in both trait data and phylogeny
common_sp <- intersect(names(trait), phy$tip.label)

trait <- trait[common_sp]
phy_sub <- drop.tip(phy,
                    setdiff(phy$tip.label, common_sp))

# Ensure order matches
trait <- trait[phy_sub$tip.label]

# Blomberg's K with randomization test
K_result <- phylosig(
  tree = phy_sub,
  x = trait,
  method = "K",
  test = TRUE,
  nsim = 999
)



### Phylosignal in endset fruit drop

trait_df <- aggregate(end_doy ~ species,
                      data = data_fruit,
                      FUN = mean,
                      na.rm = TRUE)

trait <- trait_df$end_doy
names(trait) <- trait_df$species

common_sp <- intersect(names(trait), phy$tip.label)

trait <- trait[common_sp]
phy_sub <- drop.tip(phy,
                    setdiff(phy$tip.label, common_sp))

trait <- trait[phy_sub$tip.label]

phylosig(
  phy_sub,
  trait,
  method = "K",
  test = TRUE,
  nsim = 999
)



## signal in fruit drop duration
trait_df <- aggregate(duration_doy ~ species,
                      data = data_fruit,
                      FUN = mean,
                      na.rm = TRUE)

trait <- trait_df$duration_doy
names(trait) <- trait_df$species

common_sp <- intersect(names(trait), phy$tip.label)

trait <- trait[common_sp]
phy_sub <- drop.tip(phy,
                    setdiff(phy$tip.label, common_sp))

trait <- trait[phy_sub$tip.label]

phylosig(
  phy_sub,
  trait,
  method = "K",
  test = TRUE,
  nsim = 999
)





##### Phylogeny for fruit phenology duration




data_fruit<-read.csv("/Users/kowiyou/Desktop/ACDB/students/PhD students/2026/phd/Kandolo/paper 1/data_f/fruit_final.csv", header=T)
attach(data_fruit)
names(data_fruit)


# 1. assembling the phylogeny

# Prepare species list

sp.list <- data_fruit %>%
  select(species, genus, family) %>%
  distinct()

# Replace spaces with underscores in species names
sp.list$species <- gsub(" ", "_", sp.list$species)

# Fix obvious typo
data_fruit$species[data_fruit$species == "Randia_aruminata"] <-
  "Randia_acuminata"

# Fix Grewia/Microcos spelling issue
data_fruit$species[data_fruit$species %in%
                     c("Grewia_maleocarpoides",
                       "Microcos_malacocarpoides")] <-
  "Grewia_malacocarpoides"

# Load V.PhyloMaker data


data(GBOTB.extended)
data(nodes.info.1)

# update family names in my data to match the megatree

sp.list$family[sp.list$genus == "Aptandra"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Cordia"] <- "Boraginaceae"
sp.list$family[sp.list$genus == "Diogoa"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Heisteria"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Octoknema"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Okoubaka"] <- "Santalaceae"
sp.list$family[sp.list$genus == "Ongokea"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Strombosia"] <- "Olacaceae"
sp.list$family[sp.list$genus == "Strombosiopsis"] <- "Olacaceae"


# replace NA with the correct genus name in the genus columns for 2 species

sp.list$genus[sp.list$species == "Chrysophyllum_claussensii"] <- "Chrysophyllum"

sp.list$genus[sp.list$species == "Randia_acuminata"] <- "Randia"
sp.list$genus[sp.list$species == "Vernonia_conferta"] <- "Vernonia"




# Build phylogeny

phylo_result <- phylo.maker(
  sp.list = sp.list,
  tree = GBOTB.extended,
  nodes = nodes.info.1,
  scenarios = "S3"
)

# Extract the phylogenetic tree
phy <- phylo_result$scenario.3


# Plot tree

plot(
  phy,
  cex = 0.5,
  no.margin = TRUE
)

# Match dataset to tree

data_fruit$species_tree <- gsub(" ", "_", data_fruit$species)

df_phy <- data_fruit %>%
  filter(species_tree %in% phy$tip.label)


# Lets Create phylogenetic covariance matrix for the modelling 


phy_cov <- vcv(phy, corr = TRUE)


# Save outputs


write.tree(phy, file = "tree_congo_basin_fruitproduction")



##
# SIGNAL TEST IN FRUCTIFICATION ONSET, ENDSET AND DURATION

library(ape)
library(phytools)

# Function to calculate Blomberg's K
calc_K <- function(data, trait_col, phy){
  
  # Species means
  trait_df <- aggregate(
    data[[trait_col]],
    by = list(species = data$species),
    FUN = mean,
    na.rm = TRUE
  )
  
  names(trait_df) <- c("species", "trait")
  
  # Species present in both data and phylogeny
  common_sp <- intersect(trait_df$species,
                         phy$tip.label)
  
  trait_df <- trait_df[
    trait_df$species %in% common_sp,
  ]
  
  phy_sub <- drop.tip(
    phy,
    setdiff(phy$tip.label,
            common_sp)
  )
  
  # Create named vector
  trait <- trait_df$trait
  names(trait) <- trait_df$species
  
  # Reorder to match tree
  trait <- trait[phy_sub$tip.label]
  
  # Blomberg's K
  phylosig(
    tree = phy_sub,
    x = trait,
    method = "K",
    test = TRUE,
    nsim = 999
  )
}


# Fruit production onset


K_onset <- calc_K(
  data = data_fruit2,
  trait_col = "onset_doy",
  phy = phy
)

K_onset


# Fruit production end date


K_end <- calc_K(
  data = data_fruit2,
  trait_col = "end_doy",
  phy = phy
)

K_end


# Fruit production duration


K_duration <- calc_K(
  data = data_fruit2,
  trait_col = "duration_doy",
  phy = phy
)

K_duration


# Summary table


results <- data.frame(
  Trait = c("Onset", "End date", "Duration"),
  K = c(K_onset$K,
        K_end$K,
        K_duration$K),
  P = c(K_onset$P,
        K_end$P,
        K_duration$P)
)

results


# PART 2 MODELLING OF PHENOLOGY


# A- FRUIT DROP DATA

### Dataset 

file.choose()
data_fruitdrop<-read.csv("/Users/kowiyou/Desktop/ACDB/students/PhD students/2026/phd/Kandolo/paper 1/data_f/fruit_drop_final.csv", header=T)
names(data_fruitdrop)
str(data_fruitdrop)
# reconstruct the phylogeny using phylomaker

library(V.PhyloMaker)
library(dplyr)
library(ape)
library(ggplot2)
library(patchwork)
library(phytools)
library(glmmTMB)

data_fruitdrop$temp_z <- scale(data_fruitdrop$mean_temp)
data_fruitdrop$precip_z <- scale(data_fruitdrop$total_precip)

data_fruitdrop$temp_z <- as.numeric(data_fruitdrop$temp_z)
data_fruitdrop$precip_z <- as.numeric(data_fruitdrop$precip_z)

class(data_fruitdrop$temp_z)
class(data_fruitdrop$precip_z)

data_fruitdrop$site <- factor(data_fruitdrop$site)

levels(data_fruitdrop$site)


newdat$family  <- unique(data_fruitdrop$family)[1]
newdat$genus   <- na.omit(unique(data_fruitdrop$genus))[1]
newdat$species <- unique(data_fruitdrop$species)[1]
newdat$year    <- unique(data_fruitdrop$year)[1]



# modeling fruit drop phenology

# 1. onset fruit drop
mod1a <- glmmTMB(
  onset_doy ~ temp_z * precip_z *site +
    (1|family/genus/species) +
    (1|year),
  family = gaussian(),
  data = data_fruitdrop
)

mod1b<-update(mod1a, ~.-temp_z:precip_z:site)
mod1c<-update(mod1b, ~.-temp_z:precip_z)
mod1d<-update(mod1c, ~.-precip_z:site)
mod1e<-update(mod1d, ~.-temp_z)
mod1f<-update(mod1e, ~.-precip_z)
mod1g<-update(mod1f, ~.-site)

AIC(mod1a,mod1b,mod1c,mod1d,mod1e,mod1f,mod1g)



summary(mod1a)





## NOW Lets plot the 3 way interactions

library(glmmTMB)
library(ggplot2)


# 1. Standardise predictors

data_fruitdrop$temp_z <- as.numeric(scale(data_fruitdrop$mean_temp))
data_fruitdrop$precip_z <- as.numeric(scale(data_fruitdrop$total_precip))

# Ensure site is factor
data_fruitdrop$site <- factor(data_fruitdrop$site)


# 2. Fit model

mod1a <- glmmTMB(
  onset_doy ~ temp_z * precip_z * site +
    (1 | family/genus/species) +
    (1 | year),
  family = gaussian(),
  data = data_fruitdrop
)


# 3. Create prediction grid

temp_seq <- seq(min(data_fruitdrop$temp_z),
                max(data_fruitdrop$temp_z),
                length.out = 100)

precip_seq <- seq(min(data_fruitdrop$precip_z),
                  max(data_fruitdrop$precip_z),
                  length.out = 100)

newdat <- expand.grid(
  temp_z = temp_seq,
  precip_z = precip_seq,
  site = levels(data_fruitdrop$site)
)

# Add random-effect placeholders
newdat$family  <- data_fruitdrop$family[1]
newdat$genus   <- na.omit(data_fruitdrop$genus)[1]
newdat$species <- data_fruitdrop$species[1]
newdat$year    <- data_fruitdrop$year[1]


# 4. Predictions

newdat$pred <- predict(
  mod1a,
  newdata = newdat,
  type = "response",
  re.form = NA
)


# 5. Plot all figures

library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)
library(glmmTMB)


# FIXED EFFECTS TABLE (Panel C data prep)

coefs <- summary(mod1a)$coefficients$cond

coef_df <- data.frame(
  term = rownames(coefs),
  estimate = coefs[,1],
  se = coefs[,2]
) %>%
  filter(term != "(Intercept)")


# PANEL A: Observed data trends

pB <- ggplot(data_fruitdrop,
             aes(x = temp_z, y = onset_doy, color = site)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic(base_size = 13) +
  labs(
    x = "Temperature (standardised)",
    y = "Observed fruit drop onset (DOY)",
    title = "B: Observed climate–phenology relationship"
  )


# PANEL B: Model prediction surface

pC <- ggplot(newdat,
             aes(x = temp_z, y = precip_z, z = pred)) +
  geom_contour_filled(bins = 12) +
  facet_wrap(~site, nrow = 1) +
  scale_fill_viridis_d(option = "C", name = "Onset DOY") +
  theme_classic(base_size = 13) +
  labs(
    x = "Temperature (standardised)",
    y = "Precipitation (standardised)",
    title = "C: Climate interaction surface"
  ) +
  theme(strip.text = element_text(face = "bold"))


# PANEL C: Standardised effect sizes

pA <- ggplot(coef_df,
             aes(x = reorder(term, estimate),
                 y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = estimate - se,
                    ymax = estimate + se),
                width = 0.2) +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Standardised effect size (β ± SE)",
    title = "A: Climate effect sizes on fruit drop onset"
  )


# COMBINE FIGURE

final_fig <- (pA | pB) / pC +
  plot_layout(heights = c(2, 1))

final_fig



ggsave(
  "Figure1_fruit_phenology_synthesis.png",
  final_fig,
  width = 12,
  height = 9,
  dpi = 600
)




# 2. endset fruit drop
mod2a <- glmmTMB(
  end_doy ~ temp_z * precip_z *site +
    (1|family/genus/species) +
    (1|year),
  family = gaussian(),
  data = data_fruitdrop
)

summary(mod2a)



# overall temp effect on end drop
a<-ggplot(data=data_fruitdrop, aes(x=mean_temp, y=end_doy)) +
  geom_point(na.rm = TRUE, size=4) +
  geom_smooth(method='lm', linewidth=2.5) +
  labs(title="a)", x='temperature (oC)', y='end fruit drop phenology (DOY)') +
  scale_color_manual(values=c("#0073C2FF", "#EFC000FF")) +
  theme_bw(14) +
  theme(legend.position.inside = c(0.83, 0.2),
        legend.key = element_rect(fill="white", color="grey50"),
        legend.background = element_blank()) +   # remove legend background'
  guides(color=guide_legend("site"))  # add guide properties by aesthetic

#a

# effects of temp on end drop by site since mod2a shows temp effect is site dependent
a1<-ggplot(data=data_fruitdrop, aes(x=mean_temp, y=end_doy, group=site, color=site)) +
  geom_point(na.rm = TRUE, size=4) +
  geom_smooth(method='lm', linewidth=2.5) +
  labs(title="b)", x='temperature (oC)', y='') +
  scale_color_manual(values=c("#0073C2FF", "#EFC000FF")) +
  theme_bw(14) +
  theme(legend.position.inside = c(0.83, 0.2),
        legend.key = element_rect(fill="white", color="grey50"),
        legend.background = element_blank()) +   # remove legend background'
  guides(color=guide_legend("site"))  # add guide properties by aesthetic

# a1

a+a1



# 3. duration fruit drop
mod3a <- glmmTMB(
  duration_doy ~ temp_z * precip_z *site +
    (1|family/genus/species) +
    (1|year),
  family = gaussian(),
  data = data_fruitdrop
)


summary(mod3a)




# FIGURE DURATION FRUIT DROP


library(glmmTMB)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)


# 1. STANDARDISE CLIMATE VARIABLES

# 2. FIT MODEL (DURATION)

mod3a <- glmmTMB(
  duration_doy ~ temp_z * precip_z * site +
    (1 | family/genus/species) +
    (1 | year),
  family = gaussian(),
  data = data_fruitdrop
)


# 3. PREDICTION GRID

temp_seq <- seq(min(data_fruitdrop$temp_z),
                max(data_fruitdrop$temp_z),
                length.out = 100)

precip_seq <- seq(min(data_fruitdrop$precip_z),
                  max(data_fruitdrop$precip_z),
                  length.out = 100)

newdat <- expand.grid(
  temp_z = temp_seq,
  precip_z = precip_seq,
  site = levels(data_fruitdrop$site)
)

# Random-effect placeholders
newdat$family  <- data_fruitdrop$family[1]
newdat$genus   <- na.omit(data_fruitdrop$genus)[1]
newdat$species <- data_fruitdrop$species[1]
newdat$year    <- data_fruitdrop$year[1]


# 4. PREDICTIONS

newdat$pred <- predict(
  mod3a,
  newdata = newdat,
  type = "response",
  re.form = NA
)

# 5. PANEL A: OBSERVED DATA
pBd <- ggplot(data_fruitdrop,
              aes(x = temp_z, y = duration_doy, color = site)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic(base_size = 13) +
  labs(
    x = "Temperature (standardised)",
    y = "Observed fruit-drop duration (days)",
    title = "B: relationship between climate and fruit drop duration "
  )


# 6. PANEL B: MODEL SURFACE

pCd <- ggplot(newdat,
              aes(x = temp_z, y = precip_z, z = pred)) +
  geom_contour_filled(bins = 12) +
  facet_wrap(~site, nrow = 1) +
  scale_fill_viridis_d(option = "C", name = "Duration (days)") +
  theme_classic(base_size = 13) +
  labs(
    x = "Temperature (standardised)",
    y = "Precipitation (standardised)",
    title = "C: Climate interaction surface (duration fruit drop)"
  ) +
  theme(strip.text = element_text(face = "bold"))


# 7. PANEL C: EFFECT SIZES

coefs <- summary(mod3a)$coefficients$cond

coef_df <- data.frame(
  term = rownames(coefs),
  estimate = coefs[,1],
  se = coefs[,2]
) %>%
  filter(term != "(Intercept)")

pAd <- ggplot(coef_df,
              aes(x = reorder(term, estimate),
                  y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = estimate - se,
                    ymax = estimate + se),
                width = 0.2) +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Standardised effect sizes (β ± SE)",
    title = "A: Climate effect sizes on fruit drop duration"
  )


# 8. COMBINE FIGURE

final_fig_duration <- (pAd | pBd) / pCd +
  plot_layout(heights = c(2, 1))

final_fig_duration

# 9. SAVE FIGURE

ggsave(
  "Figure_duration_fruitdrop_synthesis.png",
  final_fig_duration,
  width = 12,
  height = 9,
  dpi = 600
)






# B- modeling fructification phenology



data_fruit<-read.csv("/Users/kowiyou/Desktop/ACDB/students/PhD students/2026/phd/Kandolo/paper 1/data_f/fruit_final.csv", header=T)
attach(data_fruit)
names(data_fruit)


data_fruit$temp_z <- scale(data_fruit$mean_temp)
data_fruit$precip_z <- scale(data_fruit$total_precip)

# modeling fructification phenology

# 1. onset fructification phenology
mod1b <- glmmTMB(
  onset_doy ~ temp_z * precip_z *site +
    (1|family/genus/species) +
    (1|year),
  family = gaussian(),
  data = data_fruit
)

summary(mod1b)

# FIGURES ONSET FRUCTIFICATION



# 1. STANDARDISE CLIMATE VARIABLES

data_fruit$temp_z <- as.numeric(scale(data_fruit$mean_temp))
data_fruit$precip_z <- as.numeric(scale(data_fruit$total_precip))
data_fruit$site <- factor(data_fruit$site)


# 2. FIT MODEL

mod1b <- glmmTMB(
  onset_doy ~ temp_z * precip_z * site +
    (1 | family/genus/species) +
    (1 | year),
  family = gaussian(),
  data = data_fruit
)


# 3. CREATE PREDICTION GRID

temp_seq <- seq(min(data_fruit$temp_z),
                max(data_fruit$temp_z),
                length.out = 100)

precip_seq <- seq(min(data_fruit$precip_z),
                  max(data_fruit$precip_z),
                  length.out = 100)

newdat <- expand.grid(
  temp_z = temp_seq,
  precip_z = precip_seq,
  site = levels(data_fruit$site)
)

# Add placeholders for random effects
newdat$family  <- data_fruit$family[1]
newdat$genus   <- na.omit(data_fruit$genus)[1]
newdat$species <- data_fruit$species[1]
newdat$year    <- data_fruit$year[1]


# 4. PREDICT

newdat$pred <- predict(
  mod1b,
  newdata = newdat,
  type = "response",
  re.form = NA
)


# 5. PANEL A: OBSERVED DATA

pBf <- ggplot(data_fruit,
              aes(x = temp_z, y = onset_doy, color = site)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic(base_size = 13) +
  labs(
    x = "Temperature (standardised)",
    y = "Observed fructification onset (DOY)",
    title = "B: Observed climate–fructification relationship"
  )


# 6. PANEL B: MODEL SURFACE

pCf <- ggplot(newdat,
              aes(x = temp_z, y = precip_z, z = pred)) +
  geom_contour_filled(bins = 12) +
  facet_wrap(~site, nrow = 1) +
  scale_fill_viridis_d(option = "C", name = "Onset DOY") +
  theme_classic(base_size = 13) +
  labs(
    x = "Temperature (standardised)",
    y = "Precipitation (standardised)",
    title = "C: Climate interaction surface (fructification onset)"
  ) +
  theme(strip.text = element_text(face = "bold"))


# 7. PANEL C: EFFECT SIZES

coefs <- summary(mod1b)$coefficients$cond

coef_df <- data.frame(
  term = rownames(coefs),
  estimate = coefs[,1],
  se = coefs[,2]
) %>%
  filter(term != "(Intercept)")

pAf <- ggplot(coef_df,
              aes(x = reorder(term, estimate),
                  y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = estimate - se,
                    ymax = estimate + se),
                width = 0.2) +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Standardised effect size (β ± SE)",
    title = "A: Climate effect sizes (fructification onset)"
  )


# 8. COMBINE FIGURE

final_fig_fruct <- (pAf | pBf) / pCf +
  plot_layout(heights = c(2, 1))

final_fig_fruct


# 9. SAVE FIGURE

ggsave(
  "Figure_fructification_onset_synthesis.png",
  final_fig_fruct,
  width = 12,
  height = 9,
  dpi = 600
)



# 2. endset fructification phenology


mod2b <- glmmTMB(
  end_doy ~ temp_z * precip_z *site +
    (1|family/genus/species) +
    (1|year),
  family = gaussian(),
  data = data_fruit
)

summary(mod2b)




# overall temp effect on end fructification
b<-ggplot(data=data_fruit, aes(x=mean_temp, y=end_doy)) +
  geom_point(na.rm = TRUE, size=4) +
  geom_smooth(method='lm', linewidth=2.5) +
  labs(title="c)", x='temperature (oC)', y='end fructification phenology (DOY)') +
  scale_color_manual(values=c("#0073C2FF", "#EFC000FF")) +
  theme_bw(14) +
  theme(legend.position.inside = c(0.83, 0.2),
        legend.key = element_rect(fill="white", color="grey50"),
        legend.background = element_blank()) +   # remove legend background'
  guides(color=guide_legend("site"))  # add guide properties by aesthetic



# effects of temp on end fructification by site since mod2b shows temp effect is site dependent
b1<-ggplot(data=data_fruit, aes(x=mean_temp, y=end_doy, group=site, color=site)) +
  geom_point(na.rm = TRUE, size=4) +
  geom_smooth(method='lm', linewidth=2.5) +
  labs(title="d)", x='temperature (oC)', y='') +
  scale_color_manual(values=c("#0073C2FF", "#EFC000FF")) +
  theme_bw(14) +
  theme(legend.position.inside = c(0.83, 0.2),
        legend.key = element_rect(fill="white", color="grey50"),
        legend.background = element_blank()) +   # remove legend background'
  guides(color=guide_legend("site"))  # add guide properties by aesthetic

# a1

b+b1


gridExtra::grid.arrange(a,a1,b,b1, ncol=2)



# 3. duration fructification phenology

mod3b <- glmmTMB(
  duration_doy ~ temp_z * precip_z *site +
    (1|family/genus/species) +
    (1|year),
  family = gaussian(),
  data = data_fruit
)

summary(mod3b)

# FIGURES fructification duration
# 



# 1. STANDARDISE VARIABLES

data_fruit$temp_z <- as.numeric(scale(data_fruit$mean_temp))
data_fruit$precip_z <- as.numeric(scale(data_fruit$total_precip))
data_fruit$site <- factor(data_fruit$site)


# 2. FIT MODEL

mod3b <- glmmTMB(
  duration_doy ~ temp_z * precip_z * site +
    (1 | family/genus/species) +
    (1 | year),
  family = gaussian(),
  data = data_fruit
)


# 3. PREDICTION GRID

temp_seq <- seq(min(data_fruit$temp_z),
                max(data_fruit$temp_z),
                length.out = 100)

precip_seq <- seq(min(data_fruit$precip_z),
                  max(data_fruit$precip_z),
                  length.out = 100)

newdat <- expand.grid(
  temp_z = temp_seq,
  precip_z = precip_seq,
  site = levels(data_fruit$site)
)

# Random-effect placeholders
newdat$family  <- data_fruit$family[1]
newdat$genus   <- na.omit(data_fruit$genus)[1]
newdat$species <- data_fruit$species[1]
newdat$year    <- data_fruit$year[1]


# 4. PREDICTIONS

newdat$pred <- predict(
  mod3b,
  newdata = newdat,
  type = "response",
  re.form = NA
)


# 5. PANEL A: OBSERVED DATA

pBfd <- ggplot(data_fruit,
               aes(x = temp_z, y = duration_doy, color = site)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic(base_size = 13) +
  labs(
    x = "Temperature (standardised)",
    y = "Observed fructification duration (days)",
    title = "B: Observed climate–duration relationship"
  )


# 6. PANEL B: MODEL SURFACE

pCfd <- ggplot(newdat,
               aes(x = temp_z, y = precip_z, z = pred)) +
  geom_contour_filled(bins = 12) +
  facet_wrap(~site, nrow = 1) +
  scale_fill_viridis_d(option = "C", name = "Duration (days)") +
  theme_classic(base_size = 13) +
  labs(
    x = "Temperature (standardised)",
    y = "Precipitation (standardised)",
    title = "C: Climate interaction surface (fructification duration, not significant)"
  ) +
  theme(strip.text = element_text(face = "bold"))


# 7. PANEL C: EFFECT SIZES

coef_df <- data.frame(
  term = rownames(coefs),
  estimate = coefs[,1],
  se = coefs[,2]
)

coef_df <- coef_df %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  )


pAfd <- ggplot(coef_df,
               aes(x = reorder(term, estimate),
                   y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)+
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Standardised effect size (β ± SE)",
    title = "A: Climate effect sizes (fructification duration)"
  )



# 8. COMBINE FIGURE

final_fig_fruct_dur <- (pAfd | pBfd) / pCfd +
  plot_layout(heights = c(2, 1))

final_fig_fruct_dur


# 9. SAVE FIGURE

ggsave(
  "Figure_fructification_duration_synthesis.png",
  final_fig_fruct_dur,
  width = 12,
  height = 9,
  dpi = 600
)

# COMBINING ALL RESULTS IN A FRAMEWORK
library(DiagrammeR)

grViz("

digraph FruitPhenology {

graph [
  layout = dot,
  rankdir = TB,
  fontsize = 24,
  labelloc = t,
  label = 'Climate-driven shifts in tropical fruit phenology'
]

node [
  shape = box,
  style = 'rounded,filled',
  fontname = Helvetica,
  fontsize = 18,
  color = black
]


# Climate drivers


Climate [
label = 'Climate change\n\nTemperature ↑\nPrecipitation ↑',
fillcolor = '#DCEAF7'
]


# Fructification


F_onset [
label = 'Fructification onset\nEarlier onset',
fillcolor = '#E8F5E9'
]

F_end [
label = 'Fructification end\nSlightly earlier',
fillcolor = '#E8F5E9'
]

F_duration [
label = 'Fructification duration\nLonger duration',
fillcolor = '#C8E6C9'
]


# Fruit drop


D_onset [
label = 'Fruit-drop onset\nEarlier onset',
fillcolor = '#FFF3E0'
]

D_end [
label = 'Fruit-drop end\nSlightly earlier',
fillcolor = '#FFF3E0'
]

D_duration [
label = 'Fruit-drop duration\nLonger duration',
fillcolor = '#FFE0B2'
]


# Ecological consequence


Availability [
label = 'Earlier and more prolonged\nfruit availability window',
fillcolor = '#FFF9C4'
]

Ecology [
label = 'Potential consequences\n\n• Frugivore resource dynamics\n• Seed dispersal timing\n• Plant–animal interactions\n• Forest regeneration',
fillcolor = '#F3E5F5'
]


# Climate → Fructification


Climate -> F_onset [
label='Temp: −125.5***\nPrecip: −69.6***',
color='forestgreen',
penwidth=3
]

Climate -> F_end [
label='Temp: −44.0***\nPrecip: NS',
color='darkgreen',
penwidth=2
]

Climate -> F_duration [
label='Temp: +82.5***\nPrecip: +47.1**',
color='forestgreen',
penwidth=3
]


# Climate → Fruit Drop


Climate -> D_onset [
label='Temp: −112.4***\nPrecip: −55.1***',
color='darkorange3',
penwidth=3
]

Climate -> D_end [
label='Temp: −37.7**\nPrecip: NS',
color='darkorange3',
penwidth=2
]

Climate -> D_duration [
label='Temp: +73.9***\nPrecip: +59.0***',
color='darkorange3',
penwidth=3
]


# Internal pathways


F_onset -> F_end [arrowhead=normal]
F_end -> F_duration [arrowhead=normal]

D_onset -> D_end [arrowhead=normal]
D_end -> D_duration [arrowhead=normal]


# Consequences


F_duration -> Availability [penwidth=2]
D_duration -> Availability [penwidth=2]

Availability -> Ecology [penwidth=3]

}
")


# EXPORT FIGURE

library(DiagrammeRsvg)
library(rsvg)

library(DiagrammeR)

g <- grViz("
digraph {
Climate [label='Climate']
Phenology [label='Phenology']
Climate -> Phenology
}
")

g

svg <- export_svg(g)

writeLines(svg, "Phenology_framework.svg")

rsvg_pdf(
  "Phenology_framework.svg",
  "Phenology_framework.pdf"
)

rsvg_png(
  "Phenology_framework.svg",
  "Phenology_framework.png",
  width = 5000,
  height = 4000
)