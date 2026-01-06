############################################
# 1. Chargement des librairies nécessaires
############################################

# Packages
library(readxl)
library(R2jags)
library(dplyr)
library(tidyr)
library(tidyverse)
library(MCMCvis)

############################################
# 2. Importation et préparation des données
############################################

# Importation du fichier de pêche (données terrain)
peche <- read_excel("datasets/Pêche_IOTA_LAPITXURI_21oct2025.xlsx", 
                                              sheet = "data")
View(peche)

# Suppression d'une colonne inutile (coquille dans les données)
peche <- peche[,-5]

# Conversion de la variable poids en numérique
peche$poids <- as.numeric(peche$poids)

############################################
# 3. Sélection des données d’intérêt
############################################

# Création d'un sous-jeu de données contenant uniquement les truites (TRF)
trout <- subset(peche, espece == "TRF")

colnames(trout)[which(colnames(trout) == 'numero_poisson')] <- 'num_poisson'
colnames(trout)[which(colnames(trout) == 'numero_passage')] <- 'num_passage'

############################################
# 4. Construction des données de capture
############################################

# Nombre de poissons capturés par passage et par station
capture <- aggregate(
  num_poisson ~ num_passage + station,
  FUN = length,
  data = trout
)
colnames(capture) <- c("passage", "station", "C")

# Extraction des surfaces échantillonnées par station
surfaces <- aggregate(surface ~ station, FUN = unique, data = trout)
colnames(surfaces) <- c("station", "S")

# Fusion des données de capture et de surface
capture_with_surface <- merge(capture, surfaces, by = "station")

# Mise en forme large : une colonne par passage
capture_wide <- capture_with_surface %>%
  pivot_wider(
    names_from = passage,
    values_from = C
  )

# Suppression de la colonne surface pour créer la matrice de captures
capture_matrix <- capture_wide %>%
  select(-S)

# Conversion des colonnes de captures en numérique
capture_matrix[, 2:3] <- lapply(capture_matrix[, 2:3], as.numeric)

# Création de la matrice C (stations × passages)
C_matrix <- as.matrix(capture_matrix[, -1])
rownames(C_matrix) <- capture_matrix$station

# Vecteur des surfaces par station
S_vector <- capture_wide$S

############################################
# 5. Données transmises à JAGS
############################################

# Cas du modèle 1 : probabilité de capture fixée
data <- list(
  C = C_matrix,
  S = S_vector,
  pi = 0.5
)


############################################
# 6. MODÈLE 1 : pi fixée
############################################

m1 <- function() {
  
  # Boucle sur les stations (1 = amont, 2 = aval)
  for (k in 1:2) {
    
    # Processus écologique : abondance latente
    N[k] ~ dpois(d[k] * S[k])
    
    # Processus d'observation : captures successives
    C[k,1] ~ dbin(pi, N[k])
    R1[k] <- N[k] - C[k,1]
    C[k,2] ~ dbin(pi, R1[k])
    
    # Prior sur la densité
    d[k] ~ dgamma(1, 0.001)
  }
  
  # Quantité dérivée : différence de densité
  diff_d <- d[2] - d[1]
  prob_diff <- step(diff_d)
}

# Valeurs initiales
inits <- function() {
  list(d = rgamma(2, 1, 0.001))
}

# Paramètres surveillés
parameters <- c("N", "d", "diff_d", "prob_diff")

# Exécution du modèle
j1 <- jags(
  data = data,
  parameters.to.save = parameters,
  model.file = m1,
  n.chains = 3,
  inits = inits,
  n.iter = 50000,
  n.burnin = 1000,
  n.thin = 10
)
print(j1)


############################################
# 7. MODÈLE 2 : pi estimée (commune)
############################################

data2 <- list(
  C = C_matrix,
  S = S_vector
)

m2 <- function() {
  
  # Prior sur la probabilité de capture
  pi ~ dbeta(20,20)
  
  for (k in 1:2) {
    
    N[k] ~ dpois(d[k] * S[k])
    
    C[k,1] ~ dbin(pi, N[k])
    R1[k] <- N[k] - C[k,1]
    C[k,2] ~ dbin(pi, R1[k])
    
    d[k] ~ dgamma(1, 0.001)
  }
  
  diff_d <- d[2] - d[1]
  prob_diff <- step(diff_d)
}

parameters2 <- c("pi", "N", "d", "diff_d", "prob_diff")

j2 <- jags(
  data = data2,
  parameters.to.save = parameters2,
  model.file = m2,
  n.chains = 3,
  inits = inits,
  n.iter = 50000,
  n.burnin = 1000,
  n.thin = 10
)
print(j2)


############################################
# 8. MODÈLE 3 : pi différente entre passages
############################################

m3 <- function() {
  
  # Probabilités de capture par passage
  pi1 ~ dbeta(20,20)
  pi2 ~ dbeta(20,20)
  
  for (k in 1:2) {
    
    N[k] ~ dpois(d[k] * S[k])
    
    C[k,1] ~ dbin(pi1, N[k])
    R1[k] <- N[k] - C[k,1]
    C[k,2] ~ dbin(pi2, R1[k])
    
    d[k] ~ dgamma(1, 0.001)
  }
  
  # Quantités dérivées
  diff_d <- d[2] - d[1]
  prob_d <- step(diff_d)
  
  diff_pi <- pi1 - pi2
  prob_pi <- step(diff_pi)
}

parameters3 <- c("pi1", "pi2", "diff_pi", "prob_pi", "N", "d", "diff_d")

j3 <- jags(
  data = data2,
  parameters.to.save = parameters3,
  model.file = m3,
  n.chains = 3,
  inits = inits,
  n.iter = 50000,
  n.burnin = 1000,
  n.thin = 10
)
print(j3)


############################################
# 9. Sélection du modèle
############################################

# Comparaison des DIC
j1$BUGSoutput$DIC
j2$BUGSoutput$DIC
j3$BUGSoutput$DIC


############################################
# 10. Diagnostic de convergence
############################################

# Résumé convergence (Rhat, n.eff)
print(j3)

# Résumé numérique avec IC
MCMCvis::MCMCsummary(
  j3,
  params = c("pi1", "pi2", "d", "N", "diff_d", "diff_pi"),
  probs = c(0.025, 0.975)
)

# Traces MCMC
MCMCvis::MCMCtrace(j3, params = c("d"), pdf = FALSE)
MCMCvis::MCMCtrace(j3, params = c("N"), pdf = FALSE)
MCMCvis::MCMCtrace(j3, params = c("pi1", "pi2"), pdf = FALSE)


############################################
# 11. Visualisation des postérieurs
############################################

# Densités postérieures (MCMCvis)
MCMCvis::MCMCplot(
  j3,
  params = c('d'),
  col = c("#1F78B4", "#E31A1C"),
  rank = TRUE
)

# Extraction des valeurs a posteriori
post_d <- j3$BUGSoutput$sims.list$d

df_d <- tibble(
  d = c(post_d[,1], post_d[,2]),
  station = rep(c("Amont", "Aval"), each = nrow(post_d))
)

# Densités a posteriori avec ggplot2
ggplot(df_d, aes(x = d, fill = station, color = station)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = c("Amont" = "#1F78B4", "Aval" = "#E31A1C")) +
  scale_color_manual(values = c("Amont" = "#1F78B4", "Aval" = "#E31A1C")) +
  labs(
    x = "Densité (individus / m²)",
    y = "Densité de probabilité",
    title = "Distributions a posteriori de la densité",
    fill = "Station",
    color = "Station"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    legend.position = "top"
  )
  
