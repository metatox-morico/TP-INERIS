# Exemple d'analyse dose réponse avec le modele de Hill

## Chargement des packages
# install.packages("drc")
library(drc)

## Importation des donnees
#data_croissance <- data.frame(Concentration = rep(c(0, 0.067, 0.18, 0.48, 1.3, 3.5), 
#                                                  times = c(6, 3, 3, 3, 3, 3)),
#                              Taux_Croissance = c(1.748, 1.745, 1.655, 1.757, 1.786, 1.578, 1.727, 1.667, 1.739, 1.599, 1.491, 1.728, 1.147, 1.280, 1.423, 0.492, 0.820, 0.558, -0.144, 0.267, 0.175))

data_croissance <- read.table("data_croissance.txt", header = TRUE, sep = "\t")
data_croissance

## Modelisation dose reponse
?drm
getMeanFunctions()
res_LL_croissance <- drm(Taux_Croissance ~ Concentration, fct = LL.3(), data = data_croissance)
plot(res_LL_croissance, type = "all", broken = TRUE)
modelFit(res_LL_croissance)

res_LL4_croissance <- drm(Taux_Croissance ~ Concentration, fct = LL.4(), data = data_croissance)
plot(res_LL4_croissance, type = "all", broken = TRUE)

modelFit(res_LL4_croissance)

## Choix de modele
mselect(res_LL_croissance, fctList = list(LL.4(), LL.5()), nested = TRUE)

## Calcul des doses effectives 10%
?ED
ED(res_LL_croissance, respLev = 10, interval = "delta")
ED(res_LL_croissance, respLev = 0.1, interval = "delta")


res_LL2_croissance <- drm(Taux_Croissance ~ Concentration, fct = LL2.3(), data = data_croissance)
ED(res_LL2_croissance, respLev = 10, interval = "delta")
ED_res <- ED(res_LL2_croissance, respLev = 10, interval = "delta")
names(ED_res)
exp(ED_res[1])
exp(ED_res[3])
exp(ED_res[4])
