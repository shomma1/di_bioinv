# data integration

# ==============================================================================
# date

set.seed(777)
date <- format(Sys.time(), "%y%m%d-%H%M")
stop()

# ==============================================================================
# libraries

library(tidyverse)
library(INLA)

# ==============================================================================
# input

# dataset <- read.csv("./0_data/data_240715.csv")
dataset <- read.csv("./0_data/data.csv")
# dataset <- dataset %>%
#   filter(!port %in% c("OFT", "KRE", "HRR"))

model_list <- read.csv("./0_data/model_info.csv")

# ==============================================================================
# create "data" for analysis 
# =============================================================================

# number of interest regions
n_port <- n <- nrow(dataset)

## data processing for INLA

data <- list()
sp_list <- c("invicta", "geminata")
for (i in 1:2) {
  sp <- sp_list[[i]]
  Y_sur <- dataset[,paste0("survey_", sp)] # port survey
  Y_inc <- dataset[,paste0("incident_", sp)] # incident
  Y <- matrix(c(Y_sur, rep(NA, n), rep(NA, n), Y_inc), ncol = 2) # joint

  # explanatory variable: number of import container
  n_container <- dataset[,paste0("import_", sp)] # for separate
  n_container_joint <- matrix(
    c(n_container, rep(NA, n), rep(NA, n), n_container), ncol = 2) #for joint

  # explanatory variable: temperature
  min_temp_joint <- matrix(
    c(dataset$min_temp, rep(NA, n), rep(NA, n), dataset$min_temp), ncol = 2)
  min_temp <- dataset$min_temp # for separate
  
  data[[sp]] <- list(
    Y = Y, Y_sur = Y_sur, Y_inc = Y_inc, 
    I1 = c(rep(1,n),rep(NA,n)),
    I2 = c(rep(NA,n),rep(1,n)),
    id_joint = 1:(2*n) ,
    id_sep = 1:n,
    id_data_source = c(rep(1,n), rep(2,n)),
    n_survey = dataset$n_survey,
    n_container_joint = n_container_joint, 
    n_container = n_container,
    # n_container_offset = c(rep(NA, n), n_container),
    min_temp_joint = min_temp_joint,
    min_temp = min_temp)
}

data[["joint.sp"]] <- list(
  Y = c(dataset$survey_invicta, dataset$survey_geminata),
  I1 = c(rep(1,n),rep(NA,n)),
  I2 = c(rep(NA,n),rep(1,n)),
  n_container_i = c(dataset$import_invicta, rep(NA,n)),
  n_container_g = c(rep(NA,n),dataset$import_geminata),
  id_data_source = c(rep(1,n), rep(2,n)),
  id_port = rep(1:n,2),
  id_joint = 1:(n_port*2),
  n_survey = rep(dataset$n_survey,2),
  min_temp_i = c(dataset$min_temp, rep(NA, n_port)),
  min_temp_g = c(rep(NA, n_port), dataset$min_temp)
)

# ==============================================================================
# create W

# spatial weight matrix
W_org  <- as.matrix(read.csv("./0_data/matrix_container.csv")[,-1])
rownames(W_org) <- colnames(W_org)
W_org <- W_org[unique(dataset$port), unique(dataset$port)]

W <- W_org/max(W_org)
W <- round(W, 3)

### quick check ###
if ( any(rowSums(W) == 0) |
     any(W != t(W))) { 
  print("check W")
  stop()
}

# ==============================================================================
# define rgeneric model

# read iMCAR function
source("./1_script/functions/PCAR_Leroux.R") 
source("./1_script/functions/PMCAR_Leroux.R")

# Spacial random effect
# Univariate Leroux
leroux_car <- inla.rgeneric.define(PCAR_Leroux, W=W)

# Multivariate Leroux
k <- 2
mpcar_leroux <- inla.rgeneric.define(PMCAR_Leroux, W=W, k = k)

###============================================================================
# define "run.inla" and set family

run.inla <- function(formula, d, family) {
  res <- inla(
    formula = formula,
    data = d,
    family = family,
    Ntrials = n_survey,
    # E = E,
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
    verbose = T)
  return(res)
}

# family 
fam.joint.dt <- c("binomial", "poisson")
fam.sur <- "binomial"
fam.po <- "poisson"

###=============================================================================
# Bayesian inference
###=============================================================================

res.list <- list()
# ==============================================================================
# ===== Separate model: Survey 
# ==============================================================================

# glm
formula.fixed <- Y_sur ~
  log(n_container+1) + min_temp
  #n_container + min_temp

name <- model_list$model.name[[1]]
res.list[[name]] <- run.inla(formula.fixed, data[["invicta"]], fam.sur)
name <- model_list$model.name[[5]]
res.list[[name]] <- run.inla(formula.fixed, data[["geminata"]], fam.sur)

# Leroux
formula.pcar <- update(
  formula.fixed, . ~ . +
    f(id_sep, model = leroux_car)
)

name <- model_list$model.name[[2]]
res.list[[name]] <- run.inla(formula.pcar, data[["invicta"]], fam.sur)
name <- model_list$model.name[[6]]
res.list[[name]] <- run.inla(formula.pcar, data[["geminata"]], fam.sur)

# ======================================================
# ===== Separate model: Presence-Only 
# ======================================================

# glm
formula.fixed <- update(
  formula.fixed, Y_inc ~ . )

name <- model_list$model.name[[3]]
res.list[[name]] <- run.inla(formula.fixed, data[["invicta"]], fam.po) 
name <- model_list$model.name[[7]]
res.list[[name]] <- run.inla(formula.fixed, data[["geminata"]], fam.po) 

# pcar
formula.pcar <- update(
  formula.pcar, Y_inc ~ . )

name <- model_list$model.name[[4]]
res.list[[name]] <- run.inla(formula.pcar, data[["invicta"]], fam.po)
name <- model_list$model.name[[8]]
res.list[[name]] <- run.inla(formula.pcar, data[["geminata"]], fam.po)

# ==============================================================================
# ===== Joint model: datatype
# ==============================================================================

# fix
formula.fixed <- Y ~ 
  -1 + I1 + I2 +
  log(n_container_joint+1) + min_temp_joint

# Leroux CAR model
formula.mpcar <- update(
  formula.fixed, . ~ . +
    f(id_joint, model = mpcar_leroux)
)
name <- model_list$model.name[[9]]
res.list[[name]] <- run.inla(formula.mpcar, data[["invicta"]],  fam.joint.dt)
name <- model_list$model.name[[10]]
res.list[[name]] <- run.inla(formula.mpcar, data[["geminata"]],  fam.joint.dt)

# ==============================================================================
# ===== Joint model: species
# ==============================================================================

# fix
formula.fixed <- Y ~ 
  -1 + I1 + I2 +
  log(n_container_i+1) + log(n_container_g+1) +
  # n_container_i + n_container_g +
  min_temp_i + min_temp_g

# Leroux CAR
formula.mpcar <- update(
  formula.fixed, . ~ . +
  f(id_joint, model = mpcar_leroux)
)

name <- model_list$model.name[[11]]
res.list[[name]] <- run.inla(formula.mpcar, data[["joint.sp"]],  fam.sur)

###============================================================================
# Result
###============================================================================

# reload
# date <- "240715-0034"
# file <- paste0(
#   "./2_output/output_Data/res_inla_", date, ".RDS")
# res.list <- readRDS(file)

# Information Criterion ---------------------------------------------------------

ic <- cbind(model_list, DIC = NA, WAIC = NA, 
            local.DIC1 = NA, local.DIC2 = NA, 
            local.WAIC1 = NA, local.WAIC2 = NA)
for (i in 1:nrow(model_list)) {
  res <- res.list[[ model_list[[i, "model.name"]] ]]
  type <- model_list[[i, "type"]]
  
  tmp <- c(DIC = res$dic$dic, WAIC = res$waic$waic)
  ic[i, c("DIC", "WAIC")] <- tmp
  
  if ( type == "joint") {
    tmp <- rbind(
      DIC = tapply(res$dic$local.dic, res$.args$data$id_data_source, sum), 
      WAIC = tapply(res$waic$local.waic, res$.args$data$id_data_source, sum))
    ic[i, c("local.DIC1", "local.DIC2", "local.WAIC1", "local.WAIC2")] <- 
      c(tmp["DIC", ], tmp["WAIC", ])
    
  }
}

# file <- paste0(
#   "./2_output/output_Data/ic_", date, ".csv")
# write.csv(ic, file, row.names = F)

# Information Criterion ---------------------------------------------------------
# WAIC

ic_table <- model_list %>% 
  filter(species %in% c("invicta", "Joint-sp") & source != "PO") 

ic_table <- model_list %>%
  filter(species %in% c("geminata", "Joint-sp") & source != "PO") %>%
  rbind(ic_table, .) %>%
  mutate(species_2 = rep(c("invicta", "geminata"), each = 4))

for (i in 1:nrow(ic_table)) {
  
  model <- ic_table[i, ]
  
  if (model$source == "Joint-dt") {
    ic_table[i, "Survey"] <- 
      ic[ic$model.name == model$model.name, "local.WAIC1"]
    ic_table[i, "PO"] <- 
      ic[ic$model.name == model$model.name, "local.WAIC2"]
  } else if (model$species == "Joint-sp") {
    
    # "Joint-sp"
    if (model$species_2 == "invicta") { # invicta
      ic_table[i, "Survey"] <- 
        ic[ic$model.name == model$model.name, "local.WAIC1"]
    } else if (model$species_2 == "geminata") { # geminata
      ic_table[i, "Survey"] <- 
        ic[ic$model.name == model$model.name, "local.WAIC2"]
    }
    
  } else {
    ic_table[i, "Survey"] <- ic[ic$model.name == ic_table$model.name[i], "WAIC"]
    
    ic_table[i, "PO"] <- 
      ic[( ic$source == "PO" & 
             ic$model == model$model &
             ic$species == model$species_2) , "WAIC"]
    
  }
}
ic_table

file <- paste0(
  "./2_output/output_Data/ic_table_WAIC", date, ".csv")
write.csv(ic_table, file, row.names = F)

# Information Criterion ---------------------------------------------------------
# DIC

ic_table <- model_list %>% 
  filter(species %in% c("invicta", "Joint-sp") & source != "PO") 

ic_table <- model_list %>%
  filter(species %in% c("geminata", "Joint-sp") & source != "PO") %>%
  rbind(ic_table, .) %>%
  mutate(species_2 = rep(c("invicta", "geminata"), each = 4))

for (i in 1:nrow(ic_table)) {
  
  model <- ic_table[i, ]
  
  if (model$source == "Joint-dt") {
    ic_table[i, "Survey"] <- 
      ic[ic$model.name == model$model.name, "local.DIC1"]
    ic_table[i, "PO"] <- 
      ic[ic$model.name == model$model.name, "local.DIC2"]
  } else if (model$species == "Joint-sp") {
    
    # "Joint-sp"
    if (model$species_2 == "invicta") { # invicta
      ic_table[i, "Survey"] <- 
        ic[ic$model.name == model$model.name, "local.DIC1"]
    } else if (model$species_2 == "geminata") { # geminata
      ic_table[i, "Survey"] <- 
        ic[ic$model.name == model$model.name, "local.DIC2"]
    }
  
  } else {
    ic_table[i, "Survey"] <- ic[ic$model.name == ic_table$model.name[i], "DIC"]
  
    ic_table[i, "PO"] <- 
        ic[( ic$source == "PO" & 
             ic$model == model$model &
             ic$species == model$species_2) , "DIC"]

  }
}
ic_table

file <- paste0(
  "./2_output/output_Data/ic_table_", date, ".csv")
write.csv(ic_table, file, row.names = F)

# ==============================================================================

source("./1_script/functions/hyperparams.R")
hyper <- trans.hyper(model_list)

margs <- hyper$marges
df.hyper <- hyper$df.hyper

write.csv(df.hyper, 
          paste0("2_output/output_data/hyper/hyper_",date,".csv"))

# ###============================================================================
# # save
# ###============================================================================
# # save inla output
# 
file <- paste0(
  "./2_output/output_Data/res_inla_", date, ".RDS")
saveRDS(res.list, file)

# # summary ==========================
file <- paste0(
  "./2_output/output_Data/summary_", date, ".RDS")
sink(file, append = F)
print(paste(date, "analysis.R"))
print("=================================================================")
print(model_list)
print("=================================================================")
# 
lapply(res.list, function (x) {summary(x)})
sink()

###
###
print(date)
print(df.hyper[6:12])
print(ic_table)

if (!exists("result")) { result <- list() }
result[[date]] <- ic_table
#======================
#### End of Script ####
#======================