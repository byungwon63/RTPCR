# set working directory
setwd("./Simulation")
# install the functions for simulation
source("./Simulation_functions.r")

# Simulation
# Original p-values (large sample normal approximation)

## 1. null: error distribution - normal; sample size - 12
setting <- data.frame(Rep = 1000,
                      n_control = 3,
                      n_treat = 3,
                      n_repeat = 2,
                      beta0 = 0,
                      beta1 = 0,
                      sigmab = 0.5,
                      sigmae = 0.5,
                      err.dist = "Normal")

set.seed(1)
Simulation_normal_null_n12 <- Simulation_test_stat_distribution(setting)
Simulation_normal_null_n12_perm <- Simulation_permutation_test(setting)
  
write.csv(Simulation_normal_null_n12, "Simul_Normal_null_n12.csv")
write.csv(Simulation_normal_null_n12_perm, "Simul_Normal_null_n12_perm.csv")


## 2. null: error distribution - normal; sample size - 24
setting2 <- data.frame(Rep = 1000,
                      n_control = 6,
                      n_treat = 6,
                      n_repeat = 2,
                      beta0 = 0,
                      beta1 = 0,
                      sigmab = 0.5,
                      sigmae = 0.5,
                      err.dist = "Normal")

#Simulation_normal_null_n24 <- Simulation_test_stat_distribution(setting2)
Simulation_normal_null_n24_perm <- Simulation_permutation_test(setting2)

#write.csv(Simulation_normal_null_n24, "Simul_Normal_null_n24.csv")
write.csv(Simulation_normal_null_n24_perm, "Simul_Normal_null_n24_perm.csv")


## 3. null: error distribution - normal; sample size - 48
setting3 <- data.frame(Rep = 1000,
                      n_control = 12,
                      n_treat = 12,
                      n_repeat = 2,
                      beta0 = 0,
                      beta1 = 0,
                      sigmab = 0.5,
                      sigmae = 0.5,
                      err.dist = "Normal")

#Simulation_normal_null_n48 <- Simulation_test_stat_distribution(setting3)
Simulation_normal_null_n48_perm <- Simulation_permutation_test(setting3)

#write.csv(Simulation_normal_null_n48, "Simul_Normal_null_n48.csv")
write.csv(Simulation_normal_null_n48_perm, "Simul_Normal_null_n48_perm.csv")


## 4. null: error distribution - Laplace; sample size - 12
setting4 <- data.frame(Rep = 1000,
                      n_control = 3,
                      n_treat = 3,
                      n_repeat = 2,
                      beta0 = 0,
                      beta1 = 0,
                      sigmab = 0.5,
                      sigmae = 1.2,
                      err.dist = "Laplace")

#Simulation_laplace_null_n12 <- Simulation_test_stat_distribution(setting4)
Simulation_laplace_null_n12_perm <- Simulation_permutation_test(setting4)

#write.csv(Simulation_laplace_null_n12, "Simul_Laplace_null_n12.csv")
write.csv(Simulation_laplace_null_n12_perm, "Simul_Laplace_null_n12_perm.csv")


## 5. null: error distribution - Laplace; sample size - 24
setting5 <- data.frame(Rep = 1000,
                      n_control = 6,
                      n_treat = 6,
                      n_repeat = 2,
                      beta0 = 0,
                      beta1 = 0,
                      sigmab = 0.5,
                      sigmae = 1.2,
                      err.dist = "Laplace")

#Simulation_laplace_null_n24 <- Simulation_test_stat_distribution(setting5)
Simulation_laplace_null_n24_perm <- Simulation_permutation_test(setting5)

#write.csv(Simulation_laplace_null_n24, "Simul_Laplace_null_n24.csv")
write.csv(Simulation_laplace_null_n24_perm, "Simul_Laplace_null_n24_perm.csv")


## 6. null: error distribution - Laplace; sample size - 48
setting6 <- data.frame(Rep = 1000,
                      n_control = 12,
                      n_treat = 12,
                      n_repeat = 2,
                      beta0 = 0,
                      beta1 = 0,
                      sigmab = 0.5,
                      sigmae = 1.2,
                      err.dist = "Laplace")

#Simulation_laplace_null_n48 <- Simulation_test_stat_distribution(setting6)
Simulation_laplace_null_n48_perm <- Simulation_permutation_test(setting6)

#write.csv(Simulation_laplace_null_n48, "Simul_Laplace_null_n48.csv")
write.csv(Simulation_laplace_null_n48_perm, "Simul_Laplace_null_n48_perm.csv")

# Summary
colMeans(Simulation_normal_null_n12 <= 0.05)
colMeans(Simulation_normal_null_n24 <= 0.05)
colMeans(Simulation_normal_null_n48 <= 0.05)
colMeans(Simulation_laplace_null_n12 <= 0.05)
colMeans(Simulation_laplace_null_n24 <= 0.05)
colMeans(Simulation_laplace_null_n48 <= 0.05)

colMeans(Simulation_normal_null_n12_perm <= 0.05)
colMeans(Simulation_normal_null_n24_perm <= 0.05)
colMeans(Simulation_normal_null_n48_perm <= 0.05)
colMeans(Simulation_laplace_null_n12_perm <= 0.05)
colMeans(Simulation_laplace_null_n24_perm <= 0.05)
colMeans(Simulation_laplace_null_n48_perm <= 0.05)


