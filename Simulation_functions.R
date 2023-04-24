library(gee)
library(gtools)
library(tidyverse)
library(geesmv)
library(extraDistr)
#library(OOmisc)
library(lme4)

## Function for generating simulation data
##### Requires 
##### 1. sample size: n1, n2 for treatment and control groups
##### 2. nrepeat: number of technical replication
##### 3. beta0, beta1: true fixed effect parameters; 
#####                  beta1 is essentially required for true mean difference
##### 4. sigmab: scale parameter for random effect (effect of technical replication)
##### 5. sigmae: scale parameter for error 
#####            if "Normal", standard deviation
#####            if "Laplace", scale parameter (variance = 2*sigmae^2)
##### 6. err.dist: distribution of error
#####              - Normal (default)
#####              - Laplace <- requires "extraDistr" (or "OOmisc") package 
#####                           for random data generation

data.generation <- function(n1, n2, nrepeat, beta0, beta1, sigmab, sigmae, err.dist = "Normal"){
  n = n1 + n2
  N = n * nrepeat
  x = rbind( matrix(1,nrow = n1, ncol = nrepeat) , 
             matrix(0,nrow = n2, ncol = nrepeat))
  y = beta0 + beta1 * x + 
    rnorm(n, sd = sigmab) %>% outer(rep(1, nrepeat), "*")
  if (err.dist == "Normal"){
    y = y + rnorm(N, sd = sigmae) %>% matrix(nrow = n)
  }else if (err.dist == "Laplace"){
    y = y + rlaplace(N, sigmae) %>% matrix(nrow = n)
  }else{
    print("The specified error distribution is not appropriate for this function!")
    break;
  }
  colnames(y) = c(paste0("Rep", 1:nrepeat))
  rownames(y) = 1:n
  treatment = c( rep(1,n1), rep(2,n2) ) 
  y.dat = data.frame( y = as.vector(y),
                      treatment = rep(treatment, nrepeat),
                      sample = rep(rownames(y), times = nrepeat))
  return(Y = list(y, treatment, y.dat))
}

t.ind.stat = function(y, perm.id = 1:nrow(y), treatment, nrepeat){
  t.test(as.vector(y[perm.id,]) ~ 
           rep(treatment, nrepeat), equal.var = T, alternative = "greater")$statistic  
  
}

t.avg.stat = function(y, perm.id = 1:nrow(y), treatment){
  t.test(rowMeans(y[perm.id,]) ~ treatment, 
         equal.var = T, alternative = "greater")$statistic  
}

md.stat = function(y, perm.id = 1:nrow(y), treatment){
  # treatment 1 - treatment 2
  y.perm = y[perm.id,]
  mean( y.perm[treatment == 1]) - mean( y.perm[treatment == 2])    
}

lme.stat = function(y.dat, y, perm.id = 1:nrow(y), md.stat.){
  y.dat$y <- as.vector(y[perm.id,])
  lme.fit <- lmer(y ~ treatment + (1 | sample), data = y.dat) %>% summary()
  lme.stat <- lme.fit$coefficients[2,3]
  if (md.stat. > 0){
    lme.stat <- abs(lme.stat)
  }else{
    lme.stat <- -abs(lme.stat)
  }
}

wilcoxon.stat = function(y, perm.id = 1:nrow(y), treatment, nrepeat){
  wilcox.test(as.vector(y[perm.id,]) ~ 
                rep(treatment, nrepeat))$statistic
}

GEE.stat = function(y.dat, y, perm.id = 1:nrow(y), md.stat.){
  y.dat$y <- as.vector(y[perm.id,])
  gee.fit <- gee(y ~ treatment,
                 id = sample,
                 data = y.dat,
                 family = gaussian,
                 corstr = "exchangeable") %>%
    summary()
  gee.stat <- gee.fit$coefficients[2,5]
  
  # Mancl and DeRouen (2001)
  gee.stat.md <- gee.fit$coefficients[2,1] / 
    sqrt(GEE.var.md(y ~ treatment,
                    id = "sample",
                    family = gaussian,
                    data = y.dat,
                    corstr = "exchangeable")$cov.beta[2])
  # Fay and Graubard (2001) 
  gee.stat.fg <- gee.fit$coefficients[2,1] / 
    sqrt(GEE.var.fg(y ~ treatment,
                    id = "sample",
                    family = gaussian,
                    data = y.dat,
                    corstr = "exchangeable")$cov.beta[2])
  # Morel et al. (2003)
  gee.stat.mbn <- gee.fit$coefficients[2,1] / 
    sqrt(GEE.var.mbn(y ~ treatment,
                     id = "sample",
                     family = gaussian,
                     data = y.dat,
                     corstr = "exchangeable")$cov.beta[2])
  if (md.stat. > 0){
    GEE.stat <- abs(c(gee.stat, gee.stat.md, gee.stat.fg, gee.stat.mbn))
  }else{
    GEE.stat <- -abs(c(gee.stat, gee.stat.md, gee.stat.fg, gee.stat.mbn))
  }
  return(GEE.stat)
}

t.ind.pvalue = function(y, treatment, nrepeat){
  t.test(as.vector(y) ~ 
           rep(treatment, nrepeat), equal.var = T, alternative = "greater")$p.value
  
}

t.avg.pvalue = function(y, treatment){
  t.test(rowMeans(y) ~ treatment, 
         equal.var = T, alternative = "greater")$p.value  
}

lme.pvalue = function(y.dat, n1, n2, nrepeat, md.stat.){
  lme.fit <- lmer(y ~ treatment + (1 | sample), data = y.dat) %>% summary()
  if (md.stat. > 0){
    lme.stat <- abs(lme.fit$coefficients[2,3])
  }else{
    lme.stat <- -abs(lme.fit$coefficients[2,3])
  }
  lme.pvalue <- 1 - pt(lme.stat, df = 2*(n1 + n2)*nrepeat - (n1+n2) - nrepeat + 1)
}

wilcoxon.pvalue = function(y, treatment, nrepeat){
  wilcox.test(as.vector(y) ~ 
                rep(treatment, nrepeat), alternative = "greater")$p.value
}

GEE.pvalue = function(y.dat, md.stat.){
  gee.fit <- gee(y ~ treatment,
                 id = sample,
                 data = y.dat,
                 family = gaussian,
                 corstr = "exchangeable") %>%
    summary()
  if (md.stat. > 0){
    gee.stat <- abs(gee.fit$coefficients[2,5])
  }else{
    gee.stat <- -abs(gee.fit$coefficients[2,5])
  }
  gee.pvalue <- 1 - pnorm(gee.stat)
  
  # Mancl and DeRouen (2001)
  gee.stat.md <- gee.fit$coefficients[2,1] / 
    sqrt(GEE.var.md(y ~ treatment,
                    id = "sample",
                    family = gaussian,
                    data = y.dat,
                    corstr = "exchangeable")$cov.beta[2])
  if (md.stat. > 0){
    gee.stat.md <- abs(gee.stat.md)
  }else{
    gee.stat.md <- -abs(gee.stat.md)
  }
  gee.pvalue.md <- 1 - pnorm(gee.stat.md)
  # Fay and Graubard (2001) 
  gee.stat.fg <- gee.fit$coefficients[2,1] / 
    sqrt(GEE.var.fg(y ~ treatment,
                    id = "sample",
                    family = gaussian,
                    data = y.dat,
                    corstr = "exchangeable")$cov.beta[2])
  if (md.stat. > 0){
    gee.stat.fg <- abs(gee.stat.fg)
  }else{
    gee.stat.fg <- -abs(gee.stat.fg)
  }
  gee.pvalue.fg <- 1 - pnorm(gee.stat.fg)
  # Morel et al. (2003)
  gee.stat.mbn <- gee.fit$coefficients[2,1] / 
    sqrt(GEE.var.mbn(y ~ treatment,
                     id = "sample",
                     family = gaussian,
                     data = y.dat,
                     corstr = "exchangeable")$cov.beta[2])
  if (md.stat. > 0){
    gee.stat.mbn <- abs(gee.stat.mbn)
  }else{
    gee.stat.mbn <- -abs(gee.stat.mbn)
  }
  gee.pvalue.mbn <- 1 - pnorm(gee.stat.mbn)
  
  return(GEE.pvalue = c(gee.pvalue, gee.pvalue.md, gee.pvalue.fg, gee.pvalue.mbn))
}

Simulation_test_stat_distribution <- function(setting){
  attach(setting)
  pvalues.out <- matrix(0, ncol = 8, nrow = Rep)
  for (rr in 1:Rep){
    y.list <- data.generation(n_control, n_treat, n_repeat, 
                              beta0, beta1, 
                              sigmab, sigmae, 
                              err.dist)
    y <- y.list[[1]]
    y.dat <- y.list[[3]]
    treatment = y.list[[2]]
    md.stat_ <- md.stat(y, perm.id = 1:nrow(y), treatment)
    sink("worksave1.txt")
    test.pvalue = c( t.ind.pvalue(y, treatment = treatment, nrepeat = n_repeat), 
                     t.avg.pvalue(y, treatment = treatment), 
                     wilcoxon.pvalue(y, treatment = treatment, nrepeat = n_repeat),
                     lme.pvalue(y.dat, n1 = n_control, n2 = n_treat, nrepeat = n_repeat, md.stat_),
                     suppressMessages(GEE.pvalue(y.dat, md.stat_)) )
    sink()
    names(test.pvalue) = c("t_ind", "t_avg", "Wilcoxon", "Mixed",
                           "gee.wald", "gee.wald.md", "gee.wald.fg", "gee.wald.mbn")
    pvalues.out[rr,] <- test.pvalue
    
    if (rr %% 20 == 0){
      print(paste0("this is the ", rr, "th repetition"))
    }
  }
  detach(setting)
  colnames(pvalues.out) = c("t_ind", "t_avg", "Wilcoxon", "Mixed",
                            "gee.wald", "gee.wald.md", "gee.wald.fg", "gee.wald.mbn")
  return(pvalues.out)
}

Simulation_permutation_test <- function(setting){
  attach(setting)
  # permutation index
  perm.ids = gtools::combinations(n, n_control)
  if (nrow(perm.ids) > 2000){
    perm.ids <- perm.ids[sample(nrow(perm.ids), 2000, replace = FALSE),]
  }
  allperms_ids = apply(perm.ids, 1,  function(x){
    c(x, (1:n)[-x])
  }) %>% t()
  
  pvalues.out <- matrix(0, ncol = 9, nrow = Rep)
  for (rr in 1:Rep){
    y.list <- data.generation(n_control, n_treat, n_repeat, 
                              beta0, beta1, 
                              sigmab, sigmae,
                              err.dist)
    y <- y.list[[1]]
    y.dat <- y.list[[3]]
    treatment = y.list[[2]]
    md.stat_ <- md.stat(y, perm.id = 1:nrow(y), treatment)
    sink("worksave2.txt")
    test.stat = c( t.ind.stat(y, treatment = treatment, nrepeat = n_repeat), 
                   t.avg.stat(y, treatment = treatment), 
                   md.stat(y, treatment = treatment), 
                   wilcoxon.stat(y, treatment = treatment, nrepeat = n_repeat),
                   lme.stat(y.dat, y, md.stat. = md.stat_),
                   suppressMessages(GEE.stat(y.dat,y, md.stat. = md.stat_)) )
    names(test.stat) = c("t_ind", "t_avg", "md", "Wilcoxon", "Mixed",
                         "gee.wald", "gee.wald.md", "gee.wald.fg", "gee.wald.mbn")
    # permutation p-value
    perm.stat = apply(allperms_ids, 1, function(index){
      c( t.ind.stat(y, index, treatment = treatment, nrepeat = n_repeat),
         t.avg.stat(y, index, treatment = treatment), 
         md.stat(y, index, treatment = treatment), 
         wilcoxon.stat(y, index, treatment = treatment, nrepeat = n_repeat),
         lme.stat(y.dat, y, index, md.stat. = md.stat(y, index, treatment = treatment)),
         suppressMessages(GEE.stat(y.dat, y, index, md.stat. = md.stat(y, index, treatment = treatment))) )
    }) %>% t()
    sink()
    pvalues = vector(length = length(test.stat))
    for (j in 1:length(test.stat)){
      pvalues[j] = mean(perm.stat[,j] >= test.stat[j]) # the larger the extreme (one-sided test)
    }
    
    pvalues.out[rr,] <- pvalues
    
    if (rr %% 20 == 0){
      print(paste0("this is the ", rr, "th repetition"))
    }
  }
  detach(setting)
  colnames(pvalues.out) = c("t_ind", "t_avg", "md", "Wilcoxon", "Mixed",
                            "gee.wald", "gee.wald.md", "gee.wald.fg", "gee.wald.mbn")
  return(pvalues.out)
}
