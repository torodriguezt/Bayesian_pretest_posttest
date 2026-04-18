Sys.setenv("R_MAX_NUM_DLL"=9999999)
mtime<-Sys.time()
library(ALA)
library(dplyr)
library(purrr)
library(hypergeo)
library(ggplot2)
library(parallel)
library(gtools)
library(GA)
library(future)
#library(brms)
library(rstan)
#library(cmdstanr)
library(furrr)
library(MCMCpack)
options(mc.cores = 10)
datos<-tvsfp
attach(datos)
datos=mutate(datos,binTHKS=ifelse(THKS >= 3,1,0))

#modelo <- cmdstan_model('BBpost3.stan')
seed <- 123
n_bootstrap <- 1000

###Datos de 4 colegios, uno para cada combinación de tratamientos

##school and tv
# ##yy
datos=mutate(datos1,binTHKS=ifelse(THKS >= 3,1,0))%>% filter(school=="404")

datos_yy <- datos %>%
  filter(school.based=="yes" & tv.based == "yes")

X_yy=datos_yy %>% group_by(stage) %>% summarise(Bin=sum(binTHKS))
n_yy=datos_yy %>% group_by(stage) %>% summarise(n=n())

##school - no tv
##yn
datos=mutate(datos1,binTHKS=ifelse(THKS >= 3,1,0))%>% filter(school=="408")

datos_yn <- datos %>%
  filter(school.based=="yes" & tv.based == "no")

X_yn=datos_yn %>% group_by(stage) %>% summarise(Bin=sum(binTHKS))
n_yn=datos_yn %>% group_by(stage) %>% summarise(n=n())

##no school - tv
# ##ny
datos=mutate(datos1,binTHKS=ifelse(THKS >= 3,1,0))%>% filter(school=="508")

datos=mutate(datos1,binTHKS=ifelse(THKS >= 3,1,0))%>% filter(school=="409")

datos_ny <- datos1 %>% 
  filter(school.based=="no" & tv.based == "yes")

X_ny=datos_ny %>% group_by(stage) %>% summarise(Bin=sum(binTHKS))
n_ny=datos_ny %>% group_by(stage) %>% summarise(n=n())

##no school - no tv
# ##nn
datos=mutate(datos1,binTHKS=ifelse(THKS >= 3,1,0))%>% filter(school=="409")

datos_nn <- datos %>% 
  filter(school.based=="no" & tv.based == "no")

X_nn=datos_nn %>% group_by(stage) %>% summarise(Bin=sum(binTHKS))
n_nn=datos_nn %>% group_by(stage) %>% summarise(n=n())

##school and tv
n1=n_yy$n[1]
n2=n_yy$n[2]
# x1=X_yy$Bin[1]
# x2=X_yy$Bin[2]

##school - no tv
n1=n_yn$n[1]
n2=n_yn$n[2]
# x1=X_yn$Bin[1]
# x2=X_yn$Bin[2]

##no school -  tv
n1=n_ny$n[1]
n2=n_ny$n[2]
# x1=X_ny$Bin[1]
# x2=X_ny$Bin[2]

##no school - no tv
n1=n_nn$n[1]
n2=n_nn$n[2]
# x1=X_nn$Bin[1]
# x2=X_nn$Bin[2]

##Posterior bajo H
densBB_H_fucntor<-function(x1, x2, n1, n2, u, l, a, a0, a1, a2, i, use_log=FALSE) {
  
  densBB_H <- function(theta1) {
    log_f<-lgamma(n1+a)+lgamma(n2+a)-
      lgamma(x1[i]+a1)-lgamma(n1-x1[i]+a-a1)-
      lgamma(x2[i]+a2)-lgamma(n2-x2[i]+a-a2)-
      log(genhypergeo(U=u, L=l, check_mod=T, z=1))+
      (a1+x1[i]-1)*log(theta1)+(a2+a0+(n1-x1[i])-1)*log(1-theta1)+
      (a2+x2[i]-1)*log(theta1)+(a1+a0+(n2-x2[i])-1)*log(1-theta1)-
      (a)*log(1-theta1*theta1) 
    
    fun <- log_f
    if(!use_log) {
      fun<-exp(log_f)
    }
    return(fun)
    
  }
  
  return(densBB_H)
}

##Posterior
densBB_functor<-function(x1, x2, n1, n2, u, l, a, a0, a1, a2, i, use_log=FALSE) {
  densBB <- function(theta1, theta2) {
    log_f<-lgamma(n1+a)+lgamma(n2+a)-
      lgamma(x1[i]+a1)-lgamma(n1-x1[i]+a-a1)-
      lgamma(x2[i]+a2)-lgamma(n2-x2[i]+a-a2)-
      log(genhypergeo(U=u, L=l, check_mod=T, z=1))+
      (a1+x1[i]-1)*log(theta1)+(a2+a0+(n1-x1[i])-1)*log(1-theta1)+
      (a2+x2[i]-1)*log(theta2)+(a1+a0+(n2-x2[i])-1)*log(1-theta2)-
      (a)*log(1-theta1*theta2) 
    
    fun <- log_f
    if(!use_log) {
      fun<-exp(log_f)
    }
    
    return(fun)
    
  }
  return(densBB)
}

##Error tipo II
##Generación de valores de la apriori para aproximar la integral MC

##Prior no Informativa
alphas_opt=c(0.8373879, 0.8410984, 0.8053298)
a1<-alphas_opt[2]
a2<-alphas_opt[3]
a0<-alphas_opt[1]
a<-a1+a2+a0

alpha2<-alphas_opt[2]
alpha3<-alphas_opt[3]
alpha1<-alphas_opt[1]
alphas=alphas_opt

f2 <- function(X1,X2,theta1,theta2) {exp(dbinom(X1,n1,theta1,log=T)+dbinom(X2,n2,theta2,log=T))}

stan_model <- rstan::stan_model(file = 'BBpost3.stan')


priorB2<- function(theta) {
  
  if( any(theta > 1) | any(theta < 0) ) {
    return(-Inf)
  }
  
  a<-alphas[1]+alphas[2]+alphas[3]
  log_f<-lgamma(a)-
    lgamma(alphas[1])-lgamma(alphas[2])-
    lgamma(alphas[3])+
    (alphas[2]-1)*log(theta[1])+(alphas[3]+alphas[1]-1)*log(1-theta[1])+
    (alphas[3]-1)*log(theta[2])+(alphas[2]+alphas[1]-1)*log(1-theta[2])-
    (a)*log(1-theta[1]*theta[2])
  return(log_f)
}

theta.samp <- MCMCmetrop1R(priorB2, theta.init=c(0.3,0.3),
                           thin=10, mcmc=200000, burnin=5000,
                           # tune=c(1, 1),
                           verbose=500, logfun=TRUE,
                           force.samp=T,
                           optim.lower=c(0.01, 0.01),
                           optim.upper = c(0.99, 0.99),
                           optim.method = "L-BFGS-B")

# theta.samp[, 1] |> hist()
# theta.samp[, 2] |> hist()
# 
# theta.samp[, 1]|> acf(lag.max=1000)

Theta11<-theta.samp[, 1][seq(20,20000,40)]
Theta22<-theta.samp[, 2][seq(20,20000,40)]
N2=length(Theta11)
N1=length(Theta22)

integ11<-function(i,Theta11,Theta22,N1,k,j){
  x1<-rbinom(N1,n1,Theta11[j])
  x2<-rbinom(N1,n2,Theta22[j])
  #f2 <- function(X1,X2,theta1,theta2) {dbinom(X1,n1,theta1)*dbinom(X2,n2,theta2)}
  f_ver<-f2(x1[i],x2[i],Theta11[j], Theta22[j])
  
  
  ##########
  ####Biv beta binom posterior
  #####
  if(((n1+a)+(n2+a)) <= (a+(x1[i]+a1)+(x2[i]+a2))) {
    # browser()
    return(NA_real_)
    # return(0)
  }
  #####
  u<-c(a,x1[i]+a1,x2[i]+a2)
  l<-c(n1+a,n2+a)
  
  densBB_H <- densBB_H_fucntor(x1, x2, n1, n2, u, l, a, a0, a1, a2, i)
  # browser()
  ###MSUP
  GA_ev <- ga(type = "real-valued", 
              fitness = function(x) {
                densBB_H_fucntor(x1, x2, n1, n2, u, l, a, a0, a1, a2, i, TRUE)(x[1])
              },
              lower = c(0), upper = c(0.99), 
              popSize = 50, maxiter = 1000, run = 100,
              monitor = FALSE)
  
  #Supremo bajo H
  #summary(GA_ev)
  #
  # supremo_ga_ev <- summary(GA_ev)$solution[1,] |> as.numeric() 
  supremo_ga_ev <- GA_ev@solution[1,] |> as.numeric() 
  # 
  # if( N1 ==  200) {
  #   browser()
  # }
  
  # if(is.na(supremo_ga_ev) |> any()) {
  #   browser()
  # }
  
  ###########
  n=as.integer(c(n1,n2))
  X=c(x1[i],x2[i])
  
  P=length(X)
  
  stan_data<-list(P=P, X=X,n=n,alpha1=alpha1,alpha2=alpha2, alpha3=alpha3)
  
  ##compile model
  #model<-stan_model('BBpost3.stan')
  ##Simulación de theta1 y theta2 de la posteriori
  ##pass data to stan and run model
  #fit<- sampling(model, list(P=P, X=X,n=n,alpha1=alpha1,alpha2=alpha2, alpha3=alpha3), iter=200, chains=4)
  # fit<-stan(file = 'BBpost3.stan',
  #           data = stan_data,
  #           iter = 10000,
  #           # cores = 4,
  #           cores = 1,
  #           verbose = FALSE,
  #           open_progress = FALSE,
  #           refresh = FALSE)
  # #diagnose
  # #print(fit)
  # 
  # #graph
  # params<- rstan::extract(fit)
  
  # suppressMessages({
  
  #   suppressWarnings({
  #     fit <- modelo$sample(
  #       data = stan_data,
  #       # seed = 123,
  #       iter_warmup = 500,
  #       iter_sampling = 500,
  #       chains = 4,
  #       parallel_chains = 4,
  #       refresh = 0 # print update every 500 iters
  #     )})})
  
  # #Thetas<- cbind(params$Theta[,1],params$Theta[,2])
  # Thetas<- cbind(fit$draws("Theta")[,,1] |> as.vector() ,
  #                fit$draws("Theta")[,,2] |> as.vector())
  
  # Variational Bayes in RStan
  sink(nullfile())
  stan_vb <- rstan::vb(object = stan_model, data = stan_data, seed = seed,
                       output_samples = n_bootstrap)
  sink()
  #stan_vb_sample <- rstan::extract(stan_vb)$beta
  
  # Thetas<- cbind(params$Theta[,1],params$Theta[,2])
  Thetas<- cbind(extract(stan_vb)$Theta[,1] |> as.vector() ,
                 extract(stan_vb)$Theta[,2]|> as.vector())
  ####Biv beta binom posterior
  densBB<- densBB_functor(x1, x2, n1, n2, u, l, a, a0, a1, a2, i)
  
  #Cálculo de la evidencia ev
  ev10 <- mean(apply(Thetas,1,
                     function(t){I(densBB(t[1],t[2])<densBB_H(supremo_ga_ev))
                     }), na.rm=T)
  
  Ind = I(ev10>k)
  if(is.na(Ind)  ) {
    Ind=FALSE
    return(NA_real_)
  }
  
  
  error2<-(f_ver*Ind)
  return(error2)
}

bet_integ2<-function(j,k){
  
  i<-seq(1,N1)
  # suppressMessages({
  # 
  #   suppressWarnings({
  MC11<-try({mean(sapply(i, integ11,Theta11=Theta11,Theta22=Theta22,N1=N1,k=k, j=j))})
  
  if(!is.numeric(MC11)) {
    return(NA_real_)
  }
  #   })
  # 
  # })
  
  return(MC11)
}


#m_grid_k <- round(seq(0.0000, 0.0909, length.out=10),4)
#m_grid_k <- round(seq(0.1010, 0.1919, length.out=10),4)
#m_grid_k <- round(seq(0.2020, 0.2929, length.out=10),4)
# m_grid_k <- round(seq(0.3030, 0.3939, length.out=10),4)
# m_grid_k <- round(seq(0.4040, 0.4949, length.out=10),4)
#m_grid_k <- round(seq(0.5051, 0.5960, length.out=10),4)
#m_grid_k <- round(seq(0.6061, 0.6970, length.out=10),4)
# m_grid_k <- round(seq(0.7071, 0.7980, length.out=10),4)
# m_grid_k <- round(seq(0.8081, 0.8990, length.out=10),4)
m_grid_k <- round(seq(0.9091, 1.0000, length.out=10),4)
#cores=availableCores()

j<-seq(1,N2)
cores=50
plan(multisession, workers = cores)
beta_result_integral100 <- future_map_dbl(
  m_grid_k,
  ~mean(sapply(j,bet_integ2,k=.x),na.rm=TRUE),
  .options = furrr_options(seed=TRUE)
  
)


#save(beta_result_integral10,file="beta_k_N1000k10.RData")
#save(beta_result_integral20,file="beta_k_N1000k20.RData")
#save(beta_result_integral30,file="beta_k_N1000k30.RData")
#save(beta_result_integral40,file="beta_k_N1000k40.RData")
#save(beta_result_integral50,file="beta_k_N1000k50.RData")
#save(beta_result_integral60,file="beta_k_N1000k60.RData")
#save(beta_result_integral70,file="beta_k_N1000k70.RData")
#save(beta_result_integral80,file="beta_k_N1000k80.RData")
# save(beta_result_integral90,file="betak_N1000k90.RData")
save(beta_result_integral100,file="beta_k_N1000k100.RData")
Sys.time()-mtime
