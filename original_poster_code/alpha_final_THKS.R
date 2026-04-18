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
#library(MCMCpack)
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

##Error tipo I
# N2=10
# N1=10

# N2=200
# N1=200

N2=500
N1=500

##Generación de variables función auxiliar para aproximar la integral MC
#set.seed(1144)
Theta1<-rbeta(N2,1,1)
#Theta2<-rbeta(N2,1,1)

#f2 <- function(X1,X2,theta1) {dbinom(X1,n1,theta1)*dbinom(X2,n2,theta1)}
f2 <- function(X1,X2,theta1) {exp(dbinom(X1,n1,theta1,log=T)+dbinom(X2,n2,theta1,log=T))}
##Prior no Informativa
alphas_opt=c(0.8373879, 0.8410984, 0.8053298)
a1<-alphas_opt[2]
a2<-alphas_opt[3]
a0<-alphas_opt[1]
a<-a1+a2+a0

alpha2<-alphas_opt[2]
alpha3<-alphas_opt[3]
alpha1<-alphas_opt[1]


stan_model <- rstan::stan_model(file = 'BBpost3.stan', verbose = TRUE)

integ1<-function(i,Theta1,N1,k,j){
  
  # 
  # if( N1 ==  200) {
  #   browser()
  # }
  
  x1<-rbinom(N1,n1,Theta1[j])
  x2<-rbinom(N1,n2,Theta1[j])
  # esto esta raro
  #f <- function(theta1) {dbeta(theta1,1,1)*dbeta(theta1,1,1)}
  # f <- function(theta1,theta2) {dunif(theta1,0,1)*dunif(theta2,0,1)}
  #f_aux<-f(Theta1[j], Theta1[j])
  #f2 <- function(X1,X2,theta1) {dbinom(X1,n1,theta1)*dbinom(X2,n2,theta1)}
  f_ver<-f2(x1[i],x2[i],Theta1[j])
  
  
  
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
  
  # if( N1 ==  200) {
  #   browser()
  # }
  
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
  #           # iter = 10000, 
  #           iter = 500, 
  #           #cores = 20,
  #           cores = 4,
  #           verbose = FALSE,
  #           open_progress = FALSE,
  #           refresh = FALSE)
  #           
  #           
  #diagnose
  #print(fit)
  #
  #
  # browser()
  # 
  # suppressMessages({
  
  # suppressWarnings({
  # fit <- modelo$sample(
  #   data = stan_data,
  #   # seed = 123,
  #   iter_warmup = 500,
  #   iter_sampling = 500,
  #   chains = 4,
  #   parallel_chains = 4,
  #   refresh = 0 # print update every 500 iters
  # )})})
  
  
  # Variational Bayes in RStan
  sink(nullfile())
  stan_vb <- rstan::vb(object = stan_model, data = stan_data, seed = seed,
                       output_samples = n_bootstrap)
  sink()
  #stan_vb_sample <- rstan::extract(stan_vb)$beta
  
  # Thetas<- cbind(params$Theta[,1],params$Theta[,2])
  Thetas<- cbind(extract(stan_vb)$Theta[,1] |> as.vector() ,
                 extract(stan_vb)$Theta[,2]|> as.vector())
  
  # # Thetas<- cbind(params$Theta[,1],params$Theta[,2])
  # Thetas<- cbind(fit$draws("Theta")[,,1] |> as.vector() ,
  #                fit$draws("Theta")[,,2] |> as.vector())
  
  # if(is.na(Thetas) |> any()) {
  #   browser()
  # }
  
  ####Biv beta binom posterior
  densBB<- densBB_functor(x1, x2, n1, n2, u, l, a, a0, a1, a2, i)
  
  
  
  #Cálculo de la evidencia ev
  ev10 <- mean(apply(Thetas,1,
                     function(t){I(densBB(t[1],t[2])<densBB_H(supremo_ga_ev))
                     }), na.rm=T)
  
  
  # if(is.na(ev10) |> any()) {
  #   browser()
  # }
  Ind = I(ev10<=k)
  if(is.na(Ind)  ) {
    Ind=FALSE
    return(NA_real_)
  }
  
  # if(is.na(f_ver)) {
  #   f_ver = NA_real_
  # }
  
  #error1<-(f_ver*I(ev10<=k))/f_aux
  error1<-(f_ver*Ind)
  
  # if( N1 ==  200) {
  #   browser()
  # }
  
  return(error1)
}

alf_integ2<-function(j,k){
  
  # browser()
  
  i<-seq(1,N1)
  # suppressMessages({
  
  # suppressWarnings({
  #MC2<-sapply(i, integ1,Theta1=Theta1,N1=N1,k=k, j=j)
  MC1<-try({mean(sapply(i, integ1,Theta1=Theta1,N1=N1,k=k, j=j),na.rm=TRUE)})
  
  ###Ojo!
  # if(!is.numeric(MC1)) {
  #   return(0)
  # }
  
  # if(is.na(MC1)) {
  #   return(0)
  # }
  ####
  
  if(!is.numeric(MC1)) {
    return(NA_real_)
  }
  
  # browser()
  
  # })
  
  # })
  
  return(MC1)
}

#m_grid_k <- round(seq(0.0000, 0.0909, length.out=10),4)
#m_grid_k <- round(seq(0.1010, 0.1919, length.out=10),4)
#m_grid_k <- round(seq(0.2020, 0.2929, length.out=10),4)
#m_grid_k <- round(seq(0.3030, 0.3939, length.out=10),4)
# m_grid_k <- round(seq(0.4040, 0.4949, length.out=10),4)
#m_grid_k <- round(seq(0.5051, 0.5960, length.out=10),4)
# m_grid_k <- round(seq(0.6061, 0.6970, length.out=10),4)
# m_grid_k <- round(seq(0.7071, 0.7980, length.out=10),4)
# m_grid_k <- round(seq(0.8081, 0.8990, length.out=10),4)
m_grid_k <- round(seq(0.9091, 1.0000, length.out=10),4)


j<-seq(1,N2)
#j<-200
cores=10
plan(multisession, workers = cores)
alpha_result_integral100<- future_map_dbl(
  m_grid_k,
  ~mean(sapply(j,alf_integ2,k=.x),na.rm=TRUE),
  .options = furrr_options(seed=TRUE)
  
)


# result_integral <- purrr::map_dbl(
#   m_grid_k,
#   ~mean(sapply(j,alf_integ2,k=.x))
# 
# )

#print(result_integral)

#save(alpha_result_integral,file="alpha_k_N1000k10.RData")
#save(alpha_result_integral20,file="alpha_k_N1000k20.RData")
#save(alpha_result_integral30,file="alpha_k_N1000k30.RData")
#save(alpha_result_integral40,file="alpha_k_N1000k40.RData")
#save(alpha_result_integral,file="alpha_k_N1000k50.RData")
#save(alpha_result_integral60,file="alpha_k_N1000k60.RData")
#save(alpha_result_integral70#,file="alpha_k_N1000k70.RData")
# save(alpha_result_integral,file="alpha_k_N1000k80.RData")
# save(alpha_result_integral,file="alpha_k_N1000k90.RData")
save(alpha_result_integral100,file="alpha_k_N1000k100.RData")
Sys.time()-mtime