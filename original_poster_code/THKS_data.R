library(ALA)
library(dplyr)
library(hypergeo)
library(ggplot2)
datos1<-tvsfp
attach(datos1)

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


##Prior no Informativa: Hiperparámetros
alphas_opt=c(0.8373879, 0.8410984, 0.8053298)
a1<-alphas_opt[2]
a2<-alphas_opt[3]
a0<-alphas_opt[1]
a<-a1+a2+a0

####Gráficos a posteriori, media y moda a posteriori
##school and tv
n1=n_yy$n[1]
n2=n_yy$n[2]
x1=X_yy$Bin[1]
x2=X_yy$Bin[2]

##school - no tv
n1=n_yn$n[1]
n2=n_yn$n[2]
x1=X_yn$Bin[1]
x2=X_yn$Bin[2]

##no school -  tv
n1=n_ny$n[1]
n2=n_ny$n[2]
x1=X_ny$Bin[1]
x2=X_ny$Bin[2]

##no school - no tv
n1=n_nn$n[1]
n2=n_nn$n[2]
x1=X_nn$Bin[1]
x2=X_nn$Bin[2]


#####
(n1+a)+(n2+a)>a+(x1+a1)+(x2+a2)
#####
u<-c(a,x1+a1,x2+a2)
l<-c(n1+a,n2+a)


####Biv beta binom posterior
densBB<-function(theta1, theta2) {
  log_f<-lgamma(n1+a)+lgamma(n2+a)-
    lgamma(x1+a1)-lgamma(n1-x1+a-a1)-
    lgamma(x2+a2)-lgamma(n2-x2+a-a2)-
    log(genhypergeo(U=u, L=l, check_mod=T, z=1))+
    (a1+x1-1)*log(theta1)+(a2+a0+(n1-x1)-1)*log(1-theta1)+
    (a2+x2-1)*log(theta2)+(a1+a0+(n2-x2)-1)*log(1-theta2)-
    (a)*log(1-theta1*theta2)
  fun<-exp(log_f)
  return(fun)
}

x     <- seq(0.01, 0.8, 0.01) 
y     <- seq(0.01, 0.8, 0.01)
z     <- outer(x, y, densBB)

persp(x, y, z, theta = -30, phi = 25,
      shade = 0.75, col = "gold", expand = 0.5, r = 2,
      ltheta = 25, ticktype = "detailed", xlab="theta_1",
      ylab="theta_2", zlab="")


###Moda
GA <- ga(type = "real-valued",
         fitness =  function(x) densBB(x[1], x[2]),
         lower = c(0, 0), upper = c(0.99, 0.99),
         popSize = 50, maxiter = 1000, run = 100)
summary(GA)

##Media

esp_th1<-function(theta1, theta2) {theta1*densBB(theta1, theta2)}
esp_th2<-function(theta1, theta2) {theta2*densBB(theta1, theta2)}

library(pracma)

est_mean_th1<-integral2(esp_th1,0,1,0,1)$Q
est_mean_th2<-integral2(esp_th2,0,1,0,1)$Q
est_mean_th1
est_mean_th2

##Cálculo probabilidad a posteriori de H: Probabilidad theta1> theta2

##school and tv
n1=n_yy$n[1]
n2=n_yy$n[2]
x1=X_yy$Bin[1]
x2=X_yy$Bin[2]

##school - no tv
n1=n_yn$n[1]
n2=n_yn$n[2]
x1=X_yn$Bin[1]
x2=X_yn$Bin[2]

##no school - tv
n1=n_ny$n[1]
n2=n_ny$n[2]
x1=X_ny$Bin[1]
x2=X_ny$Bin[2]

##no school - no tv
n1=n_nn$n[1]
n2=n_nn$n[2]
x1=X_nn$Bin[1]
x2=X_nn$Bin[2]

###########
n=as.integer(c(n1,n2))
X=c(x1,x2)
alphas_opt=c(0.8373879, 0.8410984, 0.8053298)
alpha2<-alphas_opt[2]
alpha3<-alphas_opt[3]
alpha1<-alphas_opt[1]
P=length(X)

stan_data<-list(P=P, X=X,n=n,alpha1=alpha1,alpha2=alpha2, alpha3=alpha3)

##compile model
library(rstan)

#model<-stan_model('BBpost3.stan')

##Simulación de theta1 y theta2 de la posteriori
##pass data to stan and run model
#fit<- sampling(model, list(P=P, X=X,n=n,alpha1=alpha1,alpha2=alpha2, alpha3=alpha3), iter=200, chains=4)
fit<-stan(file = 'BBpost3.stan', data = stan_data, iter = 10000, cores = 4)
#diagnose
print(fit)

##Valores para theta1 y theta2 simulados vía MCMC de la posteriori
params<- rstan::extract(fit)
# hist(params$Theta[,1])
# hist(params$Theta[,2])

##Gráficos de las muestras MCMC
library("bayesplot")
library("rstanarm")
library("ggplot2")

posterior <- as.matrix(fit)

plot_title <- ggtitle("Posterior distributions",
                      "with medians  and 80% intervals")
mcmc_areas(posterior,
           pars = c("Theta[1]", "Theta[2]"),
           prob = 0.8) + plot_title

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(posterior,  pars = c("Theta[1]", "Theta[2]"), n_warmup = 300,
                facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 15)

thetas<- data.frame(theta1_est=params$Theta[,1],theta2_est=params$Theta[,2])
##P(theta1>theta2)
probab<-mean(I(thetas$theta1_est>thetas$theta2_est))
probab

##P(theta1<=theta2)
probab<-mean(I(thetas$theta1_est<=thetas$theta2_est))
probab

##########Cálculo de el e-value FBST
####Biv beta binom posterior
#####
(n1+a)+(n2+a)>a+(x1+a1)+(x2+a2)
#####
u<-c(a,x1+a1,x2+a2)
l<-c(n1+a,n2+a)

##Posterior bajo H
densBB_H<-function(theta1) {
  log_f<-lgamma(n1+a)+lgamma(n2+a)-
    lgamma(x1+a1)-lgamma(n1-x1+a-a1)-
    lgamma(x2+a2)-lgamma(n2-x2+a-a2)-
    log(genhypergeo(U=u, L=l, check_mod=T, z=1))+
    (a1+x1-1)*log(theta1)+(a2+a0+(n1-x1)-1)*log(1-theta1)+
    (a2+x2-1)*log(theta1)+(a1+a0+(n2-x2)-1)*log(1-theta1)-
    (a)*log(1-theta1*theta1) 
  fun<-exp(log_f)
  return(fun)
}

library(GA)

###MSUP: Supremo bajo H
GA_ev <- ga(type = "real-valued", 
            fitness =  function(x) densBB_H(x[1]),
            lower = c(0), upper = c(0.99), 
            popSize = 50, maxiter = 1000, run = 100)
#Supremo bajo H
summary(GA_ev)

# lower <- c(0)
# upper <- c(0.8)
# opt3<-DEoptim(densBB_H,lower,upper,control=DEoptim.control(NP = 70,itermax = 400,trace = FALSE))
# out=as.vector(opt3$optim$bestmem)

Thetas<- cbind(params$Theta[,1],params$Theta[,2])
####Biv beta binom posterior
densBB<-function(theta1, theta2) {
  log_f<-lgamma(n1+a)+lgamma(n2+a)-
    lgamma(x1+a1)-lgamma(n1-x1+a-a1)-
    lgamma(x2+a2)-lgamma(n2-x2+a-a2)-
    log(genhypergeo(U=u, L=l, check_mod=T, z=1))+
    (a1+x1-1)*log(theta1)+(a2+a0+(n1-x1)-1)*log(1-theta1)+
    (a2+x2-1)*log(theta2)+(a1+a0+(n2-x2)-1)*log(1-theta2)-
    (a)*log(1-theta1*theta2) 
  fun<-exp(log_f)
  return(fun)
}

#Cálculo de la evidencia ev
ev10 <- mean(apply(Thetas,1,function(t){I(densBB(t[1],t[2])<densBB_H(as.numeric(summary(GA_ev)$solution[1,])))}))
ev10


