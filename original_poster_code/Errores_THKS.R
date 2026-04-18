##no school -  tv

#alpha_500=c(alpha_result_integral10,alpha_result_integral20,alpha_result_integral30,alpha_result_integral40,alpha_result_integral50,alpha_result_integral60,
#          alpha_result_integral70,alpha_result_integral80,alpha_result_integral90,alpha_result_integral100)

k=round(seq(0.0000, 1.0000, length.out=100),4)
f1=spline(k,alpha_500, n = 1000, method = "fmm",
         xmin = min(k), xmax = max(k))
alpha_500_s=smooth.spline(f1$x,f1$y, df=6)$y

#beta_500=c(beta_result_integral10,beta_result_integral20,beta_result_integral30,beta_result_integral40,beta_result_integral50,beta_result_integral60,
#            beta_result_integral70,beta_result_integral80,beta_result_integral90,beta_result_integral100)

k=round(seq(0.0000, 1.0000, length.out=100),4)
f2=spline(k,beta_500, n = 1000, method = "fmm",
          xmin = min(k), xmax = max(k))
beta_500_s=smooth.spline(f2$x,f2$y, df=5)$y

k1=f1$x
sumerrors_s<-alpha_500_s+beta_500_s
kop_s <- min(k1[sumerrors_s==min(sumerrors_s)])
kop_s 

sumerrors<-alpha_500+beta_500
kop <- min(k[sumerrors==min(sumerrors)])
kop 

library(ggplot2)
labels =  c(expression(alpha), expression(beta),expression(alpha+beta))
df1 = data.frame(k1, cat_erros=rep(c("alpha", "betha.mean","somaerros"), each=length(k1)), erros=c(alpha_500_s,beta_500_s,sumerrors_s))
#df1 = data.frame(k, cat_erros=rep(c("alpha", "betha.mean","somaerros"), each=length(k)), erros=c(alpha_500,beta_500,sumerrors))
ggplot(df1,aes(x=k1, y=erros,group=cat_erros)) +
  xlab('k')+ ylab('')+
  theme(legend.title = element_blank())+ 
  geom_line(aes(linetype=cat_erros,color=cat_erros))+
  scale_linetype_manual(values =c("dashed","twodash","solid"),breaks=c("alpha", "betha.mean","somaerros"), labels=labels)+
  scale_colour_manual(values=c("darkcyan", "#E69F00", "gray46"),breaks=c("alpha", "betha.mean","somaerros"), labels=labels)

#######no school - yes tv

##no school -  tv

#alpha_500_ny=c(alpha_result_integral10,alpha_result_integral20,alpha_result_integral30,alpha_result_integral40,alpha_result_integral50,alpha_result_integral60,
#            alpha_result_integral70,alpha_result_integral80,alpha_result_integral90,alpha_result_integral100)

k=round(seq(0.0000, 1.0000, length.out=100),4)
f1=spline(k,alpha_500_ny, n = 1000, method = "fmm",
          xmin = min(k), xmax = max(k))
alpha_500_ny_s=smooth.spline(f1$x,f1$y, df=6)$y

#beta_500_ny=c(beta_result_integral10,beta_result_integral20,beta_result_integral30,beta_result_integral40,beta_result_integral50,beta_result_integral60,
#           beta_result_integral70,beta_result_integral80,beta_result_integral90,beta_result_integral100)

k=round(seq(0.0000, 1.0000, length.out=100),4)
f2=spline(k,beta_500_ny, n = 1000, method = "fmm",
          xmin = min(k), xmax = max(k))
beta_500_ny_s=smooth.spline(f2$x,f2$y, df=5)$y

k1=f1$x
sumerrors_ny_s<-alpha_500_ny_s+beta_500_ny_s
kop_s <- min(k1[sumerrors_ny_s==min(sumerrors_ny_s)])
kop_s 

sumerrors_ny<-alpha_500_ny+beta_500_ny
kop <- min(k[sumerrors_ny==min(sumerrors_ny)])
kop 

library(ggplot2)
labels =  c(expression(alpha), expression(beta),expression(alpha+beta))
df1 = data.frame(k1, cat_erros=rep(c("alpha", "betha.mean","somaerros"), each=length(k1)), erros=c(alpha_500_ny_s,beta_500_ny_s,sumerrors_ny_s))
df1 = data.frame(k, cat_erros=rep(c("alpha", "betha.mean","somaerros"), each=length(k)), erros=c(alpha_500,beta_500,sumerrors))
ggplot(df1,aes(x=k1, y=erros,group=cat_erros)) +
  xlab('k')+ ylab('')+
  theme(legend.title = element_blank())+ 
  geom_line(aes(linetype=cat_erros,color=cat_erros))+
  scale_linetype_manual(values =c("dashed","twodash","solid"),breaks=c("alpha", "betha.mean","somaerros"), labels=labels)+
  scale_colour_manual(values=c("darkcyan", "#E69F00", "gray46"),breaks=c("alpha", "betha.mean","somaerros"), labels=labels)

#######yes school - yes tv

## school -  tv

#alpha_500_yy=c(alpha_result_integral10,alpha_result_integral20,alpha_result_integral30,alpha_result_integral40,alpha_result_integral50,alpha_result_integral60,
#            alpha_result_integral70,alpha_result_integral80,alpha_result_integral90,alpha_result_integral100)

k=round(seq(0.0000, 1.0000, length.out=100),4)
f1=spline(k,alpha_500_yy, n = 1000, method = "fmm",
          xmin = min(k), xmax = max(k))
alpha_500_yy_s=smooth.spline(f1$x,f1$y, df=6)$y

#beta_500_yy=c(beta_result_integral10,beta_result_integral20,beta_result_integral30,beta_result_integral40,beta_result_integral50,beta_result_integral60,
#           beta_result_integral70,beta_result_integral80,beta_result_integral90,beta_result_integral100)

k=round(seq(0.0000, 1.0000, length.out=100),4)
f2=spline(k,beta_500_yy, n = 1000, method = "fmm",
          xmin = min(k), xmax = max(k))
beta_500_yy_s=smooth.spline(f2$x,f2$y, df=5)$y

k1=f1$x
sumerrors_yy_s<-alpha_500_yy_s+beta_500_yy_s
kop_s <- min(k1[sumerrors_yy_s==min(sumerrors_yy_s)])
kop_s 

sumerrors_yy<-alpha_500_yy+beta_500_yy
kop <- min(k[sumerrors_yy==min(sumerrors_yy)])
kop 

library(ggplot2)
labels =  c(expression(alpha), expression(beta),expression(alpha+beta))
df1 = data.frame(k1, cat_erros=rep(c("alpha", "betha.mean","somaerros"), each=length(k1)), erros=c(alpha_500_yy_s,beta_500_yy_s,sumerrors_yy_s))
#df1 = data.frame(k, cat_erros=rep(c("alpha", "betha.mean","somaerros"), each=length(k)), erros=c(alpha_500_yy,beta_500_yy,sumerrors_yy))
ggplot(df1,aes(x=k1, y=erros,group=cat_erros)) +
  xlab('k')+ ylab('')+
  theme(legend.title = element_blank())+ 
  geom_line(aes(linetype=cat_erros,color=cat_erros))+
  scale_linetype_manual(values =c("dashed","twodash","solid"),breaks=c("alpha", "betha.mean","somaerros"), labels=labels)+
  scale_colour_manual(values=c("darkcyan", "#E69F00", "gray46"),breaks=c("alpha", "betha.mean","somaerros"), labels=labels)

#######yes school - no tv

## school -  no tv

alpha_500_yn=c(alpha_result_integral10,alpha_result_integral20,alpha_result_integral30,alpha_result_integral40,alpha_result_integral50,alpha_result_integral60,
            alpha_result_integral70,alpha_result_integral80,alpha_result_integral90,alpha_result_integral100)

k=round(seq(0.0000, 1.0000, length.out=100),4)
f1=spline(k,alpha_500_yn, n = 1000, method = "fmm",
          xmin = min(k), xmax = max(k))
alpha_500_yn_s=smooth.spline(f1$x,f1$y, df=6)$y

beta_500_yn=c(beta_result_integral10,beta_result_integral20,beta_result_integral30,beta_result_integral40,beta_result_integral50,beta_result_integral60,
           beta_result_integral70,beta_result_integral80,beta_result_integral90,beta_result_integral100)

k=round(seq(0.0000, 1.0000, length.out=100),4)
f2=spline(k,beta_500_yn, n = 1000, method = "fmm",
          xmin = min(k), xmax = max(k))
beta_500_yn_s=smooth.spline(f2$x,f2$y, df=5)$y

k1=f1$x
sumerrors_yn_s<-alpha_500_yn_s+beta_500_yn_s
kop_s <- min(k1[sumerrors_yn_s==min(sumerrors_yn_s)])
kop_s 

sumerrors_yn<-alpha_500_yn+beta_500_yn
kop <- min(k[sumerrors_yn==min(sumerrors_yn)])
kop 

library(ggplot2)
labels =  c(expression(alpha), expression(beta),expression(alpha+beta))
df1 = data.frame(k1, cat_erros=rep(c("alpha", "betha.mean","somaerros"), each=length(k1)), erros=c(alpha_500_yn_s,beta_500_yn_s,sumerrors_yn_s))
#df1 = data.frame(k, cat_erros=rep(c("alpha", "betha.mean","somaerros"), each=length(k)), erros=c(alpha_500,beta_500,sumerrors))
ggplot(df1,aes(x=k1, y=erros,group=cat_erros)) +
  xlab('k')+ ylab('')+
  theme(legend.title = element_blank())+ 
  geom_line(aes(linetype=cat_erros,color=cat_erros))+
  scale_linetype_manual(values =c("dashed","twodash","solid"),breaks=c("alpha", "betha.mean","somaerros"), labels=labels)+
  scale_colour_manual(values=c("darkcyan", "#E69F00", "gray46"),breaks=c("alpha", "betha.mean","somaerros"), labels=labels)
