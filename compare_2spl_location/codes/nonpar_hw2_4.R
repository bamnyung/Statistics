###### Minyoung Bae ##
###### diff var ######
d_var<-4
###### A.type1err ######
m<-5; n<-5;
mu_x<-0; mu_y<-0; var_x<-1; var_y<-d_var

err1 = c(0,0,0,0,0)
for(i in 1:100){
  spl_x<-rnorm(n=m,mu_x, sqrt(var_x))
  spl_y<-rnorm(n=n,mu_y, sqrt(var_y))
  if(t.test(spl_x, spl_y, var.equal=TRUE)$p.value < 0.05) err1[1]=err1[1]+1
  if(t.test(spl_x, spl_y, var.equal=FALSE)$p.value < 0.05)err1[2]=err1[2]+1
  if(wilcox.test(spl_x, spl_y)$p.value < 0.05) err1[3]=err1[3]+1
  if(perm_eqvar(spl_x,spl_y,1000) < 0.05) err1[4]=err1[4]+1
  if(perm_uneqvar(spl_x,spl_y,1000) < 0.05) err1[5]=err1[5]+1
}
err1_A_dvar = err1/100

###### A.power ######
power_A_dvar<-function(d,d_var){
  m<-5; n<-5;
  mu_x<-0; mu_y<-d; var_x<-1; var_y<-d_var
  
  err2 = c(0,0,0,0,0) #err2=1-power=p(H0|true par) 
  for(i in 1:100){
    spl_x<-rnorm(m, mu_x, sqrt(var_x))
    spl_y<-rnorm(n, mu_y, sqrt(var_y))
    #t_equal
    if(t.test(spl_y,spl_x,var.equal=TRUE)$p.value>0.05)err2[1]=err2[1]+1
    #t_unequal
    if(t.test(spl_y,spl_x,var.equal=FALSE)$p.value>0.05)err2[2]=err2[2]+1
    #wilcox
    if(wilcox.test(spl_y,spl_x,mu=0)$p.value>0.05) err2[3]=err2[3]+1
    #perm_t_equal
    if(perm_eqvar(spl_y,spl_x,1000)>0.05) err2[4]=err2[4]+1
    #perm_t_unequal
    if(perm_uneqvar(spl_y,spl_x,1000)>0.05) err2[5]=err2[5]+1
  }
  return(1-err2/100)
}
d1<-power_A_dvar(1,d_var); d2<-power_A_dvar(2,d_var); d3<-power_A_dvar(3,d_var)
pow<-data.frame(d1=d1,d2=d2,d3=d3)
plot(c(1,2,3),pow[1,],xlab="delta",ylab="power",type="l",col="orange",
     main="case A\":X~N & Y~N",ylim=c(0,1))
lines(c(1,2,3),pow[2,],type="l",col="red")
lines(c(1,2,3),pow[3,],type="l",col="blue")
lines(c(1,2,3),pow[4,],type="l",col="green",lwd=3)
lines(c(1,2,3),pow[5,],type="l",col="black")
legend("topleft",c("t_homo","t_hete","wilcox","perm_t_homo","perm_t_hete"),text.col=c("orange","red","blue","green","black"))

###### B.type1err ######
m<-5; n<-5
var_x<-1; var_y<-d_var; mu_x<-2; mu_y<-2
alpha<-4/d_var; beta<-d_var/2

err1_B = c(0,0,0,0,0)
for(i in 1:100){
  spl_x<-rnorm(n=m,mu_x, sqrt(var_x))
  spl_y<-rgamma(n,shape=alpha, scale=beta)
  if(t.test(spl_x, spl_y, var.equal=TRUE)$p.value < 0.05) err1_B[1]=err1_B[1]+1
  if(t.test(spl_x, spl_y, var.equal=FALSE)$p.value < 0.05)err1_B[2]=err1_B[2]+1
  if(wilcox.test(spl_x, spl_y)$p.value < 0.05) err1_B[3]=err1_B[3]+1
  if(perm_eqvar(spl_x,spl_y,1000) < 0.05) err1_B[4]=err1_B[4]+1
  if(perm_uneqvar(spl_x,spl_y,1000) < 0.05) err1_B[5]=err1_B[5]+1
}
err1_B_dvar = err1_B/100

###### B.power ######
power_B_dvar<-function(d,d_var){
  m<-5; n<-5;
  mu_x<-0; mu_y<-d; var_x<-1; var_y<-d_var
  alpha<-d^2/d_var; beta<-d_var/d
  
  err2 = c(0,0,0,0,0) #err2=1-power=p(H0|true par) 
  for(i in 1:100){
    spl_x<-rnorm(m, mu_x, sqrt(var_x))
    spl_y<-rgamma(n, shape=alpha, scale=beta)
    #t_equal
    if(t.test(spl_y,spl_x,var.equal=TRUE)$p.value>0.05)err2[1]=err2[1]+1
    #t_unequal
    if(t.test(spl_y,spl_x,var.equal=FALSE)$p.value>0.05)err2[2]=err2[2]+1
    #wilcox
    if(wilcox.test(spl_y,spl_x,mu=0)$p.value>0.05) err2[3]=err2[3]+1
    #perm_t_equal
    if(perm_eqvar(spl_y,spl_x,1000)>0.05) err2[4]=err2[4]+1
    #perm_t_unequal
    if(perm_uneqvar(spl_y,spl_x,1000)>0.05) err2[5]=err2[5]+1
  }
  return(1-err2/100)
}
d1<-power_B_dvar(1,d_var); d2<-power_B_dvar(2,d_var); d3<-power_B_dvar(3,d_var)
pow<-data.frame(d1=d1,d2=d2,d3=d3)
plot(c(1,2,3),pow[1,],xlab="delta",ylab="power",type="l",col="orange",
     main="case B\":X~N & Y~Gamma",ylim=c(0,1))
lines(c(1,2,3),pow[2,],type="l",col="red")
lines(c(1,2,3),pow[3,],type="l",col="blue")
lines(c(1,2,3),pow[4,],type="l",col="green")
lines(c(1,2,3),pow[5,],type="l",col="black")
legend("bottomright",c("t_homo","t_hete","wilcox","perm_t_homo","perm_t_hete"),text.col=c("orange","red","blue","green","black"))

###### C.type1err ######
m<-5; n<-5; 
mu_x<-1; mu_y<-1; var_x<-1; var_y<-d_var
alpha_x<-1; beta_x<-1; alpha_y<-1/d_var; beta_y<-d_var

err1_C = c(0,0,0,0,0)
for(i in 1:100){
  spl_x<-rgamma(n=m,shape=alpha_x, scale=beta_x)
  spl_y<-rgamma(n=n,shape=alpha_y, scale=beta_y)
  if(t.test(spl_x, spl_y, var.equal=TRUE)$p.value < 0.05) err1_C[1]=err1_C[1]+1
  if(t.test(spl_x, spl_y, var.equal=FALSE)$p.value < 0.05)err1_C[2]=err1_C[2]+1
  if(wilcox.test(spl_x, spl_y)$p.value < 0.05) err1_C[3]=err1_C[3]+1
  if(perm_eqvar(spl_x,spl_y,1000) < 0.05) err1_C[4]=err1_C[4]+1
  if(perm_uneqvar(spl_x,spl_y,1000) < 0.05) err1_C[5]=err1_C[5]+1
}
err1_C_dvar = err1_C/100

###### C.power ######
power_C_dvar<-function(d,d_var){
  m<-5; n<-5;
  mu_x<-1; mu_y<-1+d; var_x<-1; var_y<-d_var
  alpha_x<-1; beta_x<-1; alpha_y<-(1+d)^2/d_var; beta_y<-d_var/(1+d)
  
  err2 = c(0,0,0,0,0) #err2=1-power=p(H0|true par) 
  for(i in 1:100){
    spl_x<-rgamma(m, shape=alpha_x, scale=beta_x)
    spl_y<-rgamma(n, shape=alpha_y, scale=beta_y)
    #t_equal
    if(t.test(spl_y,spl_x,var.equal=TRUE)$p.value>0.05)err2[1]=err2[1]+1
    #t_unequal
    if(t.test(spl_y,spl_x,var.equal=FALSE)$p.value>0.05)err2[2]=err2[2]+1
    #wilcox
    if(wilcox.test(spl_y,spl_x,mu=0)$p.value>=0.05) err2[3]=err2[3]+1
    #perm_t_equal
    if(perm_eqvar(spl_y,spl_x,1000)>=0.05) err2[4]=err2[4]+1
    #perm_t_unequal
    if(perm_uneqvar(spl_y,spl_x,1000)>=0.05) err2[5]=err2[5]+1
  }
  return(1-err2/100)
}
d1<-power_C_dvar(1,d_var); d2<-power_C_dvar(2,d_var); d3<-power_C_dvar(3,d_var)
pow<-data.frame(d1=d1,d2=d2,d3=d3)
plot(c(1,2,3),pow[1,],xlab="delta",ylab="power",type="l",col="orange"
     ,main="case C\":X~Gamma & Y~Gamma",ylim=c(0,1))
lines(c(1,2,3),pow[2,],type="l",col="red")
lines(c(1,2,3),pow[3,],type="l",col="blue")
lines(c(1,2,3),pow[4,],type="l",col="green")
lines(c(1,2,3),pow[5,],type="l",col="black")
legend("topleft",c("t_homo","t_hete","wilcox","perm_t_homo","perm_t_hete"),text.col=c("orange","red","blue","green","black"))
