###### Minyoung Bae #####################
###### equal var & without outlier ######
#rm(list=ls())
m<-5; n<-5; d<-0
mu_x<-0; mu_y<-d; var_x<-1; var_y<-1
#1.1
spl_x<-rnorm(n=m,mu_x, sqrt(var_x))
spl_y<-rnorm(n=n,mu_y, sqrt(var_y))
t.test(spl_x, spl_y, var.equal=TRUE)
#1.2
t.test(spl_x, spl_y, var.equal=FALSE)
#1.3
wilcox.test(spl_x, spl_y)
#1.4
perm_eqvar<-function(x,y,rep){
  cnt = 0 
  t0 = t.test(y,x,var.equal=TRUE)$statistic
  for(i in 1:rep){
    label<-sample(1:(length(x)+length(y)), length(x))
    x_i = c(x,y)[label]
    y_i = c(x,y)[-label]
    t_i = t.test(y_i,x_i,var.equal=TRUE)$statistic
    if(abs(t_i)>abs(t0)) cnt=cnt+1
  }
  return(cnt/rep)
}
perm_eqvar(spl_x,spl_y,1000) #p-value
#1.5
perm_uneqvar<-function(x,y,rep){
  cnt = 0 
  t0 = t.test(y,x,var.equal=FALSE)$statistic
  for(i in 1:rep){
    label<-sample(1:(length(x)+length(y)), length(x))
    x_i = c(x,y)[label]
    y_i = c(x,y)[-label]
    t_i = t.test(y_i,x_i,var.equal=FALSE)$statistic
    if(abs(t_i)>abs(t0)) cnt=cnt+1
  }
  return(cnt/rep)
}
perm_uneqvar(spl_x,spl_y,1000) #p-value

###### 2 ######
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
err1_A_ = err1/100

###### 3 ######
power_A<-function(d){
  m<-5; n<-5
  mu_x<-0; mu_y<-d; var_x<-1; var_y<-1

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
d1<-power_A(1)
###### 4 ######
d2<-power_A(2)
###### 5 ######
d3<-power_A(3)
###### 6 ######
pow<-data.frame(d1=d1,d2=d2,d3=d3)
plot(c(1,2,3),pow[2,],xlab="delta",ylab="power",type="l",col="red"
     ,main="case A:X~N & Y~N",ylim=c(0.2,1))
lines(c(1,2,3),pow[1,],type="l",col="orange",lw=2.5)
lines(c(1,2,3),pow[3,],type="l",col="blue")
lines(c(1,2,3),pow[4,],type="l",col="green")
lines(c(1,2,3),pow[5,],type="l",col="black")
legend("bottomright",c("t_homo","t_hete","wilcox","perm_t_homo","perm_t_hete"),text.col=c("orange","red","blue","green","black"))

###### B.1 ######
m<-5; n<-5; 
alpha<-4; beta<-1/sqrt(alpha); 
var_x=var_y=1; mu_x<-2; mu_y<-alpha*beta

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
err1_B_ = err1_B/100

###### B.2 ######
power_B<-function(d){
  m<-5; n<-5;
  mu_x<-0; mu_y<-d; var_x<-1; var_y<-1
  alpha<-d^2; beta<-1/sqrt(alpha)

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
d1<-power_B(1); d2<-power_B(2); d3<-power_B(3)
pow<-data.frame(d1=d1,d2=d2,d3=d3)
plot(c(1,2,3),pow[1,],xlab="delta",ylab="power",type="l",col="orange",
     main="case B:X~N & Y~Gamma",ylim=c(0,1))
lines(c(1,2,3),pow[2,],type="l",col="red")
lines(c(1,2,3),pow[3,],type="l",col="blue")
lines(c(1,2,3),pow[4,],type="l",col="green")
lines(c(1,2,3),pow[5,],type="l",col="black")
legend("bottomright",c("t_homo","t_hete","wilcox","perm_t_homo","perm_t_hete"),text.col=c("orange","red","blue","green","black"))

###### C.1 ######
m<-5; n<-5; 
var_x=var_y=1; mu_x<-2; mu_y<-2
alpha_x<-4; beta_x<-1/2; alpha_y<-mu_y^2; beta_y<-1/mu_y

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
err1_C_ = err1_C/100

###### C.2 ######
power_C<-function(d){
  m<-5; n<-5;
  var_x=var_y=1; mu_x<-2; mu_y<-2+d
  alpha_x<-4; beta_x<-1/2; alpha_y<-mu_y^2; beta_y<-1/mu_y
  

  err2 = c(0,0,0,0,0) #err2=1-power=p(H0|true par) 
  for(i in 1:100){
    spl_x<-rgamma(m, shape=alpha_x, scale=beta_x)
    spl_y<-rgamma(n, shape=alpha_y, scale=beta_y)
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
d1<-power_C(1); d2<-power_C(2); d3<-power_C(3)
pow<-data.frame(d1=d1,d2=d2,d3=d3)
plot(c(1,2,3),pow[1,],xlab="delta",ylab="power",type="l",col="orange",
     main="case C:X~Gamma & Y~Gamma",ylim=c(0,1))
lines(c(1,2,3),pow[2,],type="l",col="red")
lines(c(1,2,3),pow[3,],type="l",col="blue")
lines(c(1,2,3),pow[4,],type="l",col="green",lw=2.5)
lines(c(1,2,3),pow[5,],type="l",col="black")
legend("bottomright",c("t_homo","t_hete","wilcox","perm_t_homo","perm_t_hete"),text.col=c("orange","red","blue","green","black"))
