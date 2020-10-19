library(MASS)

rm(list = ls()) 


q<-1000;p0<-500

library(parallel)
library(foreach)
library(iterators)
library(doParallel)

######Register a cluster first
cores<-detectCores(logical=F);cores
cl = makeCluster(cores)
registerDoParallel(cl,cores=cores)
chunk.size<-q/cores;chunk.size



####### Parallel computing
res2.p <- foreach(k=1:cores, .combine="cbind") %dopar% {
  library(MASS)
  N<-1000;r=0.3;v=499;C=6;lambda=exp(-C*N^(-2/3+r))

  Betahat<-matrix(NA, p0,chunk.size) 
  alphahat<-matrix(NA, p0,chunk.size)
  BetahatPackage<-matrix(NA, p0,chunk.size)
  alphahatPackage<-matrix(NA, p0,chunk.size) 
 
  
  Betanew<-matrix(NA, p0,chunk.size)   
  alphanew<-matrix(NA, p0,chunk.size) 
    
  Finalj<-matrix(NA, p0,chunk.size) 
  rate<-matrix(NA, p0,chunk.size)
  rate1<-matrix(NA, p0,chunk.size)
  rate2<-matrix(NA, p0,chunk.size)
  criteria.matrix<-matrix(NA, p0,chunk.size)

for (n in 1:chunk.size){

###Generate sequential variable####

x1 = rnorm(N,0);x1[1:10]

Beta<-matrix(NA,1,N);#Beta
alpha<-matrix(NA,1,N);#alpha
x<-matrix(NA,1,N);#x
y<-matrix(NA,1,N);#y
      
###Generate the first one####
 s=1
 x[s]=x1[s];x[s] 
 Beta[s] = 3+1.2*cos(s*2*pi/N-pi/4); Beta[s]
 z<-Beta[s]*x[s]
 pr = 1/(1+exp(-z)) 
 y1 = rbinom(1,1,pr)
 df1 = data.frame(y=y1,x=x[s])
 y[s]<-df1$y;y[s]
      
 ###Generate the rest one####
      for (s in 2:N){
  x[s]=x1[s];x[s] 
  Beta[s] = 3+1.2*cos(s*2*pi/N-pi/4); Beta[s]
  alpha[s] = 2+1.5*cos(s*2.5*pi/N-5/6*pi); alpha[s]
  z<-Beta[s]*x[s]+alpha[s]*y[s-1]
  pr = 1/(1+exp(-z)) 
  y1 = rbinom(1,1,pr)
  df1 = data.frame(y=y1,x=x[s])
  y[s]<-df1$y;y[s]
}

 df0<-cbind.data.frame(y=matrix(y),x=matrix(x))



########Estimation Variable#######################
for (m in 1:p0){
     df<-df0[(m:(m+v)),]
####################My own Estimation############################
########3.Preparation Matrix####
########3-1.Predictor Varibale#(X)###
   X1<-df$x

#######3-2.Response Varibale(Y)####
   Y<-as.numeric(as.matrix(df$y))
   typeY=2;
      

 #####3-4 sequential data analysis####
 ####3-4-1 Sequential Predictor Variable#####
    X1<-as.matrix(X1[-1]);dim(X1)
    Y1<-as.numeric(Y[-(v+1)]);Y1;length(Y1)
    X<-cbind.data.frame(x1=X1,x2=Y1);dim(X);X[1:5,]
   
      
 ####3-4-2 Sequential Response Variable#####
     Y<-Y[-1];length(Y)
     mydata<-cbind.data.frame(x1=X1,x2=Y1,y=Y)
     #res1.p<-cbind.data.frame(x1=c(X$x1,NA),x2=c(X$x2,NA),y=c(Y,NA),xo=df$x,yo=df$y);dim(res1.p)
     #write.csv(res1.p,"testdata1.csv",row.names=FALSE)
     

      X<-as.matrix(X)
#####3-4-3test with package############
      #mylogit <- glm(factor(y)~0+x1+x2, data = mydata, family = "binomial")
      #theta0<-summary(mylogit)
      #BetahatPackage[m,n]<-theta0$coefficients[1,1]
      #alphahatPackage[m,n]<-theta0$coefficients[2,1]


  ################################
  #####4.Softmax Matrix####


#4.Softmax Matrix####
YHatMulti <-
  function(w, X){
    class_inner_prod <-exp(rowSums(X %*%t(w)));#exp(X[1:5,]%*%t(w));
    #class_inner_prod[1:5,]
    denom <- (1 + (class_inner_prod))^-1
    #rowSums(class_inner_prod)[1:5]
    #denom[1:5]
    g<-class_inner_prod * denom
    #g[1:5]
    return(g)
  }

#######5.Initial 	Vector####
 initial_theta<-c(1,1)
 initial_theta<-matrix(initial_theta,nrow =1, ncol = ncol(X),byrow=TRUE)
      
#######6. Edtimation#########################
A<-matrix(0,2,2)
theta<-initial_theta

#############################
tol=0.001;MaxStep=100;
   for(j in 1:MaxStep){
    p<-YHatMulti(theta,X)
    z=length(p);
    A<-matrix(apply(apply(cbind(p,X,((1+m-1):(m+z-1))),1,function(x){-x[1]*(1-x[1])*x[2:3]%*%t(x[2:3])*lambda^(((m+v-1)-x[4])/(N-1)^r)}),1,sum),2,2)
    B<--apply(t(sweep(t(X),2,Y-p,'*'))*(lambda^(((z+m-1)-(1+m-1):(z+m-1))/1000^r)),2,sum)+t(A%*%(as.matrix(as.vector(t(theta)))))
    theta1<-theta###
    theta <- try({ t(ginv(A)%*%t(B))},silent = TRUE)  ###
    criteria=0
    if ('try-error' %in% class(theta)) {criteria=1
    criteria.matrix[m,n]<-1
    break}
    theta<- matrix(data = theta, nrow = 1, ncol = ncol(X), byrow = TRUE)
        if (min(abs(theta1-theta)) < tol) {
          print(paste("iteration fina step:",j))
          print(paste("theta is:",theta))
          break }
        if (j>MaxStep-1){
          print(paste("theta is:",theta))
          print('Too many iterations in method')  
        }
      }
      if ( criteria==1) next
      theta
      
      
      
      Finalj[m,n]<-j;
      
      Betahat[m,n]<-theta[1,1];
      alphahat[m,n]<-theta[1,2];




###################train#####################

x=x1[(m+v+1)]
######################true##################
#Calculation of denominator for probability calculation
Betanew[m,n]= 3+1.2*cos((m+v+1)*2*pi/N-pi/4);  
alphanew[m,n]=2+1.5*cos((m+v+1)*2.5*pi/N-5/6*pi); 
z0 =Betanew[m,n]*x+alphanew[m,n]*df0$y[m+v]
pr0 = 1/(1+exp(-z0)) 
y0 = ifelse(pr0>=0.5,1,0)
dfM0 = data.frame(y=y0,x=x)
ytrue<-df0$y[m+v+1]

#####################estimator#############
z1<-Betahat[m,n]*x+alphahat[m,n]*df0$y[m+v]
pr1 = 1/(1+exp(-z1)) 
y1 = ifelse(pr1>=0.5,1,0)
dfM1 = data.frame(y=y1,x=x)
 
rate[m,n]=sum((ytrue-dfM0$y)!=0)   
rate1[m,n]=sum((ytrue-dfM1$y)!=0) 
rate2[m,n]=sum((dfM0$y-dfM1$y)!=0)   
 }
}


mylist=rbind.data.frame(Betahat=Betahat,alphahat=alphahat,Betanew=Betanew,alphanew=alphanew,criteria=criteria.matrix,Finalj=Finalj,rate=rate,rate1=rate1,rate2=rate2)
#mylist=rbind(Betahat=Betahat,BetahatPackage=BetahatPackage)
  return(mylist)
}



save(res2.p, file = "/home/rnie/simulation/multinormRegressionexp(C6)r3-1000.RData")
      
      

      


