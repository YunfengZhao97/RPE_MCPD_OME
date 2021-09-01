rm(list=ls())

library(BB)
PLRM <- function (x,beta)
{
    # ##############################################
    # INPUTS
    #x: vector of explanatory variables , dim --1xk
    #beta : vector of unknown parameters : dxk
    # OUTPUT
    # pisol: vector of probabilities of all categories , dim --1x(d+1)
    # ##############################################
    k=length (x)
    d=length (beta )/k
    betat= matrix(beta ,k,d)
    pisol=exp(x%*% betat)/(1+sum( exp(x%*% betat)))
    pisol=c(pisol ,1-sum( pisol))
    return( pisol)
}

MLEZZ <- function (X,Y,n)
{
    s1=0
    s2=0
    s3=0
    s4=0
    s5=0
    s6=0
    library(BB)

    fun <- function(x) { 
        for (i in 1:N)
        {
            s1=s1+(Y[i,1]-exp(x[1]*X[i,1]+x[2]*X[i,2])/( 1+exp(x[1]*X[i,1]+x[2]*X[i,2]) +exp(x[3]*X[i,1]+x[4]*X[i,2]) ) )*X[i,1]
            s2=s2+(Y[i,1]-exp(x[1]*X[i,1]+x[2]*X[i,2])/( 1+exp(x[1]*X[i,1]+x[2]*X[i,2]) +exp(x[3]*X[i,1]+x[4]*X[i,2]) ) )*X[i,2]

            s3=s3+(Y[i,2]-exp(x[3]*X[i,1]+x[4]*X[i,2])/( 1+exp(x[1]*X[i,1]+x[2]*X[i,2]) +exp(x[3]*X[i,1]+x[4]*X[i,2]) ) )*X[i,1]
            s4=s4+(Y[i,2]-exp(x[3]*X[i,1]+x[4]*X[i,2])/( 1+exp(x[1]*X[i,1]+x[2]*X[i,2]) +exp(x[3]*X[i,1]+x[4]*X[i,2]) ) )*X[i,2]

       } 
        f <- numeric(length(x))                   
        f[1] <-  s1
        f[2] <-  s2
        f[3] <-  s3
        f[4] <-  s4 
       
        f 
    } 
    startx <- c(0,0,0,0)
    result = dfsane(startx,fun,control=list(maxit=2500,trace = FALSE))
    theta = result$par

    return(theta)
}



DISTANCE4 <- function (W,Y,n)
{
    s1=0
    s2=0
    library(BB)

    fun <- function(x) { 
        for (i in 1:N)
        {W
            bW1=(x[1]*W[i,1]+x[2]*W[i,2])
            bW2=(x[3]*W[i,1]+x[4]*W[i,2])
            b1=t(c(x[1],x[2]))
            b2=t(c(x[3],x[4]))


            s1=s1+(Y[i,1]-1)*exp(bW1/2-b1%*%Suu%*%t(b1)/8)[1]*(W[i,]-Suu%*%t(b1)/2)+Y[i,1]*exp(-bW1/2-b1%*%Suu%*%t(b1)/8)[1]*(W[i,]+Suu%*%t(b1)/2)+Y[i,1]*exp(bW2-bW1/2-(b2-b1/2)%*%Suu%*%t(b2-b1/2)/2)[1]*(W[i,]-Suu%*%t(b2-b1/2))
            s2=s2+(Y[i,2]-1)*exp(bW2/2-b2%*%Suu%*%t(b2)/8)[1]*(W[i,]-Suu%*%t(b2)/2)+Y[i,2]*exp(-bW2/2-b2%*%Suu%*%t(b2)/8)[1]*(W[i,]+Suu%*%t(b2)/2)+Y[i,2]*exp(bW1-bW2/2-(b1-b2/2)%*%Suu%*%t(b1-b2/2)/2)[1]*(W[i,]-Suu%*%t(b1-b2/2))         
         
        } 
        f <- numeric(length(x))                   
        f[1] <-  s1[1]
        f[2] <-  s1[2]
      
        f[3] <-  s2[1]
        f[4] <-  s2[2]
   
        f 
    } 
    startx <- c(0,0,0,0)
    result = dfsane(startx,fun,control=list(maxit=2500,trace = FALSE))
    theta = result$par
    return(theta)
}

Samples =1000 

beta0=c(0.4,-0.3,0.5,-0.1) #true value of the variable vector

N=500 # samples sizes considered

n= vector(,N)+1
PIS= matrix(0,N,3)
PIS1=PIS
PIS2=PIS
Y=matrix (0,N,3)

# definition of parameters which vary for each sample
# generation fo the matrix X

KK=1000
MM1=matrix(0,KK,4)
MM2=matrix(0,KK,4)
library(tcltk)
pb <- tkProgressBar("进度","已完成 %",  0, 100)

for(M in 1:KK)
{
     x0= vector(,N)+1
    x1= (rchisq(N,2))/sqrt(2)
    X= cbind(x0,x1)
    for (i in 1:N)
    {
        PIS[i ,]= PLRM (X[i,], beta0)
        Y[i ,]= rmultinom ( 1,n[i],PIS[i ,])
    }

    w1=x1+rnorm(N,0,1)
    w0=x0
    W= cbind(w0,w1)

    Suu=matrix(c(0,0,0,1),2,2)


    MM1[M,]=DISTANCE4(W,Y,n)
    MM2[M,]=MLEZZ(W,Y,n)

    info <- sprintf("已完成 %d%%", round(M*100/KK))
    setTkProgressBar(pb, M*100/KK, sprintf("进度 (%s)", info), info)

}
close(pb)#关闭进度条

#结果
me1=apply(MM1,2,mean)
se1=apply(MM1,2,sd)    
bias1=me1-beta0    
rmse1=sqrt(apply((t(MM1)-beta0)^2,1,mean))
MM1.all=rbind(beta0,me1,bias1,se1,rmse1)
round(MM1.all,4)

#MLE结果
me2=apply(MM2,2,mean)
se2=apply(MM2,2,sd)    
bias2=me2-beta0    
rmse2=sqrt(apply((t(MM2)-beta0)^2,1,mean))
MM2.all=rbind(beta0,me2,bias2,se2,rmse2)
round(MM2.all,4)





