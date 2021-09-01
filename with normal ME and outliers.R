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
        for (i in 1:N[l])
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





OutlyLDISTANCE <- function (W,Y,lambda,n)
{
    s1=0
    s2=0
    library(BB)

    fun <- function(x) { 
        for (i in 1:N)
        {
            bW1=(x[1]*W[i,1]+x[2]*W[i,2])
            bW2=(x[3]*W[i,1]+x[4]*W[i,2])

            b1=t(c(x[1],x[2]))
            b2=t(c(x[3],x[4]))      

            s1=s1+( 
                 (1-Y[i,1])*( exp(bW2*3/2+bW1/2-(3*b2+b1)%*%Suu%*%t(3*b2+b1)/8) )[1]*(W[i,]-Suu%*%t(3*b2+b1)/2)
               + (1-Y[i,1])*( exp(bW2/2+bW1/2-(b2+b1)%*%Suu%*%t(b2+b1)/8) )[1]*(W[i,]-Suu%*%t(b2+b1)/2)
               - (1-Y[i,2])*( exp(bW2*(3/2+lambda)+(1/2-lambda)*bW1-((3+2*lambda)*b2+(1-2*lambda)*b1)%*%Suu%*%t((3+2*lambda)*b2+(1- 2*lambda)*b1)/8) )[1]*(W[i,]-Suu%*%t((3+2*lambda)*b2+(1-2*lambda)*b1)/2)
               - (1-Y[i,3])*( exp(bW2/2+(1/2-lambda)*bW1-(b2+(1-2*lambda)*b1)%*%Suu%*%t(b2+(1-2*lambda)*b1)/8) )[1]*(W[i,]-Suu%*%t (b2+(1-2*lambda)*b1)/2)
               + Y[i,2]*( 
                            exp((1/2+lambda)*bW2+(3/2-lambda)*bW1-((2*lambda+1)*b2+(3-2*lambda)*b1)%*%Suu%*%t((2*lambda+1)*b2+(3- 2*lambda)*b1)/8)[1]*(W[i,]-Suu%*%t((2*lambda+1)*b2+(3-2*lambda)*b1)/2)
                           +exp((1/2+lambda)*bW2+(1/2-lambda)*bW1-((2*lambda+1)*b2+(1-2*lambda)*b1)%*%Suu%*%t((2*lambda+1)*b2+(1- 2*lambda)*b1)/8)[1]*(W[i,]-Suu%*%t((2*lambda+1)*b2+(1-2*lambda)*b1)/2)
                        )
               + Y[i,3]*( 
                            exp(bW2/2+(3/2-lambda)*bW1-(b2+(3-2*lambda)*b1)%*%Suu%*%t(b2+(3-2*lambda)*b1)/8)[1]*(W[i,]-Suu%*%t(b2+(3 -2*lambda)*b1)/2)
                           +exp(bW2*3/2+(1/2-lambda)*bW1-(3*b2+(1-2*lambda)*b1)%*%Suu%*%t(3*b2+(1-2*lambda)*b1)/8)[1]*(W[i,]-Suu%*%t (3*b2+(1-2*lambda)*b1)/2) 
                        )
               - Y[i,1]*( 
                            exp(bW2*3/2-bW1/2-(3*b2-b1)%*%Suu%*%t(3*b2-b1)/8)[1]*(W[i,]-Suu%*%t(3*b2-b1)/2)
                           +exp(bW2*5/2-bW1/2-(5*b2-b1)%*%Suu%*%t(5*b2-b1)/8)[1]*(W[i,]-Suu%*%t(5*b2-b1)/2)
                           +exp(bW2/2-bW1/2-(b2-b1)%*%Suu%*%t(b2-b1)/8)[1]*(W[i,]-Suu%*%t(b2-b1)/2)
                           +exp(bW2*3/2-bW1/2-(3*b2-b1)%*%Suu%*%t(3*b2-b1)/8)[1]*(W[i,]-Suu%*%t(3*b2-b1)/2) 
                        )
                )
           
            s2=s2+( 
                 (1-Y[i,2])*( exp(bW1*3/2+bW2/2-(3*b1+b2)%*%Suu%*%t(3*b1+b2)/8) )[1]*(W[i,]-Suu%*%t(3*b1+b2)/2)
               + (1-Y[i,2])*( exp(bW1/2+bW2/2-(b1+b2)%*%Suu%*%t(b1+b2)/8) )[1]*(W[i,]-Suu%*%t(b1+b2)/2)
               - (1-Y[i,1])*( exp(bW1*(3/2+lambda)+(1/2-lambda)*bW2-((3+2*lambda)*b1+(1-2*lambda)*b2)%*%Suu%*%t((3+2*lambda)*b1+(1- 2*lambda)*b2)/8) )[1]*(W[i,]-Suu%*%t((3+2*lambda)*b1+(1-2*lambda)*b2)/2)
               - (1-Y[i,3])*( exp(bW1/2+(1/2-lambda)*bW2-(b1+(1-2*lambda)*b2)%*%Suu%*%t(b1+(1-2*lambda)*b2)/8) )[1]*(W[i,]-Suu%*%t (b1+(1-2*lambda)*b2)/2)
               + Y[i,1]*( 
                            exp((1/2+lambda)*bW1+(3/2-lambda)*bW2-((2*lambda+1)*b1+(3-2*lambda)*b2)%*%Suu%*%t((2*lambda+1)*b1+(3- 2*lambda)*b2)/8)[1]*(W[i,]-Suu%*%t((2*lambda+1)*b1+(3-2*lambda)*b2)/2)
                           +exp((1/2+lambda)*bW1+(1/2-lambda)*bW2-((2*lambda+1)*b1+(1-2*lambda)*b2)%*%Suu%*%t((2*lambda+1)*b1+(1- 2*lambda)*b2)/8)[1]*(W[i,]-Suu%*%t((2*lambda+1)*b1+(1-2*lambda)*b2)/2)
                        )
               + Y[i,3]*( 
                            exp(bW1/2+(3/2-lambda)*bW2-(b1+(3-2*lambda)*b2)%*%Suu%*%t(b1+(3-2*lambda)*b2)/8)[1]*(W[i,]-Suu%*%t(b1+(3 -2*lambda)*b2)/2)
                           +exp(bW1*3/2+(1/2-lambda)*bW2-(3*b1+(1-2*lambda)*b2)%*%Suu%*%t(3*b1+(1-2*lambda)*b2)/8)[1]*(W[i,]-Suu%*%t (3*b1+(1-2*lambda)*b2)/2) 
                        )
               - Y[i,2]*( 
                            exp(bW1*3/2-bW2/2-(3*b1-b2)%*%Suu%*%t(3*b1-b2)/8)[1]*(W[i,]-Suu%*%t(3*b1-b2)/2)
                           +exp(bW1*5/2-bW2/2-(5*b1-b2)%*%Suu%*%t(5*b1-b2)/8)[1]*(W[i,]-Suu%*%t(5*b1-b2)/2)
                           +exp(bW1/2-bW2/2-(b1-b2)%*%Suu%*%t(b1-b2)/8)[1]*(W[i,]-Suu%*%t(b1-b2)/2)
                           +exp(bW1*3/2-bW2/2-(3*b1-b2)%*%Suu%*%t(3*b1-b2)/8)[1]*(W[i,]-Suu%*%t(3*b1-b2)/2) 
                        )
                )


        } 
        f <- numeric(length(x))                   
        f[1] <-  (lambda+1)/N[l]^(lambda+1)*s1[1]
        f[2] <-  (lambda+1)/N[l]^(lambda+1)*s1[2]

        f[3] <-  (lambda+1)/N[l]^(lambda+1)*s2[1]
        f[4] <-  (lambda+1)/N[l]^(lambda+1)*s2[2] 
        
        f 
    } 
    startx <- c(0,0,0,0)
    result = dfsane(startx,fun,control=list(maxit=2500,trace = FALSE))
    theta = result$par

    return(theta)
}  

Samples =1000 # number of Samples in the simulation

beta0=c(-0.4,-0.1,-0.5,0.2) #true value of the variable vector

N=c(100,200,500,1000) # Samples sizes considered
lambda=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9) # lambdas considered for the MDPDEs
out=0.7 # percentage of pure data

l=4 #第一个数据量，100个数据

n= vector(,N[l])+1  #全为1的向量
PIS= matrix(0,N[l],3) 
PIS1=PIS 
PIS2=PIS 
Y=matrix (0,N[l],3) 

MM00=matrix(0,Samples,4) #MLE,不处理测量误差

MM5=matrix(0,Samples,4) #lambda[5]


# definition of parameters which vary for each sample
# generation fo the matrix X

library(tcltk)
pb <- tkProgressBar("进度","已完成 %",  0, 100)

for (M in 1: Samples ) 
{
    x0= vector(,N[l])+1
    x1= rnorm(N[l],-1.5,2)
    X= cbind(x0,x1)
    for (i in 1:N[l])
    {
        PIS[i ,]= PLRM (X[i,], beta0)
        Y[i ,]= rmultinom ( 1,n[i],PIS[i ,])
    }

    #Y2是生成有异常值的，异常比例为out%,在这里为95%
    Y2=Y
    for (i in floor(out*N[l]):floor(0.85*N[l]))
    {
        aux=c(Y[i,3],Y[i,1],Y[i,2])
        Y2[i ,]= aux
    }

    w1=x1+rnorm(N[l],0,1)
    w0=x0
    W= cbind(w0,w1)

    Suu=matrix(c(0,0,0,1),2,2)


    MM00[M,]=MLEZZ(W,Y2,n)

    MM5[M,]=OutlyLDISTANCE(W,Y2,lambda[7],n)




    info <- sprintf("已完成 %d%%", round(M*100/Samples))
    setTkProgressBar(pb, M*100/Samples, sprintf("进度 (%s)", info), info)

}
close(pb)#关闭进度条

#结果MLE,不处理测量误差
me00=apply(MM00,2,mean)
se00=apply(MM00,2,sd)    
bias00=me00-beta0    
rmse00=sqrt(apply((t(MM00)-beta0)^2,1,mean))
MM00.all=rbind(beta0,me00,bias00,se00,rmse00)
round(MM00.all,4)



#结果lambda[5]
me5=apply(MM5,2,mean)
se5=apply(MM5,2,sd)    
bias5=me5-beta0    
rmse5=sqrt(apply((t(MM5)-beta0)^2,1,mean))
MM5.all=rbind(beta0,me5,bias5,se5,rmse5)
round(MM5.all,4)





 