rm(list=ls())

library(mvtnorm) # install.packages("mvtnorm")

income <- log(ts(read.csv("https://research.stlouisfed.org/fred2/data/DPIC96.csv"),start=1947,frequency=4)[,2])
cons <- log(ts(read.csv("https://research.stlouisfed.org/fred2/data/PCECC96.csv"),start=1947,frequency=4)[,2])
data <- na.omit(cbind(income,cons))

n <- ncol(data)
names <- dimnames(data)[[2]]
data <- cbind(data,lag(data,-1))
names <- c(names,paste(names[1:n],"l1",sep="."))
data <- na.omit(diff(data))
data <- cbind(data,1)
names <- c(names,"const")
dimnames(data)[[2]] <- names

y <- t(data[,1:n])
x <- t(data[,-(1:n)])
n.vars <- n*nrow(x)
t <- ncol(y)

Z <- array(NA,dim=c(n,n.vars,t))
for (i in 1:t){
  Z[,,i] <- kronecker(t(x[,i]),diag(n))
}

# Priors
Vprior.i <- diag(0,n.vars)
bprior <- matrix(0,n.vars)
Hprior <- diag(0.0001,n)
Qprior <- diag(0.0001,n.vars)

# Initialise matrices
b.filter <- matrix(0,n.vars,t+1)
b <- matrix(0,n.vars,t+1)
H <- diag(1,n)
H.i <- solve(H)
Q <- diag(1,n.vars)
y0 <- y*0
y.dk <- y*0

# Storage matrices
iter <- 5000 # Total number of iterations
burnin <- 1000 # Burn-in periods
refresh <- seq(0,iter,round(iter/100)) # Refreshment rat of the progress bar
drawsb <- array(NA,dim=c(n.vars,t+1,iter))
drawsOmega <- matrix(NA,n+n.vars,iter)

# Gibbs sampler
# Create progress bar
pb <- txtProgressBar(width = 75,style=3)
for (s in 1:iter){
  
  # Initialise the Kalman filter
  for (i in 1:t){
    y0[,i] <- y[,i]-Z[,,i]%*%b.filter[,i+1]
  }
  ZHZ <- matrix(0,n.vars,n.vars)
  ZHy <- matrix(0,n.vars,1)
  for (i in 1:t){
    ZHZ <- ZHZ + t(Z[,,i])%*%H.i%*%Z[,,i]
    ZHy <- ZHy + t(Z[,,i])%*%H.i%*%matrix(y[,i],n)
  }
  Vpost <- solve(Vprior.i + ZHZ,tol=1e-26)
  bpost <- Vpost%*%(Vprior.i%*%bprior + ZHy)
  b0 <- t(rmvnorm(1,bpost,Vpost))

  for (i in 1:t) {
    y.dk[,i] <- y[,i]-Z[,,i]%*%b0
  }
  # Draw a sample after FFBS
  b.filter <- dk.alg2(k=n,m=n.vars,y=y.dk,Z=Z,H=H,Q=Q,TT=t,a1=0,P1=diag(1,n.vars))
  b <- matrix(b0,n.vars,t+1) + b.filter
  
  # Draw coefficient covariance matrix
  DSum <- matrix(0,n.vars,n.vars)
  for (i in 1:t){
    res <- matrix(b[,i+1] - b[,i],n.vars)
    DSum <- DSum + res%*%t(res)
  }
  Q <- solve(rWishart(1,t+n.vars+2,solve(Qprior + DSum))[,,1])
  
  # Draw error covariance matrix
  DSum <- matrix(0,n,n)
  for (i in 1:t){
    res <- y[,i] - Z[,,i]%*%b[,i+1]
    DSum <- DSum + res%*%t(res)
  }
  H.i <- rWishart(1,t+n+2,solve(Hprior + DSum))[,,1]
  H <- solve(H.i)
  
  # Store draws
  drawsb[,,s] <- b
  drawsOmega[,s] <- sqrt(matrix(c(diag(H),diag(Q)),ncol=1))
  
  # Update progress bar
  if (is.element(s,refresh)) setTxtProgressBar(pb,s/iter)
}

plot.ts(t(drawsOmega))

# OLS
ols <- y%*%t(x)%*%solve(x%*%t(x))
# Transformation to get median and HPD intervals from draws
B <- array(NA,dim=c(t,4,n.vars))
alpha <- .9
for (i in 1:n.vars){
  B[,1,i] <- matrix(ols,n.vars)[i,1]
  for (j in 1:t){
    B[j,2,i] <- median(drawsb[i,j,(burnin+1):iter])
    temp <- drawsb[i,j,(burnin+1):iter][order(drawsb[i,j,(burnin+1):iter])]
    B[j,3,i] <- temp[floor((iter-burnin)*alpha)] # HPDI
    B[j,4,i] <- temp[ceiling((iter-burnin)*(1-alpha))] # HPDI
  }
}

# Plot
par(mfcol=c(2,2))
nam <- matrix(c("Lag Income -> Income","Lag Income -> Consumption",
                "Lag Consumption -> Income","Lag Consumption -> Consumption"))
for (i in 1:4) {series <- ts(B[,,i],start=tsp(data)[1],frequency = tsp(data)[3])
plot.ts(series,plot.type = "single",col=c(3,4,2,2),#ylim=c(-1,1),
        ylab=nam[i,1]);
abline(h=0)};par(mfcol=c(1,1))
