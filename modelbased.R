# Maximum entropy design

library(plgp)
maxent <- function(n, m, theta=0.1, g=0.01, T=100000) 
{  
  if(length(theta) == 1) theta <- rep(theta, m)
  X <- matrix(runif(n*m), ncol=m)
  K <- covar.sep(X, d=theta, g=g) 
  ldetK <- determinant(K, logarithm=TRUE)$modulus
  
  for(t in 1:T) {
    row <- sample(1:n, 1)
    xold <- X[row,] 
    X[row,] <- runif(m)
    Kprime <- covar.sep(X, d=theta, g=g)
    ldetKprime <- determinant(Kprime, logarithm=TRUE)$modulus
    if(ldetKprime > ldetK) { ldetK <- ldetKprime  
    } else { X[row,] <- xold }
  }
  return(X)
}

X <- maxent(25, 2)

plot(X, xlab="x1", ylab="x2")

X <- maxent(25, 2, theta=c(0.1, 0.5))

plot(X, xlab="x1", ylab="x2")

X <- maxent(25, 3)

Is <- as.list(as.data.frame(combn(ncol(X),2)))
par(mfrow=c(1,length(Is)))
for(i in Is) {
  plot(X[,i], xlim=c(0,1), ylim=c(0,1), type="n",
       xlab=paste0("x", i[1]), ylab=paste0("x", i[2]))
  text(X[,i], labels=1:nrow(X)) 
}

# IMSPE designs

library(hetGP)
imspe.criteria <- function(X, theta, g, ...) 
{
  IMSPE(X, theta=theta, Lambda=diag(g, nrow(X)), covtype="Gaussian",
        mult=rep(1, nrow(X)), nu=1) 
}

imspe <- function(n, m, theta=0.1, g=0.01, T=100000, ...) 
{  
  if(length(theta) == 1) theta <- rep(theta, m)
  X <- matrix(runif(n*m), ncol=m)
  I <- imspe.criteria(X, theta, g, ...)
  
  for(t in 1:T) {
    row <- sample(1:n, 1)
    xold <- X[row,] 
    X[row,] <- runif(m)
    Iprime <- imspe.criteria(X, theta, g, ...)
    if(Iprime < I) { I <- Iprime 
    } else { X[row,] <- xold }
  }
  return(X)
}

X <- imspe(25, 2)

plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))

X <- imspe(25, 3)

Is <- as.list(as.data.frame(combn(ncol(X),2)))
par(mfrow=c(1,length(Is)))
for(i in Is) {
  plot(X[,i], xlim=c(0,1), ylim=c(0,1), type="n",
       xlab=paste0("x", i[1]), ylab=paste0("x", i[2]))
  text(X[,i], labels=1:nrow(X)) 
}

imspe.criteria <- function(X, theta, g, Xref)
{
  K <- covar.sep(X, d=theta, g=g)
  Ki <- solve(K)
  KXref <- covar.sep(X, Xref, d=theta, g=0)
  return(mean(1 + g - diag(t(KXref) %*% Ki %*% KXref)))
}

g <- expand.grid(seq(0,1,length=10), seq(0,1, length=10))
X <- imspe(25, 2, Xref=g)

plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))
points(g, pch=20, cex=0.25, col="gray")

Xref <- rmvnorm(100, mean=c(0.25, 0.25), 
                sigma=0.005*rbind(c(2, 0.35), c(0.35, 0.1)))
Xref <- rbind(Xref, 
              rmvnorm(100, mean=c(0.25, 0.25), 
                      sigma=0.005*rbind(c(0.1, -0.35), c(-0.35, 2))))
X <- imspe(25, 2, Xref=Xref)

plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))
points(Xref, pch=20, cex=0.25, col="gray")

# Sequential design

# Active learning Mackey
set.seed(1)
library(lhs)
ninit <- 12 # initial sample size
X <- randomLHS(ninit, 2) #initial LHD
f <- function(X, sd=0.01) #true function
{
  X[,1] <- (X[,1] - 0.5)*6 + 1
  X[,2] <- (X[,2] - 0.5)*6 + 1
  y <- X[,1] * exp(-X[,1]^2 - X[,2]^2) + rnorm(nrow(X), sd=sd)
}
y <- f(X)

plot(X, xlim=c(0,1), ylim=c(0,1), xlab="x1", ylab="x2")
points(X, pch=20)

library(laGP)
g <- garg(list(mle=TRUE, max=1), y)
d <- darg(list(mle=TRUE, max=0.25), X)
gpi <- newGP(X, y, d=d$start, g=g$start, dK=TRUE)
mle <- jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), d$ab, g$ab) #Fitting an isotropic GP using laGP

x1 <- x2 <- seq(0, 1, length=100)
XX <- expand.grid(x1, x2)
yytrue <- f(XX, sd=0)

rmse <- sqrt(mean((yytrue - predGP(gpi, XX, lite=TRUE)$mean)^2)) #Estimate of RMSE

obj.alm <- function(x, gpi) 
  - sqrt(predGP(gpi, matrix(x, nrow=1), lite=TRUE)$s2) #negative since optim function minimizes in R

mymaximin <- function(n, m, T=100000, Xorig=NULL) # maximin design function from spacefilling.R
{
  X <- matrix(runif(n*m), ncol=m)     ## initial design
  d <- distance(X)
  d <- d[upper.tri(d)]
  md <- min(d)
  if(!is.null(Xorig)) {               ## new code
    md2 <- min(distance(X, Xorig))
    if(md2 < md) md <- md2
  }
  
  for(t in 1:T) {
    row <- sample(1:n, 1)
    xold <- X[row,]                   ## random row selection
    X[row,] <- runif(m)               ## random new row
    d <- distance(X)
    d <- d[upper.tri(d)]
    mdprime <- min(d)
    if(!is.null(Xorig)) {             ## new code
      mdprime2 <- min(distance(X, Xorig))
      if(mdprime2 < mdprime) mdprime <- mdprime2
    }
    if(mdprime > md) { md <- mdprime  ## accept
    } else { X[row,] <- xold }        ## reject
  }
  
  return(X)
}

xnp1.search <- function(X, gpi, obj=obj.alm, ...)
{
  start <- mymaximin(nrow(X), 2, T=100*nrow(X), Xorig=X) #start with initial point from maximin design
  xnew <- matrix(NA, nrow=nrow(start), ncol=ncol(X) + 1)
  for(i in 1:nrow(start)) { #optimize
    out <- optim(start[i,], obj, method="L-BFGS-B", lower=0, 
                 upper=1, gpi=gpi, ...)
    xnew[i,] <- c(out$par, -out$value)
  }
  solns <- data.frame(cbind(start, xnew))
  names(solns) <- c("s1", "s2", "x1", "x2", "val")
  return(solns)
}

# Visualize optimized points from different initializations
solns <- xnp1.search(X, gpi)
plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))
arrows(solns$s1, solns$s2, solns$x1, solns$x2, length=0.1)
m <- which.max(solns$val)
prog <- solns$val[m]
points(solns$x1[m], solns$x2[m], col=2, pch=20, cex=3)

xnew <- as.matrix(solns[m, 3:4]) #add the maximum variance point in (red dot)
X <- rbind(X, xnew)
y <- c(y, f(xnew))

updateGP(gpi, xnew, y[length(y)]) #update GP fit
mle <- rbind(mle, jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), 
                         d$ab, g$ab))
rmse <- c(rmse, sqrt(mean((yytrue - predGP(gpi, XX, lite=TRUE)$mean)^2))) #updating RMSE

solns <- xnp1.search(X, gpi) #getting another point and adding this in
m <- which.max(solns$val)
prog <- c(prog, solns$val[m])
xnew <- as.matrix(solns[m, 3:4])
X <- rbind(X, xnew)
y <- c(y, f(xnew))
updateGP(gpi, xnew, y[length(y)])
mle <- rbind(mle, jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), 
                         d$ab, g$ab))
p <- predGP(gpi, XX, lite=TRUE)
rmse <- c(rmse, sqrt(mean((yytrue - p$mean)^2))) #updating RMSE

plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))
arrows(solns$s1, solns$s2, solns$x1, solns$x2, length=0.1)
m <- which.max(solns$val)
points(solns$x1[m], solns$x2[m], col=2, pch=20, cex = 3)

for(i in nrow(X):24) { # run this sequentially to a total of 25 points
  solns <- xnp1.search(X, gpi)
  m <- which.max(solns$val)
  prog <- c(prog, solns$val[m])
  xnew <- as.matrix(solns[m, 3:4])
  X <- rbind(X, xnew)
  y <- c(y, f(xnew))
  updateGP(gpi, xnew, y[length(y)])
  mle <- rbind(mle, jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), 
                           d$ab, g$ab))
  p <- predGP(gpi, XX, lite=TRUE)
  rmse <- c(rmse, sqrt(mean((yytrue - p$mean)^2)))
}

mle[seq(1,nrow(mle),by=2),] #summarizes hyperparameter updates in the sequential procedure


par(mfrow=c(1,2)) #visualizing sequential points
cols <- heat.colors(128)
image(x1, x2, matrix(p$mean, ncol=length(x1)), col=cols, main="mean")
text(X, labels=1:nrow(X), cex=0.75)
image(x1, x2, matrix(sqrt(p$s2), ncol=length(x1)), col=cols, main="sd")
text(X, labels=1:nrow(X), cex=0.75)

plot((ninit+1):nrow(X), prog, xlab="n") #and its ALM criterion progression

d <- darg(list(mle=TRUE), X) #let's add 100 more points
for(i in nrow(X):99) {
  solns <- xnp1.search(X, gpi)
  m <- which.max(solns$val)
  prog <- c(prog, solns$val[m])
  xnew <- as.matrix(solns[m, 3:4])
  X <- rbind(X, xnew)
  y <- c(y, f(xnew))
  updateGP(gpi, xnew, y[length(y)])
  mle <- rbind(mle, jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), 
                           d$ab, g$ab))
  p <- predGP(gpi, XX, lite=TRUE)
  rmse <- c(rmse, sqrt(mean((yytrue - p$mean)^2)))
}

mle[seq(14,nrow(mle), by=10),]


par(mfrow=c(1,2))
plot((ninit+1):nrow(X), prog, xlab="n: design size", ylab="ALM progress")
plot(ninit:nrow(X), rmse, xlab="n: design size", ylab="OOS RMSE")

par(mfrow=c(1,2))
image(x1, x2, matrix(p$mean, ncol=length(x1)), col=cols, main="mean")
text(X, labels=1:nrow(X), cex=0.75)
image(x1, x2, matrix(sqrt(p$s2), ncol=length(x1)), col=cols, main="sd")
text(X, labels=1:nrow(X), cex=0.75)

ninit <- 100
X2 <- randomLHS(ninit, 2)
y2 <- f(X2)
gpi <- newGP(X2, y2, d=d$start, g=g$start, dK=TRUE)
mle <- jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), d$ab, g$ab) #Fitting an isotropic GP using laGP
yytrue <- f(XX, sd=0)
rmse2 <- sqrt(mean((yytrue - predGP(gpi, XX, lite=TRUE)$mean)^2))
rmse2


# Active learning Cohn
set.seed(1)
library(lhs)
ninit <- 12 # initial sample size

obj.alc <- function(x, gpi, Xref) 
  - sqrt(alcGP(gpi, matrix(x, nrow=1), Xref))
X <- randomLHS(ninit, 2) #initial LHD
f <- function(X, sd=0.01) #true function
{
  X[,1] <- (X[,1] - 0.5)*6 + 1
  X[,2] <- (X[,2] - 0.5)*6 + 1
  y <- X[,1] * exp(-X[,1]^2 - X[,2]^2) + rnorm(nrow(X), sd=sd)
}
y <- f(X)

library(laGP)
g <- garg(list(mle=TRUE, max=1), y)
d <- darg(list(mle=TRUE, max=0.25), X)
gpi <- newGP(X, y, d=d$start, g=g$start, dK=TRUE)
mle <- jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), d$ab, g$ab)
x1 <- x2 <- seq(0, 1, length=100)
XX <- expand.grid(x1, x2)
yytrue <- f(XX, sd=0)
p <- predGP(gpi, XX, lite=TRUE)
rmse.alc <- sqrt(mean((yytrue - p$mean)^2))

Xref <- randomLHS(100, 2)
solns <- xnp1.search(X, gpi, obj=obj.alc, Xref=Xref)
m <- which.max(solns$val)
xnew <- as.matrix(solns[m, 3:4])
prog.alc <- solns$val[m]

plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))
arrows(solns$s1, solns$s2, solns$x1, solns$x2, length=0.1)
points(solns$x1[m], solns$x2[m], col=2, pch=20)
points(Xref, cex=0.25, pch=20, col="gray")

X <- rbind(X, xnew)
y <- c(y, f(xnew))
updateGP(gpi, xnew, y[length(y)])
mle <- rbind(mle, jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), 
                         d$ab, g$ab))
p <- predGP(gpi, XX, lite=TRUE)
rmse.alc <- c(rmse.alc, sqrt(mean((yytrue - p$mean)^2)))

d <- darg(list(mle=TRUE), X)
for(i in nrow(X):99) {
  Xref <- randomLHS(100, 2)
  solns <- xnp1.search(X, gpi, obj=obj.alc, Xref=Xref)
  m <- which.max(solns$val)
  prog.alc <- c(prog.alc, solns$val[m])
  xnew <- as.matrix(solns[m, 3:4])
  X <- rbind(X, xnew)
  y <- c(y, f(xnew))
  updateGP(gpi, xnew, y[length(y)])
  mle <- rbind(mle, jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), 
                           d$ab, g$ab))
  p <- predGP(gpi, XX, lite=TRUE)
  rmse.alc <- c(rmse.alc, sqrt(mean((yytrue - p$mean)^2)))
}

par(mfrow=c(1,2))
plot((ninit+1):nrow(X), prog.alc, xlab="n: design size", 
     ylab="ALC progress")
plot(ninit:nrow(X), rmse, xlab="n: design size", ylab="OOS RMSE")
points(ninit:nrow(X), rmse.alc, col=2, pch=20)
legend("topright", c("alm", "alc"), pch=c(21,20), col=1:2)

par(mfrow=c(1,2))
image(x1, x2, matrix(p$mean, ncol=length(x1)), col=cols, main="mean")
text(X, labels=1:nrow(X), cex=0.75)
image(x1, x2, matrix(sqrt(p$s2), ncol=length(x1)), col=cols, main="sd")
text(X, labels=1:nrow(X), cex=0.75)
