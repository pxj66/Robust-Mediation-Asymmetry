################################################################################
#Optimal tuning constant for the Asymmetric Huber and Tukey Losses
#-------------------------------------------------------------------------------
# Skew-Normal distribution, alpha = 1
#-------------------------------------------------------------------------------
## Huber
library(sn)
alpha = -6
A <- function(k1, k2, alpha){
  integrate(function(x) x^2 * dsn(x, alpha = alpha), -k1, k2)$value + 
    k1^2 * integrate(function(x) dsn(x, alpha = alpha), -20, -k1)$value +
    k2^2 * integrate(function(x) dsn(x, alpha = alpha), k2, 20)$value
}

B <- function(k1, k2, alpha){
  integrate(function(x) dsn(x, alpha = alpha), -k1, k2)$value
}

k1 <- seq(0.1, 3, 0.1)
k2 <- k1
K <- expand.grid(k1 = k1, k2 = k2)
tau.value <- mapply(function(k1, k2) A(k1, k2, alpha) / B(k1, k2, alpha)^2, k1 = K$k1, k2 = K$k2)

Res.asy <- data.frame(k1 = K$k1, k2 = K$k2, tau = tau.value)
Res.sy <- Res.asy[Res.asy$k1==Res.asy$k2, ]
Res.sy[which.min(Res.sy$tau), ]

K[which.min(tau.value), ]
tau.value[which.min(tau.value)]

## Tukey
A <- function(k1, k2, alpha){
  integrate(function(x) (x - 2 * x^3 / k1^2 + x^5 / k1^4)^2 * dsn(x, alpha = alpha), -k1, 0)$value +
    integrate(function(x) (x - 2 * x^3 / k2^2 + x^5 / k2^4)^2 * dsn(x, alpha = alpha), 0, k2)$value
}

B <- function(k1, k2, alpha){
  integrate(function(x) (1 - 6 * x^2 / k1^2 + 5 * x^4 / k1^4) * dsn(x, alpha = alpha), -k1, 0)$value + 
    integrate(function(x) (1 - 6 * x^2 / k2^2 + 5 * x^4 / k2^4) * dsn(x, alpha = alpha), 0, k2)$value
}

k1 <- seq(0.1, 6, 0.1)
k2 <- k1
K <- expand.grid(k1 = k1, k2 = k2)
tau.value <- mapply(function(k1, k2) A(k1, k2, alpha) / B(k1, k2, alpha)^2, k1 = K$k1, k2 = K$k2)

Res.asy <- data.frame(k1 = K$k1, k2 = K$k2, tau = tau.value)
Res.sy <- Res.asy[Res.asy$k1==Res.asy$k2, ]
Res.sy[which.min(Res.sy$tau), ]

K[which.min(tau.value), ]
tau.value[which.min(tau.value)]



#-------------------------------------------------------------------------------
# t distribution with df = 2 
#-------------------------------------------------------------------------------
## Huber
df = 100

A <- function(k1, k2, df){
  integrate(function(x) x^2 * dt(x, df = df), -k1, k2)$value + 
    k1^2 * integrate(function(x) dt(x, df = df), -20, -k1)$value +
    k2^2 * integrate(function(x) dt(x, df = df), k2, 20)$value
}

B <- function(k1, k2, alpha){
  integrate(function(x) dt(x, df = df), -k1, k2)$value
}

k1 <- seq(0.1, 3, 0.1)
k2 <- k1
K <- expand.grid(k1 = k1, k2 = k2)
tau.value <- mapply(function(k1, k2) A(k1, k2, df) / B(k1, k2, df)^2, k1 = K$k1, k2 = K$k2)
K[which.min(tau.value), ]
tau.value[which.min(tau.value)]

Res.asy <- data.frame(k1 = K$k1, k2 = K$k2, tau = tau.value)
Res.sy <- Res.asy[Res.asy$k1==Res.asy$k2, ]
Res.sy[which.min(Res.sy$tau), ]


## Tukey
A <- function(k1, k2, df){
  integrate(function(x) (x - 2 * x^3 / k1^2 + x^5 / k1^4)^2 * dt(x, df = df), -k1, 0)$value +
    integrate(function(x) (x - 2 * x^3 / k2^2 + x^5 / k2^4)^2 * dt(x, df = df), 0, k2)$value
}

B <- function(k1, k2, df){
  integrate(function(x) (1 - 6 * x^2 / k1^2 + 5 * x^4 / k1^4) * dt(x, df = df), -k1, 0)$value + 
    integrate(function(x) (1 - 6 * x^2 / k2^2 + 5 * x^4 / k2^4) * dt(x, df = df), 0, k2)$value
}

k1 <- seq(0.1, 6, 0.1)
k2 <- k1
K <- expand.grid(k1 = k1, k2 = k2)
tau.value <- mapply(function(k1, k2) A(k1, k2, df) / B(k1, k2, df)^2, k1 = K$k1, k2 = K$k2)
K[which.min(tau.value), ]
tau.value[which.min(tau.value)]

Res.asy <- data.frame(k1 = K$k1, k2 = K$k2, tau = tau.value)
Res.sy <- Res.asy[Res.asy$k1==Res.asy$k2, ]
Res.sy[which.min(Res.sy$tau), ]


#-------------------------------------------------------------------------------
# normal
#-------------------------------------------------------------------------------
## Huber
A <- function(k1, k2){
  integrate(function(x) x^2 * dnorm(x), -k1, k2)$value + 
    k1^2 * integrate(function(x) dnorm(x), -20, -k1)$value +
    k2^2 * integrate(function(x) dnorm(x), k2, 20)$value
}

B <- function(k1, k2){
  integrate(function(x) dnorm(x), -k1, k2)$value
}

k1 <- seq(0.1, 3, 0.1)
k2 <- k1
K <- expand.grid(k1 = k1, k2 = k2)
tau.value <- mapply(function(k1, k2) A(k1, k2) / B(k1, k2)^2, k1 = K$k1, k2 = K$k2)
K[which.min(tau.value), ]
tau.value[which.min(tau.value)]

Res.asy <- data.frame(k1 = K$k1, k2 = K$k2, tau = tau.value)
Res.sy <- Res.asy[Res.asy$k1==Res.asy$k2, ]
Res.sy[which.min(Res.sy$tau), ]

## Tukey
A <- function(k1, k2){
    integrate(function(x) (x - 2 * x^3 / k1^2 + x^5 / k1^4)^2 * dnorm(x), -k1, 0)$value +
    integrate(function(x) (x - 2 * x^3 / k2^2 + x^5 / k2^4)^2 * dnorm(x), 0, k2)$value
}

B <- function(k1, k2){
  integrate(function(x) (1 - 6 * x^2 / k1^2 + 5 * x^4 / k1^4) * dnorm(x), -k1, 0)$value + 
    integrate(function(x) (1 - 6 * x^2 / k2^2 + 5 * x^4 / k2^4) * dnorm(x), 0, k2)$value
}

k1 <- seq(0.1, 6, 0.1)
k2 <- k1
K <- expand.grid(k1 = k1, k2 = k2)
tau.value <- mapply(function(k1, k2) A(k1, k2) / B(k1, k2)^2, k1 = K$k1, k2 = K$k2)
K[which.min(tau.value), ]
tau.value[which.min(tau.value)]

Res.asy <- data.frame(k1 = K$k1, k2 = K$k2, tau = tau.value)
Res.sy <- Res.asy[Res.asy$k1==Res.asy$k2, ]
Res.sy[which.min(Res.sy$tau), ]

#-------------------------------------------------------------------------------
# Laplace(0,1)
## Huber
library(VGAM)
A <- function(k1, k2){
  integrate(function(x) x^2 * dlaplace(x), -k1, k2)$value + 
    k1^2 * integrate(function(x) dlaplace(x), -20, -k1)$value +
    k2^2 * integrate(function(x) dlaplace(x), k2, 20)$value
}

B <- function(k1, k2){
  integrate(function(x) dlaplace(x), -k1, k2)$value
}

k1 <- seq(0.1, 3, 0.1)
k2 <- k1
K <- expand.grid(k1 = k1, k2 = k2)
tau.value <- mapply(function(k1, k2) A(k1, k2) / B(k1, k2)^2, k1 = K$k1, k2 = K$k2)
K[which.min(tau.value), ]
tau.value[which.min(tau.value)]

Res.asy <- data.frame(k1 = K$k1, k2 = K$k2, tau = tau.value)
Res.sy <- Res.asy[Res.asy$k1==Res.asy$k2, ]
Res.sy[which.min(Res.sy$tau), ]

## Tukey
A <- function(k1, k2){
  integrate(function(x) (x - 2 * x^3 / k1^2 + x^5 / k1^4)^2 * dlaplace(x), -k1, 0)$value +
    integrate(function(x) (x - 2 * x^3 / k2^2 + x^5 / k2^4)^2 * dlaplace(x), 0, k2)$value
}

B <- function(k1, k2){
  integrate(function(x) (1 - 6 * x^2 / k1^2 + 5 * x^4 / k1^4) * dlaplace(x), -k1, 0)$value + 
    integrate(function(x) (1 - 6 * x^2 / k2^2 + 5 * x^4 / k2^4) * dlaplace(x), 0, k2)$value
}

k1 <- seq(0.1, 6, 0.1)
k2 <- k1
K <- expand.grid(k1 = k1, k2 = k2)
tau.value <- mapply(function(k1, k2) A(k1, k2) / B(k1, k2)^2, k1 = K$k1, k2 = K$k2)
K[which.min(tau.value), ]
tau.value[which.min(tau.value)]

Res.asy <- data.frame(k1 = K$k1, k2 = K$k2, tau = tau.value)
Res.sy <- Res.asy[Res.asy$k1==Res.asy$k2, ]
Res.sy[which.min(Res.sy$tau), ]
