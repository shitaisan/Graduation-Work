library(stats4) #библиотека для поиска ОМП

# считываем данные
load(url("http://statweb.stanford.edu/~ckirby/brad/LSI/datasets-and-programs/data/prostatedata.RData"))

# указываем, какие данные будем обрабатывать
data <- prostatedata

# количество тестирований
N <- nrow(data)

# объемы выборок
n1 <- table(colnames(data))[1]
n2 <- table(colnames(data))[2]

# считает p-значения и статистику по критерию Уэлча
t.welch <- apply(data, 1, function(x) 
  return(t.test(x[(n1+1):(n1+n2)], x[1:n1])$statistic))

pvals.welch <- 2-2*pnorm(abs(t.welch))

# считает статистику и p-значения по критерию Уилкоксона
m <- n1*n2/2
s <- sqrt(n1*n2*(n1+n2+1)/12)
t.wilc <- apply(data, 1, function(x) 
  (wilcox.test(x[(n1+1):(n1+n2)], x[1:n1])$statistic - m)/s)

pvals.wilc <- 2-2*pnorm(abs(t.wilc))

t <- t.welch
pvals <- pvals.welch


# d-апостериорная процедура, основанная на p-value 
#======================================================

# функция плотности модели
f <- function(x, theta) {
  theta[1]+(1-theta[1])*dbeta(x, theta[2], theta[3])}

# логарифмическая функция максимального правдоподобия 
Lf <- function(theta, x) {
  fs <- f(x, theta)
  fs[fs==Inf] <- 1e100
  fs[fs==0] <- 1e-100
  return(-sum(log(fs)))
}

# функция распределения смешанной бета-модели
F <- function(x, theta) {
  theta[1]*x+(1-theta[1])*pbeta(x, theta[2], theta[3])}

# d-риск 1-го рода
R1.pvals <- function(c, theta) {
  c*theta[1]/F(c, theta)}

# d-риск 2-го рода
R0.pvals <- function(c, theta) {
  (1-theta[1])*(1-pbeta(c, theta[2], theta[3]))/(1-F(c, theta))}



# ищем оценки максимльного правдоподобия 
theta.pvals <- mle(function(p0=0.9, a1=0.5, b1=1)
  Lf(c(p0, a1, b1), x=pvals),
  lower = c(0, 0, 0), 
  upper = c(1, 10, 10),
  method = 'L-BFGS-B')@coef

# задаем ограничение на d-риск 1-го рода
alpha <- 0.1
# ищем уровень значимости из ограничения на d-риск 1-го рода
if (theta.pvals[1] >= alpha){
  q <- tryCatch({
    q <- uniroot(function(x) R1.pvals(x, theta.pvals)-alpha, c(0.000001, 1), tol = 1e-10)$root},
    error = function(e) {
      q <- 0
      return (q)})
} else {q <- 1}

# считаем количество отвергнутых гипотез
R.pvals <- sum(pvals<=q)

# процедура Бенджамини-Хохберга на уровне alpha
R.bh <- sum(p.adjust(pvals, 'BH')<=0.1)

# инициализируем значения неизвестных параметров модели оценками МП
pr0 <- theta.pvals[1]
a <- theta.pvals[2]
b <- theta.pvals[3]

# функция для стохастического моделирования
modeling <- function(x, alpha=0.1){
  set.seed(x)
  hypo <- sort(sample(c('=', '!='), N, replace = T,
                      prob = c(pr0, 1-pr0)))
  pvals <- c(rbeta(sum(hypo=='!='), a, b), 
             runif(sum(hypo=='=')))
  
  names(pvals) <- hypo
  
  tryCatch({
    theta <- mle(function(p0=0.5, a=0.5, b=2)
      Lf(c(p0, a, b), x=pvals),
      lower = c(1e5, 1e5, 1e5), 
      upper = c(1, 10, 10),
      method = 'L-BFGS-B')@coef
    
    if (theta[1] >= alpha){
      q <- tryCatch({
        q <- uniroot(function(x) R1.pvals(x, theta)-alpha, c(0.000001, 1), tol = 1e-10)$root},
        error = function(e) {
          q <- 0
          return (q)})
    } else {q <- 1}
    V <<- sum(pvals<=q & names(pvals) %in% c('='))
    R <<- sum(pvals<=q)
    
  }, error = function(e) {
    V <<- R <<- NA
  })
  
  bh <- p.adjust(pvals, 'BH')
  return (data.frame(V = V, R = R,
                     bh.V = sum(bh<=alpha & names(bh) %in% c('=')), bh.R = sum(bh<=alpha)))
}

#======================================================

# d-апостериорная процедура, основанная на статистике 
#======================================================
# фукнция плотности модели
g <- function(x, theta){
  theta[1]*dnorm(x)+(1-theta[1])*dnorm(x, theta[2], theta[3])
}

#Функция распределения модели
G <- function(x, theta){
  theta[1]*pnorm(x)+(1-theta[1])*pnorm(x, theta[2], theta[3])
}

# логарифмическая функция правдоподобия
Lg <- function(theta, x) {
  gs <- g(x, theta)
  gs[gs==Inf] <- 1e100
  gs[gs==0] <- 1e-100
  return(-sum(log(gs)))
}

# d-апостериорный риск 1-го рода
R1.stat <- function(c, theta) {
  theta[1]*2*(1-pnorm(c))/(1-G(c, theta)+G(-c, theta))}

# находим оценки максимального правдоподобия
theta.stat <- mle(function(p0=0.8, mu=0, sd=1)
  Lg(c(p0, mu, sd), x=t),
  lower = c(0, -10, 0), 
  upper = c(1, 10, 10),
  method = 'L-BFGS-B')@coef

# находим критическую константу из ограничения на d-риск 1-го рода
alpha <- 0.1
if (theta.stat[1] >= alpha){
  crit <- tryCatch({
    crit <- uniroot(function(x) R1.stat(x, theta.stat)-alpha, c(-10, 10), tol = 1e-10)$root},
    error = function(e) {
      crit <- -Inf
      return (crit)})
} else {crit <- Inf}

R.stat <- sum(abs(t)>crit)

#======================================================
