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

# считает статистику и p-значения по критерию Уэлча 
t.welch <- apply(data, 1, function(x) 
  return(t.test(x[(n1+1):(n1+n2)], x[1:n1])$statistic))

pvals.welch <- 1-pnorm(t.welch)

# считает статистику и p-значения по критерию Уилкоксона
m <- n1*n2/2
s <- sqrt(n1*n2*(n1+n2+1)/12)
t.wilc <- apply(data, 1, function(x) 
  (wilcox.test(x[(n1+1):(n1+n2)], x[1:n1])$statistic - m)/s)

pvals.wilc <- 1-pnorm(t.wilc)


# выбираем критерий, для которого нужно провести процедуру
t <- t.welch
pvals <- pvals.welch


# d-апостериорная процедура, основанная на p-value 
#======================================================

# функция плотности модели
f <- function(x, theta) {
  theta[1]+theta[2]*dbeta(x, theta[3], theta[4])+(1-theta[1]-theta[2])*dbeta(x, theta[5], theta[6])}

# логарифмическая функция максимального правдоподобия 
Lf <- function(theta, x) {
  if (theta[1]+theta[2]>1) return (1e10)
  fs <- f(x, theta)
  fs[fs==Inf] <- 1e10
  return(-sum(log(fs)))
}

# функция распределения смешанной бета-модели
F <- function(x, theta) {
  theta[1]*x+theta[2]*pbeta(x, theta[3], theta[4])+(1-theta[1]-theta[2])*pbeta(x, theta[5], theta[6])}

# d-риск 1-го рода
R1.pvals <- function(c, theta) {
  (c*theta[1]+(1-theta[1]-theta[2])*pbeta(c, theta[5], theta[6]))/F(c, theta)}

# d-риск 2-го рода
R0.pvals <- function(c, theta) {
  theta[2]*(1-pbeta(c, theta[3], theta[4]))/(1-F(c, theta))}


# ищем оценки максимльного правдоподобия 
theta.pvals <- mle(function(p0=0.9, p1=0.05, a1=0.5, b1=5, a2=5, b2=0.5)
  Lf(c(p0, p1, a1, b1, a2, b2), x=pvals),
  lower = c(1e-10, 1e-10, 1e-10, 1, 1, 1e-10), 
  upper = c(1, 1, 1, 10, 10, 1),
  method = 'L-BFGS-B')@coef

# задаем ограничение на d-риск 1-го рода
alpha <- 0.1
# ищем уровень значимости из ограничения на d-риск 1-го рода
if (sum(theta.pvals[1:2]) >= alpha){
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
p0 <- theta[1]
p1 <- theta[2]
a1 <- theta[3]
b1 <- theta[4]
a2 <- theta[5]
b2 <- theta[6]

# функция для стохастического моделирования
modeling <- function(x, alpha=0.1){
  set.seed(x)
  hypo <- sort(sample(c('=', '>', '<'), N, replace = T,
                      prob = c(p0, p1, 1-p0-p1)))
  pvals <- c( rbeta(sum(hypo=='<'), a1, b1), 
              runif(sum(hypo=='=')),
              rbeta(sum(hypo=='>'), a2, b2))
  
  names(pvals) <- hypo
  
  tryCatch({
    theta.pvals <- mle(function(p0=0.5, p1=0.2, a1=0.5, b1=5, a2=5, b2=0.5)
      Lf(c(p0, p1, a1, b1, a2, b2), x=pvals),
      lower = c(0, 0, 0, 1, 1, 0), 
      upper = c(1, 1, 1, 10, 10, 1),
      method = 'L-BFGS-B')@coef
    
    if (sum(theta.pvals[1:2]) >= alpha){
      q <- tryCatch({
        q <- uniroot(function(x) R1.pvals(x, theta.pvals)-alpha, c(0.000001, 1), tol = 1e-10)$root},
        error = function(e) {
          q <- 0
          return (q)})
    } else {q <- 1}
    V <<- sum(pvals<=q & names(pvals) %in% c('=', '<'))
    R <<- sum(pvals<=q)
    
  }, error = function(e) {
    V <<- R <<- NA
  })
  
  bh <- p.adjust(pvals, 'BH')
  return (data.frame(V = V, R = R,
                     bh.V = sum(bh<=alpha & names(bh) %in% c('=', '<')), bh.R = sum(bh<=alpha)))
}

#======================================================

# d-апостериорная процедура, основанная на статистике 
#======================================================
g <- function(x, theta){
  theta[1]*dnorm(x)+theta[2]*dnorm(x, theta[3], theta[4])+(1-theta[1]-theta[2])*dnorm(x, theta[5], theta[6])
}

G <- function(x, theta){
  theta[1]*pnorm(x)+theta[2]*pnorm(x, theta[3], theta[4])+(1-theta[1]-theta[2])*pnorm(x, theta[5], theta[6])
}

Lg <- function(theta, x) {
  if (theta[1]+theta[2]>1) return(1e10)
  return(-sum(log(g(x, theta))))
}

R1.stat <- function(c, theta) {
  (theta[1]*(1-pnorm(c))+(1-theta[1]-theta[2])*(1-pnorm(c, theta[5], theta[6])))/(1-G(c, theta))}


# находим оценки максмиального правдоподобия
theta.stat <- mle(function(p0=0.9, p1=0.05, mu1=1.5, sd1=0.5, mu2=-1.5, sd2=1)
  Lg(c(p0, p1, mu1, sd1, mu2, sd2), x=t),
  lower = c(0, 0, 0, 0, -10, 0), 
  upper = c(1, 1, 10, 10, 0, 10),
  method = 'L-BFGS-B')@coef


# находим критическую константу из ограничения на d-риск 1-го рода
alpha <- 0.1
if (sum(theta.stat[1:2]) >= alpha){
  crit <- tryCatch({
    crit <- uniroot(function(x) R1.stat(x, theta.stat)-alpha, c(-5, 5), tol = 1e-10)$root},
    error = function(e) {
      crit <- -Inf
      return (crit)})
} else {crit <- Inf}

R.stat <- sum(t>crit)

#======================================================
