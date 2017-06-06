
#' @title Bimodal Skew Symmetric Normal Distribution
#' @description These functions define the Bimodal Skew Symmetric Normal Distribution. This is a four parameter distribution and can be used to fit a GAMLSS model.The functions dBSSN, pBSSN, qBSSN and rBSSN define the probability distribution function, the cumulative distribution function,  the inverse cumulative distribution functions and  the random generation for the Bimodal Skew Symmetric Normal Distribution; respectively.
#' @param x,q Vector of quantiles
#' @param mu Vector of location parameter values
#' @param sigma Vector of scale parameter values
#' @param nu Vector of nu parameter values
#' @param tau Vector of bimodality parameter values
#' @param mu.link Defines the mu.link, with identity link as the default for the mu parameter
#' @param sigma.link Defines the sigma.link, with log link as the deafult for the sigma parameter
#' @param nu.link Defines the nu.link, with identity link as the default for the nu parameter
#' @param tau.link Defines the tau.link, with log link as the default for the tau parameter
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param p Vector of probabilities
#' @param n number of observations; if length(n) > 1, the length is taken to be the number required
#' @details The probability density function of the BSSN distribution is given by \deqn{f_Y(y|\mu, \sigma, \nu, \tau)= c[\tau + (y-\nu)^(2)]e^{-\sigma(y-\mu)^(2)}} for \eqn{-\infty < y < \infty}, where \eqn{c = 2\sigma^(3/2) / \gamma \sqrt\pi},  \eqn{\gamma= 1 + 2 \sigma \theta},  \eqn{\theta = \tau + \delta^{2}},  \eqn{\delta= \nu - \mu}. \eqn{-\infty <\mu <\infty} and  \eqn{-\infty < \nu < \infty} are location parameters and \eqn{\sigma >0} and \eqn{\tau \geq 0 } denote the scale and bimodality parameters respectively.
#' @references Hassan, M. Y. and  El-Bassiouni M. Y. (2015). Bimodal skew-symmetric normal distribution,\emph{Communications in Statistics-Theory and Methods}, \bold{45}, part 5, pp 1527--1541.
#'
#' Hossain, A.Rigby, R. A. Stasinopoulos D. M. and Enea, M. A flexible approach for modelling proportion response variable:LGD, \emph{31st International workshop for Statistical Modelling Society},\bold{1}, pp 127--132.
#'
#' Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), \emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.
#' @examples op<-par(mfrow=c(3,3))
#' curve(dBSSN(x,  mu=1, sigma=0.1, nu=1, tau=1),-12, 12, ylab="f(x)", main="BSSN")
#' curve(dBSSN(x,  mu=1, sigma=0.1, nu=1, tau=5),-12, 12,ylab="f(x)", main="BSSN")
#' curve(dBSSN(x,  mu=1, sigma=0.1, nu=1, tau=10),-12, 12, ylab="f(x)", main="BSSN")
#' curve(dBSSN(x,  mu=1, sigma=0.1, nu=1, tau=20),-12, 12, ylab="f(x)", main="BSSN")
#' curve(dBSSN(x,  mu=1, sigma=0.1, nu=0, tau=4),-12, 12, ylab="f(x)", main="BSSN")
#' curve(dBSSN(x,  mu=-1, sigma=0.1, nu=0, tau=3),-12, 12, ylab="f(x)", main="BSSN")
#' curve(dBSSN(x,  mu=1, sigma=0.1, nu=2, tau=0),-12, 12, ylab="f(x)", main="BSSN")
#' curve(dBSSN(x,  mu=-1, sigma=0.1, nu=-2, tau=0),-12, 12, ylab="f(x)", main="BSSN")
#' curve(dBSSN(x,  mu=-1, sigma=0.1, nu=-3, tau=0.8),-12, 12, ylab="f(x)", main="BSSN")
#' par(op)
#' @return NULL
#' @importFrom  gamlss.dist checklink
#' @importFrom  stats uniroot runif
#' @importFrom  gamlss.dist dNO pNO
#' @rdname BSSN
#' @export
BSSN <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
        mstats <- checklink(   "mu.link", "Bimodal skew-symmetric normal", substitute(mu.link),
                               c("inverse", "log", "identity", "own"))
        dstats <- checklink("sigma.link", "Bimodal skew-symmetric normal", substitute(sigma.link),
                            c("inverse", "log", "identity", "own"))
        vstats <- checklink(   "nu.link", "Bimodal skew-symmetric normal", substitute(nu.link),
                               c("inverse", "log", "identity", "own"))
        tstats <- checklink(  "tau.link", "Bimodal skew-symmetric normal", substitute(tau.link),
                              c("inverse", "log", "identity", "own"))
        structure(
                list(family = c("BSSN", "Bimodal skew-symmetric normal"),
                     parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE),
                     nopar = 4,
                     type = "Continuous",
                     mu.link = as.character(substitute(mu.link)),
                     sigma.link = as.character(substitute(sigma.link)),
                     nu.link = as.character(substitute(nu.link)),
                     tau.link = as.character(substitute(tau.link)),
                     mu.linkfun = mstats$linkfun,
                     sigma.linkfun = dstats$linkfun,
                     nu.linkfun = vstats$linkfun,
                     tau.linkfun = tstats$linkfun,
                     mu.linkinv = mstats$linkinv,
                     sigma.linkinv = dstats$linkinv,
                     nu.linkinv = vstats$linkinv,
                     tau.linkinv = tstats$linkinv,
                     mu.dr = mstats$mu.eta,
                     sigma.dr = dstats$mu.eta,
                     nu.dr = vstats$mu.eta,
                     tau.dr = tstats$mu.eta,
                     dldm = function(y, mu, sigma, nu, tau)
                     {     lambda <- nu-mu
                           theta <- tau + (lambda)^2
                           gamma <- 1 + 2*sigma*theta
                           c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                           dldm<- 4*sigma*lambda/gamma + 2*sigma*(y-mu)

                           dldm
                     },
                     d2ldm2 = function(y, mu, sigma, nu, tau)
                     {   lambda <- nu-mu
                         theta <- tau + (lambda)^2
                         gamma <- 1 + 2*sigma*theta
                         c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                         dldm<- 4*sigma*lambda/gamma + 2*sigma*(y-mu)
                             d2ldm2<- - (dldm) * (dldm)

                             d2ldm2
                     },


                     dldd= function(y,mu, sigma, nu, tau)
                     {      lambda <- nu-mu
                            theta <- tau + (lambda)^2
                            gamma <- 1 + 2*sigma*theta
                            c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                            dldd <- 3/(2*sigma) - 2*theta/gamma - (y-mu)^2
                            dldd
                     },

                     d2ldd2 = function(y, mu, sigma, nu, tau)
                     {       lambda <- nu-mu
                             theta <- tau + (lambda)^2
                             gamma <- 1 + 2*sigma*theta
                             c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                             dldd <- 3/(2*sigma) - 2*theta/gamma - (y-mu)^2
                             d2ldd2<- -(dldd) * (dldd)

                             d2ldd2
                     },

                     dldv= function(y,mu, sigma, nu, tau)
                     {       lambda <- nu-mu
                             theta <- tau + (lambda)^2
                             gamma <- 1 + 2*sigma*theta
                             c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                             dldv <- -4*sigma*lambda/gamma - 2*(y-nu)/(tau + (y-nu)^2)
                             dldv
                     },

                     d2ldv2 = function(y, mu, sigma, nu, tau)
                     {       lambda <- nu-mu
                             theta <- tau + (lambda)^2
                             gamma <- 1 + 2*sigma*theta
                             c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                             dldv <- -4*sigma*lambda/gamma - 2*(y-nu)/(tau + (y-nu)^2)
                             d2ldv2<- -(dldv) * (dldv)

                             d2ldv2
                     },

                     dldt= function(y,mu, sigma, nu, tau)
                     {       lambda <- nu-mu
                             theta <- tau + (lambda)^2
                             gamma <- 1 + 2*sigma*theta
                             c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                             dldt <- -2*sigma/gamma + 1/(tau + (y-nu)^2)
                             dldt
                     },

                     d2ldt2 = function(y, mu, sigma, nu, tau)
                     {       lambda <- nu-mu
                             theta <- tau + (lambda)^2
                             gamma <- 1 + 2*sigma*theta
                             c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                             dldt <- -2*sigma/gamma + 1/(tau + (y-nu)^2)
                             d2ldt2<- -(dldt) * (dldt)

                             d2ldt2
                     },

                     d2ldmdd = function(y, mu, sigma, nu, tau)
                     {       lambda <- nu-mu
                             theta <- tau + (lambda)^2
                             gamma <- 1 + 2*sigma*theta
                             c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                             dldm<- 4*sigma*lambda/gamma + 2*sigma*(y-mu)
                             dldd <- 3/(2*sigma) - 2*theta/gamma - (y-mu)^2
                             d2ldmdd<- -(dldm) * (dldd)

                             d2ldmdd
                     },

                     d2ldmdv = function(y, mu, sigma, nu, tau)
                     {       lambda <- nu-mu
                             theta <- tau + (lambda)^2
                             gamma <- 1 + 2*sigma*theta
                             c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                             dldm<- 4*sigma*lambda/gamma + 2*sigma*(y-mu)
                             dldv <- -4*sigma*lambda/gamma - 2*(y-nu)/(tau + (y-nu)^2)

                             d2ldmdv<- -(dldm) * (dldv)

                             d2ldmdv
                     },

                     d2ldmdt = function(y, mu, sigma, nu, tau)
                     {       lambda <- nu-mu
                             theta <- tau + (lambda)^2
                             gamma <- 1 + 2*sigma*theta
                             c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                             dldm<- 4*sigma*lambda/gamma + 2*sigma*(y-mu)
                             dldt <- -2*sigma/gamma + 1/(tau + (y-nu)^2)

                             d2ldmdt<- -(dldm) * (dldt)

                             d2ldmdt
                     },

                     d2ldddv = function(y, mu, sigma, nu, tau)
                     {       lambda <- nu-mu
                             theta <- tau + (lambda)^2
                             gamma <- 1 + 2*sigma*theta
                             c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                             dldd <- 3/(2*sigma) - 2*theta/gamma - (y-mu)^2
                             dldv <- -4*sigma*lambda/gamma - 2*(y-nu)/(tau + (y-nu)^2)
                             d2ldddv<- -(dldd) * (dldv)

                             d2ldddv
                     },

                     d2ldddt = function(y, mu, sigma, nu, tau)
                     {       lambda <- nu-mu
                             theta <- tau + (lambda)^2
                             gamma <- 1 + 2*sigma*theta
                             c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                             dldd <- 3/(2*sigma) - 2*theta/gamma - (y-mu)^2
                             dldt <- -2*sigma/gamma + 1/(tau + (y-nu)^2)
                             d2ldddt<- -(dldd) * (dldt)

                             d2ldddt
                     },

                     d2ldvdt = function(y, mu, sigma, nu, tau)
                     {       lambda <- nu-mu
                             theta <- tau + (lambda)^2
                             gamma <- 1 + 2*sigma*theta
                             c <- 2 *(sigma)^(3/2)/(gamma*(pi)^(1/2))
                             dldv <- -4*sigma*lambda/gamma - 2*(y-nu)/(tau + (y-nu)^2)
                             dldt <- -2*sigma/gamma + 1/(tau + (y-nu)^2)
                             d2ldvdt<- -(dldv) * (dldt)

                             d2ldvdt
                     },

                     G.dev.incr  = function(y,mu,sigma,nu,tau,...)
                     {
                             -2*dBSSN(y,mu,sigma,nu,tau,log=TRUE)
                     } ,
                     rqres = expression(
                             rqres(pfun="pBSSN", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
                     mu.initial = expression(mu <- (y+mean(y))/2),
                     sigma.initial = expression(sigma <- rep(0.1, length(y))),
                     nu.initial = expression(nu <- rep(1, length(y))),
                     tau.initial = expression(tau <-rep(1, length(y))),
                     mu.valid = function(mu) TRUE,
                     sigma.valid = function(sigma) all(sigma > 0),
                     nu.valid = function(nu) TRUE,
                     tau.valid = function(tau) all(tau > 0),
                     y.valid = function(y)  TRUE
                ),
                class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
#' @rdname BSSN
#' @export
dBSSN <- function(x, mu = 0, sigma = 1, nu = 1, tau = .5, log = FALSE)
{
        if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
        if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))
        lambda <- nu-mu
        theta <- tau + (lambda)^2
        gamma1 <- 1 + 2*sigma*theta
        c <- 2 *(sigma)^(3/2)/(gamma1*(pi)^(1/2))
        d <- c*(tau+(x-nu)^2)*exp(-sigma*(x-mu)^2)
        loglik <- log(c)+log(tau+(x-nu)^2)- sigma*(x-mu)^2
        if(log==FALSE) ft  <- exp(loglik) else ft <- loglik
        ft
}
#-----------------------------------------------------------------
# pfun needs to check
#' @rdname BSSN
#' @export
pBSSN <- function(q, mu = 0, sigma = 1, nu = 1, tau = .5, lower.tail = TRUE, log.p = FALSE,log=T)
{
        if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))
       # lambda <- nu-mu
      #  theta <- tau + (lambda)^2
      #  gamma <- 1 + 2*sigma*theta
      #      c <- (2 *(sigma)^(3/2))/(gamma*((pi)^(1/2)))
           ax <- (q + mu-(2*nu))/(1 + (2*sigma)*(tau+(nu-mu)^2))
        #  dNx <- dNO(q,0,1)
        #  pNx <- pNO(q,0,1)
            p <- pNO(q, mu=mu, sigma=sqrt(1/(2*sigma))) - ax * dNO(q, mu=mu, sigma=sqrt(1/(2*sigma)))
        if(lower.tail==TRUE) p  <- p else  p <- 1-p
        if(log.p==FALSE) p  <- p else  p <- log(p)
        p
}
#-----------------------------------------------------------------
# needs to get q function
#' @rdname BSSN
#' @export
qBSSN <-  function(p, mu = 0, sigma = 1, nu = 1, tau = .5, lower.tail = TRUE, log.p = FALSE)
{
        #---functions--------------------------------------------
        h1 <- function(q)
        {
                pBSSN(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i]) - p[i]
        }
        h <- function(q)
        {
                pBSSN(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i])
        }
        #-----------------------------------------------------------------
        #if (any(mu <= 0))  stop(paste("mu must be positive", "\n", ""))
        if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
        if (log.p==TRUE) p <- exp(p) else p <- p
        if (lower.tail==TRUE) p <- p else p <- 1-p
        if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))
        lp <-  max(length(p),length(mu),length(sigma),length(nu), length(tau))
        p <- rep(p, length = lp)
        sigma <- rep(sigma, length = lp)
        mu <- rep(mu, length = lp)
        nu <- rep(nu, length = lp)
        tau <- rep(tau, length = lp)
        q <- rep(0,lp)
        for (i in  seq(along=p))
        {
                if (h(mu[i])<p[i])
                {
                        interval <- c(mu[i], mu[i]+sigma[i])
                        j <-2
                        while (h(interval[2]) < p[i])
                        {interval[2]<- mu[i]+j*sigma[i]
                         j<-j+1
                        }
                }
                else
                {
                        interval <-  c(mu[i]-sigma[i], mu[i])
                        j <-2
                        while (h(interval[1]) > p[i])
                        {interval[1]<- mu[i]-j*sigma[i]
                         j<-j+1
                        }
                }
                q[i] <- uniroot(h1, interval)$root
                #interval <- c(.Machine$double.xmin, 20)
        }
        q
}
#-----------------------------------------------------------------
#' @rdname BSSN
#' @export
rBSSN <- function(n, mu=0, sigma=1, nu=1, tau=.5)
{
        if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
        if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))
        n <- ceiling(n)
        p <- runif(n)
        r <- qBSSN(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
        r
}
#-----------------------------------------------------------------
