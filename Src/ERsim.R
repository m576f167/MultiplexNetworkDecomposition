library(sna)

aggnet  <- function(w, g1=ga, g2=gb, g3=gc) { # agg net given weights w {{{
  w <- w/sum(w)
  agg <-  w[1]*g1+w[2]*g2+w[3]*g3
} #}}}

tracefxn <- function(agg, stepprob=NULL) { #generate trace of variable length{{{
  # given (weighted) agg network
  s <- sample(1:ncol(agg),1)
  steps <- sample(3:25, 1, prob=stepprob)
  vsteps <- matrix(NA, nr=1, nc=25)
  vsteps[1]=s
  for (i in 2:steps) {
    cur <- vsteps[i-1]
    vsteps[i] <- sample(1:ncol(agg), 1, prob=agg[cur,])
  }
  vsteps
} #}}}

llktraces <- function(w, trace, g1, g2, g3) { # compute llk given w of trace {{{
# something to think about this weights each trace, not each transition equally
  eagg  <- aggnet(w)
  prob <- matrix(NA, nr=1, nc=24)
  for (i in 2:sum(is.finite(trace))) {
    prob[i-1] <- eagg[trace[i-1], trace[i]]
  }
  prob[prob==0] <- .001
  sum(log(prob), na.rm=T)
}# }}}

llkcalc <- function(w, data) { # compute llk given w and data {{{
  sum(sapply(1:nrow(data), function(x) llktraces(w, data[x,])))
} #}}}

mcmc <- function(data, iters, std) { # run mcmc iters times on data {{{
  x0 <- runif(3)
  x0 <- x0/sum(x0)
  lk_old <- llkcalc(x0, data)
  res <- matrix(NA, nc=8, nrow=iters)
  colnames(res) <- c('x1', 'x2', 'x3', 'tllk', 'x1star', 'x2star', 'x3star', 
                     'reject')
  x_old <- x0
  res[1,1:3] <- x0
  res[1,4] <- lk_old
  for (i in 2:iters) {
    x_old <- res[i-1,1:3] # pull the last used weights
    lk_old <- res[i-1,4]  # pull the llk of those weights
    x_star <- x_old+rnorm(3, mean=0, sd=std) # new weights
    x_star[x_star<0]=0 # prevent a weight from being negative (note this messes
                       # symmetry which we should fix later)
    x_star <-  x_star/sum(x_star) # renormalize weights, also may mess up sym.
    res[i,5:7] <- x_star # store that so we can make sure we're doing things right
    lk_star <- llkcalc(x_star, data) # compute llk of xstar
    if (lk_star > lk_old) { # always accept if new proposal is better
      res[i,1:3] <- x_star
      res[i,4] <- lk_star
      res[i,8] <- 0
    } else { # if new proposal isn't better sometimes we should still take it
      a <- lk_old # helps with overflow and underflow, shift mean of log
      st <- exp(lk_star-a)/(exp(lk_star-a)+1) 
      #prop of star / sum(prop of star+prop of old) but prop of old is exp(0)=1
      ol <- 1/(exp(lk_star)+1) # proportion of old/ total proportion
      alpha <- st/ol # ratio of the two proportions
      u <- runif(1)
      if (u <= alpha) {# accept x_star
        res[i,1:3] <- x_star
        res[i,4] <- lk_star
        res[i,8] <- 0
      } else { # reject x_star
        res[i,1:3] <- x_old
        res[i,4] <- lk_old
        res[i,8] <- 1
      }
    }
  }
  res <- res[-1,]
  res <- as.data.frame(res)
} #}}}

weights <- c(.5,.2, .3)
ga <- rgraph(100, tprob=.2, mode='graph')
gb <- rgraph(100, tprob=.2, mode='graph')
gc <- rgraph(100, tprob=.2, mode='graph')

agg <- aggnet(weights)

data <- t(sapply(1:100, function(x) tracefxn(agg)))


tt <- mcmc(data,500,.015) #std chosen so acceptance prob is around .44 
                          # (Gelman, Robers, Gilks 1996) Efficient Metropolis...


c1 <- mcmc(data,5000,.015)
c1 <- c1[500:4999,] # remove burn in
plot(c1$x1, type='l') # check for mixing (should look hairballish)
apply(c1[,1:3],2,mean) # parameter estimate

# now that I have found the right std, run for a long time and for multiple
# chains
mres <- sapply(1:10, function(x) {res <- mcmc(data,500,.015);
                                  apply(res[200:499,1:3],2,mean)})
apply(mres,1,mean)
