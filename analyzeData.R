# analyze data

# fit the models for the nuisance parameters
fit_models <- function(DAT){
  rxb <- glm(rxb ~ rxb_lag + I(V_lag==1) + I(V_lag==2),
                  data=DAT, family=binomial(logit))
  
  S <- glm(S ~ day + rxb + I(V==1) + I(V==2) + U,
           data=DAT, family=binomial(logit))
  return(list(rxb=rxb,S=S))
}

# compute the g-computation integral numerically
g_comp <- function(BSLN,MODELS){
  n <- BSLN[,.N]
  
  rxb <- matrix(NA,nrow=n,ncol=50)
  S   <- matrix(NA,nrow=n,ncol=50)
  
  rxb[,1] <- BSLN[,rxb]
  S[,1] <- rep(0,n)
  
  for(t in 2:50){
    # note that `rep(1,n)` is used where V should be
    ###### rxb
    X <- model.matrix(
      rep(2,n) ~ rxb[,t-1] + rep(1,n)
    )
    rxb[,t] <- rbinom(n,1,expit(X %*% na.omit(MODELS[["rxb"]]$coefficients)))

    ###### S
    X <- model.matrix(
      rep(2,n) ~ rep(t,n) + rxb[,t] + rep(1,n) + rep(u_star,n)
    )
    S[,t] <- rbinom(n,1,expit(X %*% na.omit(MODELS[["S"]]$coefficients)))
  }
  survival <- 1 - colMeans(S)
  return(survival)
}

# compute the restricted mean difference
analyze <- function(DATA,BAND,NUMSIM){
  # fit models
  b_mods <- fit_models(DAT = DATA[U < u_star])
  
  g_mods <- fit_models(DAT = DATA[U >= u_star])
  
  # perform g computation
  bsln <- DATA[day==1][sample(.N, NUMSIM, replace=TRUE,
                              prob = dnorm(DATA[day==1]$U, mean=u_star, sd=BAND))] 
  
  brand   <- g_comp(bsln[U < u_star],b_mods)
  generic <- g_comp(bsln[U >= u_star],g_mods)
  
  # get restricted mean difference
  rmdiff <- sum(brand - generic)
  
  return(rmdiff)
}

# compute a 95% CI using bootstrap
bootstrap <- function(DATA,R,BAND,NUMSIM){
  boots <- numeric(0)
  setkey(DATA,id)
  for(r in 1:R){
    ids_resampled <- sample(unique(DATA$id), length(unique(DATA$id)), replace = TRUE)
    data_resampled <- DATA[J(ids_resampled), allow.cartesian=TRUE]
    boots[r] <- analyze(DATA = data_resampled,BAND,NUMSIM)
  }
  return(quantile(boots,probs=c(.025,.975)))
}

