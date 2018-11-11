# analyze data

# this function fits the models for the nuisance parameters
fit_models <- function(DAT){
  rxb <- rms::orm(rxb ~ day + rxb_lag + oop_lag + V_lag + cumcost_lag,
                  data=DAT)
  
  oop1 <- glm(zero ~ day + rxb + oop_lag + V_lag + cumcost_lag,
              data=DAT, family=binomial(logit))
  
  oop2 <- lm(oop ~ day + rxb + oop_lag + V_lag + cumcost_lag,
             data=DAT[zero > 0])
  
  S <- glm(S ~ day + rxb + oop + V + cumcost + U,
           data=DAT, family=binomial(logit))
  return(list(rxb=rxb,oop1=oop1,oop2=oop2,S=S))
}

# this function computes the g-computation integral numerically
g_comp <- function(BSLN,MODELS,V_lag1,V_lag2,V,COST){
  n <- BSLN[,.N]
  
  rxb <- matrix(NA,nrow=n,ncol=270)
  oop <- matrix(NA,nrow=n,ncol=270)
  S   <- matrix(NA,nrow=n,ncol=270)
  
  rxb[,1] <- BSLN[,rxb]
  oop[,1] <- BSLN[,oop]
  S[,1] <- rep(0,n)
  
  for(t in 2:270){
    ###### rxb
    df <- data.frame(
      day=t,rxb_lag=rxb[,t-1],oop_lag=oop[,t-1],V_lag=V_lag1,cumcost_lag=COST[t-1]
    )
    rxb_probs <- predict(MODELS$rxb,newdata = df,type = "fitted.ind",codes = TRUE)
    rxb[,t] <- apply(rxb_probs,1,function(x){
      base::sample(c(paste(0:4)),1,prob=x)
    })
    
    ###### oop1
    zero <- rbinom(n, 1, 
                   prob = predict(MODELS$oop1,
                                  newdata = data.frame(
                                    day=t,
                                    rxb=rxb[,t],
                                    oop_lag=oop[,t-1],
                                    V_lag=V_lag1,
                                    cumcost_lag=COST[t-1]
                                  ), type = "response")
    )
    
    ###### oop2
    
    oop[,t] <- rnorm(n, mean = predict(MODELS$oop2,
                                       newdata = data.frame(
                                         day=t,
                                         rxb=rxb[,t],
                                         oop_lag=oop[,t-1],
                                         V_lag=V_lag2,
                                         cumcost_lag=COST[t-1]
                                       )),
                     sd = sqrt(summary(MODELS$oop2)$sigma)
    )
    
    oop[,t] <- (zero == 1)*oop[,t]
    
    ###### S
    S[,t] <- rbinom(n, 1, 
                    prob = predict(MODELS$S,
                                   newdata = data.frame(
                                     day=t,
                                     rxb=rxb[,t],
                                     oop=oop[,t],
                                     V=V,
                                     cumcost=COST[t],
                                     U=u_star
                                   ), type = "response")
    )
  }
  survival <- 1 - colMeans(S)
  return(survival)
}

# this wraps `fit_models` and `g_comp` into a function that
# computes the restricted mean difference of the 
# counterfactual survival curves
analyze <- function(DATA,BAND,NUMSIM,d,p){
  # fit models
  b_mods <- fit_models(DAT = DATA[U < u_star])
  
  g_mods <- fit_models(DAT = DATA[U >= u_star])
  
  # perform g computation
  bsln <- DATA[day==1][sample(.N, NUMSIM, replace=TRUE,
                              prob = dnorm(DATA[day==1]$U, mean=u_star, sd=BAND))] 
  
  cost <- sort(rep(seq(from = p, to = 270/d * p, by = p),d))
  
  b_V      <- factor("1", levels = na.omit(unique(DATA[U < u_star,V])))
  b_V_lag1 <- factor("1", levels = na.omit(unique(DATA[U < u_star,V_lag])))
  b_V_lag2 <- factor("1", levels = na.omit(unique(DATA[U < u_star & zero==1,V_lag])))
  
  g_V      <- factor("2", levels = na.omit(unique(DATA[U >= u_star,V])))
  g_V_lag1 <- factor("2", levels = na.omit(unique(DATA[U >= u_star,V_lag])))
  g_V_lag2 <- factor("2", levels = na.omit(unique(DATA[U >= u_star & zero==1,V_lag])))
  
  
  brand   <- g_comp(bsln[U < u_star], 
                    b_mods,
                    V_lag1 = b_V_lag1,
                    V_lag2 = b_V_lag2,
                    V = b_V,
                    COST = cost)
  generic   <- g_comp(bsln[U >= u_star], 
                      g_mods,
                      V_lag1 = g_V_lag1,
                      V_lag2 = g_V_lag2,
                      V = g_V,
                      COST = cost)
  
  # get restricted mean difference
  rmdiff <- sum(brand - generic)
  
  return(rmdiff)
}

# this function resamples the data by re-sampling whole patient histories
resample <- function(DATA){
  
  
  return(data_resampled)
}

# this function computes a 95% CI using bootstrap
bootstrap <- function(DATA, R){
  
  
  return(c(lower,upper))
}

