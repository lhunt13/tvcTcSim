# this function simulates counterfactual data
# with U=u* and continuous adherance defined by p and d
sim_cf <- function(N,ssU,pR1,ssO,u_star,venform,p,d){
  data <- data.frame(id=NA,day=NA,U=NA,V=NA,D=NA,C=NA,cumcost=NA,rxb=NA,oop=NA,S=NA)
  U1 <- rep(u_star,N)
  
  #R1
  R1 <- sample(0:4, prob=pR1, size=N, replace=TRUE)
  
  #O1
  O1 <- rbinom(N, size=1, prob=0.89)
  n <- sum(O1)
  O11 <- predict(ssO, runif(n))$y
  O1 <- c(rep(0, N-n), O11)
  O1 <- sample(O1, size=N, replace=FALSE)
  
  # determine V1
  V1 <- rep(venform,N)
  
  # simulate S1
  vars <- cbind(1,1,R1==1, R1==2, R1==3, R1==4, O1, 
                V1==1,V1==2,V1==3,V1==4,V1==5,V1==6,V1==7, 0, U1) 
  S1 <- rbinom(N, size=1, prob=expit(c(vars %*% par[[1]])))
  
  for(i in 1:N){
    V <- V.day <- D <- D.day <- numeric(0)
    C <- C.day <- cc <- cc.day <- numeric(0)
    R <- O <- S <- numeric(0)
    
    V[1] <- V1[i]
    V.day[1] <- V1[i]
    O[1] <- O1[i]
    R[1] <- R1[i]
    S[1] <- S1[i]
    
    D[1] <- d
    D.day[1] <- D[1]
    C[1] <- p
    cc[1] <- C[1]
    C.day[1] <- C[1]
    cc.day <- cc[1]
    
    k <- 1
    tk <- 1
    repeat{
      if(S[tk]==1 | tk+1 > 270){break}
      t <- tk+1
      repeat{
        if(t > 270){break}
        # simulate Rt
        predRt <- c(1, t, R[t-1]==1:4, O[t-1], V[k]==1:7, cc[k])
        #phij <- 1/(1+exp(c(t(par[[2]]) %*% predRt)))
        phij <- expit(c(t(par[[2]]) %*% predRt))
        #pRt <- c(phij[1], diff(phij))
        pRt <- c(1-phij[1],phij[1]-phij[2],phij[2]-phij[3],phij[3]-phij[4],phij[4])
        R[t] <- sample(0:4, size=1, prob=pRt)
        
        # simulate Ot
        predpOt <- c(1, tk, R[t]==1:4, O[t-1], V[k]==1:7, cc[k])
        pOt <- rbinom(1, size=1, prob=expit(sum(par[[3]]*predpOt)))
        O.temp <- 0 #use temp because need to preserve last value of Ot
        if(pOt==1) {
          predOt <- c(1, tk, R[t]==1:4, O[t-1], V[k]==1:7, cc[k])
          O.temp <- rnorm(1, mean=sum(par[[4]]*predOt), sd=sqrt(0.013))
        }
        O[t] <- O.temp
        
        if(t < tk+D[k]){
          # record day level trt values having impact on the above sims of O,R,S
          V.day[t] <- V[k]
          D.day[t] <- D[k]
          C.day[t] <- 0 # 0 because no Rx was filled on any day inside the fill period
          cc.day[t] <- cc[k]
          
          predSt <- c(1,t,R[t]==1:4, O[t], V[k]==1:7, cc[k], U1[i]) 
          S[t] <- rbinom(1, size=1, prob=expit(c(predSt %*% par[[1]])))
          
          # if failure occurs here, just exit this repeat loop
          if(S[t]==1){break}
          
        } #else, still need to get S[t], but after getting V[k+1]
        
        if(t >= tk+D[k]){break}
        else{t <- t+1}
      }
      if(sum(S) > 0){break}
      tk <- t  # compute t(k) for next value of k
      if(tk > 270){break}
      k <- k+1
      
      # simulate treatment
      V[k] <- venform
      D[k] <- d
      C[k] <- p
      cc[k] <- sum(C)
      
      predSt <- c(1,tk,R[tk]==1:4, O[tk], V[k]==1:7, cc[k], U1[i]) 
      S[tk] <- rbinom(1, size=1, prob=expit(c(predSt %*% par[[1]])))
      
      V.day[tk] <- V[k]
      D.day[tk] <- D[k]
      C.day[tk] <- C[k]
      cc.day[tk] <- cc[k]
      
      if(S[tk]==1){break}
    }
    data <- rbind(data,
                  cbind(id=i,day=1:length(R),
                        U=U1[i],V=V.day,D=D.day,C=C.day,cumcost=cc.day,
                        rxb=R,oop=O,S=S))
  }
  return(data[-1,])
}

# this function computes the true value of the target parameter
get_truth <- function(n,p,d){
  brand <- sim_cf(n,ssU,pR1,ssO,u_star,1,p,d)
  generic <- sim_cf(n,ssU,pR1,ssO,u_star,2,p,d)

  setDT(brand)
  setDT(generic)
  
  truth <- sum(brand[,.(.N/n),.(day)]$V1 - 
                 generic[,.(.N/n),.(day)]$V1)
  
  return(truth)
}





