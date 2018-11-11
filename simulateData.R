# simulate data
sim_obs <- function(N,ssU,pR1,ssO,u_star){
  data <- data.frame(id=NA,day=NA,U=NA,V=NA,D=NA,C=NA,cumcost=NA,rxb=NA,oop=NA,S=NA)
  
  # simulate U, O1, R1 from empirical
  U1 <- predict(ssU, runif(N))$y
  
  #R1
  R1 <- sample(0:4, prob=pR1, size=N, replace=TRUE)
  
  #O1
  O1 <- rbinom(N, size=1, prob=0.89)
  n <- sum(O1)
  O11 <- predict(ssO, runif(n))$y
  O1 <- c(rep(0, N-n), O11)
  O1 <- sample(O1, size=N, replace=FALSE)
  
  # determine V1
  V1 <- ifelse(U1 < u_star, 1, 2)
  
  # simulate S1
  vars <- cbind(1,1,R1==1, R1==2, R1==3, R1==4, O1, 
                V1==1,V1==2,V1==3,V1==4,V1==5,V1==6,V1==7, 0, U1) 
  S1 <- rbinom(N, size=1, prob=expit(c(vars %*% par[[1]])))
  
  # generate data for each of the 1 to N subjects
  for(i in 1:N){
    V <- V.day <- D <- D.day <- numeric(0)
    C <- C.day <- cc <- cc.day <- numeric(0)
    R <- O <- S <- numeric(0)
    
    V[1] <- V1[i]
    V.day[1] <- V1[i]
    O[1] <- O1[i]
    R[1] <- R1[i]
    S[1] <- S1[i]
    
    predDk <- c(1, 1, R[1]==0:3, O[1], V[1]==2:7)
    phik <- exp(c(t(par[[7]]) %*% predDk))
    pDk1 <- 1/(1 + sum(phik))
    pDk <- pDk1*phik
    D[1] <- sample(c(1,15,30,60,90),size=1,prob=c(pDk1,pDk))
    D.day[1] <- D[1]
    
    predpCk <- c(1, 1, R[1]==0:3, O[1], V[1]==2:7, D[1]==c(15,30,60,90))
    pCk <- rbinom(1, size=1, prob=expit(sum(par[[8]]*predpCk)))
    C[1] <- 0
    if(pCk==1) {
      predCk <- c(1, 1, R[1]==0:3, O[1], V[1]==2:7, D[1]==c(15,30,60,90))
      C[1] <- rnorm(1, mean=sum(par[[9]]*predCk), sd=sqrt(13.67))
      #C[1] <- 0.5*exp(-C[1])*(exp(2*C[1])-1)
    } 
    cc[1] <- C[1]
    
    C.day[1] <- C[1]
    cc.day <- cc[1]
    
    
    k <- 1
    tk <- 1
    ### make rest of rows for subject i
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
        
        # simulate St, only if V[k] is actually current treatment
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
      
      ### simulate next values of Vk,Dk,Ck,cc_k
      ### but don't even bother if S==1 already
      if(sum(S) > 0){break}
      tk <- t  # compute t(k) for next value of k
      if(tk > 270){break}
      k <- k+1
      
      # simulate Vk
      if(U1[i] + tk - 1 < u_star){
        predVk <- c(1, tk, R[tk]==0:3, O[tk])
        phik <- exp(c(t(par[[5]]) %*% predVk))
        pVk0 <- 1/(1 + sum(phik))
        pVk <- pVk0*phik
        V[k] <- sample(c(0,1,3,7),size=1,prob=c(pVk0,pVk))
      }
      else{
        predVk <- c(1, tk, R[tk]==0:3, O[tk])
        phik <- exp(c(t(par[[6]]) %*% predVk))
        pVk0 <- 1/(1 + sum(phik))
        pVk <- pVk0*phik
        V[k] <- sample(c(0,1,2,3,4,5,6,7),size=1,prob=c(pVk0,pVk))
      }
      
      # simulate D and C
      if (V[k] == 0) {
        D[k] <- 1
        C[k] <- 0
      }
      else{
        # simulate D[k]
        predDk <- c(1, tk, R[tk]==0:3, O[tk], V[k]==2:7)
        phik <- exp(c(t(par[[7]]) %*% predDk))
        pDk1 <- 1/(1 + sum(phik))
        pDk <- pDk1*phik
        D[k] <- sample(c(1,15,30,60,90),size=1,prob=c(pDk1,pDk))
        
        # simulate Ck
        predpCk <- c(1, tk, R[tk]==0:3, O[tk], V[k]==2:7, D[k]==c(15,30,60,90))
        pCk <- rbinom(1, size=1, prob=expit(sum(par[[8]]*predpCk)))
        C[k] <- 0
        if(pCk==1) {
          predCk <- c(1, tk, R[tk]==0:3, O[tk], V[k]==2:7, D[k]==c(15,30,60,90))
          C[k] <- rnorm(1, mean=sum(par[[9]]*predCk), sd=sqrt(13.67))
          #C[k] <- 0.5*exp(-C[k])*(exp(2*C[k])-1)
        } 
      }
      
      # compute cumcost_k from Ck
      cc[k] <-  sum(C)
      
      # simulate S[tk]
      predSt <- c(1,tk,R[tk]==1:4, O[tk], V[k]==1:7, cc[k], U1[i]) 
      S[tk] <- rbinom(1, size=1, prob=expit(c(predSt %*% par[[1]])))
      
      V.day[tk] <- V[k]
      D.day[tk] <- D[k]
      C.day[tk] <- C[k]
      cc.day[tk] <- cc[k]
      
      if(S[tk]==1){break}
    }
    # assemble data for patient i
    #(id=NA,day=NA,U=NA,V=NA,D=NA,C=NA,cumcost=NA,rxb=NA,oop=NA,S=NA)
    data <- rbind(data,
                  cbind(id=i,day=1:length(R),
                        U=U1[i],V=V.day,D=D.day,C=C.day,cumcost=cc.day,
                        rxb=R,oop=O,S=S))
    #print(c(i,length(V.day),length(D.day),length(C.day),length(cc.day),
    #      length(R),length(O),length(S)))
    
  }
  setDT(data)
  data[,zero := as.numeric(C > 0)]
  
  lagcols <- c("V","cumcost","rxb","oop","zero")
  data[,(paste0(lagcols,"_lag")) := shift(.SD),id,.SDcols=lagcols]
  
  data[,rxb := as.character(rxb)]
  data[,rxb_lag := as.character(rxb_lag)]
  data[,V := factor(as.character(V),levels = paste(0:7))]
  data[,V_lag := factor(as.character(V_lag),levels = paste(0:7))]
  
  return(data[-1])
}
