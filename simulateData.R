# simulate data
sim_obs <- function(N,ssU,pR1,u_star){
  data <- data.frame(id=NA,U=NA,day=NA,V=NA,rxb=NA,S=NA)
  
  # simulate U, O1, R1 from empirical
  U1 <- predict(ssU, runif(N))$y
  
  #R1
  R1 <- rbinom(N, size=1, prob=pR1)
  
  # determine V1
  V1 <- ifelse(U1 < u_star, 1, 2)
  
  # simulate S1
  vars <- cbind(1,1,V1==1,V1==2,R1,U1)
  S1 <- rbinom(N, size=1, prob=expit(c(vars %*% par[[1]])))
  
  # generate data for each of the 1 to N subjects
  for(i in 1:N){
    V <- R <- S <- numeric(0)
    V[1] <- V1[i]
    R[1] <- R1[i]
    S[1] <- S1[i]
    
    t <- 2
    ### make rest of rows for subject i
    repeat{
      # simulate R(t)|V(t-1),R(t-1)
      vars <- c(1,V[t-1]==1,V[t-1]==2,R[t-1])
      R[t] <- rbinom(1,1,expit(c(vars %*% par[[2]])))
      
      # simulate V(t)|R(t)
      vars <- c(1,R[t])
      if(U1[i] < u_star){
        V[t] <- rbinom(1,1,expit(c(vars %*% par[[3]])))
      }
      else{
        V[t] <- 2*rbinom(1,1,expit(c(vars %*% par[[4]])))
      }
      
      # simulate S(t)|day,V(t),R(t),U
      vars <- c(1,t,V[t]==1,V[t]==2,R[t],U1[i])
      S[t] <- rbinom(1,1,expit(c(vars %*% par[[1]])))
      
      if(S[t]==1 | t+1 > 50){break}
      t <- t+1
    }
    # assemble data for patient i
    #(id,U,day,V,rxb=NA,S=NA)
    data <- rbind(data,
                  cbind(id=i,
                        U=U1[i],
                        day=1:length(R),
                        V=V,
                        rxb=R,
                        S=S))
  }
  setDT(data)
  
  lagcols <- c("V","rxb")
  data[,(paste0(lagcols,"_lag")) := shift(.SD),id,.SDcols=lagcols]
  
  return(data[-1])
}





