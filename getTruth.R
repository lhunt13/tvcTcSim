# this function simulates counterfactual data
# with U=u* and continuous adherance

sim_cf <- function(N,ssU,pR1,u_star,venform){
  data <- data.frame(id=NA,U=NA,day=NA,V=NA,rxb=NA,S=NA)
  
  # simulate U, O1, R1 from empirical
  U1 <- rep(u_star,N)
  
  #R1
  R1 <- rbinom(N, size=1, prob=pR1)
  
  # determine V1
  V1 <- rep(venform,N)
  
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
      V[t] <- venform
      
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
  
  return(data[-1])
}

# this function computes the true value of the target parameter
get_truth <- function(n,ssU,pR1,u_star){
  brand <- sim_cf(n,ssU,pR1,u_star,1)
  generic <- sim_cf(n,ssU,pR1,u_star,2)

  setDT(brand)
  setDT(generic)
  
  truth <- sum(brand[,.(.N/n),.(day)]$V1 - 
                 generic[,.(.N/n),.(day)]$V1)
  
  return(truth)
}






