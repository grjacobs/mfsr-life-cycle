
model {
  # Randomize unknown independent variables
  x.sst[32] ~ dnorm(0,3)
  x.upapr[32] ~ dnorm(0,3)
  for (u in 1:6){x.ph[u] ~ dnorm(0,1/sqrt(3))}

  #####################
  # Smolt production  # 
  #####################
  
  # per-capita fecundity * egg survival rate * parr survival rate
    log.b ~ dgamma(64,14)
    sigma.b.seg ~ dt(0, pow(2.5,-2), 1)T(0,)
    sigma.b.yr ~ dt(0, pow(2.5,-2), 1)T(0,)
    sigma.b.seg.yr ~ dt(0, pow(2.5,-2), 1)T(0,)
    tau.b.seg <- pow(sigma.b.seg, -2)
    tau.b.yr <- pow(sigma.b.yr, -2)
    tau.b.seg.yr <- pow(sigma.b.seg.yr, -2)
      for (t in 1:Y) {
        alpha.b.yr[t] ~ dnorm(0, tau.b.yr)
      }
    for (i in 1:s){
        alpha.b.seg[i] ~ dnorm(0, tau.b.seg)
        b.seg[i] <- exp(log.b + alpha.b.seg[i])
      for (t in 1:Y) {
        alpha.b.seg.yr[i,t] ~ dnorm(0, tau.b.seg.yr)
        b[i,t] <- exp(log.b + 
          alpha.b.seg[i] + alpha.b.yr[t] + alpha.b.seg.yr[i,t])
      }
    }
  
  ###################
  # Smolt transport # 
  ###################
  
    for (u in 1:2) {
      beta.del[u] ~ dnorm(0, 1/sqrt(2))
    }
    sigma.del ~ dt(0, pow(2.5,-2), 1)T(0,) 
    tau.del <- pow(sigma.del, -2)
    r.del ~ dt(0, pow(2.5,-2), 1)T(0,) 
    for (t in 1:Y){
      alpha.del[t] ~ dnorm(0, tau.del)
      logit(del[t]) <- beta.del[1] + beta.del[2]*T.sched[t] + alpha.del[t]
      y.del[t] ~ dbeta(r.del*del[t], r.del*(1-del[t]))
    }
    
  #################
  # Age at return #
  #################
  
  # Probability of returning at age 4
    mu.rho3 ~ dnorm(0, 1/sqrt(2)) #weakly informative prior centered at 50%
    sigma.rho3 ~ dt(0, pow(2.5,-2), 1)T(0,) 
    tau.rho3 <- pow(sigma.rho3, -2)
    logit(rho3[1]) <- mu.rho3 # initial value
    for(t in 2:(Y+2)){ # need to add 2 to Y to predict last 2 cohorts' survival
      alpha.rho3[t-1] ~ dnorm(0, tau.rho3)
      logit(rho3[t]) <- logit(rho3[t-1]) + alpha.rho3[t-1]
    }
  # Probability of returning at age 5 
    mu.rho4 ~ dnorm(2, 1/sqrt(2)) #weakly informative prior centered at 88%
    logit(rho4) <- mu.rho4 
    
  ######################
  # Survival Constants #
  ######################
  
  # Constants
  phi.B <- 0.98 # Barge smolt survival
  phi.A <- 0.80 # Annual survival at sea
  
  #####################
  # In-river survival #
  #####################
  
    # In-river group only
      beta.H[1] ~ dnorm(0, 1) #prec=1/sqrt(1))
      for (u in 2:3) {
        beta.H[u] ~ dnorm(0, 1/sqrt(2))
      }
      sigma.H ~ dt(0, pow(2.5,-2), 1)T(0,) 
      tau.H <- pow(sigma.H, -2)
      r.H ~ dt(0, pow(2.5,-2), 1)T(0,) 
      for (c in 1:Y){
        alpha.H[c] ~ dnorm(0, tau.H)
        logit(phi.H[c]) <- beta.H[1] + 
          beta.H[2]*x.wtt[c] + 
          beta.H[3]*x.ph[c] + 
          alpha.H[c]
        logit(phi.H.pred[c]) <- beta.H[1] + 
          beta.H[2]*x.wtt[c] + 
          beta.H[3]*x.ph[c]
        y.H[c] ~ dbeta(r.H*phi.H[c], r.H*(1-phi.H[c]))
      }
      
  ########################
  # Early ocean survival #
  ########################

    # In-river group: all covariates + random year effect
      for (u in 1:6){
        beta.O[u] ~ dnorm(0, 1/sqrt(2))
      }
      sigma.O ~ dt(0, pow(2.5, -2), 1)T(0,)
      tau.O <- pow(sigma.O, -2)
      for (c in 1:Y){
        alpha.O[c] ~ dnorm(0, tau.O)
        logit(phi.O[c]) <- beta.O[1] + 
          beta.O[2]*x.sst[c] + 
          beta.O[3]*x.upapr[c] + 
          beta.O[4]*x.pdomjj[c] + 
          beta.O[5]*x.wtt[c] + 
          beta.O[6]*x.ph[c] + 
          alpha.O[c]
        logit(phi.O.pred[c]) <- beta.O[1] + 
          beta.O[2]*x.sst[c] + 
          beta.O[3]*x.upapr[c] + 
          beta.O[4]*x.pdomjj[c] + 
          beta.O[5]*x.wtt[c] + 
          beta.O[6]*x.ph[c] 
      }

    # Transported group: all covariates + random year effect 
      for (u in 1:4){
        beta.T[u] ~ dnorm(0, 1/sqrt(2))
      }
      sigma.T ~ dt(0, pow(2.5, -2), 1)T(0,)
      tau.T <- pow(sigma.T, -2)
      for (c in 1:Y){
        alpha.T[c] ~ dnorm(0, tau.T)
        logit(phi.T[c]) <- beta.T[1] + 
          beta.T[2]*x.sst[c] + 
          beta.T[3]*x.upapr[c] + 
          beta.T[4]*x.pdomjj[c] +
          alpha.T[c] 
        logit(phi.T.pred[c]) <- beta.T[1] + 
          beta.T[2]*x.sst[c] + 
          beta.T[3]*x.upapr[c] + 
          beta.T[4]*x.pdomjj[c]
      }
      
  #########################    
  # Adult return survival #
  #########################
  
    beta.R ~ dnorm(0, 1/sqrt(2))
    logit(phi.R) <- beta.R 
 
  ###################################
  # PIT tag observation probability #
  ###################################
  
    # Objects in this section are indexed by k=recapture occasion and 
    # c=migration year cohort. 
    
    # In-river group
      for (k in 1:(n.occ-1)){
        mean.p[k] ~ dnorm(0, 1/sqrt(2))
        sigma.p[k] ~ dt(0, pow(2.5, -2), 1)T(0,)
        tau.p[k] <- pow(sigma.p[k], -2)
        for (c in 1:6){ #set p to 0 prior to 1994
          p[k,c] <- 0
        } #c
        for (c in 7:Y){
          alpha.p[k,c] ~ dnorm(0, tau.p[k])
          logit(p[k,c]) <- mean.p[k] + alpha.p[k,c]
        } #c
      } #k
    # Transported group (assume no difference in adult obs prob)
      for (c in 1:Y){
        p.T[1,c] <- p[2,c]
        p.T[2,c] <- p[3,c]
      } #c
  
  #############################################    
  # CJS multinomial formulation using m-array #
  #############################################
  
    # Note on indexes in this section only: j=release occasion, k=recapture 
    # occation, and c=migration year cohort. In the manuscript and in the rest 
    # of the model code, release occasion is indexed by j and recapture occasion
    # is indexed by k. Migration year is indexed by c elsewhere in the model
  
  	# Cohort survival vectors for multinomial CJS
  		for (c in 1:(Y)){
  			# For convenience, cohort whole ocean survival
  			Smarine[c] <- phi.A*phi.A*phi.A*(1-rho3[c+2])*(1-rho4)+
  			  phi.A*phi.A*(1-rho3[c+2])*rho4+phi.A*rho3[c+2]
  			# Survival vectors
  			# Conditional estimates wrt random year effects 
  			pi.H[1,c] <- phi.H[c]
  			pi.H[2,c] <- phi.O[c]*Smarine[c]
  			pi.H[3,c] <- phi.R
  			pi.T[1,c] <- phi.B*phi.T[c]*Smarine[c]
  			pi.T[2,c] <- phi.R
  			# Unconditional estimates wrt random year effects 
  			pi.H.pred[1,c] <- phi.H.pred[c]
  			pi.H.pred[2,c] <- phi.O.pred[c]*Smarine[c]
  			pi.H.pred[3,c] <- phi.R
  			pi.T.pred[1,c] <- phi.B*phi.T.pred[c]*Smarine[c]
  			pi.T.pred[2,c] <- phi.R
  		}
  		
    # In-river group
      for (j in 1:(n.occ-1)){
        for (c in 7:Y){
          marr[j,1:n.occ,c] ~ dmulti(pr[j, ,c], rel[j,c])
        } #c
      } #j
      # Define the cell probabilities of the m-array
      # Main diagonal
      for (j in 1:(n.occ-1)){
        for (c in 7:Y){
          q[j,c] <- 1-p[j,(c)]                # Probability of non-recapture
          pr[j,j,c] <- pi.H[j,(c)]*p[j,(c)]
          # Above main diagonal
          for (k in (j+1):(n.occ-1)){
            pr[j,k,c] <- prod(pi.H[j:k,(c)])*prod(q[j:(k-1),c])*p[k,(c)]
          } #k
          # Below main diagonal
          for (k in 1:(j-1)){
            pr[j,k,c] <- 0
          } #k
        } #c
      } #j
      # Last column: probability of non-recapture
      for (j in 1:(n.occ-1)){
        for (c in 7:Y){
          pr[j,n.occ,c] <- 1-sum(pr[j,1:(n.occ-1),c])
        } #c
      } #j
      
    # Transported group 
      for (j in 1:(n.occT-1)){
        for (c in 7:Y){
          marr.T[j,1:n.occT,c] ~ dmulti(pr.T[j, ,c], rel.T[j,c])
        } #c
      } #j
      # Main diagonal
      for (j in 1:(n.occT-1)){
        for (c in 7:Y){
          q.T[j,c] <- 1-p.T[j,(c)]                # Probability of non-recapture
          pr.T[j,j,c] <- pi.T[j,(c)]*p.T[j,(c)]
          # Above main diagonal
          for (k in (j+1):(n.occT-1)){
            pr.T[j,k,c] <- prod(pi.T[j:k,(c)])*prod(q.T[j:(k-1),c])*p.T[k,(c)]
          } #k
          # Below main diagonal
          for (k in 1:(j-1)){
            pr.T[j,k,c] <- 0
          } #k
        } #c
      } #j
      # Last column: probability of non-recapture
      for (j in 1:(n.occT-1)){
        for (c in 7:Y){
          pr.T[j,n.occT,c] <- 1-sum(pr.T[j,1:(n.occT-1),c])
        } #c
      } #j

    # In-river group m-array count discrepancy - Freeman-Tukey
      # Compute fit statistics for observed data
      for (c in 7:Y){
        for (j in 1:(n.occ-1)){
          for (k in 1:n.occ){
            expmarr[j,k,c] <- rel[j,c]*pr[j,k,c]
            E.org[j,k,c] <- pow((pow(marr[j,k,c], 0.5)-pow(expmarr[j,k,c], 0.5)), 2)
          } #k
        } #j
        # Generate replicate data and compute fit stats from them
        for (j in 1:(n.occ-1)){
          marr.new[j,1:n.occ,c] ~ dmulti(pr[j, ,c], rel[j,c])
          for (k in 1:n.occ){
            E.new[j,k,c] <- pow((pow(marr.new[j,k,c], 0.5)-pow(expmarr[j,k,c], 0.5)), 2)
          } #k
        } #j
      } #c
      fit <- sum(E.org[,,7:Y])
      fit.new <- sum(E.new[,,7:Y])
    # Transported group m-array count discrepancy - Freeman-Tukey
      # Compute fit statistics for observed data
      for (c in 7:Y){
        for (j in 1:(n.occT-1)){
          for (k in 1:n.occT){
            expmarr.T[j,k,c] <- rel.T[j,c]*pr.T[j,k,c]
            E.org.T[j,k,c] <- pow((pow(marr.T[j,k,c], 0.5)-pow(expmarr.T[j,k,c], 0.5)), 2)
          } #k
        } #j
        # Generate replicate data and compute fit stats from them
        for (j in 1:(n.occT-1)){
          marr.new.T[j,1:n.occT,c] ~ dmulti(pr.T[j, ,c], rel.T[j,c])
          for (k in 1:n.occT){
            E.new.T[j,k,c] <- pow((pow(marr.new.T[j,k,c], 0.5)-pow(expmarr.T[j,k,c], 0.5)), 2)
          } #k
        } #j
      } #c
      fit.T <- sum(E.org.T[,,7:Y])
      fit.new.T <- sum(E.new.T[,,7:Y])
      
  ######################
  # Population process #
  ######################
  
  # Derived SAR statistics
    for (c in 1:(Y)){
      # conditional estimates of SAR (random year effects in estimate)
      sar.H[c] <- prod(pi.H[,c]) # LGR-LGA SAR indexed by migration year
      sar.T[c] <- prod(pi.T[,c])
      phi.sar[c] <- (1-del[c])*sar.H[c] + del[c]*sar.T[c]
      # unconditional estimates of SAR (random year effects into uncertainty)
      sar.H.pred[c] <- prod(pi.H.pred[,c]) # LGR-LGA SAR indexed by migration year
      sar.T.pred[c] <- prod(pi.T.pred[,c])
      phi.sar.pred[c] <- (1-del[c])*sar.H.pred[c] + del[c]*sar.T.pred[c]
    } #c
  
  #  Hindcast
    for (i in 1:s){ #for each basin
      for (t in 1:(a.max+1)){ #1:7
        logratio[i,t] ~ dnorm(logr[i], sigma.logr[i]) 
        log.S[i,t] <- log.index[t] + logratio[i,t]
        S[i,t] <- exp(log.S[i,t])
      }
    }
    for (i in 1:s){ 
    
      N0[i,2] <- S[i,1] * b[i,1]
    
      N0[i,3] <- S[i,2] * b[i,2]
      N1[i,3] <- N0[i,2]
    
      N0[i,4] <- S[i,3] * b[i,3]
      N1[i,4] <- N0[i,3]
      N2[i,4] <- N1[i,3] * ((phi.H[3]*phi.O[3]*(1-del[3]))+(phi.B*phi.T[3]*del[3]))
    
      N0[i,5] <- S[i,4] * b[i,4]
      N1[i,5] <- N0[i,4]
      N2[i,5] <- N1[i,4] * ((phi.H[4]*phi.O[4]*(1-del[4]))+(phi.B*phi.T[4]*del[4]))
      N3[i,5] <- N2[i,4] * phi.A
    
      N0[i,6] <- S[i,5] * b[i,5]
      N1[i,6] <- N0[i,5]
      N2[i,6] <- N1[i,5] * ((phi.H[5]*phi.O[5]*(1-del[5]))+(phi.B*phi.T[5]*del[5]))
      N3[i,6] <- N2[i,5] * phi.A
      N4[i,6] <- N3[i,5] * phi.A * (1-rho3[5])
      R4[i,6] <- N3[i,5] * phi.A * rho3[5] * phi.R
      
      N0[i,7] <- S[i,6] * b[i,6]
      N1[i,7] <- N0[i,6]
      N2[i,7] <- N1[i,6] * ((phi.H[6]*phi.O[6]*(1-del[6]))+(phi.B*phi.T[6]*del[6]))
      N3[i,7] <- N2[i,6] * phi.A
      N4[i,7] <- N3[i,6] * phi.A * (1-rho3[6])
      N5[i,7] <- N4[i,6] * phi.A * (1-rho4)
      R4[i,7] <- N3[i,6] * phi.A * rho3[6] * phi.R
      R5[i,7] <- N4[i,6] * phi.A * rho4 * phi.R
    }
  # Population dynamics
  # note, the transitions [.]->Bon->Redd->S occur within the same year
  for (i in 1:s){
    for (t in (a.max+1):(Y-1)){
      # results of state transition probabilities
      N0[i,t+1] <- S[i,t] * b[i,t]
      N1[i,t+1] <- N0[i,t] #b predicts N1, this is just passage of time
      N2[i,t+1] <- N1[i,t] * ((phi.H[t]*phi.O[t]*(1-del[t]))+(phi.B*phi.T[t]*del[t]))
      N3[i,t+1] <- N2[i,t] * phi.A
      N4[i,t+1] <- N3[i,t] * phi.A * (1-rho3[t])
      N5[i,t+1] <- N4[i,t] * phi.A * (1-rho4)
      R4[i,t+1] <- N3[i,t] * phi.A * rho3[t] * phi.R
      R5[i,t+1] <- N4[i,t] * phi.A * rho4 * phi.R
      R6[i,t+1] <- N5[i,t] * phi.A * phi.R
      S[i,t+1] <- R4[i,t+1] + R5[i,t+1] + R6[i,t+1] 
    }
  }
  # Lambda, discrete time population growth rate and fitness measure
  # Note: Lambda is indexed by BROODYEAR COHORT, not migration year cohort
    for (c in 1:(Y-2)){
      for (i in 1:s){
        R.S[i,c] <- b[i,c]*phi.sar[c+2]
        R[i,c] <- R.S[i,c]*S[i,c] 
        R.S.pred[i,c] <- b[i,c]*phi.sar.pred[c+2]
        R.pred[i,c] <- R.S.pred[i,c]*S[i,c]
      }
      Smfsr[c] <- sum(S[,c])
      Rmfsr[c] <- sum(R[,c])
      Rmfsr.pred[c] <- sum(R.pred[,c])
      R.Smfsr[c] <- Rmfsr[c]/Smfsr[c] # Conditional estimate
      R.Smfsr.pred[c] <- Rmfsr.pred[c]/Smfsr[c] # Unconditional estimate
    }
  
  # Expected proportion-at-age by year and MPG 'population'
  for (t in (a.max+2):Y){
    hat.p.age[1,t,1] <- sum(R4[inds1,t]) / sum(S[inds1,t]) # age-4
    hat.p.age[1,t,2] <- sum(R5[inds1,t]) / sum(S[inds1,t]) # age-5
    hat.p.age[1,t,3] <- sum(R6[inds1,t]) / sum(S[inds1,t]) # age-6
    
    hat.p.age[2,t,1] <- sum(R4[inds2,t]) / sum(S[inds2,t]) # age-4
    hat.p.age[2,t,2] <- sum(R5[inds2,t]) / sum(S[inds2,t]) # age-5
    hat.p.age[2,t,3] <- sum(R6[inds2,t]) / sum(S[inds2,t]) # age-6
    
    hat.p.age[3,t,1] <- sum(R4[inds3,t]) / sum(S[inds3,t]) # age-4
    hat.p.age[3,t,2] <- sum(R5[inds3,t]) / sum(S[inds3,t]) # age-5
    hat.p.age[3,t,3] <- sum(R6[inds3,t]) / sum(S[inds3,t]) # age-6
    
    hat.p.age[4,t,1] <- sum(R4[inds4,t]) / sum(S[inds4,t]) # age-4
    hat.p.age[4,t,2] <- sum(R5[inds4,t]) / sum(S[inds4,t]) # age-5
    hat.p.age[4,t,3] <- sum(R6[inds4,t]) / sum(S[inds4,t]) # age-6
    
    hat.p.age[5,t,1] <- sum(R4[inds5,t]) / sum(S[inds5,t]) # age-4
    hat.p.age[5,t,2] <- sum(R5[inds5,t]) / sum(S[inds5,t]) # age-5
    hat.p.age[5,t,3] <- sum(R6[inds5,t]) / sum(S[inds5,t]) # age-6
    
    hat.p.age[6,t,1] <- sum(R4[inds6,t]) / sum(S[inds6,t]) # age-4
    hat.p.age[6,t,2] <- sum(R5[inds6,t]) / sum(S[inds6,t]) # age-5
    hat.p.age[6,t,3] <- sum(R6[inds6,t]) / sum(S[inds6,t]) # age-6
    
    hat.p.age[7,t,1] <- sum(R4[inds7,t]) / sum(S[inds7,t]) # age-4
    hat.p.age[7,t,2] <- sum(R5[inds7,t]) / sum(S[inds7,t]) # age-5
    hat.p.age[7,t,3] <- sum(R6[inds7,t]) / sum(S[inds7,t]) # age-6
    
    hat.p.age[8,t,1] <- sum(R4[inds8,t]) / sum(S[inds8,t]) # age-4
    hat.p.age[8,t,2] <- sum(R5[inds8,t]) / sum(S[inds8,t]) # age-5
    hat.p.age[8,t,3] <- sum(R6[inds8,t]) / sum(S[inds8,t]) # age-6
  }
  # Age distribution of returning females
    for (i in age.i){
      for (t in age.t[i,(1:age.nt[i])]){
        y.age[i,t, ] ~ dmulti(hat.p.age[i,t, ], age.draws[i,t])
      }
    }

  #########################
  # Redd count likelihood #
  #########################
  
  # multiplicative extra-Poisson variance 
  sigma.pois ~ dt(0, pow(2.5, -2), 1)T(0,)
  tau.pois <- pow(sigma.pois, -2)
  # likelihood and discrepancy
  for (i in 1:s){
    for (t in 1:Y){
    
      # State-space with over-dispersed poisson observation process
      y[i,t] ~ dpois( eta[i,t] )
      eta[i,t] <- S[i,t] * exp.eps[i,t] # log(eta) = log(S) + eps
      exp.eps[i,t] <- exp(eps[i,t])
      eps[i,t] ~ dnorm(0, tau.pois) # eps ~ N(0, sigma.pois)
      Tobs0[i,t] <- (sqrt(y[i,t]) - sqrt(eta[i,t]))^2  
      
      # Simulated data for state-space version
      ySim[i,t] ~ dpois(S[i,t] * exp.simeps[i,t] )
      exp.simeps[i,t] <- exp(simeps[i,t])
      simeps[i,t] ~ dnorm(0, tau.pois)
      Tsim0[i,t] <- (sqrt(ySim[i,t]) - sqrt(S[i,t] * exp.simeps[i,t]))^2  
      
    }
    Tobsi[i] <- sum(Tobs0[i,(a.max+2):Yy])
    Tsimi[i] <- sum(Tsim0[i,(a.max+2):Yy])
  }
  # GOF assessment
  Tobs <- sum(Tobsi[])
  Tsim <- sum(Tsimi[])
  
}  #model

