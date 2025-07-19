resample <- function(samples,weight)
{
  
  pick <- replicate(length(weight),which(rmultinom(1,prob=weight,size=1)==1))
  
  return(samples[pick,])
  
}

calc_beta<-function(R0,TE,TI,TT,eta,chi)
{
  u = 16.0/1000.0
  A0 = 1/(u*(R0-1))
  
  if(TT == 0.0)
  {
    return(((1/TE + u)*((1/TI) + u) * R0) * (TE));
  }
  else if(eta == 1.0 & chi == 1)
  {
    # Mean Age of Primary Infection
    Ap = A0*364.0;
    #Rate of Loss of Immunity
    an = 1.0/(TT*364.0);
    #Ratio of Primary to Secondary Infectious Periods
    np = 1*(TI/TI);
    np = 1.0;
    
    beta = ((1/TE + u)*((1/TI) + u) * R0) * (TE);
    
    Y = 1 - (an/(an + u))*(1 - Ap*u);
    
    return(beta = beta * Y)
  }
  else
  {
    #Mean Age of Primary Infection
    Ap = A0*364.0;
    #Rate of Loss of Immunity
    an = 1.0/(364.0*TT);
    #Ratio of Primary to Secondary Infectious Periods
    np = neta*(TI/TI);
    
    l1star = (1.0/Ap) - u;
    
    epsilon = 0.0;
    
    a0 = (chi - eta/chi)*(((1 - epsilon)*np*chi - 1)*(an/(an+u))+1);
    
    a1 = (u / l1star) + (1 - epsilon)*np*chi*(an/(an+u)) + ((an/(an+u)) - 1)*(chi - 2*(eta/chi));
    
    a2 = - (u / l1star) + ((an/(an+u)) - 1)*(eta/chi);
    
    Y = 1;
    
    roots <- polyroot(c(a0,a1,a2))
    
    if(length(roots) == 0){return(NA)
    }else if(numroots == 1){
      if(roots[1] < 0.0){return(NA)
      }else
      {
        Y = roots[1];
      }
    }else if(length(roots) == 2){
      if(roots[1] > 0.0)
      {
        Y = roots[1];
      }else if(roots[2] > 0.0)
      {
        Y = roots[2];
      }else
      {
        return(NA)
      }
    }
    
    beta = ((1/TE + u)*((1/TI) + u) * R0) * (TE);
    return(beta = beta * Y)
    
  }
}

mk_terms <- function(beta, alpha)
{
  s = 273.0/364.0
  bh = (beta*(1 + 2*alpha*(1-s))) 
  bl = (beta*(1 - 2*alpha*s)) 
  
  terms = numeric(364)
  
  terms[c(252:299)] = bh
  terms[c(300:307)] = bl
  terms[c(308:355)] = bh
  terms[c(356:364)] = bl
  terms[c(1:6)] = bl
  terms[c(7:99)] = bh
  terms[c(100:115)] = bl
  terms[c(116:199)] = bh
  terms[c(200:251)] = bl
  
  return(terms)
  
}

mk_vacc_coverage <- function(coverage,pattern,alpha)
{
  # Proportion of high coverage areas
  s = sum(pattern)/length(pattern)
  vh = (coverage*(1 + 2*alpha*(1-s))) 
  vl = (coverage*(1 - 2*alpha*s)) 

  vac_coverage <- numeric(length(pattern))
  vac_coverage[which(pattern)] = vh
  vac_coverage[which(!pattern)] = vl

  return(vac_coverage)
  
}


detSEIT <-function(x,plist)
{
  
  u = plist[1]
  beta = plist[2];
  gE = plist[3];
  gI = plist[4];
  gT = plist[5];
  gEA = plist[6];
  gIA = plist[7];
  neta = plist[8];
  eta = plist[9];
  chi = plist[10];
  
  s = x[1]
  e = x[2]
  i = x[3]
  t = x[4]
  sa = x[5]
  ea = x[6]
  ia = x[7]
  
  f0 = u - beta*s*(i+neta*chi*ia) - u * s;
  f1 = beta*s*(i+neta*chi*ia) - u * e - (gE) * e;
  f2 = (gE) * e - (gI) * i - u * i;
  f3 = (gI) * i + (gIA) * ia - (gT) * t - u * t;
  f4 = (gT) * t - beta*sa*(chi*i+neta*eta*ia) - u * sa;
  f5 =  beta*sa*(chi*i+neta*eta*ia) - (gEA) * ea - u * ea;
  f6 = (gEA) * ea - (gIA) * ia - u*ia;
  
  return(c(f0=f0,f1=f1,f2=f2,f3=f3,f4=f4,f5=f5,f6=f6))
  
}

SEITInit <- function(beta,v,TE,TI,TT,R0)
{
  require(rootSolve)
  gE = 1/TE
  gI = 1/TI
  p=c(v,beta,1/TE,1/TI,1/TT,1/TE,1/TI,1,0.75,0.5)
  S0 =  ((1.0/ R0));  
  SA0 = 2 * S0 
  E0 = (v*(v + gI)/(beta*gE))*(R0-1);
  I0 = ((v/beta)*(R0-1));
  EA0 = E0;
  IA0 = I0;
  T0 = 1 - S0 - SA0 - 2*E0 - 2*I0;
  init <- multiroot(detSEIT, c(S0, E0, I0,T0,SA0,EA0,IA0),plist=p,positive=TRUE)
  return(init$root)
}


fadeout_statistics <- function(x)
{
  y<-rle(x);
  
  max_fade <- max(y$lengths[y$values==0])
  mean_fadeL <- mean(y$lengths[y$values==0])
  no_fade  <- sum(y$values==0)
  return(c(max_fade=max_fade,no_fade=no_fade,mean_fadeL=mean_fadeL))
}


# Function runs one simulation of the fitted E&W measles SEITmeta model
# For a single vaccine coverage (proportion between 0-1)
# Returns tibble of statistics related to prevalence of affected districts, etc

vacc_scenario_init <- function(vaccine_coverage=0, import_rate=1.0/(1000*364.0))
{
# Note vacc_coverage was missing as an argument for this function, what value 
#  did you use to generate your initial condition 
# 0.5?


aito <- apply(resample(samples,weight),2,median)

theta = aito[1]
alpha = aito[2]
R0 = aito[3]
sigma = aito[4]
u = 16.0
Burnin = 40
beta = calc_beta(R0,TE,TI,TT,1,1)

terms <- mk_terms(beta,aito[2])
mu <- apply(z$weekly_births,2,median)/(7*z$N)

init<-SEITInit(beta,median(mu),TE,TI,TT*364,R0);

S = round(init[1]*z$N);
E = round(init[2]*z$N);
I = round(init[3]*z$N);
T = round(init[4]*z$N);
SA = round(init[5]*z$N);
EA = round(init[6]*z$N);
IA = round(init[7]*z$N);

Opportunities <- (sapply(1:length(z$N),function(i){max(z$s[i,which(z$r[i,] < aito[5])])}))
Opportunities <- ((Opportunities)/mean(Opportunities))^aito[6]
Opportunities[which(Opportunities==0)] = min(Opportunities[which(Opportunities!=0)])

init = tibble(S=S,E=E,I=I,T=T,SA=SA,EA=EA,IA=IA)
nji = ((1-theta*Opportunities)*diag(z$N) + theta*Opportunities * z$W * z$N)
neff = rowSums(nji)
outputC<-SEITgmeta_ext(z$N,    # Vector of population size of each patch in model (length 982)
                   init$S, # Vector of initial susceptibles in each patch (length 982)
                   init$E, # Vector of initial exposed in each patch (length 982)
                   init$I, # Vector of initial infected in each patch (length 982)
                   init$T, # Vector of initial recovered in each patch (length 982)
                   init$SA, # Vector of initial Susceptible-Again (length 982)
                   init$EA, # Vector of initial exposed again (length 982)
                   init$IA, # Vector of initial infected again (length 982)
                   terms, # Vector of 364 days (transmission rate on every day of the year calculated by mk_terms) 
                   rep(import_rate,length(S)), # lambda_ext Rate of infectious imports used during the burn-in period
                   1.0, # Chi (Not in model set to 1)
                   1.0, # neta (Not in model set to 1)
                   1.0, # eta (Not in model set to 1)
                   as.matrix(z$weekly_births), # Matrix of crude birth rates: 521 weeks for each of 982 districts
                   z$W, # Commuter matrix (982 x 982) 
                   mu*(1-vaccine_coverage), # Vector of per capita birth rates for each location (length 982)
                   TE, # Exposed Period
                   kE, # Shape parameter for exposed period
                   TI, # Infectious Period
                   kI,# Shape parameter for infectious period
                   TE, # Exposed (Again) period
                   kE, # Shape parameter for Exposed (Again) period
                   TI, # Infectious (Again) period
                   kI, # Shape parameter for Infectious (Again) period
                   1/(TT*364.0), #Rate of loss of immunity
                   theta*Opportunities, # Movement probability (per capita)
                   1, # Time step dt (days) DO NOT CHANGE
                   sigma, # Environmental stochasticity term
                   10*52, # Time to simulate for (number of weeks)
                   41*52, # Burn in period (number of weeks)
                   TRUE, # If FALSE use weekly births (as above), if TRUE use per capita birth rate mu)
                   0) # FOI_Function chooses between different commuter approximations. DO NOT CHANGE.

# Peel off the final states and use as initial conditions
colnames(outputC) <- c('time','patch','cases','casesA','S','E','I','T','SA','EA','IA')
outputC <- as_tibble(outputC)
init = outputC %>% filter(time==max(time)) %>% select(S,E,I,T,SA,EA,IA)

return(init)
}


vacc_scenario_ext <- function(vacc_coverage, init, import_rate)
{
  aito <- apply(resample(samples,weight),2,median)
  
  theta = aito[1]
  alpha = aito[2]
  R0 = aito[3]
  sigma = aito[4]
  u = 16.0
  Burnin = 40
  beta = calc_beta(R0,TE,TI,TT,1,1)
  
  terms <- mk_terms(beta,aito[2])
  mu <- apply(z$weekly_births,2,median)/(7*z$N)
  
  Opportunities <- (sapply(1:length(z$N),function(i){max(z$s[i,which(z$r[i,] < aito[5])])}))
  Opportunities <- ((Opportunities)/mean(Opportunities))^aito[6]
  Opportunities[which(Opportunities==0)] = min(Opportunities[which(Opportunities!=0)])
  
  nji = ((1-theta*Opportunities)*diag(z$N) + theta*Opportunities * z$W * z$N)
  neff = rowSums(nji)
  outputC<-outputC<-SEITgmeta_ext(z$N,    # Vector of population size of each patch in model (length 982)
                                  init$S, # Vector of initial susceptibles in each patch (length 982)
                                  init$E, # Vector of initial exposed in each patch (length 982)
                                  init$I, # Vector of initial infected in each patch (length 982)
                                  init$T, # Vector of initial recovered in each patch (length 982)
                                  init$SA, # Vector of initial Susceptible-Again (length 982)
                                  init$EA, # Vector of initial exposed again (length 982)
                                  init$IA, # Vector of initial infected again (length 982)
                                  terms, # Vector of 364 days (transmission rate on every day of the year calculated by mk_terms) 
                                  import_rate, # lambda_ext Rate of infectious imports used during the burn-in period
                                  1.0, # Chi (Not in model set to 1)
                                  1.0, # neta (Not in model set to 1)
                                  1.0, # eta (Not in model set to 1)
                                  as.matrix(z$weekly_births), # Matrix of crude birth rates: 521 weeks for each of 982 districts
                                  z$W, # Commuter matrix (982 x 982) 
                                  mu*(1-vacc_coverage), # Vector of per capita birth rates for each location (length 982)
                                  TE, # Exposed Period
                                  kE, # Shape parameter for exposed period
                                  TI, # Infectious Period
                                  kI,# Shape parameter for infectious period
                                  TE, # Exposed (Again) period
                                  kE, # Shape parameter for Exposed (Again) period
                                  TI, # Infectious (Again) period
                                  kI, # Shape parameter for Infectious (Again) period
                                  1/(TT*364.0), #Rate of loss of immunity
                                  theta*Opportunities, # Movement probability (per capita)
                                  1, # Time step dt (days) DO NOT CHANGE
                                  sigma, # Environmental stochasticity term
                                  10*52, # Time to simulate for (number of weeks)
                                  0, # Burn in period (number of weeks)
                                  TRUE, # If FALSE use weekly births (as above), if TRUE use per capita birth rate mu)
                                  0) # FOI_Function chooses between different commuter approximations. DO NOT CHANGE.

  colnames(outputC) <- c('time','patch','cases','casesA','S','E','I','T','SA','EA','IA')
  outputC <- as_tibble(outputC)
  cases <- as.matrix(outputC %>% 
                       pivot_wider(id_cols=time,names_from=patch,values_from=cases) %>% 
                       dplyr::select(-time))
  
  
  max_extent <- numeric(dim(cases)[1])
  
  for(t in 1:dim(cases)[1])
  {
    max_extent[t]<-max(dm[which(cases[t,]>0),which(cases[t,]>0)])
  }
  
  # Calculate number of districts with measles activity for each time point
 
  x <- rowSums(cases>0)
  # Runlength encoding for time series of switches between <5 districts with
  # activity and >= 5 districts with activity
  # Use to identify first outbreak with large spatial extent after drop in 
  # vaccination rate
  y <- rle(x>=5)
  
  # If there are more than 2 outbreaks with >= 5 districts with activity
  if(sum(x >= 5) > 2 & length(y$values) > 1) 
    {
    return(tibble(
      prevalence_infected_districts = max(rowSums(cases > 0)),
      time_to_first_outbreak = rle(x >= 5)$lengths[1],
      peak_first_outbreak = max(x[1:(rle(x >= 5)$lengths[1] + rle(x >= 5)$lengths[2])]),
      maximum_extent = max(max_extent)
    ))
  } else {
    return(tibble(
      prevalence_infected_districts = 0,
      time_to_first_outbreak = Inf,
      peak_first_outbreak = 0,
      maximum_extent = max(max_extent)
    ))
  }
  

}

vacc_scenario_extC <- function(vacc_coverage, init, import_rate)
{
  aito <- apply(resample(samples,weight),2,median)
  
  theta = aito[1]
  alpha = aito[2]
  R0 = aito[3]
  sigma = aito[4]
  u = 16.0
  Burnin = 40
  beta = calc_beta(R0,TE,TI,TT,1,1)
  
  terms <- mk_terms(beta,aito[2])
  mu <- apply(z$weekly_births,2,median)/(7*z$N)
  
  Opportunities <- (sapply(1:length(z$N),function(i){max(z$s[i,which(z$r[i,] < aito[5])])}))
  Opportunities <- ((Opportunities)/mean(Opportunities))^aito[6]
  Opportunities[which(Opportunities==0)] = min(Opportunities[which(Opportunities!=0)])
  
  nji = ((1-theta*Opportunities)*diag(z$N) + theta*Opportunities * z$W * z$N)
  neff = rowSums(nji)
  outputC<-outputC<-SEITgmeta_ext(z$N,    # Vector of population size of each patch in model (length 982)
                                  init$S, # Vector of initial susceptibles in each patch (length 982)
                                  init$E, # Vector of initial exposed in each patch (length 982)
                                  init$I, # Vector of initial infected in each patch (length 982)
                                  init$T, # Vector of initial recovered in each patch (length 982)
                                  init$SA, # Vector of initial Susceptible-Again (length 982)
                                  init$EA, # Vector of initial exposed again (length 982)
                                  init$IA, # Vector of initial infected again (length 982)
                                  terms, # Vector of 364 days (transmission rate on every day of the year calculated by mk_terms) 
                                  import_rate, # lambda_ext Rate of infectious imports used during the burn-in period
                                  1.0, # Chi (Not in model set to 1)
                                  1.0, # neta (Not in model set to 1)
                                  1.0, # eta (Not in model set to 1)
                                  as.matrix(z$weekly_births), # Matrix of crude birth rates: 521 weeks for each of 982 districts
                                  z$W, # Commuter matrix (982 x 982) 
                                  mu*(1-vacc_coverage), # Vector of per capita birth rates for each location (length 982)
                                  TE, # Exposed Period
                                  kE, # Shape parameter for exposed period
                                  TI, # Infectious Period
                                  kI,# Shape parameter for infectious period
                                  TE, # Exposed (Again) period
                                  kE, # Shape parameter for Exposed (Again) period
                                  TI, # Infectious (Again) period
                                  kI, # Shape parameter for Infectious (Again) period
                                  1/(TT*364.0), #Rate of loss of immunity
                                  theta*Opportunities, # Movement probability (per capita)
                                  1, # Time step dt (days) DO NOT CHANGE
                                  sigma, # Environmental stochasticity term
                                  10*52, # Time to simulate for (number of weeks)
                                  0, # Burn in period (number of weeks)
                                  TRUE, # If FALSE use weekly births (as above), if TRUE use per capita birth rate mu)
                                  0) # FOI_Function chooses between different commuter approximations. DO NOT CHANGE.
  
  colnames(outputC) <- c('time','patch','cases','casesA','S','E','I','T','SA','EA','IA')
  outputC <- as_tibble(outputC)
  cases <- as.matrix(outputC %>% 
                       pivot_wider(id_cols=time,names_from=patch,values_from=cases) %>% 
                       dplyr::select(-time))
  
  
  return(cases)

}


vacc_scenario_extT <- function(vacc_coverage, init, import_rate)
{
  aito <- apply(resample(samples,weight),2,median)
  
  theta = aito[1]
  alpha = aito[2]
  R0 = aito[3]
  sigma = aito[4]
  u = 16.0
  Burnin = 40
  beta = calc_beta(R0,TE,TI,TT,1,1)
  
  terms <- mk_terms(beta,aito[2])
  mu <- apply(z$weekly_births,2,median)/(7*z$N)
  
  Opportunities <- (sapply(1:length(z$N),function(i){max(z$s[i,which(z$r[i,] < aito[5])])}))
  Opportunities <- ((Opportunities)/mean(Opportunities))^aito[6]
  Opportunities[which(Opportunities==0)] = min(Opportunities[which(Opportunities!=0)])
  
  nji = ((1-theta*Opportunities)*diag(z$N) + theta*Opportunities * z$W * z$N)
  neff = rowSums(nji)
  outputC<-outputC<-SEITgmeta_ext(z$N,    # Vector of population size of each patch in model (length 982)
                                  init$S, # Vector of initial susceptibles in each patch (length 982)
                                  init$E, # Vector of initial exposed in each patch (length 982)
                                  init$I, # Vector of initial infected in each patch (length 982)
                                  init$T, # Vector of initial recovered in each patch (length 982)
                                  init$SA, # Vector of initial Susceptible-Again (length 982)
                                  init$EA, # Vector of initial exposed again (length 982)
                                  init$IA, # Vector of initial infected again (length 982)
                                  terms, # Vector of 364 days (transmission rate on every day of the year calculated by mk_terms) 
                                  import_rate, # lambda_ext Rate of infectious imports used during the burn-in period
                                  1.0, # Chi (Not in model set to 1)
                                  1.0, # neta (Not in model set to 1)
                                  1.0, # eta (Not in model set to 1)
                                  as.matrix(z$weekly_births), # Matrix of crude birth rates: 521 weeks for each of 982 districts
                                  z$W, # Commuter matrix (982 x 982) 
                                  mu*(1-vacc_coverage), # Vector of per capita birth rates for each location (length 982)
                                  TE, # Exposed Period
                                  kE, # Shape parameter for exposed period
                                  TI, # Infectious Period
                                  kI,# Shape parameter for infectious period
                                  TE, # Exposed (Again) period
                                  kE, # Shape parameter for Exposed (Again) period
                                  TI, # Infectious (Again) period
                                  kI, # Shape parameter for Infectious (Again) period
                                  1/(TT*364.0), #Rate of loss of immunity
                                  theta*Opportunities, # Movement probability (per capita)
                                  1, # Time step dt (days) DO NOT CHANGE
                                  sigma, # Environmental stochasticity term
                                  10*52, # Time to simulate for (number of weeks)
                                  0, # Burn in period (number of weeks)
                                  TRUE, # If FALSE use weekly births (as above), if TRUE use per capita birth rate mu)
                                  0) # FOI_Function chooses between different commuter approximations. DO NOT CHANGE.
  
  colnames(outputC) <- c('time','patch','cases','casesA','S','E','I','T','SA','EA','IA')
  outputC <- as_tibble(outputC)
  cases <- as.matrix(outputC %>% 
                       pivot_wider(id_cols=time,names_from=patch,values_from=cases) %>% 
                       dplyr::select(-time))
  
  
  return(outputC)
  
}
