cjs_taille1 <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      p[i,t_bis] <- 0
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      s[i,t] <- size[i,t]
      logit(p[i,t]) <- a_p*10^-2*s[i,t] + b_p #a2_p*10^-4*(s[i,t]^2) + 
      logit(phi[i,t]) <- a_phi*10^-2*s[i,t] + b_phi # + a2_phi*10^-4*size[i,t]^2
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  a_p ~ dunif(-1000,1000)
  a_phi ~ dunif(-1000,1000)
  #a2_p ~ dunif(-1000,1000)
  #a2_phi ~ dunif(-10,0)
  b_phi ~ dunif(-1000,1000)
  b_p ~ dunif(-1000,1000)
 
  #Priors
  sd <- 12 #~ dunif(0,100)
  K <- 0.352 #~ dunif(-1000,1000)
  Linf <- 187 #~ dunif(0,500)
  #Likelihood
  for (i in 1:nind){
    for (t in (f[i]+1):noccas){
      size[i,t] ~ dnorm(size[i,t-1] + K*(Linf-size[i,t-1]), 1/(sd^2))
    } # t
  } # i
}

inits <- function(){
  list(a_p = runif(1,5,13), a_phi = runif(1,1,20),
       #a2_p = runif(1,-7,0), #a2_phi = runif(1,-7,0),
       b_p = runif(1,-10,5), b_phi = runif(1,-20,10),
       z = zi)}#,
       #sd = runif(1,0,100),
       #g =runif(1,-100,100), 
       #Linf = runif(1,0,500)
parameters = c("a2_p","a_p","a_phi","b_p","b_phi","K","sd","Linf")

CJS_taille1 <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = cjs_taille1,
                             n.chains = 3,
                             n.iter = 1000)

taille <- function(){
  #Priors
  sd ~ dunif(0,100)
  for (tr in 1:4){
    K[tr] ~ dunif(-1000,1000)
    Linf[tr] ~ dunif(0,500)
  }
  
  #Likelihood
  for (i in 1:nind){
    for (t in (f[i]+1):noccas){
      size[i,t] ~ dnorm(size[i,t-1] + K[Treatment[i]]*(Linf[Treatment[i]]-size[i,t-1]), 1/(sd^2))
    } # t
  } # i
}

inits = function(){
  list(sd = runif(1,0,100), K = runif(4,-1000,1000), Linf = runif(4, 0,500))
}

parameters = c("K","Linf","sd")

Taille <- jags.parallel(data = jags.data,
               inits = inits,
               parameters.to.save = parameters,
               model.file = taille,
               n.iter = 1000,
               n.chains = 2)


#####################################
#####################################
s <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from = Size)

tr <- BDD_a %>% filter(Year != "2022") %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
mooved <- tr[duplicated(tr$Tag_id),]$Tag_id
tr <- tr %>% filter ( !duplicated(tr$Tag_id))
tr <- tr %>% mutate (Treatment = ifelse(Tag_id %in% mooved, NA, Treatment))

la <- BDD_a %>% group_by(Tag_id, Lake) %>% summarize() %>% ungroup()
mooved <- la[duplicated(la$Tag_id),]$Tag_id
la <- la %>% filter ( !duplicated(la$Tag_id))
la <- la %>% mutate (Lake = ifelse(Tag_id %in% mooved, NA, Lake))

s <- merge(s,tr)
s <- merge(s,la)

s <- s %>% filter(!is.na(Lake))

#set.seed(30)
#s <- s[sample(1:nrow(s), 200),]

size_break <- c(0,120,180,300)
s_group <- data.frame(apply(s[2:7], 2, 
                            function(x) discretize(x,method = "fixed",
                                                   breaks = size_break, 
                                                   labels= c(1:(length(size_break)-1)))))
s_group [] <- apply(s_group, 2, as.numeric)

fs <- apply(s_group, 1, function(x) first(na.omit(x)))

s_group <- data.frame(apply(s_group,2, function(x) ifelse(is.na(x),4,x)))

get.first <- function(x) min(which(x!=4))
f <- apply(s_group, 1, get.first)


#z.inits <- function(ch){
#  state <- ch
#  state[] <- 1
#  get.first <- function(x) min(which(x!=0))
#  f <- apply(ch, 1, get.first)
#  for (i in 1:nrow(ch)){
#    state[i,1:(f[i])] <- NA
#  }
#  return(state)
#}

zinit <- s_group
for (i in 1:nrow(s_group)) {
  for (j in 1:ncol(s_group)) {
    if (j > f[i] & s_group[i,j]==4) {zinit[i,j] <- zinit[i,j-1]}
  }
}
for (i in 1:nrow(s_group)) {
  for (j in 1:ncol(s_group)) {    
    if (j <= f[i]) {zinit[i,j] <- NA}
  }
}
for (i in 1:nrow(s_group)) {
  if (f[i]<4){
  for (j in (f[i]+2):ncol(zinit)) {
    if (zinit[i,j]<zinit[i,j-1]) {zinit[i,j] <- zinit[i,j-1]}
  }
}}
zinit <- as.matrix(zinit)

jags.data <- list(y = s_group,
                  zi = zinit,
                  f = f,
                  fs = fs,
                  nind = nrow(s_group),
                  noccas = ncol(s_group),
                  ni = 1000)

multievent <- function(){
  
  # -------------------------------------------------
  # Parameters:
  # phiB: survival probability state B
  # phiNB: survival probability state NB
  # psiBNB: transition probability from B to NB
  # psiNBB: transition probability from NB to B
  # pB: recapture probability B
  # pNB: recapture probability NB
  # piB prob. of being in initial state breeder
  # -------------------------------------------------
  # States (z):
  # 1 alive B
  # 2 alive NB
  # 3 dead
  # Observations (y):  
  # 1 = non-detected
  # 2 = seen and ascertained as breeder
  # 3 = seen and ascertained as non-breeder
  # 4 = not ascertained
  # -------------------------------------------------
  
  # priors
  phi1 ~ dunif(0,1)
  phi2 ~ dunif(0,1)
  phi3 ~ dunif(0,1)
  
  psi12 ~ dunif(0,1)
  psi23 ~ dunif(0,1)
  
  psi13 ~ dunif(0,1)
  psi32 ~ dunif(0,1)
  psi21 ~ dunif(0,1)
  psi31 ~ dunif(0,1)
  
  p1 ~ dunif(0,1)
  p2 ~ dunif(0,1)
  p3 ~ dunif(0,1)
  
  erreur ~ dunif(0,1)
  
  # probabilities of state z(t+1) given z(t)
  gamma[1,1] <- phi1 * (1-psi12-psi13)
  gamma[1,2] <- phi1 * psi12
  gamma[1,3] <- phi1 * psi13
  gamma[1,4] <- (1-phi1)
  gamma[2,1] <- phi2 * psi21
  gamma[2,2] <- phi2 * (1-psi23-psi21)
  gamma[2,3] <- phi2 * psi23
  gamma[2,4] <- (1-phi2)
  gamma[3,1] <- phi3 * psi31
  gamma[3,2] <- phi3 * psi32
  gamma[3,3] <- phi3 * (1-psi31-psi32)
  gamma[3,4] <- (1-phi3)
  gamma[4,1] <- 0
  gamma[4,2] <- 0
  gamma[4,3] <- 0
  gamma[4,4] <- 1
  
  # probabilities of y(t) given z(t)
  omega[1,1] <- p1
  omega[1,2] <- 0
  omega[1,3] <- 0
  omega[1,4] <- 1-p1
  omega[2,1] <- p2 * erreur
  omega[2,2] <- p2 * (1-erreur)
  omega[2,3] <- 0
  omega[2,4] <- 1-p2
  omega[3,1] <- p3 * erreur
  omega[3,2] <- p3 * erreur
  omega[3,3] <- p3 * (1- 2*erreur)
  omega[3,4] <- 1-p3
  omega[4,1] <- 0
  omega[4,2] <- 0
  omega[4,3] <- 0
  omega[4,4] <- 1
  
  # likelihood 
  for (i in 1:nind){
    # latent state at first capture
    z[i,f[i]] <- fs[i]
    for (t in (f[i]+1):noccas){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:4])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:4])
    }
  }
}

inits = function(){
  list(phi1 = runif(1,0,1),phi2 = runif(1,0,1),phi3 = runif(1,0,1),
       p1 = runif(1,0,1),p2 = runif(1,0,1),p3 = runif(1,0,1),
       psi12 = runif(1,0,0.5),psi23 = runif(1,0,0.5),
       psi32 = runif(1,0,0.5),psi21 = runif(1,0,0.5), psi31 = dunif(1,0,0.5), psi13 = runif(1,0,0.5),
       erreur = runif(1,0,0.4),
       z = zi)}

parameters = c("phi1","phi2","phi3",
               "p1","p2","p3",
               "psi12","psi21",
               "psi32","psi21","psi31","psi13",
               "erreur")

Model_multi <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = multievent,
                             n.chains = 2,
                             n.iter = ni)



# probabilities of state z(t+1) given z(t)
gamma[1,1] <- phi1 * (1-psi12-psi13)
gamma[1,2] <- phi1 * psi12
gamma[1,3] <- phi1 * psi13
gamma[1,4] <- (1-phi1)
gamma[2,1] <- phi2 * psi21
gamma[2,2] <- phi2 * (1-psi23-psi21)
gamma[2,3] <- phi2 * psi23
gamma[2,4] <- (1-phi2)
gamma[3,1] <- phi3 * psi31
gamma[3,2] <- phi3 * psi32
gamma[3,3] <- phi3 * (1-psi31-psi32)
gamma[3,4] <- (1-phi3)
gamma[4,1] <- 0
gamma[4,2] <- 0
gamma[4,3] <- 0
gamma[4,4] <- 1

# probabilities of y(t) given z(t)
omega[1,1] <- p1
omega[1,2] <- 0
omega[1,3] <- 0
omega[1,4] <- 1-p1
omega[2,1] <- 0
omega[2,2] <- p2
omega[2,3] <- 0
omega[2,4] <- 1-p2
omega[3,1] <- p3 * erreur
omega[3,2] <- 0
omega[3,3] <- p3 * (1-erreur)
omega[3,4] <- 1-p3
omega[4,1] <- 0
omega[4,2] <- 0
omega[4,3] <- 0
omega[4,4] <- 1
