################# data
s <- taillle

tr <- BDD_a %>% filter(Year !="2022") %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
mooved <- tr[duplicated(tr$Tag_id),]$Tag_id
tr <- tr %>% filter ( !duplicated(tr$Tag_id))
tr <- tr %>% mutate (Treatment = ifelse(Tag_id %in% mooved, NA, Treatment))

la <- BDD_a %>% group_by(Tag_id, Lake) %>% summarize() %>% ungroup()
mooved <- la[duplicated(la$Tag_id),]$Tag_id
la <- la %>% filter ( !duplicated(la$Tag_id))
la <- la %>% mutate (Lake = ifelse(Tag_id %in% mooved, NA, Lake))

t <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from =Lake)
t <- merge(t,tr)
t <- merge(t,la)

t <- data.frame(t$Tag_id, apply(t[2:8],2, function(x) ifelse(is.na(x),0,1)), t$Treatment, t$Lake)
t <- na.omit(t)

set.seed(11)
t <- t[sample(1:nrow(t), 200),]

s <- s %>% filter (Tag_id %in% t$t.Tag_id)
t <- t %>% arrange(t.Tag_id)
s <- s %>% arrange(Tag_id)

CH <- as.matrix(t[2:7])

get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

z.inits <- function(ch){
  state <- ch
  state[state==0] <- 1
  get.first <- function(x) min(which(x!=0))
  f <- apply(ch, 1, get.first)
  for (i in 1:nrow(ch)){
    state[i,1:(f[i])] <- NA
  }
  return(state)
}

z.inits2 <- function(ch){
  state <- ch
  state[state==0] <- 1
  get.first <- function(x) min(which(x!=0))
  f <- apply(ch, 1, get.first)
  for (i in 1:nrow(ch)){
    #print(f[i]-1)
    if (f[i]>10){
      state[i,1:6] <- 0
    }
    if (f[i]>1 & f[i]<10){
      state[i,1:(f[i]-1)] <- 0
    }
  }
  state[,1] <- NA
  return(state)
}

jags.data <- list(y = CH,
                  size = as.matrix(s[2:7]),
                  Treatment = t$t.Treatment,
                  f = apply(CH, 1, get.first), 
                  nind = nrow(CH), 
                  noccas = ncol(CH),
                  ni = 100,
                  #zi = z.inits(CH),
                  zi2 = z.inits2(CH))

############################ Model 1 phi(~1)p(~1)n(~t)

inits <- function(){
  list(phi = runif(1,0,1), p = runif(1,0.5,1), n = runif(5,0.5,1), z = zi2)}


parameters <- c("phi", "p", "n", "N")

sa_1 <- function(){
  # Likelihood
  for (i in 1:nind){
    # State at first capture is alive
    z[i,1] <- ifelse(f[i]==1,1,0)
    for (t in 2:noccas){
      # determine the state alive/dead
      z[i,t] ~ dbern(phi* z[i,t-1] + n[t-1]*(1-max(z[i,(1:(t-1))]))) #(1-max(z[i,(1:(t-1))]))) 
      # determine the capture status
      y[i,t] ~ dbern(p * z[i,t])
    }
  }
  # Priors
  phi ~ dunif(0,1)
  p ~ dunif(0,1)
  for (t in 1:(noccas-1)){
    n[t] ~ dunif(0,1)
  }
  # Derived variable
  for (t in 1:noccas){
    N[t] <- sum(z[,t])
  }
}

SA_1 <- jags.parallel(data = jags.data,
                       inits = inits,
                       parameters.to.save = parameters,
                       model.file = sa_1,
                       n.chains = 4,
                       n.iter = ni)

save(SA_1, file = "object/SA_1.RData")

