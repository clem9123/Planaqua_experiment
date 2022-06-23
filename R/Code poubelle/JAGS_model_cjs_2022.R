################# data
s <- s2

tr <- BDD_a %>% filter(Tag_year != "2022") %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
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

#set.seed(10)
#t <- t[sample(1:nrow(t), 100),]

s <- s %>% filter (Tag_id %in% t$t.Tag_id)
t <- t %>% arrange(t.Tag_id)
s <- s %>% arrange(Tag_id)

CH <- as.matrix(t[2:8])

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

jags.data <- list(y = CH,
                  size = as.matrix(s[2:8]),
                  Treatment = t$t.Treatment,
                  f = apply(CH, 1, get.first), 
                  nind = nrow(CH), 
                  noccas = ncol(CH),
                  ni = 1000,
                  zi = z.inits(CH))

################### model p(~1)phi(~1)

inits <- function(){
  list(phi = runif(1,0,1), p = runif(1,0,1), phi_v = runif(1,0,1), z = zi)}


parameters <- c("phi", "p", "phi_v")

cjs_0_2022 <- function(){
  # Likelihood
  for (i in 1:nind){
    # State at first capture is alive
    z[i,f[i]] <- 1
    for (t in (f[i]+1):(noccas-1)){
      # determine the state alive/dead
      z[i,t] ~ dbern(phi * z[i,t-1]) 
      # determine the capture status
      y[i,t] ~ dbern(p * z[i,t]) 
    }
    z[i,noccas] ~ dbern(phi_v*z[i,noccas-1])
    y[i,noccas] ~ dbern(z[i,noccas])
  }
  # Priors
  phi ~ dunif(0,1)
  p ~ dunif(0,1)
  phi_v ~ dunif(0,1)
}


CJS_0_2022 <- jags.parallel(data = jags.data,
                       inits = inits,
                       parameters.to.save = parameters,
                       model.file = cjs_0_2022,
                       n.chains = 4,
                       n.iter = ni)

save(CJS_0_2022, file = "object/CJS_0_2022.RData")
