# n = 0.2
# phi = 0.9
# p = 0.7

simul.data <- function(nind,noccas){
  z = data.frame()
  y = data.frame()
  for (i in 1:nind){
    # State at first capture is alive
    z[i,1] <- rbinom(1,1,0.4)
    y[i,1] <- z[i,1]
    for (t in 2:noccas){
      # determine the state alive/dead
      z[i,t] <- rbinom(1,1,0.9* z[i,t-1] + 0.2*(1-max(z[i,(1:(t-1))]))) #(1-max(z[i,(1:(t-1))]))) 
      # determine the capture status
      y[i,t] <- rbinom(1,1,0.7 * z[i,t])
    }
  }
  return(y)
}

CH <- simul.data(200,6)
#CH <- CH %>% mutate(vide = V1+V2+V3+V4+V5+V6) %>% filter(vide !=0)
CH <- CH[1:6]
head(CH)
CH <- as.matrix(CH)

jags.data <- list(y = CH,
                  #size = as.matrix(s[2:7]),
                  #Treatment = t$t.Treatment,
                  f = apply(CH, 1, get.first), 
                  nind = nrow(CH), 
                  noccas = ncol(CH),
                  ni = 1000,
                  #zi = z.inits(CH),
                  zi2 = z.inits2(CH))



######################################################
#############                       ##################
#############      Simul taille     ##################
#############                       ##################
######################################################



# n = 0.2
# phi = 0.9
# p = 0.7

get.first <- function(x) min(which(x!=0))

simul.data.taille <- function(s){
  y <- data.frame(apply(s,2, function(x) ifelse(is.na(x),0,1)))
  z <- data.frame(apply(s,2, function(x) ifelse(is.na(x),0,1)))
  f <- apply(y, 1, get.first)
  for (i in 1:nrow(s)){
    # State at first capture is alive
    z[i,f[i]] <- 1
    y[i,f[i]] <- 1
    if (f[i]<6){
    for (t in (f[i]+1):ncol(s)){
      #print(s[i,t])
      phi <- invlogit(-0.001*s[i,t]+2)
      p <- invlogit(0.003*s[i,t]+2)
      # determine the state alive/dead
      z[i,t] <- rbinom(1,1,phi* z[i,t-1]) #(1-max(z[i,(1:(t-1))]))) 
      # determine the capture status
      y[i,t] <- rbinom(1,1,p * z[i,t])
    }}
  }
  return(y)
}

CH <- simul.data.taille(s)
#CH <- CH %>% mutate(vide = V1+V2+V3+V4+V5+V6) %>% filter(vide !=0)
CH <- CH[1:6]
head(CH)
CH <- as.matrix(CH)

jags.data <- list(y = CH,
                  size = s,
                  #size = as.matrix(s[2:7]),
                  #Treatment = t$t.Treatment,
                  f = apply(CH, 1, get.first), 
                  nind = nrow(CH), 
                  noccas = ncol(CH),
                  ni = 1000,
                  zi = z.inits(CH))#,
                  #zi2 = z.inits2(CH))