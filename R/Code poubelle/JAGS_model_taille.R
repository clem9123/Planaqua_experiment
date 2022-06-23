
################## data
tr <- BDD_a %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
mooved <- tr[duplicated(tr$Tag_id),]$Tag_id
tr <- tr %>% filter ( !duplicated(tr$Tag_id))
tr <- tr %>% mutate (Treatment = ifelse(Tag_id %in% mooved, NA, Treatment))

la <- BDD_a %>% group_by(Tag_id, Lake) %>% summarize() %>% ungroup()
mooved <- la[duplicated(la$Tag_id),]$Tag_id
la <- la %>% filter ( !duplicated(la$Tag_id))
la <- la %>% mutate (Lake = ifelse(Tag_id %in% mooved, NA, Lake))

s <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from = Size)
s <- s %>% merge(la) %>% merge(tr)

s <- s %>% filter (!is.na(Treatment) & !is.na(Lake))

set.seed(1234)
s <- s[sample(1:nrow(s), 500),]

get.first <- function(x) min(which(!is.na(x)))
f <- apply(s, 1, get.first)

jags.data <- list(y = s[2:7],
                  f = f,
                  Lake = s$Lake,
                  Treatment = s$Treatment,
                  Lt = data.frame(Lake = as.numeric(as.character(Lake_treatment$Lake)), 
                                  Treatment = as.numeric(as.character(Lake_treatment$Treatment))),
                  nind = nrow(s),
                  noccas = 6,
                  ni = 1000)

################## Model 1

inits <- function(){
  list(t0 = runif(nind, 2014, 2015), K = runif(1,0,10), Linf = runif(1,130,250), p = runif(1,0,1))}

parameters = c("K","Linf","p")

model_g_t0i <- function ()
{
  #Priors
  p ~ dunif(0,1)
  Linf ~ dunif (130,250)
  K ~ dunif(0,10)
  for (i in 1:nind){
    t0[i] ~ dunif(2014,f[i]+2015)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linf*(1-exp(-K*(2015+t-t0[i]))), p)
    } # t
  } # i
} # func

Model_g_t0i <- jags.parallel(data = jags.data,
                                     inits = inits,
                                     parameters.to.save = parameters,
                                     model.file = model_g_t0i,
                                     n.chains = 4,
                                     n.iter = ni)

#print(Model_g_t0i, digits = 3)
#traceplot(Model_g_t0i, ask =FALSE)
#autocorr.plot(Model_g_t0i)

save(Model_g_t0i, file = "R/object/Model_g_t0i.RData")

################## Model 2

inits <- function(){
  list(t0 = runif(nind, 2014, 2015), K = runif(4,0,10), Linf = runif(1,130,250), p = runif(1,0,1))}

parameters = c("K","Linf","p")

model_g_Ktr_t0i <- function ()
{
  #Priors
  p ~ dunif(0,1)
  Linf ~ dunif (130,250)
  for(tr in 1:4){
    K[tr] ~ dunif(0,10)
  }
  for (i in 1:nind){
    t0[i] ~ dunif(2014,f[i]+2015)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linf*(1-exp(-K[Treatment[i]]*(2015+t-t0[i]))), p)
    } # t
  } # i
} # func

Model_g_Ktr_t0i <- jags.parallel(data = jags.data,
                                     inits = inits,
                                     parameters.to.save = parameters,
                                     model.file = model_g_Ktr_t0i,
                                     n.chains = 4,
                                     n.iter = ni)

#print(Model_g_Ktr_t0i, digits = 3)
#traceplot(Model_g_Ktr_t0i, ask =FALSE)
#autocorr.plot(Model_g_Ktr_t0i)

save(Model_g_Ktr_t0i, file = "R/object/Model_g_Ktr_t0i.RData")

################## Model 3

inits <- function(){
  list(t0 = runif(nind, 2014, 2015), K = runif(1,0,10), Linf = runif(4,130,250), p = runif(1,0,1))}

parameters = c("K","Linf","p")

model_g_Ltr_t0i <- function ()
{
  #Priors
  p ~ dunif(0,1)
  K ~ dunif(0,10)
  for(tr in 1:4){
    Linf[tr] ~ dunif (130,250)
  }
  for (i in 1:nind){
    t0[i] ~ dunif(2014,f[i]+2015)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linf[Treatment[i]]*(1-exp(-K*(2015+t-t0[i]))), p)
    } # t
  } # i
} # func

Model_g_Ltr_t0i <- jags.parallel(data = jags.data,
                                     inits = inits,
                                     parameters.to.save = parameters,
                                     model.file = model_g_Ltr_t0i,
                                     n.chains = 4,
                                     n.iter = ni)

#print(Model_g_Ltr_t0i, digits = 3)
#traceplot(Model_g_Ltr_t0i, ask =FALSE)
#autocorr.plot(Model_g_Ltr_t0i)

save(Model_g_Ltr_t0i, file = "R/object/Model_g_Ltr_t0i.RData")

################## Model 4

inits <- function(){
  list(t0 = runif(nind, 2014, 2015), K = runif(4,0,10), Linf = runif(4,130,250), p = runif(1,0,1))}

parameters = c("K","Linf","p")

model_g_Ktr_Ltr_t0i <- function ()
{
  #Priors
  p ~ dunif(0,1)
  for(tr in 1:4){
    Linf[tr]~ dunif (130,250)
    K[tr] ~ dunif(0,10)
  }
  for (i in 1:nind){
    t0[i] ~ dunif(2014,f[i]+2015)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linf[Treatment[i]]*(1-exp(-K[Treatment[i]]*(2015+t-t0[i]))), p)
    } # t
  } # i
} # func

Model_g_Ktr_Ltr_t0i <- jags.parallel(data = jags.data,
                                     inits = inits,
                                     parameters.to.save = parameters,
                                     model.file = model_g_Ktr_Ltr_t0i,
                                     n.chains = 4,
                                     n.iter = ni)

#print(Model_g_Ktr_Ltr_t0i, digits = 3)
#traceplot(Model_g_Ktr_Ltr_t0i, ask =FALSE)
#autocorr.plot(Model_g_Ktr_Ltr_t0i)

save(Model_g_Ktr_Ltr_t0i, file = "R/object/Model_g_Ktr_Ltr_t0i.RData")


################## data 2016

tr <- BDD_a %>% filter (Tag_year == "2016") %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
mooved <- tr[duplicated(tr$Tag_id),]$Tag_id
tr <- tr %>% filter ( !duplicated(tr$Tag_id))
tr <- tr %>% mutate (Treatment = ifelse(Tag_id %in% mooved, NA, Treatment))

la <- BDD_a %>% group_by(Tag_id, Lake) %>% summarize() %>% ungroup()
mooved <- la[duplicated(la$Tag_id),]$Tag_id
la <- la %>% filter ( !duplicated(la$Tag_id))
la <- la %>% mutate (Lake = ifelse(Tag_id %in% mooved, NA, Lake))

s <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from = Size)
s <- s %>% merge(la) %>% merge(tr)

s <- s %>% filter (!is.na(Treatment)) %>% filter(!is.na(Lake))

#set.seed(1234)
#s <- s[sample(1:nrow(s), 100),]

get.first <- function(x) min(which(!is.na(x)))
f <- apply(s, 1, get.first)

jags.data <- list(y = s[2:7],
                  f = f,
                  Lake = s$Lake,
                  Treatment = s$Treatment,
                  nind = nrow(s),
                  noccas = 6,
                  ni = 100000,
                  Lake_treatment = Lake_treatment)

################## Model 5

inits <- function(){
  list(t0 = runif(nind, 2014, 2015), K = runif(1,0,10), Linf = runif(1,130,250), p = runif(1,0,1))}

parameters = c("K","Linf","p")

model_g_t0 <- function ()
{
  #Priors
  p ~ dunif(0,1)
  Linf ~ dunif (130,250)
  K ~ dunif(0,10)
  for (i in 1:nind){
    t0[i] ~ dunif(2014,f[i]+2015)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linf*(1-exp(-K*(2015+t-t0[i]))), p)
    } # t
  } # i
} # func

Model_g_t0i_2016 <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = model_g_t0,
                             n.chains = 4,
                             n.iter = ni)

#print(Model_g_t0i, digits = 3)
#traceplot(Model_g_t0i, ask =FALSE)
#autocorr.plot(Model_g_t0i)

save(Model_g_t0i_2016, file = "R/object/Model_g_t0i_2016.RData")

################## Model 6

inits <- function(){
  list(t0 = runif(nind, 2014, 2015), K = runif(4,0,10), Linf = runif(1,130,250), p = runif(1,0,1))}

parameters = c("K","Linf","p")

model_g_Ktr_t0 <- function ()
{
  #Priors
  p ~ dunif(0,1)
  Linf ~ dunif (130,250)
  for(tr in 1:4){
    K[tr] ~ dunif(0,10)
  }
  for (i in 1:nind){
    t0[i] ~ dunif(2014,f[i]+2015)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linf*(1-exp(-K[Treatment[i]]*(2015+t-t0[i]))), p)
    } # t
  } # i
} # func

Model_g_Ktr_t0i_2016 <- jags.parallel(data = jags.data,
                                 inits = inits,
                                 parameters.to.save = parameters,
                                 model.file = model_g_Ktr_t0,
                                 n.chains = 4,
                                 n.iter = ni)

#print(Model_g_Ktr_t0i, digits = 3)
#traceplot(Model_g_Ktr_t0i, ask =FALSE)
#autocorr.plot(Model_g_Ktr_t0i)

save(Model_g_Ktr_t0i_2016, file = "R/object/Model_g_Ktr_t0i_2016.RData")

################## Model 7

inits <- function(){
  list(t0 = runif(nind, 2014, 2015), K = runif(1,0,10), Linf = runif(4,130,250), p = runif(1,0,1))}

parameters = c("K","Linf","p")

model_g_Ltr_t0 <- function ()
{
  #Priors
  p ~ dunif(0,1)
  K ~ dunif(0,10)
  for(tr in 1:4){
    Linf[tr] ~ dunif (130,250)
  }
  for (i in 1:nind){
    t0[i] ~ dunif(2014,f[i]+2015)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linf[Treatment[i]]*(1-exp(-K*(2015+t-t0[i]))), p)
    } # t
  } # i
} # func

Model_g_Ltr_t0i_2016 <- jags.parallel(data = jags.data,
                                 inits = inits,
                                 parameters.to.save = parameters,
                                 model.file = model_g_Ltr_t0,
                                 n.chains = 4,
                                 n.iter = ni)

#print(Model_g_Ltr_t0i, digits = 3)
#traceplot(Model_g_Ltr_t0i, ask =FALSE)
#autocorr.plot(Model_g_Ltr_t0i)

save(Model_g_Ltr_t0i_2016, file = "R/object/Model_g_Ltr_t0i_2016.RData")

################## Model 8

inits <- function(){
  list(t0 = runif(nind, 2014, 2015), K = runif(4,0,10), Linf = runif(4,130,250), p = runif(1,0,1))}

parameters = c("K","Linf","p")

model_g_Ktr_Ltr_t0 <- function ()
{
  #Priors
  p ~ dunif(0,1)
  for(tr in 1:4){
    Linf[tr]~ dunif (130,250)
    K[tr] ~ dunif(0,10)
  }
  for (i in 1:nind){
    t0[i] ~ dunif(2014,f[i]+2015)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linf[Treatment[i]]*(1-exp(-K[Treatment[i]]*(2015+t-t0[i]))), p)
    } # t
  } # i
} # func

Model_g_Ktr_Ltr_t0i_2016 <- jags.parallel(data = jags.data,
                                     inits = inits,
                                     parameters.to.save = parameters,
                                     model.file = model_g_Ktr_Ltr_t0,
                                     n.chains = 4,
                                     n.iter = ni)

#print(Model_g_Ktr_Ltr_t0i, digits = 3)
#traceplot(Model_g_Ktr_Ltr_t0i, ask =FALSE)
#autocorr.plot(Model_g_Ktr_Ltr_t0i)

save(Model_g_Ktr_Ltr_t0i_2016, file = "R/object/Model_g_Ktr_Ltr_t0i_2016.RData")

#################### Model 9
# 2016 sans t0 variable

inits <- function(){
  list(t0 = runif(1, 2010, 2016), K = runif(16,0,10), Linf = runif(16,130,250), p = runif(1,0,1))}

parameters = c("K","Linf","p", "t0")

model_g_Kla_Lla <- function ()
{
  #Priors
  p ~ dunif(0,1)
  t0 ~ dunif(2010,2016)
  for(la in 1:16){
    Linf[la]~ dunif (130,250)
    K[la] ~ dunif(0,10)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linf[Lake[i]]*(1-exp(-K[Lake[i]]*(2015+t-t0))), p)
    } # t
  } # i
} # func

Model_g_Kla_Lla_2016 <- jags.parallel(data = jags.data,
                                      inits = inits,
                                      parameters.to.save = parameters,
                                      model.file = model_g_Kla_Lla,
                                      n.chains = 4,
                                      n.iter = ni)

#print(Model_g_Ktr_Ltr_t0i, digits = 3)
#traceplot(Model_g_Ktr_Ltr_t0i, ask =FALSE)
#autocorr.plot(Model_g_Ktr_Ltr_t0i)

save(Model_g_Kla_Lla_2016, file = "R/object/Model_g_Kla_Lla_2016.RData")

####################### Model 10

inits <- function(){
  list(t0 = runif(nind, 2010, 2016), K = runif(16,0,10), Linf = runif(16,130,250), p = runif(1,0,1))}

parameters = c("K","Linf","p")

model_g_Kla_Lla_t0i <- function ()
{
  #Priors
  p ~ dunif(0,1)
  for (i in 1: nind){
    t0[i] ~ dunif(2010,2016)
  }
  
  for(la in 1:16){
    Linf[la]~ dunif (130,250)
    K[la] ~ dunif(0,10)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linf[Lake[i]]*(1-exp(-K[Lake[i]]*(2015+t-t0[i]))), p)
    } # t
  } # i
} # func

Model_g_Kla_Lla_t0i_2016 <- jags.parallel(data = jags.data,
                                          inits = inits,
                                          parameters.to.save = parameters,
                                          model.file = model_g_Kla_Lla_t0i,
                                          n.chains = 4,
                                          n.iter = ni)

#print(Model_g_Ktr_Ltr_t0i, digits = 3)
#traceplot(Model_g_Ktr_Ltr_t0i, ask =FALSE)
#autocorr.plot(Model_g_Ktr_Ltr_t0i)

save(Model_g_Kla_Lla_t0i_2016, file = "R/object/Model_g_Kla_Lla_t0i_2016.RData")

################### Model 11
# 2016 sans t0 variable

inits <- function(){
  list(t0 = runif(1, 2010, 2016), K = runif(4,0,10), Linf = runif(4,130,250), p = runif(1,0,1))}

parameters = c("K","Linf","p","t0")

model_g_Ktr_Ltr <- function ()
{
  #Priors
  p ~ dunif(0,1)
  t0 ~ dunif(2010,2016)
  for(tr in 1:4){
    Linf[tr]~ dunif (130,250)
    K[tr] ~ dunif(0,10)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linf[Treatment[i]]*(1-exp(-K[Treatment[i]]*(2015+t-t0))), p)
    } # t
  } # i
} # func

Model_g_Ktr_Ltr_2016 <- jags.parallel(data = jags.data,
                                      inits = inits,
                                      parameters.to.save = parameters,
                                      model.file = model_g_Ktr_Ltr,
                                      n.chains = 4,
                                      n.iter = ni)

#print(Model_g_Ktr_Ltr_t0i, digits = 3)
#traceplot(Model_g_Ktr_Ltr_t0i, ask =FALSE)
#autocorr.plot(Model_g_Ktr_Ltr_t0i)

save(Model_g_Ktr_Ltr_2016, file = "R/object/Model_g_Ktr_Ltr_2016.RData")

#################### Model 12
# 2016 avec Lake variable et pas t0

inits <- function(){
  list(t0 = runif(1, 2010, 2016), 
       K = runif(4,0,10), 
       Kla = runif(16,0,10),
       pKla = runif(4,0,1),
       Linf = runif(4,130,250), 
       Linfla = runif(16,0,10),
       pLla = runif(4,0,1),
       p = runif(1,0,1))}

parameters = c("K","Linf","p","Kla","Linfla", "pKla", "pLla","t0")

model_g_Ktrla_Ltrla <- function ()
{
  #Priors
  p ~ dunif(0,1)
  t0 ~ dunif(2010,2016)
  for(tr in 1:4){
    Linf[tr] ~ dunif (130,250)
    K[tr] ~ dunif(0,10)
    pKla[tr] ~ dunif(0,100)
    pLla[tr] ~ dunif(0,100)
  }
  for (la in 1:16){
    Linfla[la] ~ dnorm (Linf[Lake_treatment[la,2]], pLla[Lake_treatment[la,2]])
    Kla[la] ~ dnorm (K[Lake_treatment[la,2]], pKla[Lake_treatment[la,2]])
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linfla[Lake[i]]*(1-exp(-Kla[Lake[i]]*(2015+t-t0))), p)
    } # t
  } # i
} # func

Model_g_Ktrla_Ltrla_2016 <- jags.parallel(data = jags.data,
                                          inits = inits,
                                          parameters.to.save = parameters,
                                          model.file = model_g_Ktrla_Ltrla,
                                          n.chains = 4,
                                          n.iter = ni)

#print(Model_g_Ktr_Ltr_t0i, digits = 3)
#traceplot(Model_g_Ktr_Ltr_t0i, ask =FALSE)
#autocorr.plot(Model_g_Ktr_Ltr_t0i)

save(Model_g_Ktrla_Ltrla_2016, file = "R/object/Model_g_Ktrla_Ltrla_2016.RData")


################### Model 13
# 2016 avec Lake variable et t0

inits <- function(){
  list(t0 = runif(nind, 2014, 2015), 
       K = runif(4,0,10), 
       Kla = runif(16,0,10),
       pKla = runif(4,0,1),
       Linf = runif(4,130,250), 
       Linfla = runif(16,0,10),
       pLla = runif(4,0,1),
       p = runif(1,0,1))}

parameters = c("K","Linf","p","Kla","Linfla", "pKla", "pLla")

model_g_Ktrla_Ltrla_t0 <- function ()
{
  #Priors
  p ~ dunif(0,1)
  for(tr in 1:4){
    Linf[tr]~ dunif (130,250)
    K[tr] ~ dunif(0,10)
    pKla[tr] ~ dunif(0,1)
    pLla[tr] ~ dunif(0,1)
  }
  for (la in 1:16){
    Linfla[la] ~ dnorm (Linf[Lake_treatment[la,2]], pLla[Lake_treatment[la,2]])
    Kla[la] ~ dnorm (K[Lake_treatment[la,2]], pKla[Lake_treatment[la,2]])
  }
  for (i in 1:nind){
    t0[i] ~ dunif(2014,f[i]+2015)
  }
  #Likelihood
  for (i in 1:nind){
    for (t in f[i]:noccas){
      y[i,t] ~ dnorm(Linfla[Lake[i]]*(1-exp(-Kla[Lake[i]]*(2015+t-t0[i]))), p)
    } # t
  } # i
} # func

Model_g_Ktrla_Ltrla_t0i_2016 <- jags.parallel(data = jags.data,
                                              inits = inits,
                                              parameters.to.save = parameters,
                                              model.file = model_g_Ktrla_Ltrla_t0,
                                              n.chains = 4,
                                              n.iter = ni)

#print(Model_g_Ktr_Ltr_t0i, digits = 3)
#traceplot(Model_g_Ktr_Ltr_t0i, ask =FALSE)
#autocorr.plot(Model_g_Ktr_Ltr_t0i)

save(Model_g_Ktrla_Ltrla_t0i_2016, file = "R/object/Model_g_Ktrla_Ltrla_t0i_2016.RData")

