################## data

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
                  Lake_treatment = Lake_treatment %>%
                    mutate(Lake = as.integer(Lake)) %>% 
                    mutate(Treatment = as.numeric(Treatment)) %>% select(Lake, Treatment),
                  nind = nrow(s),
                  noccas = 6,
                  ni = 100000)


#################### Model 1
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

#################### Model 2
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
    pKla[tr] ~ dunif(0,1)
    pLla[tr] ~ dunif(0,1)
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

################### Model 3
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
