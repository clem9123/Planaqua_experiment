library(R2jags)
print(CJS_0)
print(CJS_1)

rjags::load.module('dic')

recompile(CJS_0)
recompile(CJS_1)

samples <- jags.samples(CJS_0$model, c("WAIC", "deviance"), type = "mean", n.iter = 200)

samples$p_waic <- samples$WAIC
samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
CJS_0_waic <- round(c(waic = tmp[["waic"]], 
                p_waic = tmp[["p_waic"]]),1)


samples <- jags.samples(CJS_1$model, c("WAIC", "deviance"), type = "mean", n.iter = 200)

samples$p_waic <- samples$WAIC
samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
CJS_1_waic <- round(c(waic = tmp[["waic"]], 
                      p_waic = tmp[["p_waic"]]),1)

df <- data.frame(c("CJS_0", "CJS_1"), c(CJS_0[[2]]$DIC, CJS_1[[2]]$DIC), c(CJS_0_waic[1] , CJS_1_waic[1]),c(CJS_0_waic[2] , CJS_1_waic[2]))
colnames(df) <- c("model", "DIC", "WAIC", "p_WAIC")
df
