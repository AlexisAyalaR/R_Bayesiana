#Regresión Poisson

#log(E(y)) = lambda = b0 + b1x1 + ... bnxn 

library("rjags")
library("COUNT")
data("badhealth")
?badhealth

#Modelo 1 con age*badhealth
mod_string = " model {
    for (i in 1:length(numvisit)) {
        numvisit[i] ~ dpois(lam[i])
        log(lam[i]) = int + b_badh*badh[i] + b_age*age[i] + b_intx*age[i]*badh[i]
        }
    int ~ dnorm(0, 1/1e6)
    b_badh ~ dnorm(0, 1/1e4)
    b_age ~ dnorm(0, 1/1e4)
    b_intx ~ dnorm(0, 1/1e4)
} "

set.seed(92)

data_jags = as.list(badhealth)

params = c("int", "b_badh", "b_age", "b_intx")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)

update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                        variable.names=params,
                        n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))



#Modelo 2 sin age*badhealth
mod1_string = " model {
    for (i in 1:length(calls)) {
		calls[i] ~ dpois( days_active[i] * lam[i] )
		log(lam[i]) = b0 + b[1]*age[i] + b[2]*isgroup2[i]
	}
    int ~ dnorm(0, 1/1e6)
    b_badh ~ dnorm(0, 1/1e4)
    b_age ~ dnorm(0, 1/1e4)
} "

set.seed(92)

data_jags = as.list(badhealth)

params1 = c("int", "b_badh", "b_age")

mod1 = jags.model(textConnection(mod1_string), data=data_jags, n.chains=3)

update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                       variable.names=params1,
                       n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))


#DIC 

dic1 <- dic.samples(mod, 1e3)
dic2 <- dic.samples(mod1, 1e3)

dic1
dic2


ppois(21, 30)

dat = read.csv(file="Desktop/Bayesiana_R/callers.csv", header=TRUE)
head(dat)
pairs(dat)

plot(calls/days_active ~ isgroup2, data = dat)
plot(age ~ isgroup2, data = dat)

mod2_string = " model {
    for (i in 1:length(calls)) {
        calls[i] ~ dpois( lam[i] )
        log(lam[i]) = b0 + b[1]*age[i] + b[2]*isgroup2[i]
    }
    b0 ~ dnorm(0, 10^2)
    
    for (j in 1:2) {
    b[j] ~ dnorm(0, 10^2)
    }
} "


set.seed(92)

data2_jags = as.list(dat)

params2 = c("b0", "b")

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, n.chains=3)

update(mod2, 5e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

## Convergencia
par(mfrow=c(3,2))
plot(mod2_sim, ask=TRUE)

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)

coefs = colMeans(mod2_csim)

mean(mod2_csim[,2] > 0)

