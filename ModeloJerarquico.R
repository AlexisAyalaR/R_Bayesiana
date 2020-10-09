#Modelos Jerárquicos

library("rjags")
dat = read.csv(file="Desktop/Bayesiana_R/pctgrowth.csv", header=TRUE)
head(dat)


#Modelo 
mod_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dnorm(theta[grp[i]], sig2)
    }
    for (j in 1:5) {
        theta[j] ~ dnorm(mu, tau2)
    }
    mu ~ dnorm(0, 1e6)
    tau2 ~ dgamma(1/2, 1*3/2)
    sig2 ~ dgamma(2/2, 2*1/2)
} "

set.seed(111)

data_jags = as.list(dat)

params = c("theta", "mu", "tau2", "sig2")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)

update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

summary(mod_csim)

par(mfrow=c(3,2))
par(mar=c(1,1,1,1))
plot(mod_sim, ask=TRUE)


#Comparamos con un modelo anova
means_anova = tapply(dat$y, INDEX=dat$grp, FUN=mean)

means_jags = colMeans(mod_csim)[4:8]

plot(means_anova)
points(means_jags, col="red")





#Ejemplo 2
#Modelo jerárquico en regresión lineal

library("MASS")
data("OME")

dat1 = subset(OME, OME != "N/A")
dat1$OME = factor(dat1$OME) # relabel OME
dat1$ID = as.numeric(factor(dat1$ID)) # relabel ID so there are no gaps in numbers (they now go from 1 to 63)

## Original reference model and covariate matrix
mod_glm = glm(Correct/Trials ~ Age + OME + Loud + Noise, data=dat1, weights=Trials, family="binomial")
X = model.matrix(mod_glm)[,-1]

## Original model (that needs to be extended)
mod1_string = " model {
	for (i in 1:length(y)) {
		y[i] ~ dbin(phi[i], n[i])
		logit(phi[i]) = alpha[ID[i]] + b[1]*Age[i] + b[2]*OMElow[i] + b[3]*Loud[i] + b[4]*Noiseincoherent[i]
	}
	for (j in 1:max(ID)) {
	  alpha[j] ~ dnorm(mu, tau2)
	}
	
	mu ~ dnorm(0, 10^2)
	sig2 ~ dgamma(1/2, 1/2)
	tau2 = 1 / sig2 
	
	for (j in 1:4) {
		b[j] ~ dnorm(0.0, 1.0/4.0^2)
	}
	
} "

data1_jags = as.list(as.data.frame(X))
data1_jags$y = dat1$Correct
data1_jags$n = dat1$Trials
data1_jags$ID = dat1$ID

set.seed(20)

data1_jags = as.list(data1_jags)

params1 = c("alpha", "b", "mu", "tau2")

mod1 = jags.model(textConnection(mod1_string), data=data1_jags, n.chains=3)

update(mod1, 1e3)

mod_sim1 = coda.samples(model=mod1,
                       variable.names=params1,
                       n.iter=5e3)
mod_csim1 = as.mcmc(do.call(rbind, mod_sim1))

## Convergencia
par(mfrow=c(3,2))
plot(mod_sim1, ask=TRUE)

gelman.diag(mod_sim1)
autocorr.diag(mod_sim1)
autocorr.plot(mod_sim1)
effectiveSize(mod_sim1)

dic = dic.samples(mod1, 1e3)
