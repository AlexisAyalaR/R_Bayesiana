library("rjags")
library('coda')

##Modelando Bayesiana con Jags

#1.- Especificar el modelo

mod_string = "model {
  for (i in 1:n) {
    y[i] ~ dnorm(mu, 1.0/sig2)
  }
  mu ~ dt(0.0, 1.0/1.0, 1)
  sig2 = 1.0
} "

# 2.- Preparar modelo
set.seed(50)

y <- c(1.2, 1.4, -0.5, 0.3, 0.9, 
       2.3, 1.0, 0.1, 1.3, 1.9)

n <- length(y)

data_jags = list(y=y, n=n)
params = c("mu")

inits = function() {
  inits = list("mu"= 0.0)
}

mod = jags.model(textConnection(mod_string), 
                 data = data_jags, inits = inits)

#3.- Corres MCMC

update(mod, 500)

mod_sim = coda.samples(model = mod, variable.names = params, n.iter = 1000)


#4.- Post procesamiento

plot(mod_sim)
summary(mod_sim)

##### Autocorrelaciones

#Cadena
traceplot(as.mcmc(mod_sim))

#Autocorrelaciones y valores
autocorr.plot(as.mcmc(mod_sim))

autocorr.diag(as.mcmc(mod_sim))

#Effective sample size
effectiveSize(as.mcmc(mod_sim))
raftery.diag(as.mcmc(mod_sim))






##########
#Regresión Lineal

#Regresión frecuentista:

library("car")  # load the 'car' package
data("Anscombe")  # load the data set
Anscombe  # read a description of the data
head(Anscombe)  # look at the first few lines of the data
pairs(Anscombe)  # scatter plots for each pair of variables

plot(education ~ income + young + urban, data = Anscombe)

mod = lm(education ~ income + young + urban, data = Anscombe)





#Regresión Bayesiana

library("rjags")

mod_string = " model {
    for (i in 1:length(education)) {
        education[i] ~ dnorm(mu[i], prec)
        mu[i] = b0 + b[1]*income[i] + b[2]*young[i] + b[3]*urban[i]
    }
    
    b0 ~ dnorm(0.0, 1.0/1.0e6)
    for (i in 1:3) {
        b[i] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(1.0/2.0, 1.0*1500.0/2.0)
    	## Initial guess of variance based on overall
    	## variance of education variable. Uses low prior
    	## effective sample size. Technically, this is not
    	## a true 'prior', but it is not very informative.
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
} "

data_jags = as.list(Anscombe)

params =c("b0","b", "sig")

inits1 = function() {
  inits =list("b0" = rnorm(1, 0, 100), "b"=rnorm(3, 0, 100), "prec"=rgamma(1, 1, 1))
}

#Generamos 3 cadenas
mod = jags.model(textConnection(mod_string), data= data_jags, inits = inits1, n.chains = 3)

update(mod, 1000)

mod_sim = coda.samples(model = mod, variable.names = params, n.iter = 1e5)

#Juntamos las cadenas
mod_csim = do.call(rbind, mod_sim)
means = colMeans(mod_csim)

plot(mod_sim)  



#Convergencia
#Si es mucho mayor a uno, no están convergiendo a la misma distribución
gelman.diag(mod_sim)

autocorr.diag(mod_sim)
  
effectiveSize(mod_sim)




#Análisis de Residuales en lm
modRes = lm(education ~ income, data = Anscombe)

#Residuales
#Si son independientes las muestras entonces no debe de haber patrón
plot(resid(modRes))
plot(predict(modRes), resid(modRes))

#Checar si viene de una normal:
#Si viene debería ser una línea recta
qqnorm(resid(modRes))
  




#Análisis de Residuales en bayesiana
#Primero obtenemos predicciones:
X = cbind(rep(1, length(data_jags$education)), data_jags$income, 
          data_jags$young, data_jags$urban)
head(X)
yhat = drop(X %*% means[1:4])

#Obtenemos residuales
resid1 = data_jags$education - yhat

plot(resid1)
plot(yhat, resid1)
qqnorm(resid1)



#Comparación entre modelos:

#Mod: mu[i] = b0 + b[1]*income[i] + b[2]*young[i] + b[3]*urban[i]
dic.samples(mod, 1e5)

#Mod2: mu[i] = b0 + b[1]*income[i] + b[2]*young[i]
dic.samples(mod2, 1e5)

#Mod3: mu[i] = b0 + b[1]*income[i] + b[2]*young[i] + b[3]*income[i]*young[i]
dic.samples(mod3, 1e5)


#Prob de b1 (income) > 0 en mod

mean(mod_csim[,1] > 0)


###### 
#Mod2
mod_string = " model {
    for (i in 1:length(education)) {
        education[i] ~ dnorm(mu[i], prec)
        mu[i] = b0 + b[1]*income[i] + b[2]*young[i]
    }
    
    b0 ~ dnorm(0.0, 1.0/1.0e6)
    for (i in 1:2) {
        b[i] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(1.0/2.0, 1.0*1500.0/2.0)
    	## Initial guess of variance based on overall
    	## variance of education variable. Uses low prior
    	## effective sample size. Technically, this is not
    	## a true 'prior', but it is not very informative.
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
} "

data_jags = as.list(Anscombe)

params =c("b0","b", "sig")

inits1 = function() {
  inits =list("b0" = rnorm(1, 0, 100), "b"=rnorm(2, 0, 100), "prec"=rgamma(1, 1, 1))
}

mod2 = jags.model(textConnection(mod_string), data= data_jags, inits = inits1, n.chains = 3)

update(mod2, 1000)

mod_sim = coda.samples(model = mod2, variable.names = params, n.iter = 1e5)

#Juntamos las cadenas
mod_csim = do.call(rbind, mod_sim)



###### 
#Mod3

mod_string = " model {
    for (i in 1:length(education)) {
        education[i] ~ dnorm(mu[i], prec)
        mu[i] = b0 + b[1]*income[i] + b[2]*young[i] + b[3]*income[i]*young[i]
    }
    
    b0 ~ dnorm(0.0, 1.0/1.0e6)
    for (i in 1:3) {
        b[i] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(1.0/2.0, 1.0*1500.0/2.0)
    	## Initial guess of variance based on overall
    	## variance of education variable. Uses low prior
    	## effective sample size. Technically, this is not
    	## a true 'prior', but it is not very informative.
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
} "

data_jags = as.list(Anscombe)

params =c("b0","b", "sig")

inits1 = function() {
  inits =list("b0" = rnorm(1, 0, 100), "b"=rnorm(3, 0, 100), "prec"=rgamma(1, 1, 1))
}

mod3 = jags.model(textConnection(mod_string), data= data_jags, inits = inits1, n.chains = 3)

update(mod3, 1000)

mod_sim = coda.samples(model = mod3, variable.names = params, n.iter = 1e5)

#Juntamos las cadenas
mod_csim = do.call(rbind, mod_sim)



