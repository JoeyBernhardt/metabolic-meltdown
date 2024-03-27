

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

#### Learning about metabolic meltdown predictions


### the TPC

A <- 1
f <- 1

G <- function(x) (A*f*C) - R

G <- function(x, A, f=1) ((A*f*0.4*exp(0.09*x)*(1-((x-15)/(34/2))^2)) - 0.01*exp(0.15*x))

C <- function(x, A, f) A*f*exp(0.07*x)*(1-((x-15)/(34/2))^2)
R <- function(x) 0.02*exp(0.15*x)

briere <- function(x, Tmin = 5, Tmax = 32, m = 1, c = 1) c*(x*(x-Tmin)*(((Tmax-x))^(1/m)))

f <- function(x, fmax = 1, fmin = 0, a = 0.3, xmid = 15) fmax - ((fmax - fmin)/(1 + exp(-a*(x - xmid))))
f(10)

C2 <- function(x, fmax = 1, fmin = 0, a = 0.3, xmid = 15)  (fmax - ((fmax - fmin)/(1 + exp(-a*(x - xmid)))))* exp(0.07*x)*(1-((x-15)/(34/2))^2)
G2 <- function(x, A, fmax = 1, fmin = 0, a = 0.3, xmid = 15) ((A*(fmax - ((fmax - fmin)/(1 + exp(-a*(x - xmid)))))*0.4*exp(0.09*x)*(1-((x-15)/(34/2))^2)) - 0.01*exp(0.15*x))
G <- function(x, A, f=1, z = 15, w = 34, b = 0.09) ((A*f*0.4*exp(b*x)*(1-((x-z)/(w/2))^2)) - 0.01*exp(0.15*x))


briere <- function(x, Tmin, Tmax, m, c){
	res <- c*(x*(x-Tmin)*(((Tmax-x))^(1/m)))
	res
}

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 

p + 
	# stat_function(fun = C, args = list(f = 1, A = 0.3), color = "black", size = 1) +
	# stat_function(fun = C2, color = "black", size = 1) +
	# stat_function(fun = G2, args = list(A = 1), color = "pink", size = 1) +
	# stat_function(fun = briere, args = list(Tmin = 0, Tmax = 30, m = 2, c = 0.0045),  color = "pink", size = 1) +
	# stat_function(fun = briere, args = list(Tmin = 5, Tmax = 35, m = 2.2, c = 0.0040),  color = "purple", size = 1) +
	# stat_function(fun = G, args = list(A = 1), color = "purple", size = 1) +
	stat_function(fun = C, args = list(f = 0.8, A = 0.3), color = "blue", size = 1) +
	# stat_function(fun = C, args = list(f = 0.7, A = 0.3), color = "cadetblue", size = 1) +
	# stat_function(fun = C, args = list(f = 0.6, A = 0.3), color = "lightblue", size = 1) +
	# stat_function(fun = C, args = list(f = 0.5, A = 0.3), color = "lightgreen", size = 1) +
	# stat_function(fun = R, color = "orange", size = 1) +
	stat_function(fun = G, args = list(f = 1, A = 1), color = "black", size = 1) +
	stat_function(fun = G, args = list(f = 1, A = 1, z = 18, b = 0.0815), color = "grey", size = 1) +
	# stat_function(fun = G, args = list(f = 0.7, A = 1), color = "lightblue", size = 1) +
	# stat_function(fun = G, args = list(f = 0.8, A = 1), color = "blue", size = 1) +
	# stat_function(fun = G, args = list(f = 0.6, A = 1), color = "lightgreen", size = 1) +
	xlim(0, 35) + ylim(-1, 3) + 
	xlab("Temperature (°C)") + ylab("Growth rate") + geom_hline(yintercept = 0)
ggsave("figures/metabolic-meltdown-growth.png", width = 6, height = 4)


mod_briere <- function(x, Tmin, Tmax, m, c){
	res <- c*(-x + Tmin + Tmax)*(Tmax - x)*((x - Tmin)^(1/m))
	res
}


briere <- function(x, Tmin, Tmax, m, c){
	res <- c*(x*(x-Tmin)*(((Tmax-x))^(1/m)))
	res
}

briere2 <- function(x, Tmin = 0, Tmax = 22, m = 1.5,  c = 0.004){
	res <- c*(x*(x-Tmin)*(((Tmax-x))^(1/m)))
	res
}

briere3 <- function(x, Tmin = 0, Tmax = 22, m = 1.5,  c = 0.0045){
	res <- c*(x*(x-Tmin)*(((Tmax-x))^(1/m)))
	res
}

cheese <- function(x, alpha = 0.0005, Tmin = 0, Tmax = 35){
	res <- alpha*x*(x - Tmin)*(Tmax - x)
	res
}

cheese_vals <- data.frame(x =seq(from = 0, to = 35, by = 0.1), performance =  sapply(X = seq(from = 0, to = 35, by = 0.1), FUN = cheese))

Tmax <- 35
Tmin <- 0
alpha <- 0.0005
topt <- (1/3)*(Tmax + Tmin + sqrt(Tmax^2 - (Tmin *Tmax) + Tmin^2))

briere_vals <- data.frame(x =seq(from = 0, to = 33, by = 0.1), performance =  sapply(X = seq(from = 0, to = 33, by = 0.1), FUN = briere2))
briere_vals2 <- data.frame(x =seq(from = 0, to = 33, by = 0.1), performance =  sapply(X = seq(from = 0, to = 33, by = 0.1), FUN = briere3))

G3 <- function(x, A = 1, fmax = 1, fmin = 0, a = 0.3, xmid = 15) log(((A*(fmax - ((fmax - fmin)/(1 + exp(-a*(x - xmid)))))*alpha*x*(x - Tmin)*(Tmax - x)) - 0.01*exp(0.15*x)))
G4 <- function(x, A = 1, fmax = 1, fmin = 1, a = 0.3, xmid = 15) log(((A*(fmax - ((fmax - fmin)/(1 + exp(-a*(x - xmid)))))*alpha*x*(x - Tmin)*(Tmax - x)) - 0.01*exp(0.15*x)))


p + 
	stat_function(fun = G3, args = list(A = 1), color = "pink", size = 1) +
	stat_function(fun = G4, args = list(A = 1), color = "purple", size = 1) +
	geom_vline(xintercept = topt, color = "grey", linetype = "dotted") + 
	xlim(0, 35) + ylim(-0.05, 4) +  xlab("Temperature (°C)") + ylab("Fitness") +
	geom_hline(yintercept = 0)


### distribution of TPCs
p + 
	# geom_line(aes(x = x, y = performance), data = cheese_vals, size = 1, color = "cadetblue") + 
	stat_function(fun = cheese, args = list(Tmin = 0, Tmax = 35, alpha = 0.0005), color = "pink", size = 1) +
	stat_function(fun = cheese, args = list(Tmin = 2, Tmax = 37, alpha = 0.0005), color = "purple", size = 1) +
	stat_function(fun = cheese, args = list(Tmin = 0, Tmax = 35, alpha = 0.0007), color = "green", size = 1) +
	stat_function(fun = cheese, args = list(Tmin = 3, Tmax = 40, alpha = 0.0004), color = "blue", size = 1) +
	stat_function(fun = cheese, args = list(Tmin = 0, Tmax = 35, alpha = 0.0006), color = "orange", size = 1) +
	xlim(0, 45) + ylim(-1, 5) +  xlab("Temperature (°C)") + ylab("Population growth rate") +
	geom_hline(yintercept = 0)

p + 
	# geom_line(aes(x = x, y = performance), data = cheese_vals, size = 1, color = "cadetblue") + 
	stat_function(fun = G3, args = list(A = 1), color = "pink", size = 1) +
	stat_function(fun = G4, args = list(A = 1), color = "purple", size = 1) +
	geom_vline(xintercept = topt, color = "grey", linetype = "dotted") + 
	xlim(0, 35) + ylim(-0.05, 1.2) +  xlab("Temperature (°C)") + ylab("Fitness") +
	geom_hline(yintercept = 0)


# ok let’s simulate some distrubution of Tmaxes and some day time  --------

G_norm <- function(x = 15, A = 1, fmax = 1, fmin = 1, a = 0.3, xmid = 15, tmax) (((A*(fmax - ((fmax - fmin)/(1 + exp(-a*(x - xmid)))))*alpha*x*(x - Tmin)*(Tmax - x)) - 0.01*exp(0.15*x)))
G_melt <- function(x = 15, A = 1, fmax = 1, fmin = 0.2, a = 0.3, xmid = 15, tmax) (((A*(fmax - ((fmax - fmin)/(1 + exp(-a*(x - xmid)))))*alpha*x*(x - Tmin)*(Tmax - x)) - 0.01*exp(0.15*x)))

tmaxes <- rnorm(n = 10000, mean = 32, sd = 7)

ggplot() +
	geom_histogram(aes(x = tmaxes))
performance_15_melt <- data.frame(tmax = tmaxes, fitness = sapply(X = tmaxes, FUN = G_melt), type = "melt")
performance_15_norm <- data.frame(tmax = tmaxes, fitness = sapply(X = tmaxes, FUN = G_norm), type = "norm")

all_perf <- bind_rows(performance_15_melt, performance_15_norm) %>% 
	group_by(tmax, type) %>% 
	summarise(mean_fitness = mean(fitness)) %>% 
	mutate(log_mean_fitness = log(mean_fitness))

all_perf %>% 
	ggplot(aes(x = tmax, y = log_mean_fitness, color = type)) + geom_line(size = 1) +
	ylim(-.3, 1.2) + xlim(0, 40) + geom_hline(yintercept = 0) + ylab("ln(mean fitness)")

all_perf %>% 
	filter(type == "norm") %>% 
	filter(mean_fitness > 0) %>% View

omega_melt <- (27.50346 -  4.456408)
omega_norm <- 32.51401 - 4.456408

opt_melt <- 13.89434
opt_norm <- 22.08327

B_melt <- exp(0.2436506)
B_norm <- exp(1.056060)

kcrit_melt <- 7*(sqrt(2*log(B_melt)/omega_melt))
kcrit_norm <- 7*(sqrt(2*log(B_norm)/omega_norm))

## melt
## 9.939055, 18.39713

## normal
## 9.279005, 30.75217

G <- function(x, A, f=1) ((A*f*0.4*exp(0.09*x)*(1-((x-15)/(34/2))^2)) - 0.01*exp(0.15*x))
gain <- function(x, A, f=1) ((A*f*0.4*exp(0.09*x)*(1-((x-15)/(34/2))^2)))
R <- function(x) 0.01*exp(0.15*x)

p + 
	# stat_function(fun = C, args = list(f = 1, A = 0.3), color = "black", size = 1) +
	# stat_function(fun = C, args = list(f = 0.8, A = 0.3), color = "blue", size = 1) +
	# stat_function(fun = C, args = list(f = 0.7, A = 0.3), color = "cadetblue", size = 1) +
	# stat_function(fun = C, args = list(f = 0.6, A = 0.3), color = "lightblue", size = 1) +
	# stat_function(fun = C, args = list(f = 0.5, A = 0.3), color = "lightgreen", size = 1) +
	# stat_function(fun = R, color = "orange", size = 1) +
	stat_function(fun = R, color = "pink", size = 1) +
	stat_function(fun = gain, color = "purple", size = 1) +
	# stat_function(fun = G, args = list(f = 1, A = 1), color = "black", size = 1) +
	# stat_function(fun = G, args = list(f = 0.7, A = 1), color = "cadetblue", size = 1) +
	# stat_function(fun = G, args = list(f = 0.8, A = 1), color = "blue", size = 1) +
	# stat_function(fun = G, args = list(f = 0.6, A = 1), color = "lightblue", size = 1) +
	# stat_function(fun = G, args = list(f = 0.5, A = 1), color = "lightgreen", size = 1) +
	xlim(0, 35) + ylim(-0.1, 3) + 
	xlab("Temperature (°C)") + ylab("ln(fitness)") + geom_hline(yintercept = 0)
ggsave("figures/fitness_functions.png", width = 6, height = 4)

grfunc <- function(x){
	-G(x,f = 0.6, A = 1)
}
optinfo<-optim(c(x= 1),grfunc)
opt<-optinfo$par[[1]]
opt

## 17.44609, 16.09219, 15.17578, 14.00781

#### the fitness function

w_curve <- function(x, B, theta, omega) log(B*exp(-((x - theta)^2)/(2*(omega^2))))
B <- 10
theta <- 15
omega <- 35

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) 

p + 
	stat_function(fun = w_curve, args = list(B = 10, theta = 15, omega = 5), color = "purple", size = 2) +
	stat_function(fun = w_curve, args = list(B = 5, theta = 10, omega = 1), color = "orange", size = 2) +
	# stat_function(fun = w_curve(B = 10, theta = 20, omega = 20), color = "orange", size = 2) +
	ylim(-1, 5) + xlim(0, 30) + ylab("ln(fitness)") + xlab("Trait (Topt)") + coord_cartesian() +
	geom_hline(yintercept = 0)



### normal curve

p1 <- ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
	stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), size = 3, color = "red") + ylab("") +
	scale_y_continuous(breaks = NULL)
p1
ggsave("figures/normal-curve.pdf", width = 6, height = 4)
