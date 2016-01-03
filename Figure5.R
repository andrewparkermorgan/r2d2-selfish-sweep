#! /usr/bin/env Rscript

library(plyr)
library(expm)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(MASS)
library(scales)

source("themes.R")

## utility functions for meitoic drive model from Hedrick PW (1981) Evolution 35:322-332
hedrick.bound <- function(s, m) (1/2)*(1-s)*(1+2*m)
hedrick.step <- function(p, N, s, m) {
	q <- 1-p
	p + p*q*(s*(4*q-2*m-1)+2*m-1)/(2*(1-2*s*p*q))
}
hedrick.trans <- function(pp,j,N,s,m) {
	p <- pp/(2*N)
	q <- 1-p
	qq <- hedrick.step(p, N, s, m)
	#print(qq)
	prob <- dbinom(j, 2*N, qq)
	return(prob)
}
gens <- c(1,500,1000,10000,1000000,10000000,100000000,1000000000)
hedrick.iterate <- function(N, s, m) {
	n <- 2*N
	tau <- outer(0:n, 0:n, hedrick.trans, s = s, m = m, N = N)
	init <- c(0,1,rep(0,n-1))
	sapply(gens[4], function(k) (init %*% (tau%^%k))[n+1])
}

## compute fixation probabilities for varying s,m,N under Hedrick model
params <- expand.grid(N = c(10,25,50,100,250), s = seq(0, 0.3, 0.05), m = seq(0.5, 1.0, 0.1))
rez <- ddply(params, .(N,s,m), function(d) {
	y <- hedrick.iterate(d$N[1], s = d$s[1], m = d$m[1])
	data.frame(gen = gens[ seq_along(y) ], prob = y, bounded = hedrick.bound(d$s[1], d$m[1]) > 1)
}, .progress = "text")
fixprobs <- rez
## stash result for later
save(fixprobs, params, file = "hedrick.fixprobs.Rdata")

rez.summ <- ddply(fixprobs, .(N,s,m), summarize, prob = max(prob))
ggplot(rez.summ, aes(x = s, y = prob, colour = factor(m))) +
	geom_line() +
	geom_point() +
	facet_grid(N ~ .)

## Figure 5A (phase diagram)
p0 <- ggplot(rez.summ) +
	geom_tile(aes(x = s, y = m, fill = hedrick.bound(s,m) > 1), colour = "white") +
	scale_fill_manual(values = c("grey","black")) +
	guides(fill = FALSE) +
	xlab(expression(atop("","Selection coefficient ("~italic(s)~")"))) +
	ylab(expression(atop("","Transmission ratio ("~italic(m)~")"))) +
	coord_equal() +
	theme_classic() + theme(strip.background = element_blank(),
							axis.line = element_blank())

## Figure 5B (Hedrick model)
p1 <- ggplot(rez.summ) +
	geom_tile(aes(x = s, y = m, fill = prob), colour = "white") +
	scale_fill_distiller("Fixation\nprobability", palette = "RdBu") +
	scale_colour_manual(values = c(NA,"black")) +
	facet_grid(. ~ N, labeller = label_bquote(bold(bolditalic(N) == .(x)))) +
	xlab(expression(atop("","Selection coefficient ("~italic(s)~")"))) +
	ylab(expression(atop("","Transmission ratio ("~italic(m)~")"))) +
	coord_equal() +
	theme_classic() + theme(strip.background = element_blank(),
							axis.line = element_blank())

## load forward simulations from simtrd.py
fixtime <- rbind( transform(read.table("simulations/n100.t5.s80.p80.m0.txt", header = TRUE), model = "m0", freq = 0),
				  transform(read.table("simulations/n100.t5.s80.p80.m4.txt", header = TRUE), model = "m4", freq = 0.01),
				  transform(read.table("simulations/n100.t5.s80.p80.m3.txt", header = TRUE), model = "m3", freq = 0.1),
				  transform(read.table("simulations/n100.t5.s80.p80.m1.txt", header = TRUE), model = "m1", freq = 0.2),
				  transform(read.table("simulations/n100.t5.s80.p80.m2.txt", header = TRUE), model = "m2", freq = 0.5) )

logit <- function(x) log(x/(1-x))

fixprobs.sim <- ddply(fixtime, .(model, freq), function(d) {
	x <- fitdistr(d$attempts, "geometric")
	rez <- data.frame(value = x$estimate, se = x$sd)
	transform(rez, mid = logit(value), lo = logit(value-2*se), hi = logit(value+2*se))
})

## Figure 5C: forward simulation results
p2 <- ggplot(fixprobs.sim) +
	geom_hline(yintercept = logit(1/200), lty = "dashed", colour = "grey") +
	geom_pointrange(aes(x = factor(freq), y = mid, ymax = hi, ymin = lo,
						colour = freq > 0)) +
	scale_colour_manual(values = c("darkgrey","black")) +
	scale_x_discrete(labels = function(x) ifelse(x == 0, "(no modifier)",
												   sprintf("%.2f", as.numeric(as.character(x))))) +
	guides(colour = FALSE) +
	xlab("\nModifier allele frequency") + ylab("Log-odds of fixation\n") +
	theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Figure 5D: forward simulation results
p3 <- ggplot(fixtime) +
	geom_boxplot(aes(x = factor(freq), y = generations, fill = freq != 0)) +
	scale_fill_manual(values = c("grey90","grey40")) +
	scale_y_log10("Generations to fixation\n", breaks = c(10,50,100,500,1000)) +
	scale_x_discrete("\nModifier allele frequency", labels = function(x) ifelse(x == 0, "(no modifier)",
												 sprintf("%.2f", as.numeric(as.character(x))))) +
	guides(fill = FALSE) +
	theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Figure 5E: example frequency trajectory from simulation
onerun <- read.table("simulations/collapse.example.txt", header = TRUE)
onerun.m <- melt(onerun, id.var = "generation")
(p4 <- ggplot(onerun.m) +
	geom_line(aes(x = generation, y = value, colour = variable)) +
	scale_colour_manual("allele", values = c("darkgrey","darkblue")) +
	scale_y_continuous(limits = c(0,1)) +
	ylab("Allele frequency\n") + xlab("\nGeneration") +
	theme_clean("topleft") )

pp1 <- arrangeGrob(panel.label.grob("A", fontface = "bold"), p0,
				   panel.label.grob("B", fontface = "bold"), p1,
				   nrow = 1, widths = c(0.2,1,0.2,4.2))
pp2 <- arrangeGrob(panel.label.grob("C", fontface = "bold"), p2,
				   panel.label.grob("D", fontface = "bold"), p3,
				   panel.label.grob("E", fontface = "bold"), p4,
				   nrow = 1, widths = c(0.2,1,0.2,1,0.2,2))

pdf("figures/fixation.probs.pdf", width = 12, height = 6)
grid.arrange(pp1, pp2, nrow = 2, heights = c(1,1))
dev.off()
