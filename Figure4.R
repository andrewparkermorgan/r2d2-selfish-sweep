#! /usr/bin/env Rscript

library(plyr)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(scales)

load("simulations/HR.sims.Rdata")
hr.rez$fixed <- factor(hr.rez$fixed, labels = c("lost","fixed"))

hr.col <- "darkblue"
(p1 <- ggplot(subset(hr.rez, gen > -1)) +
 	geom_line(aes(x = gen, y = maf, group = run, colour = fixed), alpha = 0.5) +
 	annotate("point", x = 0, y = 0.75, colour = hr.col, fill = "white", size = 4, pch = 21) +
 	scale_colour_manual("fate", values = c("grey",hr.col)) +
 	theme_clean("topright") +
 	xlab("\ngenerations") + ylab("focal allele frequency\n"))

hr.rez.summ <- ddply(hr.rez, .(run), summarize, gen = max(gen), fixed = fixed[1])
med.fix <- median(subset(hr.rez.summ, fixed == "fixed")$gen)
(p2 <- ggplot(hr.rez.summ) +
 	stat_ecdf(aes(colour = fixed, x = gen)) +
 	geom_vline(xintercept = med.fix, colour = "grey30", lty = "dashed") +
 	annotate("text", x = med.fix, y = 0.1, label = paste0(" (",round(med.fix),")"),
 			 hjust = 0, size = 5, fontface = "italic") +
 	scale_x_log10(breaks = c(1, 10, 100, 1000)) +
 	scale_colour_manual("fate", values = c("grey","darkblue")) +
 	theme_clean("topleft") +
 	xlab("\ngenerations to fixation") + ylab("cumulative proportion\n"))

load("simulations/AIL.sims.Rdata")
ail.rez$fixed <- factor(ail.rez$fixed, labels = c("lost","fixed"))

ail.col <- muted("red")
(p3 <- ggplot(subset(ail.rez, gen > -1)) +
 	geom_line(aes(x = gen, y = maf, group = factor(batch):factor(run), colour = fixed), alpha = 0.5) +
 	annotate("point", x = 0, y = 0.50, colour = ail.col, fill = "white", size = 4, pch = 21) +
 	scale_colour_manual("fate", values = c("grey",ail.col)) +
 	theme_clean("topright") +
 	xlab("\ngenerations") + ylab("focal allele frequency\n"))

ail.rez.summ <- ddply(ail.rez, .(run, batch), summarize, gen = max(gen), fixed = fixed[1])
med.fix <- median(subset(ail.rez.summ, fixed == "fixed")$gen)
(p4 <- ggplot(ail.rez.summ) +
 	stat_ecdf(aes(colour = fixed, x = gen)) +
 	geom_vline(xintercept = med.fix, colour = "grey30", lty = "dashed") +
 	annotate("text", x = med.fix, y = 0.1, label = paste0(" (",round(med.fix),")"),
 			 hjust = 0, size = 5, fontface = "italic") +
 	scale_x_log10(breaks = c(1, 10, 100, 1000)) +
 	scale_colour_manual("fate", values = c("grey",ail.col)) +
 	theme_clean("topleft") +
 	xlab("\ngenerations to fixation") + ylab("cumulative proportion\n"))

pdf("figures/hr.simulations.pdf", width = 11, height = 5)
grid.arrange(panel.label.grob("A", "topleft", fontface = "bold"), p1,
			 panel.label.grob("B", "topleft", fontface = "bold"), p2,
			 panel.label.grob("C", "topleft", fontface = "bold"), p3,
			 panel.label.grob("D", "topleft", fontface = "bold"), p4,
			 nrow = 2, ncol = 4, widths = c(1/5,2,1/5,1))
dev.off()