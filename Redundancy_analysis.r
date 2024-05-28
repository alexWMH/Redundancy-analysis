#!/usr/bin/env Rscript 
#usage: Rscript Redundancy_analysis.r -t [transpose taxonomy table] -e [env data] -g [group]
#### library package 

library("vegan")
library("psych")
library("ggplot2")
library("ggrepel")
library("optparse")
library("dplyr")

# Get Option 
option_list <- list(
					make_option(c("-t","--table"), type = "character", help = "\ttranspose taxonomy table\n"),
					make_option(c("-e","--env"), type = "character", help = "\tenv data\n"),
					make_option(c("-g", "--group"), type = "character", help = "\tgroup\n")
)
custom_usage <- "\n\tRedundancy_analysis.r -t [transpose taxonomy table] -e [env data] -g [group]\n"
parser <- OptionParser(option_list = option_list, usage = custom_usage)
arguments <- parse_args(parser, positional_arguments = 0)
opt <- arguments$options

#if (!is.null(opt$help)) {
	##### load file 
	env <- read.csv(opt$env, sep = "\t", header = T, row.names = 1)
	table <- read.csv(opt$table, sep = "\t", header = T, row.names = 1)
	spe_helling <- decostand(table, method="hellinger")
	##### RDA 
	max <- 0
	env_data <- env[c(3:ncol(env))]
	
	##### env pair plot before vif
	png(paste0("pairs_before_vif.png"),width=2000, height=2000, res= 120)
	pairs.panels(env_data, scale=T)
	dev.off()

	while (max == 0 || max > 10) { 
		RDA <- rda(spe_helling ~ ., data = env_data, scale = T )
		vif <- vif.cca(RDA)
		print (vif)
		max <- 0
		for (i in 1:dim(env_data)[2]){
			if(vif[i]>max){
				maxindex <- i
				max <- vif[i]
			}
		}
		if(max > 10) {
			env_data <- env_data[,-c(maxindex)]
		}
	}
	print (RDA)
	
	###### cca to see envs regression 
	cca <- anova.cca(RDA, step = 1000, by = "term")   #by = "term", "axis"
	sink("cca_out.txt")
	cca
	sink()


	##### env pair plot after vif
	png(paste0("pairs_after_vif.png"),width=2000, height=2000, res= 120)
	pairs.panels(env_data, scale=T)
	dev.off()

	# preprocess
	p <- plot(RDA, type = "n", scaling=3)
	env[,opt$group] <- as.factor(env[,opt$group])
	levels(env[,opt$group]) <- sort(unique(env[,opt$group]))
	eco <- env[,opt$group]
	bg <- c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c") # 6 nice colors for our ecotypes

	# change RDA information into dataframe
	species <- as.data.frame(p$species)
	sites <- as.data.frame(p$sites)
	bp <- as.data.frame(p$biplot)
	bind  <- bind_rows(sites,species)
	
	# plot title and x y axis 
	Proportion <- round((RDA$CCA$tot.chi/RDA$tot.chi)*100, digits = 2)
	rda1 <- round(RDA$CCA$eig[1]/sum(RDA$CCA$eig)*100, digits = 2)
	rda2 <- round(RDA$CCA$eig[2]/sum(RDA$CCA$eig)*100, digits = 2)

	# plot RDA 

	ggplot(sites,aes(x= RDA1, y = RDA2)) + 
		geom_hline(yintercept = 0, linewidth=.2, linetype = "dashed") + 
		geom_vline(xintercept = 0, linewidth=.2, linetype = "dashed") + 
		geom_segment(data = bp, aes(x = 0, y = 0, xend = RDA1*2.5, yend = RDA2*2.5), arrow = arrow(length = unit(0.5,"cm")), colour = "#0868ac") +
		geom_text_repel(aes(label = rownames(sites)),box.padding = 0.3,min.segment.length = 0,segment.color = 'black') +   # aes(label = rownames(sites),fill = eco)
		geom_text_repel(data = bp*2.65,aes(label = rownames(bp)), colour = "#0868ac") +
		geom_point(aes(colour = eco),size = 3,show.legend = TRUE) +
		geom_point(data = species,aes(x = RDA1, y = RDA2), colour = "black", size = 2) + 
		geom_text_repel(data = species,aes(label = rownames(species))) +
		scale_y_continuous() + 
		theme_bw() + 
		labs(title = paste0("Constrained = ",Proportion,"%")) +
		xlab(paste0("RDA1 (",rda1,"%)")) + 
		ylab(paste0("RDA2 (",rda2, "%)")) + 
		theme(legend.position="none", plot.title=element_text(hjust=0.5, size = 30),
			  axis.title.x = element_text(size = 20),
  			  axis.title.y = element_text(size = 20)) # dont show legend : legend.position = "none"
	ggsave("my_RDA.png", width = 15, height=10)
