#load required packages
library(plyr)
library(qqman)

#set working directory
setwd("~/path/to/your/working/directory/")

#prepare chromosome.number dataframe to substitute non numerical chromosome ID with integers
#This step is required for species where chromosome IDs are non-integers
chromosome.numbers <- read.table("contig.to.chr", stringsAsFactors = F)
chromosome.numbers <- chromosome.numbers[,c(1,3)]
chromosome.numbers$V1 <- c(1:21)
chromosome.numbers <- chromosome.numbers[1:19, ]

#load tassel output
tassel.out.example <- read.table("Tassel_output_example_stats.txt",
                           sep = "\t", header = T, stringsAsFactors = F)[-1,]

#
manh.qq <- function(arg_1){
  arg_2 <- arg_1[-c(grep("NW", arg_1$Chr)),]
  arg_2$Chr <- mapvalues(arg_2$Chr, from = chromosome.numbers$V3, to = chromosome.numbers$V1)
  arg_2 <- arg_2[which(arg_2$Chr %in% c(1:19)),]
  for (trait in unique(arg_2$Trait)){
    trait.df <- subset(arg_2, arg_2$Trait == trait)
    trait.df$Pos <- as.numeric(trait.df$Pos)
    trait.df$Chr <- as.numeric(trait.df$Chr)
    trait.df <- trait.df[which(complete.cases(trait.df$p)),]
    fdr <- as.data.frame(cbind(p.adjust(trait.df$p, method = "fdr"), trait.df$p))
    fdr <- max(fdr$V2[which(fdr$V1 < 0.05)])
    par(mfrow = c(2, 1), tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))
    manhattan(trait.df, chr="Chr", bp="Pos", p="p", snp = "Marker", 
              genomewideline = -log10(0.05/(nrow(trait.df))), suggestiveline = F)
    abline(h = -log10(fdr), col = "blue")
    mtext(trait, side = 3)
    qq(trait.df$p) 
}
}

#print pdf for one trait per sheet
pdf("Manhattan_qq_example.pdf")
manh.qq(tassel.out.example)
dev.off()
