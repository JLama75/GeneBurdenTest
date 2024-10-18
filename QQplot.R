#!/usr/bin Rscript

#install.packages("qqman")
library("qqman")
library("ggplot2")
library("dplyr")
library("stringr")
library("data.table")

# Load required libraries
library(optparse)

# Define command-line options
option_list <- list(
  make_option("--file1", type="character", default=NULL, help="Path to file 1", metavar="character"),
  make_option("--file2", type="character", default=NULL, help="Path to file 2", metavar="character"),
  make_option("--file3", type="character", default=NULL, help="Path to file 3", metavar="character"),
  make_option("--file4", type="character", default=NULL, help="Path to file 4", metavar="character"),
  make_option("--trait", type="character", default=NULL, help="Trait name", metavar="character"),
  make_option("--output", type="character", default=NULL, help="path to output file", metavar="character")
)

# Parse the arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

file1=opt$file1
file2=opt$file2
file3=opt$file3
file4=opt$file4
trait=opt$trait
output=opt$output

# Print the arguments for debugging
print(paste("File 1:", opt$file1))
print(paste("File 2:", opt$file2))
print(paste("File 3:", opt$file3))
print(paste("File 4:", opt$file4))
print(paste("Trait:", opt$trait))
print(paste("Trait:", opt$output))

qqplot <- function(df, trait, Run, output){
  
  df$P <- as.numeric(df$P)
  #df <- run1[!is.na(df$P)]
  summary(df$P)
  
  name=paste0(output,"/",trait,".",Run,"_qq.jpeg")
  DIR=
  jpeg(name,
       height=1350,width=1200,res=300)
  qq(df$P)
  
  dev.off()
  
  ##GC Lambda
  
  df$CHISQ <- qchisq(df$P, df = 1, lower = F)
  #df <- df[complete.cases(df$CHISQ), ]
  lambda<- median(df$CHISQ)/qchisq(0.5,1)
  
  # Return lambda and run identifier
  return(data.frame(Run = Run, Lambda = lambda))
  
}

run1 <- data.table::fread(file1)
run2 <- data.table::fread(file2)
run3 <- data.table::fread(file3)
run4 <- data.table::fread(file4)

R1 <- qqplot(run1, trait, "Run1", output)
R2 <- qqplot(run2, trait, "Run2", output)
R3 <- qqplot(run3, trait, "Run3", output)
R4 <- qqplot(run4, trait, "Run4", output)

# Combine the results into a single dataframe
lambda_df <- rbind(R1, R2, R3, R4)
write.table(lambda_df, paste0(output, "/", trait,".lambda_values.tsv"), quote = F, row.names = F)
