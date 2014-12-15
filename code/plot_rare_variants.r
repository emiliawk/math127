# this code takes in the filtered data files obtained by running parseVcf.py and 
# it then filters for rare variants and produces another 345 plots of the minor allele frequency as a function of position
# this makes the plots for rare variants only while makePlots.r graphs both all variants and rare variants


for (chr in 4:22) {
  # paths and names ---------------------------------------------------------
  
  input <- paste0('~/data/filtered/chr.',chr,'_analysis')
  output <- '/Users/Emilka/Desktop/math127_finalProject/plots/release2/'
  data <- read.table(input,stringsAsFactors = FALSE)
  
  overall_frequency <- data[,13]
  populations <- c('AFR', 'AMR', 'EAS', 'FIN', 'NFE', 'SAS')
  
  AC_orig <- data[,c(4,5,7,8,11,12)]
  colnames(AC_orig) <- populations
  AN_orig <- data[,c(15,16,18,19,20,21)]
  colnames(AN_orig) <- populations
  
  # functions ---------------------------------------------------------------
  
  get_AC <- function(listOfNames) {
    tmp <- strsplit(listOfNames, "=")
    ac <- sapply(tmp, "[", 2)
    return(as.double(ac))
  }
  
  
  # script ------------------------------------------------------------------
  
  # extract the numbers
  n <- nrow(AC_orig)
  
  AC <- data.frame(AFR=numeric(n), AMR=numeric(n), EAS=numeric(n), FIN=numeric(n), NFE=numeric(n), SAS=numeric(n))
  AN <- data.frame(AFR=numeric(n), AMR=numeric(n), EAS=numeric(n), FIN=numeric(n), NFE=numeric(n), SAS=numeric(n))
  frequency <- data.frame(position = data$V2, AFR=numeric(n), AMR=numeric(n), EAS=numeric(n), FIN=numeric(n), NFE=numeric(n), SAS=numeric(n))
  
  for (i in 1:ncol(AC_orig)) {
    tmp <- AC_orig[,i]
    AC[,i] <- get_AC(tmp)
    tmp <- AN_orig[,i]
    AN[,i] <- get_AC(tmp)
    frequency[,i+1] <- AC[,i]/AN[,i]
  }
  overall_frequency <- get_AC(overall_frequency)
  # add overal_frequency info to data frame
  frequency$overall_freq <- overall_frequency
  
  # plots -------------------------------------------------------------------
  
  # remove missing entries just in case (not quite sure how plot handles NAs)
  frequency <- na.omit(frequency)
  position <- frequency$position
  names <- names(frequency)
  
  rare_variants <- frequency[frequency$overall_freq <= 0.01,]
  
  for (i in 2:7) {
    for (j in 2:7) {
      if (j > i) {
        filename <- paste0(output,'chr', chr, '/', names[i], '_', names[j], '.png')
        png(file = filename,width=1164, height=578)
        plot(rare_variants[,1], rare_variants[,i], type='l', lwd=3, col='red', lty=1, main = paste('Rare Variants Chromosome', chr), xlab='bp', ylab = 'Minor Allele Frequency', cex.axis=1.5,cex.lab=1.5,cex=1.5, cex.main=1.5, ylim=c(0,1))
        points(rare_variants[,1], rare_variants[,j], type='l', lwd=3, col='cornflowerblue', lty=1)
        legend('topright', c(names[i],names[j]),col=c("red","cornflowerblue"), lty=1, cex = 1, inset=0.01)
        dev.off()
      }
    }
  }
}