library(optparse)
library("stringr")
library(ggplot2)
library(gridExtra)
library("data.table")
library(scales)

### HERE GO FUNCTIONS UTILIZED FOR THIS SCRIPT ###

returnFeaturesFromFileName <- function(filePaths) {
  # This function extract information from path such as parent directory, filename or file extension
  
  path_features <- as.data.frame(str_match(string = filePaths,pattern = "(?:(?:.*\\/)*(.*)\\/)*(.*)\\.(.*)"))
  colnames(path_features) <- c("orig","parent_dir","f-name","ext")
  path_features$parent_dir <- str_replace(path_features$parent_dir,"/","")
  
  if ( all(is.na(path_features)) ) {
    print ("it seems the path is wrong.. returning NULL")
    return (NULL)
  }
  
  return (path_features)
}

returnFPKMsFromMIX2 <- function(files) {
  # This function reads MIX2 outputs from specified paths, merges the data and returns requisted data in 3 different formats: abundance, FPKM or counts
  # input: mix2 output files in .dat format
  # output: data.frame with abundance, FPKM or counts in columns as samples
  
  it <- 1
  
  for (file in files) {
    
    features <- returnFeaturesFromFileName(file)
    
    if (features$ext=="dat") {
      
      system(paste0('grep "SIRV\\|FPKM_CHN\\|ERCC" ',file, ' > metadata_MIX2table.dat'))
      FPKM_table <- read.table(file = "metadata_MIX2table.dat", sep="\t", header = T, stringsAsFactors = F)
      
      temp <- data.frame(x=FPKM_table$FPKM_CHN)
      colnames(temp) <- features$parent_dir
      row.names(temp) <- FPKM_table$tracking_ID
    }
    
    if (it==1) {
      out <- temp
    }else{
      out <- cbind(out,temp)
    }
    
    it <- it + 1
  }
  
  return (out)
}

### SECTION FOR STANDALONE (TERMINAL) USAGE ###

option_list <- list(
  make_option(c("-n","--normalization"),action = "store", type = "character", default = "nm1", 
              help = "Specifying method for normalizing counts. Options: nm1 (default) - normalizing technique based on expected abudance for SIRV, i.e. 1/68 of measured abundance for each transcript"),
  make_option(c("-o","--output"), action = "store", type = "character", default = NA,
              help = "Specificy output name for SIRV boxplots with extenstion (supported .png, .svg, .ps, .eps, .pdf, .jpeg, .tiff, .bmp, .svg)."),
  make_option(c("-l","--sample-list"), action = "store", type = "character", default = NA,
              help = "Specify path to a *.csv file containing sample paths + conditions for technical replicate distinction. If chosen this option no argument should be specified."),
  make_option(c('-a','--average-replicates'), action = "store_true", dest = "aver",
              help = "Specify path to a *.csv file containing sample paths + conditions for technical replicate distinction. If chosen this option no argument should be specified.")
)

parser <- OptionParser(option_list=option_list,usage = "%prog [options] files")

args <- parse_args(parser, positional_arguments = T)
opt <- args$options
file <- args$args

if ( length(file) == 0 & is.na(opt[["sample-list"]]) ) {
  stop("No input files paths provided. Specify files or use -l option to provide file path list. Use -h option for help.")
} else if ( !length(file) == 0 & !is.na(opt[["sample-list"]]) ) {
  stop("Provided both file list and files.. please use only files or file lists. Use -h option for help.")
} else if ( length(file) == 0 & !is.na(opt[["sample-list"]]) ) {
  file_list <- read.csv2(opt[["sample-list"]],stringsAsFactors = F)
  file <- file_list$sample
  condition <- file_list$condition
  print (paste0("Reading from file list: ",opt[["sample-list"]]))
} else if ( length(file) > 0 & is.na(opt[["sample-list"]]) ) {
  print ("Reading files from arguments...")
}

### CONCENTRATION QUANTIFICATION ###

FPKM_table <- returnFPKMsFromMIX2(file)
samples <- colnames(FPKM_table)
method <- opt$normalization

if (!is.null(opt$aver)) {
  
  condition_table <- data.frame(row.names = row.names(FPKM_table))
  
  for (cond in unique(condition)) {
    condition_selection <- FPKM_table[,samples[condition == cond]]
    condition_table <- cbind(condition_table, rowMeans(condition_selection))
  }
  
  FPKM_table <- condition_table
  colnames(FPKM_table) <- unique(condition)
  
  samples <- unique(condition)
}


print (paste0("Quantifying SIRV concentration using method: ", opt[["normalization"]]))

spikeData <- list()
spikeData[["SIRV"]][["FPKM"]] <- FPKM_table[startsWith(row.names(FPKM_table),"SIRV"), , drop=FALSE]
spikeData[["ERCC"]][["FPKM"]] <- FPKM_table[which(grepl("DQ|EF",row.names(FPKM_table))), , drop=FALSE]
out_temp <- data.frame(id=c(),FPKM_norm=c(),sample=c())


for (sample in samples) {
  if (method == "nm1") {
    FPKM_temp <- spikeData[["SIRV"]][["FPKM"]][,sample]
    expected_abundance <- sum(FPKM_temp)/length(FPKM_temp)
    out_temp <- rbind(out_temp,data.frame(
      id = row.names(spikeData[["SIRV"]][["FPKM"]]), 
      FPKM_norm = FPKM_temp/expected_abundance, 
      sample = rep(sample,length(FPKM_temp)),
      condition = rep( condition[which(samples %in% sample)] , length(FPKM_temp) )))
    spikeData[["SIRV"]][["NORM_FPKM"]] <- out_temp
  }else if(method == "nm2") {
    expected_abundance <- spikeData$SIRV$abundance[ind]/length(SIRV_genes)
    spikeData$SIRV$FPKM_prop[ind] <- expected_abundance/(1/length(ind))
  }
}

# PLOTTING THE RESULT INTO BOXPLOTS #

oldw <- getOption("warn")
options(warn = -1)

print (paste0("Creating boxplots..."))

width_box <- .5
size = 1

plot <- ggplot(spikeData[["SIRV"]][["NORM_FPKM"]], aes(x=sample,y=FPKM_norm))
plot <- plot + stat_boxplot(size = size, geom = "errorbar", width = width_box)
plot <- plot + theme_update(panel.grid.major.x = element_blank()) + theme_minimal()
plot <- plot + ggtitle(paste0("Concentration evaluation - ","SIRV"))

plot <- plot + theme(panel.background = element_rect(fill="white",color="black"),
                     panel.grid.major = element_line(color="gray",size = .3),
                     panel.grid.minor = element_line(color="gray",size = .3))

plot <- plot + theme(text = element_text(size = 30),
                     plot.margin = unit(c(.5,.5,.5,.5),"cm"),
                     plot.title = element_text(hjust = .5, margin = margin(r=0,l=0,t=0,b=0.5,unit="cm")),
                     axis.text.y = element_text(margin = margin(r=0.1,l=0,t=0,b=0,unit="cm")),
                     axis.title.y = element_text(margin = margin(r=0.5,l=0,t=0,b=0,unit="cm")),
                     axis.text.x = element_text(margin = margin(r=0,l=0.5,t=0,b=0,unit="cm")),
                     axis.title.x = element_text(margin = margin(r=0,l=0,t=0.5,b=0,unit="cm")))

plot <- plot + geom_boxplot(color = "black",fill="#63ff86",outlier.shape = NA, width = width_box) 
plot <- plot + geom_jitter(width = 0.1, color = "gray35", alpha = .5)
plot <- plot + stat_summary(geom = "crossbar", width=0.78, fatten= 2, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })
plot <- plot + stat_summary(geom = "crossbar", width=0.50, fatten= 1, color="black", fun.data = function(x){ return(c(y = median(x) - .005, ymin = median(x) - .005, ymax = median(x) - .005)) })
plot <- plot + stat_summary(geom = "crossbar", width=0.50, fatten= 1, color="black", fun.data = function(x){ return(c(y = median(x) + .005, ymin = median(x) + .005, ymax=median(x) + .005 )) })
plot <- plot + scale_y_continuous(trans = "log10",limits = c(10^(-1.666),10^(1.17)),
                                  labels = format_format(big.mark = " ", decimal.mark = ",", scientific = F),
                                  minor_breaks = 10^(seq(-1.75,1.9,.5)),
                                  breaks = round(10^(seq(-2,2,0.5)),2))
plot <- plot + ylab("normalized abundance") 
plot <- plot + stat_summary(fun.y = mean, colour = "darkred", geom = "point", size = 5, shape = 21, stroke = 2, alpha = .8)
plot <- plot + geom_hline(yintercept = 1, size = 1, color = "darkblue", alpha = .5)
plot <- plot + coord_flip() 

if (is.null(opt$aver)) {
  plot <- plot + facet_wrap(~condition, scales = "free", ncol = 1)
}

# SAVING CREATED PLOT
print (paste0("Saving graphics: ", opt[["output"]]))

ggsave(plot = plot,
       filename = opt$output,
       width = 40,
       height = 12+length(samples),
       units = "cm",
       dpi = 250, 
       scale = 1)

options(warn = oldw)
