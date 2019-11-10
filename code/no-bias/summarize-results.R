# if results.summarized is not saved correctly ----------------------------
# we need to load all individual files and summarize those.
rm(list = ls())
options(scipen=999)
x = c("data.table")
sapply(x,library,character.only=T)
results.summarized <- c() # empty vector to store results
getwd()
setwd("../idealworld/results")

# Load all individual result files 
#filelist = list.files(pattern = "res_")

# Load all individual result files 
files <- list.files()
res.files <- files[-grep("str", files, fixed=T)]
res.files <- res.files[-grep(".txt", res.files, fixed=T)]
res.files <- res.files[-grep(".out", res.files, fixed=T)]
head(res.files);last(res.files)

# Check which conditions are saved
sel <- substr(res.files, 11, 14) # select only the 11-14th character from res.files

sel <- gsub('_0.', '', sel)
sel <- gsub('_0', '', sel)
sel <- gsub('_', '', sel)
sel <- as.numeric(sel)
sel <- sort(sel)

# enter number of conditions here
cond <- c(1:1620)

# check if all conditions have been checked
length(cond) == length(sel)
sel %in% cond
which(!sel %in% cond)

# save the conditions that are missing
missing <- which(!sel %in% cond)

# load all the datafiles
datalist = lapply(res.files, fread)

for (i in seq_along(datalist)){
  
  # assign NA to string variables first, otherwise mean cannot be calculated
  # technically this isnt necessary because Ive filtered out the strings beforehand
  datalist[[i]] <- sapply(datalist[[i]], as.numeric)
  
  # take mean of all iterations per condition
  results.summarized.new <- apply(datalist[[i]],2,mean,na.rm=T)
  
  # assign it to summarized results
  results.summarized <- rbind(results.summarized,results.summarized.new)
  
}

# assign column names to results.summarized and sort by cond
colnames(results.summarized) <- colnames(datalist[[1]]) 
results.summarized <- results.summarized[order(results.summarized[,1]),]

# write summarized results to file
write.table(results.summarized, file = "../idealworld_summarized_results.dat", row.names=F, col.names=T)
