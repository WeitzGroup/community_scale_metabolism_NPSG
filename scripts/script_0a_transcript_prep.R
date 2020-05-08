## The purpose of this script is to prepare >5um fraction transcripts for
## RAIN analysis etc. Because there are a lot of CAMNT contigs, this dataset is somewhat large and cumbersome
## to deal with in memory, so this script is meant to show you the steps for how we got from camnt_genus_ko_cleaned.txt to 
## 10_2_transcript_full_diels.csv

## Next is reading initial dataset
genus_data<-fread('../raw_data/camnt_genus_ko_cleaned.txt')
## Adding contig id and taxonomic information to row headers for easier reading down the line
row_names<-paste0(genus_data$V1,'_',genus_data$V2,'_',genus_data$V3)
## pulling tiume series data
numerics<-genus_data[,4:24]
## Finding redundant contigs
camnts<-unique(genus_data$V1)
all_camnts<-genus_data$V1
signal_indices<-c()
## This is very slow
for(i in 1:length(camnts)){
  signal_indices[i]<-which(all_camnts==camnts[i])[1]
}
## Saving output so you never have to run that loop again
write(signal_indices,'../intermediate_data/signal_indices.txt')
## Reading in output to pick up where we left off 
signal_indices<-fread('../intermediate_data/signal_indices.txt')

## Formatting so we can extract just the unique contigs
signal_indices<-c(apply(signal_indices,1,c))

## Getting the timeseries for those contigs
just_contigs<-numerics[signal_indices,]

write.csv(just_contigs,quote=FALSE,row.names=FALSE,'../intermediate_data/big_fraction_transcripts.csv')

