### SCRIPT ONE
### GOALS: LOAD ALL DATASETS, PERFORM RAIN ANALYSIS, SUBSET DIEL SIGNALS FROM ALL DATASETS, LOAD FUNCTIONAL DATA FOR TRANSCRIPTS
## Start by importing libraries
library(pracma)
library(rain)
library(tidyverse)
library(vegan)
library(gplots)
library(kohonen)
library(data.table)
library(cluster)
library(BioPhysConnectoR)

## Establishing Working directory
setwd('raw_data')

### We'll organize this one dataset at a time bc each one has different processing procedures
### Dataset 1: Metabolomics Data GOAL 1: PREPROCESSING + READING DATA
raw_file<-read.csv('More_Columns_DielTargetedSamples_KM1513-trimed-2018-06-08.csv',header=TRUE)
## putting together the replicate and metabolite for each measurement
time_compound<-paste0(raw_file$replicate,',',raw_file$Compound.Name)
## pulling sample times
all_sample_times<-unique(raw_file$day.time)
max_samples<-length(all_sample_times)
total_samples<-unique(time_compound)
## Iterating to find measurements for all samples at each given time chopped into replicates
indi_series<-list(NA)
for(i in 1:length(total_samples)){
  indices<-which(time_compound==total_samples[i])
  times<-raw_file$day.time[indices]
  value<-raw_file$NormalizedPkArea[indices]
  scheme<-match(times,all_sample_times)
  finalvalue<-rep(NA,max_samples)
  finalvalue[scheme]<-value
  indi_series[[i]]<-matrix(data=finalvalue,nrow=1)
  colnames(indi_series[[i]])<-all_sample_times
  rownames(indi_series[[i]])<-gsub('.,','',total_samples[i])
}
## Turning the output into a matrix
series_mat<-do.call(rbind,indi_series)
time_order<-order(as.numeric(colnames(series_mat)))
sorted_mat<-series_mat[,time_order]
## Writing to a file
rownames(sorted_mat)<-paste0(rownames(sorted_mat),'_',rep(c(1,2,3),nrow(sorted_mat)/3))
transposed<-t(sorted_mat)
cols<-colnames(transposed)
## Finding unique metabolite names
mets<-gsub('_.','',cols)
mets_unique<-unique(mets)
## Performing detrending step using replicates
final_col<-c()
for(i in 1:length(mets_unique)){
  ## Grab the columns for that metabolite
  met_sub<-transposed[,which(mets==mets_unique[i])]
  ## Using first cruise half for consistency w/analysis
  first_half<-met_sub[1:20,]
  first_half_collapsed<-as.numeric(c(t(first_half)))
  time_points<-rep(seq(0,76,by=4),each=3)
  ## Detrending as linear regression over time
  regression<-lm(first_half_collapsed~time_points)
  stretched_values<-rep(1,length(first_half_collapsed))
  ## Accounting for missing measurements
  stretched_values[which(is.na(first_half_collapsed)==TRUE)]<-NA
  stretched_values[which(is.na(stretched_values)==FALSE)]<-regression$fitted.values
  ## Performing detrending
  detrended_data<-first_half_collapsed-stretched_values
  detrended_data<-cbind(detrended_data)
  ## Output detrended metabolite intensities
  final_col<-cbind(final_col,detrended_data)
}
colnames(final_col)<-mets_unique
## DATASET 1: METABOLOMICS GOAL 2: PERFORM RAIN ANALYSIS
met_rain_output<-rain(final_col,deltat=4,period=24,nr.series=3,na.rm=TRUE)
met_rain_ps<-met_rain_output$pVal[order(met_rain_output$pVal)]
bh_met<-1:length(met_rain_ps)*0.05/length(met_rain_ps)
bh_pass<-max(which(met_rain_ps<=bh_met))
met_sig_subset<-met_rain_output[which(met_rain_output$pVal %in% met_rain_ps[1:bh_pass]),]
## Calculating the averages for each metabolite at each time point
averages<-matrix(NA,nrow=nrow(transposed),ncol=length(mets_unique))
for(i in 1:length(mets_unique)){
  series<-which(mets==mets_unique[i])
  series<-transposed[,series]
  series<-apply(series,1,as.numeric)
  series<-colMeans(series,na.rm=TRUE)
  averages[,i]<-series
}
#Reformatting for output with correct row and column names
time_titles<-as.numeric(rownames(transposed))
format_avg<-cbind(time_titles,averages)
format_avg<-as.matrix(format_avg)
colnames(format_avg)<-c('time_titles',mets_unique)
format_avg_sig<-format_avg[,rownames(met_sig_subset)]
met_first_cruise<-format_avg_sig[1:20,]
metabolite_signals<-t(scale(detrend(met_first_cruise)))


### Writing intermediate data to outputs
write.csv(metabolite_signals,'../intermediate_data/metabolite_signals.csv',quote=FALSE)
write.csv(data.frame(met_rain_output[order(met_rain_output$pVal),],
                     diel=c(rep('YES',bh_pass),rep('NO',nrow(met_rain_output)-bh_pass))),
          '../outputs_and_summaries/metabolite_rain_out.csv',quote=FALSE)

## DATASET 2: LIPIDOMICS 
lipid_data<-read.delim('KBecker_VanMooyLab_Lipids.csv',sep=';')
##Conducting RAIN analysis on lipid data based on first cruise half
lip_subset<-lipid_data[2:nrow(lipid_data),4:23]
lip_subset<-t(lip_subset)
lip_subset<-t(apply(lip_subset,1,as.numeric))
colnames(lip_subset)<-lipid_data[2:nrow(lipid_data),1]
measurement_times<-c('07_27_0200','07_27_0600','07_27_1000','07_27_1400','07_27_1800','07_27_2200',
                     '07_28_0200','07_28_0600','07_28_1000','07_28_1400','07_28_1800','07_28_2200',
                     '07_29_0200','07_29_0600','07_29_1000','07_29_1400','07_29_1800','07_29_2200',
                     '07_30_0200','07_30_0600')
lip_detrend<-detrend(lip_subset)
##Doing RAIN
lip_rain_results<-rain(lip_detrend,deltat=4,period=24)
lip_pvals_order<-order(lip_rain_results$pVal)
lip_rain_sorted<-lip_rain_results[lip_pvals_order,]
##FDR Control
lip_q<-0.05
lip_is<-1:nrow(lip_rain_sorted)
lip_abh_metric<-lip_is*lip_q/(nrow(lip_rain_sorted))
lip_keeper_rank<-max(which(lip_rain_sorted$pVal<=lip_abh_metric))
## In case we want to look at this for supplemental purposes
plot(lip_is,lip_rain_sorted$pVal,col='red',type='p',
     xlab='Lipid by p-value rank',
     ylab='p-value',
     main='Adaptive FDR Adjustment for Lipid RAIN Analysis',
     ylim=c(0,0.1))
lines(lip_is,lip_abh_metric,lwd=2)
legend(0,0.1,fill=c('red','black'),legend=c('Analytic p-val','ABH Cutoff'))
##Retaining series which pass ABH criterion
lip_diels<-lip_subset[,rownames(lip_rain_sorted)[1:lip_keeper_rank]]
lipid_signals<-t(scale(detrend(lip_diels)))

## Writing intermediate data to output
write.csv(lipid_signals,'../intermediate_data/lipid_signals.csv',quote=FALSE)
write.csv(t(lip_diels),'../intermediate_data/diel_lipid_units.csv',quote=FALSE)
write.csv(data.frame(lip_rain_sorted,
                     diel=c(rep('YES',lip_keeper_rank),rep('NO',nrow(lip_rain_sorted)-lip_keeper_rank))),
          '../outputs_and_summaries/lipid_rain_out.csv',quote=FALSE)

## DATASET 3: <5UM TRANSCRIPTS
by_clade<-read.delim('transcripts_table_ko_byClade.txt',header=TRUE)
## Subsetting for first half of cruise
## Separating first cruise half
by_clade_first_half<-by_clade[,3:27]
rownames(by_clade_first_half)<-paste0(by_clade[,1],"_",by_clade[,2])
## Detrending
clade_dt<-detrend(t(by_clade_first_half))

## Performing rain analysis
clade_rain_results<-rain(clade_dt,deltat=4,period=24)

## BH FDR Control
clade_rain_p_order<-clade_rain_results[order(clade_rain_results$pVal),]
q_val<-0.05
abh_metric<-1:nrow(clade_rain_p_order)*(q_val/nrow(clade_rain_p_order))
max_keeper<-max(which(clade_rain_p_order$pVal<=abh_metric))


## Gathering results for clustering:
keeper_ids<-rownames(clade_rain_p_order)[1:max_keeper]
## Subsetting data
clade_diels<-by_clade_first_half[keeper_ids,]

## writing to output
write.csv(clade_diels,'../intermediate_data/clade_diel_data.csv')
write.csv(data.frame(clade_rain_p_order,
                     diel=c(rep('YES',max_keeper),rep('NO',nrow(clade_rain_p_order)-max_keeper))),
          '../outputs_and_summaries/small_transcripts_rain_out.csv',quote=FALSE)

## recording the proportion diel signals by taxon
total_clade_count<-table(gsub('_.*$','',rownames(clade_rain_results)))
diel_clade_count<-table(gsub('_.*$','',rownames(clade_diels)))
clade_prop_diel<-diel_clade_count/total_clade_count
diel_count_frame<-data.frame(clade=names(total_clade_count),
                             number_kos_analyzed=as.numeric(total_clade_count),
                             number_diel_kos=as.numeric(diel_clade_count),
                             proportion_kos_diel=as.numeric(clade_prop_diel),
                             fraction='<5um')

## Now >5um Transcripts
## Because some of the processes here are memory intensive and therefore slow, I will leave
## some of the preprocessing steps in a separate script script_0a_transcript_prep.R but also leave intermediate files
## preloaded in the directory so you don't have to run this step unless you want to 
just_contigs<-read.csv('../intermediate_data/big_fraction_transcripts.csv')
genus_data<-fread('../raw_data/camnt_genus_ko_cleaned.txt')
## Getting necessary information from the full contig table
all_camnts<-genus_data$V1
## Time series intensities
numerics<-genus_data[,4:24]


## Organizing contig information into something useful for clustering analyses:
## Getting taxonomic assignment information for contigs + resolving some spelling inconsistencies
all_taxa<-read.delim('../raw_data/mmetsp_idlist.txt',stringsAsFactors = FALSE)
all_taxa_real<-all_taxa[-which(all_taxa$PHYLUM=='Unknown' & all_taxa$GENUS=='Unknown'),] #removing no annotation
all_taxa_real$GENUS[which(all_taxa_real$GENUS=="Chrysochromulina")]<-"Chyrsochromulina" #resolving spelling mistake
all_taxa_real$GENUS[which(all_taxa_real$GENUS=="Aristerostoma")]<-"Aristerstoma" #resolving spelling
all_taxa_real$PHYLUM[which(all_taxa_real$PHYLUM=="Forminafera")]<-"Foraminifera" #resolving spelling
## Simplifying taxonomic collapse
all_taxa_real$PHYLUM[which(all_taxa_real$PHYLUM=="Pyrrophycophyta")]<-"Dinophyta"

## Making secondary lookup table for dealing with ambiguous taxonomic annotations
mmetsp_to_cam<-read.delim('../raw_data/amb_mmetsp.txt',stringsAsFactors=FALSE,sep=' ',header=FALSE)
## Identifying which contigs in our analysis have ambiguous annotations
unknowns_plus_mmetsp<-which(camnts %in% mmetsp_to_cam$V1)

## Identifying all organism annotations
unique_mmetsps<-unique(mmetsp_to_cam$V2)
## Identifying which taxonomic annotations are in the taxonomies for our contigs
new_additions<-all_taxa_real$PHYLUM[which(all_taxa_real$mmetspID %in% unique_mmetsps)]
## Adding in sarcomastigophora because it didn't make it into our taxonomy guide 
new_additions[4]<-'Sarcomastigophora'
new_additions<-cbind(new_additions,all_taxa_real$mmetspID[which(all_taxa_real$mmetspID %in% unique_mmetsps)])
## Getting genus level annotations
genus_names<-genus_data$V2
## Finding all ambiguous identifies
ambiguous<-which(genus_names %in% c('unknown',
                                    'unid',
                                    'non',
                                    'CAMNT',
                                    'unident.',
                                    'unid.',
                                    'unident.',
                                    'Unid.'))
## Identifying all contig annotations which are still ambiguous
ambiguous_camnts<-unique(all_camnts[ambiguous])

## Matching taxonomic assignment with KO putative function
ko_names<-genus_data$V3

## Now we generate a function which identifies all contigs with the same taxonomic and KO assignments and
## adds them together for a 'total phylum KO activity signal'
subset_adder<-function(text){
  # Initializing output
  additional_subset<-c()
  # Finding all the contigs w/this taxonomic annotation
  if(text %in% new_additions[,1]){
    additional_subset<-which(all_camnts %in% mmetsp_to_cam$V1[which(mmetsp_to_cam$V2 %in% new_additions[which(new_additions[,1]==text),2])])
  }
  taxon_subset<-which(genus_names %in% all_taxa_real$GENUS[which(all_taxa_real$PHYLUM==text)])
  taxon_subset<-c(taxon_subset,additional_subset)
  ## Identifying which KO annotations go with those contigs
  taxon_kos<-ko_names[taxon_subset]
  ## Finding out which KOs are represented
  unique_taxon_kos<-unique(taxon_kos)
  ## Initializing final output (taxon_ko summed signal)
  output_phylum<-rep(text,length(unique_taxon_kos))
  output_timeseries<-c()
  ## Adding iteratively per KO
  for(i in 1:length(unique_taxon_kos)){
    rows_to_use<-which(taxon_kos==unique_taxon_kos[i])
    full_rows_to_use<-taxon_subset[rows_to_use]
    summed_timeseries<-colSums(numerics[full_rows_to_use,])
    output_timeseries<-rbind(output_timeseries,summed_timeseries)
  }
  ## Printing a progress message
  print(paste('Finished long part for ',text))
  ## Carrying over rownames to output
  rownames(output_timeseries)<-paste0(output_phylum,'_',unique_taxon_kos)
  return(output_timeseries)
}

## Subsetting for included phyla:
phy_to_include<-c("Dinophyta",
                  "Chromerida",
                  "Bacillariophyta",
                  "Bicosoecida",
                  "Cercozoa",
                  "Chlorarachniophyta",
                  "Chlorophyta",
                  "Cryptophyta",
                  "Euglenozoa",
                  "Glaucophyta",
                  "Haptophyta",
                  "Ochrophyta",
                  "Rhodophyta",
                  "Sarcomastigophora")

## Implementing and reformatting for writing to output
full_result<-sapply(phy_to_include,subset_adder)
full_result<-do.call(rbind,full_result) 
write.csv(full_result,'../intermediate_data/big_transcript_signals.csv')
## Detrending for rain analysis
fr_dt<-detrend(t(full_result))
## Implementing rain analysis and BH FDR control 
phylum_rain_output<-rain(fr_dt,deltat=4,period=24)
phylum_p_order<-order(phylum_rain_output$pVal,decreasing=FALSE)
q_thresh<-0.05
is<-1:length(phylum_p_order)
abh_metric<-is*q_thresh/length(phylum_p_order)
max_keeper<-max(which(phylum_rain_output$pVal[phylum_p_order]<=abh_metric))
sig_subset<-phylum_rain_output[phylum_p_order[1:max_keeper],]


## Saving output
new_phylum_diels<-full_result[which(rownames(full_result) %in% rownames(sig_subset)),]
write.csv(new_phylum_diels,'../intermediate_data/diel_big_transcripts.csv')
write.csv(data.frame(phylum_rain_output[phylum_p_order,],
                     diel=c(rep('YES',max_keeper),rep('NO',length(phylum_p_order)-max_keeper))),
          '../outputs_and_summaries/big_transcripts_rain_out.csv',quote=FALSE)

## Writing proportion kos diel
total_phy_count<-table(gsub('_.*$','',rownames(phylum_rain_output)))
diel_phy_count<-table(gsub('_.*$','',rownames(new_phylum_diels)))
phy_prop_diel<-diel_phy_count/total_phy_count
diel_count_frame_big<-data.frame(clade=names(total_phy_count),
                             number_kos_analyzed=as.numeric(total_phy_count),
                             number_diel_kos=as.numeric(diel_phy_count),
                             proportion_kos_diel=as.numeric(phy_prop_diel),
                             fraction='>5um')
## Saving to supp file
write.csv(rbind(diel_count_frame,
                diel_count_frame_big),
          '../outputs_and_summaries/taxa_proportions_diel.csv',quote=FALSE)

## Combining transcript output for further analyses
full_phylum_dtsc<-t(scale(detrend(t(new_phylum_diels))))
small_tran_dtsc<-t(scale(detrend(t(clade_diels))))
full_data<-rbind(full_phylum_dtsc,small_tran_dtsc[,5:25])
write.csv(full_data,'../intermediate_data/full_diel_transcripts.csv',quote=FALSE)

## That's everything for this preprocessing so now we can move into clustering analysis, transcript-specific analyses, etc.
## depending on what I end up calling them all.

## Now synthesizing rescaled data for F1 

## large transcripts
taxonomic_appearances_euk<-sort(table(gsub('_.*$','',rownames(sig_subset))),decreasing=TRUE)
big_taxonomic_groups_euk<-names(taxonomic_appearances_euk[1:5])
top_hits_per_euk<-do.call(c,lapply(big_taxonomic_groups_euk,function(x) rownames(sig_subset)[grep(x,rownames(sig_subset))[1:3]]))
euk_top10_values<-fr_dt[,top_hits_per_euk]
euk_top10_rescaled<-apply(euk_top10_values,2,function(x) (x-min(x))/(max(x)-min(x)))
## small transcripts
taxonomic_appearances<-sort(table(gsub('_.*$','',keeper_ids)),decreasing=TRUE)
big_taxonomic_groups<-names(taxonomic_appearances[1:5])
top_hits_per_taxon<-do.call(c,lapply(big_taxonomic_groups,function(x) keeper_ids[grep(x,keeper_ids)[1:3]]))
clade_top10_values<-detrend(t(by_clade_first_half[top_hits_per_taxon,]))
clade_top10_rescaled<-apply(clade_top10_values,2,function(x) (x-min(x))/(max(x)-min(x)))
## lipids
lip_p_top10<-rownames(lip_rain_results[order(lip_rain_results$pVal)[1:15],])
lip_top10_values<-lip_detrend[,lip_p_top10]
lip_top10_rescaled<-apply(lip_top10_values,2,function(x) (x-min(x))/(max(x)-min(x)))
## metabolites
met_p_top10<-rownames(met_rain_output[order(met_rain_output$pVal)[1:15],])
met_top10_values<-detrend(met_first_cruise[,met_p_top10])
met_top10_rescaled<-apply(met_top10_values,2,function(x) (x-min(x))/(max(x)-min(x)))

new_timecodes<-c('0200_727',
                 '0600_727',
                 '1000_727',
                 '1400_727',
                 '1800_727',
                 '2200_727',
                 '0200_728',
                 '0600_728',
                 '1000_728',
                 '1400_728',
                 '1800_728',
                 '2200_728',
                 '0200_729',
                 '0600_729',
                 '1000_729')
top_mets<-t(met_top10_rescaled[6:20,])
colnames(top_mets)<-new_timecodes
top_lips<-t(lip_top10_rescaled[1:15,])
colnames(top_lips)<-new_timecodes
top_small_tran<-t(clade_top10_rescaled[5:25,])
top_big_tran<-t(euk_top10_rescaled)
top_alltran<-rbind(top_small_tran,top_big_tran)
all_tran_subset<-top_alltran[,2:16]
colnames(all_tran_subset)<-new_timecodes
## function we will use for ordering
id_highest_mean_rank_intro<-function(v){
  tenpm<-mean(v[c(6,12)])
  twoam<-mean(v[c(1,7,13)])
  sixam<-mean(v[c(2,8,14)])
  tenam<-mean(v[c(3,9,15)])
  twopm<-mean(v[c(4,10)])
  sixpm<-mean(v[c(5,11)])
  meanranks<-c(twoam,sixam,tenam,twopm,sixpm,tenpm)
  maxmeanrank<-which(meanranks==max(meanranks))
  if(1 %in% maxmeanrank & 6 %in% maxmeanrank){
    output<-6
  } else output<-min(maxmeanrank)
  return(output)
}
f1_ready_mat<-data.frame(rbind(top_mets,
                         top_lips,
                         all_tran_subset)) %>%
  rownames_to_column(var='signal_name') %>%
  mutate(signal_type=rep(c('M','L','P','E'),c(15,15,15,15))) %>%
  gather(key='sample_time',value='intensity',contains('X')) %>%
  mutate(sample_time=factor(sample_time,levels=paste0('X',new_timecodes)),
         sample_number=as.numeric(sample_time),
         mod_time=(sample_number%%6)+1) %>%
  group_by(signal_name) %>%
  mutate(int_rank=order(intensity,decreasing=TRUE),
         weighted_rank_time=weighted.mean(mod_time,int_rank/15),
         max_time=mod_time[which(int_rank==15)],
         new_rank_time=1/id_highest_mean_rank_intro(rank(intensity))) %>%
  group_by(new_rank_time) %>%
  arrange(signal_name,.by_group=TRUE)

adj_factor<-rep(cumsum(rep(1,60))-1,each=15)
f1_ready_mat<-data.frame(f1_ready_mat,adjusted_int=f1_ready_mat$intensity+adj_factor)

f1_panel<-ggplot(f1_ready_mat)+
  geom_rect(data=data.frame(xmin=c(1,5,11),xmax=c(2,8,14),ymin=rep(-Inf,3),ymax=rep(Inf,3)),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            alpha=0.2,col='darkgrey')+
  geom_line(aes(x=sample_number,y=adjusted_int,group=signal_name))+
  theme_bw()+
  theme(legend.position='none',
        axis.text.y=element_text(face='bold'),
        #axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        #axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,vjust=0.5,size=14),
        axis.ticks.length.x=unit(0.5,'cm'))+
  scale_x_continuous(breaks=c(2,5,8,11,14),
                     labels=c('0600 7/27',
                              '1800 7/27',
                              '0600 7/28',
                              '1800 7/28',
                              '0600 7/29'))+
  scale_y_continuous(position='right',breaks=seq(0.5,59.5,by=1),
                     labels=gsub('High-light Pro','HL Prochlorococcus',
                       gsub('E$','',gsub('_.*$','',
                                 unique(f1_ready_mat$signal_name)))))+
  #scale_y_continuous(breaks=unique(adj_factor)+0.5,
  #                   labels=f1_ready_mat$signal_type[seq(1,nrow(f1_ready_mat),by=15)])+
  #scale_color_discrete(name='Peak Time',
  #                     labels=c('2200','1800','1600','1400','1000','0600','0200'),
  #                    guide=guide_legend(nrow=7))+
  #geom_text(data=data.frame(x=rep(15.5,60),
  #                          y=unique(adj_factor)+0.5,
  #                          label=f1_ready_mat$signal_type[seq(1,nrow(f1_ready_mat),by=15)]),
  #          aes(x=x,y=y,label=label,col=label),size=2.85)+
  #scale_color_manual(values=rev(c('chocolate4','goldenrod3','seagreen')))+
  #geom_rect(data=data.frame(xmins=rep(-1.75,6),
  #                          xmaxs=rep(-0.75,6),
  #                          ymins=head(c(60-cumsum(table(1/f1_ready_mat$new_rank_time))/15,0),-1),
  #                          ymaxs=head(c(60,60-cumsum(table(1/f1_ready_mat$new_rank_time))/15),-1),
  #                          time=c('0200','0600','1000','1600','1800','2200')),
  #          aes(xmin=xmins,xmax=xmaxs,ymin=ymins,ymax=ymaxs),
  #          fill=c('blue4','navy','royalblue3','lightblue1',
  #                 'khaki1','gold','darkorange','red','salmon',
  #                 'maroon','magenta4','dodgerblue4')[c(3,5,7,9,11,1)],
  #          color='black')+
  geom_label(data=data.frame(x=rep(-0.15,60),
                            y=unique(adj_factor)+0.5,
                            label=f1_ready_mat$signal_type[seq(1,nrow(f1_ready_mat),by=15)]),
                            aes(x=x,y=y,fill=label,label=label),
             color='white',label.r=unit(0,'pt'),
             label.size=0,size=2.75,label.padding=unit(0.05,'lines'),
             fontface='bold')+
  geom_vline(xintercept=0.5,color='grey80')+
  #geom_segment(data=data.frame(x=rep(0,length(unique(adj_factor))),
  #                             xend=rep(0.5,length(unique(adj_factor))),
  #                             yend=unique(adj_factor)+0.5),
  #             aes(x=x,xend=xend,y=yend,yend=yend),color='grey80')+
  scale_fill_manual(values=c('red2','chocolate4','palegreen4','steelblue3'))+
  coord_fixed(ratio=1.1)+
  geom_point(aes(x=sample_number,y=adjusted_int),size=0.35)
  #geom_point(aes(x=1/new_rank_time,y=adjusted_int),col='red')+
  #geom_point(aes(x=(1/new_rank_time)+6,y=adjusted_int),col='red')+
  #geom_point(aes(x=(1/new_rank_time)+12,y=adjusted_int),col='red')

ggsave(f1_panel,filename='test_colorbar_and_lines.pdf',device='pdf',height=8.5,width=11,units='in')
