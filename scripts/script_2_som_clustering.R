## CLUSTERING ANALYSIS
## First, after running script_1_dataprep.R, we can conduct clustering analysis and perform clustering model selection
## for all datasets
## Loading necessary libraries
library(cluster)
library(fpc)
library(factoextra)
library(kohonen)
library(clustertend)
library(pracma)
library(gplots)
library(ggplot2)
library(tidyverse)
library(permute)

## Reading in diel data for all datasets
setwd('intermediate_data/')
lipid_signals<-read.csv('../intermediate_data/lipid_signals.csv',row.names=1)
metabolite_signals<-read.csv('../intermediate_data/metabolite_signals.csv',row.names=1)
transcript_signals<-read.csv('../intermediate_data/full_diel_transcripts.csv',row.names=1)


## Aligning all time series based off of which measurement times
## have no gaps.
## The first lipid measurement is at 0200
## the first metabolite measurement is at 0600
## and the first transcript measurement is at 2200
## let's resolve dates
## first metabolite measurement is 0600 on 7/26
## first lipid measurement is 0200 on 7/27
## first transcript measurement is 2200 on 7/26
## So for everything to line up and be happy, we need to start at lipid measurement 1
## which is metabolite measurement 6
## and gene measurement 2
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
met_subset<-metabolite_signals[,6:20]
colnames(met_subset)<-new_timecodes
lip_subset<-lipid_signals[,1:15]
colnames(lip_subset)<-new_timecodes
tran_subset<-transcript_signals[,2:16]
colnames(tran_subset)<-new_timecodes
lined_up_full_mat<-rbind(met_subset,
                         lip_subset,
                         tran_subset)
keeper_names<-rownames(lined_up_full_mat)
keeper_names<-gsub('E$','',keeper_names)
write.csv(lined_up_full_mat,'../intermediate_data/lined_up_full_mat.csv',quote=FALSE,row.names=TRUE)

## Cluster validation statistics -- hopkins statistic identifies the degree to which
## modularity in the distance matrix exceeds that of uniformly distributed distances
## Warning: because large matrix this is very slow
full_h<-clustertend::hopkins(lined_up_full_mat,1000) ## value is 0.781

## Preliminarily calculate distance matrix to ease calculation
full_dist<-dist(lined_up_full_mat)

## Calculate associated metrics fro pam clustering algorithm
pam_asws<-c()
pam_chs<-c()
for(i in 1:8){
  pam_clust<-clara(lined_up_full_mat,i+1,metric="euclidean")
  pam_asws[i]<-summary(silhouette(pam_clust,full=TRUE))$avg.width
  pam_chs[i]<-calinhara(lined_up_full_mat,pam_clust$clustering)
  print(paste('completed cluster',i))
}

## Calculate associated metrics for hc clustering algorithm
ch_metrics<-c()
asw_metrics<-c()
for(i in 1:8){
  clust<-hcut(full_dist,k=i+1,isdiss=TRUE)
  ch_metrics[i]<-calinhara(lined_up_full_mat,clust$cluster)
  asw_metrics[i]<-summary(silhouette(clust$cluster,full_dist))$avg.width
  print(paste0('completed cluster ',i))
}

## Calculate some SOMs and their associated metrics
x_grid_dim<-c(2,2,3,3,1,7,2,3,5,11,3,6)
y_grid_dim<-c(1,2,1,2,5,1,4,3,2,1,4,2)
ch_val<-c()
asw<-c()
for(i in 1:length(x_grid_dim)){
  set.seed(7583176)
  test_som<-som(as.matrix(lined_up_full_mat),grid=somgrid(x_grid_dim[i],y_grid_dim[i],'rectangular'))
  ch_val[i]<-calinhara(lined_up_full_mat,test_som$unit.classif)
  asw[i]<-summary(silhouette(test_som$unit.classif,full_dist))$avg.width
  print(paste0('completed grid ',i))
}
som_nclust<-x_grid_dim*y_grid_dim
som_ch<-rep(NA,12)
som_asw<-rep(NA,12)
som_ch[som_nclust]<-ch_val
som_asw[som_nclust]<-asw

## Visualizing preliminary ODI for all data to get sense of data structure (WARNING: THIS IS VERY SLOW)
## ODI: Ordered dissimilarity image, representation of the distance matrix with rows/columns sorted to 
## maximize modularity
auto_odi<-get_clust_tendency(lined_up_full_mat,10,
                             gradient=list(low="Darkblue",mid="White",high="Gold"))

## Based off of manually inspecting these, it appears the data could roughly correspond to
## either four or five partially overlapping clusters. We can compare the metrics for these
## amongst algorithms

## So let's compare amongst algorithms + metrics
algs<-rep(c('pam','hc','som'),each=8)
chs<-c(pam_chs,ch_metrics[1:8],som_ch[2:9])
asws<-c(pam_asws,asw_metrics[1:8],som_asw[2:9])
metric_frame<-data.frame(alg=algs,ch=chs,asw=asws,nclust=rep(2:9,3))
ch_plot<-ggplot(metric_frame,aes(x=nclust,col=alg,y=ch))+
  geom_point(size=4)+
  xlab('Number Clusters')+ylab('C-H Value')+theme_bw()+
  scale_color_manual(name='Algorithm',labels=c('HC','PAM','SOM'),
                     values=c('gold','red','navy'))+
  theme(text=element_text(size=16,face='bold'))

ggsave(plot=ch_plot,filename='../figures/ch_metrics.pdf',device='pdf')


asw_plot<-ggplot(metric_frame,aes(x=nclust,col=alg,y=asw))+geom_point(size=4)+
  xlab('Number Clusters')+ylab('Avg Silhouette Width Value')+theme_bw()+
  scale_color_manual(name='Algorithm',labels=c('HC','PAM','SOM'),
                     values=c('gold','red','navy'))+
  theme(text=element_text(size=16,face='bold'))+
  ggrepel::geom_label_repel(data=data.frame(),aes(x=4,y=0.18,label='Elbow'),
                            col='navy',nudge_y=0.1,nudge_x=0.25)
ggsave('../figures/avg_silhouette_metric.pdf',plot=asw_plot,device='pdf')

## Under the decision to compare between 4 and 5, an SOM with 4 clusters appears to be the
## most supported
## Recapitulating 3 cluster SOM for SI
set.seed(1674)
som3<-som(as.matrix(lined_up_full_mat),grid=somgrid(3,1,'hexagonal'))

## Recapitulating that clustering:
set.seed(7583176)
som4<-som(as.matrix(lined_up_full_mat),grid=somgrid(2,2,'hexagonal'))

## Making a 5 cluster SOM for SI
set.seed(37111)
som5<-som(as.matrix(lined_up_full_mat),grid=somgrid(5,1,'hexagonal'))
## Silhouette profile for comparison
somsil5<-silhouette(dist=full_dist,x=som5$unit.classif)
somsil3<-silhouette(dist=full_dist,x=som3$unit.classif)

## Calculating individual silhouette profile
somsil4<-silhouette(dist=full_dist,x=som4$unit.classif)

## Rearranging distance matrix for generating ODI
mornings<-which(som4$unit.classif==2)
nights<-which(som4$unit.classif==4)
afternoons<-which(som4$unit.classif==1)
dusks<-which(som4$unit.classif==3)

rearranged_dmat<-as.matrix(full_dist)[c(nights,dusks,afternoons,mornings),
                           c(nights,dusks,afternoons,mornings)]

## Doing the same for 5 cluster som
clust_1<-which(som5$unit.classif==4)
clust_2<-which(som5$unit.classif==3)
clust_3<-which(som5$unit.classif==2)
clust_4<-which(som5$unit.classif==1)
clust_5<-which(som5$unit.classif==5)
dmat_5<-as.matrix(full_dist)[c(clust_5,clust_4,clust_3,clust_2,clust_1),
                             c(clust_5,clust_4,clust_3,clust_2,clust_1)]


## Partitioning transcript data taxonomically for plotting purposes
signal_names<-rownames(som4$data[[1]])
t_names<-rownames(transcript_signals)
t_names<-gsub('_.*$','',t_names)
euk_names<-grep('^[a-z].*[0-9]$',t_names)
t_names[euk_names]<-gsub('[0-9].*$','',t_names[euk_names])
unique_t_names<-unique(t_names)
appears<-sapply(unique_t_names,function(x) length(which(t_names==x)))
## Focusing on taxa w/at least 100 diel signals for designations in plots, binning remaining taxa as 'other'
keeper_names<-c('Haptophyta','Bacillariophyta','Crocosphaera','Ochrophyta','Chlorophyta',
                'High-light Pro','SAR11','Dinophyta','Cryptophyta','Rhodophyta',
                'Bicosoecida','Cercozoa','Sarcomastigophora','SAR92')
keeper_names<-keeper_names[order(keeper_names)]
final_names<-rep('other',length(t_names))
final_names[which(t_names %in% keeper_names)]<-t_names[which(t_names %in% keeper_names)]
ko_names<-gsub('^.*_','',rownames(transcript_signals))
descriptors<-c(rep('metabolite',nrow(metabolite_signals)),
               rep('lipid',nrow(lipid_signals)),
               final_names)

## Generating a dataframe for the purpose of plotting
sil_comp_frame<-data.frame(id=descriptors,
                           clus4id=som4$unit.classif,
                           clus4sil=somsil4[,3],
                           full_id=c(rep('metabolite',nrow(metabolite_signals)),
                                     rep('lipid',nrow(lipid_signals)),
                                     t_names))

sil_comp_5<-data.frame(id=descriptors,
                       clus4id=som5$unit.classif,
                       clus4sil=somsil5[,3],
                       full_id=c(rep('metabolite',nrow(metabolite_signals)),
                                 rep('lipid',nrow(lipid_signals)),
                                 t_names))

sil_comp_3<-data.frame(id=descriptors,
                       clus4id=som3$unit.classif,
                       clus4sil=somsil3[,3],
                       full_id=c(rep('metabolite',nrow(metabolite_signals)),
                                 rep('lipid',nrow(lipid_signals)),
                                 t_names))

## Generating silhouette profilesfor each cluster
output_figures<-list()
output_figures_5<-list()
make_figure<-function(frame,i){
  clust1order<-frame %>% filter(clus4id==i)
  clust1order<-clust1order[order(clust1order$clus4sil),]
  output_figure<-ggplot(clust1order,aes(x=1:nrow(clust1order),
                                        y=clus4sil,
                                        col=id))+geom_point()+
    theme_bw()+xlab('Ascending Silhouette Width')+ylab('Silhouette Width')+
    scale_color_discrete(name="Signal Type",
                        labels=c('Metabolite','Lipid','Dinophyta','Other Transcript','Bacillariophyta',
                                 'Bicosoecida','Cercozoa','Chlorophyta','Cryptophyta','Haptophyta','Ochrophyta',
                                 'Rhodophyta','Sarcomastigophora','High-light Prochlorococcus','SAR92','Crocosphaera','SAR11'))+
    scale_size_continuous(guide="none")+geom_hline(yintercept=0)
  return(output_figure)
}
for(i in 1:4){
  output_figures[[i]]<-make_figure(sil_comp_frame,i)
  output_figures_5[[i]]<-make_figure(sil_comp_5,i)
}
output_figures_5[[5]]<-make_figure(sil_comp_5,5)
ggsave(filename='../figures/silhouette_morning.pdf',
       device='pdf',
       plot=output_figures[[2]])
ggsave(filename='../figures/silhouette_night.pdf',
       device='pdf',
       plot=output_figures[[4]])
ggsave(filename='../figures/silhouette_afternoon.pdf',
       device='pdf',
       plot=output_figures[[1]])
ggsave(filename='../figures/silhouette_dusk.pdf',
       device='pdf',
       plot=output_figures[[3]])
ggsave(filename='../figures/silhouette_clust1_5clust.pdf',
       device='pdf',
       plot=output_figures_5[[1]])
ggsave(filename='../figures/silhouette_clust2_5clust.pdf',
       device='pdf',
       plot=output_figures_5[[2]])
ggsave(filename='../figures/silhouette_clust3_5clust.pdf',
       device='pdf',
       plot=output_figures_5[[3]])
ggsave(filename='../figures/silhouette_clust4_5clust.pdf',
       device='pdf',
       plot=output_figures_5[[4]])
ggsave(filename='../figures/silhouette_clust5_5clust.pdf',
       device='pdf',
       plot=output_figures_5[[5]])

## Comparing maximum and minimum silhouette widths for each cluster for 4 vs 5 clusters
som_4_summary_stats<-sil_comp_frame %>% group_by(clus4id) %>% summarize(max_width=max(clus4sil),
                                                   min_width=min(clus4sil),
                                                   tot_negative=sum(clus4sil<0),
                                                   tot_sig=length(clus4sil),
                                                   mean_sil=mean(clus4sil),
                                                   med_sil=median(clus4sil),
                                                   sd_sil=sd(clus4sil),
                                                   tot_neg_normed=tot_negative/4)
som_3_summary_stats<-sil_comp_3 %>% group_by(clus4id) %>% summarize(max_width=max(clus4sil),
                                                                        min_width=min(clus4sil),
                                                                        tot_negative=sum(clus4sil<0),
                                                                        tot_sig=length(clus4sil),
                                                                        mean_sil=mean(clus4sil),
                                                                        med_sil=median(clus4sil),
                                                                        sd_sil=sd(clus4sil),
                                                                    tot_neg_normed=tot_negative/3)
som_5_summary_stats<-sil_comp_5 %>% group_by(clus4id) %>% summarize(max_width=max(clus4sil),
                                                                    min_width=min(clus4sil),
                                                                    tot_negative=sum(clus4sil<0),
                                                                    tot_sig=length(clus4sil),
                                                                    mean_sil=mean(clus4sil),
                                                                    med_sil=median(clus4sil),
                                                                    sd_sil=sd(clus4sil),
                                                                    tot_neg_normed=tot_negative/5)

full_som_summary<-cbind(rbind(som_3_summary_stats,
                              som_4_summary_stats,
                              som_5_summary_stats),SOM_version=c(rep('3 cluster',3),
                                                                 rep('4 cluster',4),
                                                                                     rep('5 cluster',5)))

## Writing to output
write.csv(full_som_summary,'../outputs_and_summaries/data_s9_som_summary.csv',quote=FALSE,row.names=FALSE)

## Summarizing clustering compositions
clust_summary<-sapply(1:4,function(x) summary(sil_comp_frame$full_id[which(sil_comp_frame$clus4id==x)]))
clust_summary<-t(rbind(1:4,clust_summary))
clust_proportion<-apply(clust_summary[,-1],2,function(x) x/sum(x))
## Writing to output
write.csv(clust_proportion,'../outputs_and_summaries/clust_proportion_summary.csv')

## Testing for equal distribution of signals amongst taxa
chisq_output<-chisq.test(clust_summary)

## Generate plotted summary
clust_preframe<-clust_summary
rownames(clust_preframe)<-clust_preframe[,1]
clust_preframe<-t(clust_preframe[,-1])
tax_assigns<-rep(rownames(clust_preframe),4)
tax_assigns[-which(tax_assigns %in% c(keeper_names,'lipid','metabolite'))]<-'other'
clust_assigns<-rep(colnames(clust_preframe),each=nrow(clust_preframe))
long_assigns<-c(clust_preframe[,1],clust_preframe[,2],clust_preframe[,3],clust_preframe[,4])
clust_plot_frame<-data.frame(tax=tax_assigns,cluster=clust_assigns,numsig=long_assigns)
clust_plot_frame$cluster<-factor(clust_plot_frame$cluster,levels=c('1','3','4','2'))

## Writing this output to table as well
write.csv(file='../outputs_and_summaries/cluster_membership_summary.csv',clust_summary)


## Organizing cluster membership in terms of signal type
lip_clust_list<-list()
met_clust_list<-list()
ko_clust_list<-list()
clade_clust_list<-list()
transcript_info_list<-list()
for(i in 1:4){
  lip_clust_list[[i]]<-rownames(lined_up_full_mat)[which(som4$unit.classif==i & sil_comp_frame$id=="lipid")]
  met_clust_list[[i]]<-rownames(lined_up_full_mat)[which(som4$unit.classif==i & sil_comp_frame$id=="metabolite")]
  ko_list_full<-ko_names[which(som4$unit.classif==i & sil_comp_frame$id!="lipid" & sil_comp_frame$id!="metabolite")]
  unique_kos<-unique(ko_list_full)
  ko_clust_list[[i]]<-sapply(unique_kos,function(x) length(which(ko_list_full==x)))
  clade_list_full<-t_names[which(som4$unit.classif==i & sil_comp_frame$id!="lipid" & sil_comp_frame$id!="metabolite")]
  unique_clades<-unique(clade_list_full)
  clade_clust_list[[i]]<-sapply(unique_clades,function(x) length(which(clade_list_full==x)))
  transcript_info_list[[i]]<-rownames(lined_up_full_mat)[which(som4$unit.classif==i & sil_comp_frame$id!="lipid" & sil_comp_frame$id!="metabolite")]
}

## Writing output
output_writer<-function(l,sil_comp_frame,signal_names){
  numbers<-do.call(rbind,lapply(l,length))
  labels<-c(rep('Afternoon',numbers[1]),
            rep('Morning',numbers[2]),
            rep('Dusk',numbers[3]),
            rep('Night',numbers[4]))
  full_names<-do.call(c,l)
  sil_width<-sil_comp_frame$clus4sil[sapply(full_names,function(x) which(signal_names==x))]
  neighbor<-somsil4[sapply(full_names,function(x) which(signal_names==x)),2]
  neighbor[which(neighbor=='3')]<-'Dusk'
  neighbor[which(neighbor=='2')]<-'Morning'
  neighbor[which(neighbor=='1')]<-'Afternoon'
  neighbor[which(neighbor=='4')]<-'Night'
  result<-data.frame(home_cluster=labels,signal=full_names,silhouette_width=sil_width,
                     neighbor_cluster=neighbor)
  return(result)
}
## Writing summary files to output
lip_spreadsheet<-output_writer(lip_clust_list,sil_comp_frame,signal_names)
met_spreadsheet<-output_writer(met_clust_list,sil_comp_frame,signal_names)
transcript_spreadsheet<-output_writer(transcript_info_list,sil_comp_frame,signal_names)
write.csv(lip_spreadsheet,'../outputs_and_summaries/lipid_cluster_details.csv')
write.csv(met_spreadsheet,'../outputs_and_summaries/met_cluster_details.csv')
write.csv(transcript_spreadsheet,'../outputs_and_summaries/tran_cluster_details.csv')


