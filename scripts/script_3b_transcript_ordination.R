## This script performs transcript-specific analyses for diel data including NMDS ordination and pathway distribution analysis
library(vegan)
setwd('../intermediate_data/')
## Reading in diel transcript data
transcript_data<-read.csv('full_diel_transcripts.csv',row.names=1)
## Reading in other data
all_types_data<-read.csv('intermediate_data/lined_up_full_mat.csv',row.names=1)
## Calculating Euclidean distance matrix
full_dmat<-dist(transcript_data)
all_types_dmat<-dist(all_types_data)
## Heads up this line will take FOREVER (like an hour) so I just ran it once and saved the results
## but the solution does appear to converge somewhat nicely (reached similar stress < initial conditions on most runs)
full_nmds<-metaMDS(full_dmat)
all_types_nmds<-metaMDS(all_types_dmat)
write.csv(full_nmds$points,'../intermediate_data/transcript_nmds_coordinates.csv')
write.csv(all_types_nmds$points,'intermediate_data/all_data_nmds_coordinates.csv')

## Reading in nmds coordinates
full_nmds<-read.csv('all_data_nmds_coordinates.csv',row.names=1)
## Additionally convert data to ranks to identify peak rank expression 
rank_mat<-t(apply(transcript_data,1,rank))
rank_mat_all<-t(apply(all_types_data,1,rank))
id_highest_mean_rank<-function(v){
  tenpm<-mean(v[c(1,7,13,19)],na.rm=TRUE)
  twoam<-mean(v[c(2,8,14,20)],na.rm=TRUE)
  sixam<-mean(v[c(3,9,15,21)],na.rm=TRUE)
  tenam<-mean(v[c(4,10,16)],na.rm=TRUE)
  twopm<-mean(v[c(5,11,17)],na.rm=TRUE)
  sixpm<-mean(v[c(6,12,18)],na.rm=TRUE)
  meanranks<-c(tenpm,twoam,sixam,tenam,twopm,sixpm)
  maxmeanrank<-which(meanranks==max(meanranks))
  if(1 %in% maxmeanrank & 6 %in% maxmeanrank){
    output<-6.5
  } else output<-mean(maxmeanrank)
  return(output)
}
mean_ranks<-apply(rank_mat,1,id_highest_mean_rank)
mean_rank_all<-apply(rank_mat_all,1,id_highest_mean_rank)
ko_annotes_all<-gsub('^.*_K','K',rownames(all_types_nmds$points))
ko_annotes_all<-ifelse(ko_annotes_all %in% rownames(met_subset),'Metabolite',
                       ifelse(ko_annotes_all %in% rownames(lip_subset),'Lipid',ko_annotes_all))
ko_annotes<-gsub('^.*_','',rownames(full_nmds))
taxa_names<-gsub('_.*$','',rownames(full_nmds))
taxa_names_all<-gsub('_K.*$','',rownames(all_types_nmds$points)) %>%
  gsub('E$','',.) %>%
  gsub(' [0-9].*$','',.) %>%
  gsub('Q.*$','Q',.)
#taxa_names_all<-ifelse(taxa_names_all %in% rownames(met_subset),'Metabolite',
                       #ifelse(taxa_names_all %in% rownames(lip_subset),'Lipid',taxa_names_all))
all_frame<-data.frame(x=all_types_nmds$points[,1],
                      y=all_types_nmds$points[,2],
                      full_ids=rownames(all_types_nmds$points),
                      taxa=taxa_names_all,kos=ko_annotes_all,time_rank=mean_rank_all) %>%
  mutate(big_class=ifelse(kos!='Metabolite',ifelse(kos!='Lipid','Transcript',as.character(kos)),
                          as.character(kos))) %>%
  mutate(new_tax=ifelse(taxa %in% keeper_tax, as.character(taxa),
                                  ifelse(big_class=='Metabolite','Metabolite',
                                         ifelse(big_class=='Transcript','Other Transcript',
                                                ifelse(taxa %in% c('PQ','SQ','UQ'),'Quinone',
                                                       ifelse(taxa=='TAG','TAG',
                                                              ifelse(taxa %in% c('19prime_hex_fuco_a',
                                                                                 '19prime_but_fuco',
                                                                                 'Fuco',
                                                                                 'Lut',
                                                                                 'Zeax','Dd_Ddc'), 'Carotenoid',
                                                                     ifelse(taxa %in% c('Chl-a1',
                                                                            'Chl-a2','Chl-b1','Chl-b2','Pheophytin_a1'),
                                                                            'Chloropigment',
                                                                            ifelse(taxa %in% c('PG','PE',
                                                                                               'PC','DGTS','DGTS_DGTA','DGTS DGTA'),
                                                                                   'Cell Membrane',
                                                                                   ifelse(taxa %in% c('SQDG','MGDG',
                                                                                                      'DGDG','DGCC'),
                                                                                          'Chloroplast Membrane',as.character(taxa)))
                                                                            ))))))))
ggplot(filter(all_frame,big_class!='Transcript'))+
  geom_point(aes(x=x,y=y,col=factor(floor(time_rank))),size=2.5)+
  ggforce::geom_circle(aes(x0=0,y0=0,r=3))+
  coord_fixed()+
  theme_bw()+
  facet_wrap(~new_tax,nrow=2)+
  scale_color_discrete(name='Peak Time',labels=c('0200','0600','1000','1400','1800','2200'))+
  theme(axis.title=element_blank(),
        legend.position='bottom',
        strip.background=element_blank())
ggplot()+
  geom_point(data=filter(all_frame,big_class=='Metabolite'),
             aes(x=x,y=y,
                 col=factor(floor(time_rank))),
             size=4)+
  ggrepel::geom_label_repel(data=filter(all_frame,big_class=='Metabolite'),
                            aes(x=x,y=y,label=full_ids))+
  coord_fixed()+
  theme_bw()+
  scale_color_discrete(name='Peak Time',labels=c('0200','0600','1000','1400','1800','2200'))+
  theme(axis.title=element_blank(),
        legend.position='bottom')+
  ggtitle('Metabolites')
ggsave(filename='figures/metabolite_ordination.pdf',device='pdf',scale=2)
nmds_frame<-data.frame(x=full_nmds[,1],
                       y=full_nmds[,2],
                       full_ids=rownames(full_nmds),
                       taxa=taxa_names,kos=ko_annotes,time_rank=mean_ranks)
appearances<-sapply(unique(taxa_names),function(x) length(grep(x,taxa_names)))
keeper_tax<-c('Haptophyta','Bacillariophyta','Crocosphaera','Ochrophyta','Chlorophyta',
              'High-light Pro','SAR11','Dinophyta','Cryptophyta','Rhodophyta',
              'Bicosoecida','Cercozoa','Sarcomastigophora','SAR92')
new_tax<-rep('other',nrow(nmds_frame))
new_tax[which(nmds_frame$taxa %in% keeper_tax)]<-as.character(nmds_frame$taxa[which(nmds_frame$taxa %in% keeper_tax)])
new_tax<-gsub('_.*$','',new_tax)
graphing_frame<-cbind(nmds_frame,new_tax)
graphing_frame<-all_frame
rownames(graphing_frame)<-graphing_frame$full_ids

## Adding pathway annotations
functional_annotations<-read.csv('../intermediate_data/pathway_annotations.csv')
rank_abundance<-functional_annotations %>%
  group_by(pathway) %>%
  count()
rank_abundance<-rank_abundance[order(rank_abundance$n,decreasing=TRUE),]
top_paths<-rank_abundance$pathway[3:32]
keeper_path_assigns<-which(functional_annotations$signal %in% rownames(graphing_frame))
org_type<-c()
multiple_annotes<-c()
k<-1
path_assign<-c()
clust_assign<-c()
for(i in 1:nrow(graphing_frame)){
  org_type[i]<-as.character(functional_annotations$taxon_group[which(functional_annotations$signal==rownames(graphing_frame)[i])[1]])
  clust_assign[i]<-as.character(functional_annotations$cluster[which(functional_annotations$signal==rownames(graphing_frame)[i])[1]])
  path_sub<-which(functional_annotations$signal==rownames(graphing_frame)[i] & functional_annotations$pathway %in% top_paths)
  path_sub<-as.character(functional_annotations$pathway[path_sub])
  if(length(path_sub)>1){
    keep_path<-path_sub[sample(1:length(path_sub),1)]
    multiple_annotes[k]<-i
    k<-k+1
  }else if(length(path_sub)==1){
    path_assign[i]<-path_sub
  }else path_assign[i]<-'other'
}
graphing_frame<-data.frame(graphing_frame,path=path_assign,
                           tax_group=org_type,
                           cluster=clust_assign,
                           stringsAsFactors=FALSE)
### Need to assign cluster to metabolites and lipids
lipid_entries<-match(graphing_frame$full_ids,lip_spreadsheet$signal)
met_entries<-match(graphing_frame$full_ids,met_spreadsheet$signal)
graphing_frame$cluster[!(is.na(lipid_entries))]<-as.character(lip_spreadsheet$home_cluster[lipid_entries[which(is.na(lipid_entries)==FALSE)]])
graphing_frame$cluster[!(is.na(met_entries))]<-as.character(met_spreadsheet$home_cluster[met_entries[which(is.na(met_entries)==FALSE)]])
graphing_frame$tax_group[is.na(graphing_frame$tax_group)]<-'Chemical'

## Rotating NMDS ordination so that 0000 hours is at the top of the ordination and 1200 hours is at the bottom
## For just transcripts we reflect about the y axis and rotate 3/48 pi
#graphing_frame_rotated<-graphing_frame %>%
  #mutate(x_rot=(-1*x*cos(3*pi/48))+(y*sin(3*pi/48)),y_rot=y*cos(3*pi/48)+(x*sin(3*pi/48)))
graphing_frame_rotated<-graphing_frame %>%
  mutate(x_rot=(x*cos(pi/16))-(y*sin(pi/16)),y_rot=(y*cos(pi/16)+(x*sin(pi/16))))
write.csv(graphing_frame_rotated,quote=FALSE,row.names=TRUE,'../intermediate_data/rotated_nmds_coordinates_413.csv')

## Summarizing pathway distrbutions across clusters
heatmap_frame<-graphing_frame_rotated %>%
  group_by(tax_group,path,cluster) %>%
  filter(tax_group %in% c('eukaryotic','prokaryotic_photoauto','prokaryotic_non_photoauto')) %>%
  filter(path %in% top_paths) %>%
  count() %>%
  group_by(tax_group,path) %>%
  mutate(percent_count=n/sum(n)) %>%
  filter(tax_group %in% c('eukaryotic','prokaryotic_photoauto','prokaryotic_non_photoauto')) %>%
  filter(path %in% top_paths,n!=0)

heatmap_frame$path<-gsub('^.*[0-9][[:space:]]','',heatmap_frame$path)
pathway_order<-subset(heatmap_frame,tax_group=='eukaryotic')
weighted_path_sum<-c()
for(i in 1:length(unique(pathway_order$path))){
  sub_path<-subset(pathway_order,path==unique(pathway_order$path)[i])
  sub_path$cluster<-factor(sub_path$cluster,levels=c('Night','Dusk','Afternoon','Morning'))
  weighted_path_sum[i]<-pracma::dot(as.numeric(sub_path$cluster),sub_path$percent_count)
}
next_path_order<-unique(pathway_order$path)[order(weighted_path_sum,decreasing=FALSE)]
next_path_order<-c(as.character(next_path_order),setdiff(unique(heatmap_frame$path),next_path_order))
heatmap_frame$path<-factor(heatmap_frame$path,levels=next_path_order)


