## Getting photosynthesis and TCA cycle genes
photo_gene_kos<-as.character(unique(filter(graphing_frame_rotated,path=='Photosynthesis')$kos))
photo_genes_scaled<-lined_up_full_mat[(gsub('^.*_','',rownames(lined_up_full_mat)) %in% photo_gene_kos),]
tca_gene_kos<-as.character(unique(filter(graphing_frame_rotated,path=='Citrate cycle (TCA cycle)')$kos))
tca_genes_scaled<-lined_up_full_mat[(gsub('^.*_','',rownames(lined_up_full_mat)) %in% tca_gene_kos),]
sugar_transporters<-as.character(unique(filter(new_transporters,General_category=='sugar')$kos))
sugar_transporters_scaled<-lined_up_full_mat[(gsub('^.*_','',rownames(lined_up_full_mat)) %in% sugar_transporters),]
heterotroph_names<-setdiff(unique(gsub('_.*$','',clade_diels$X)),c('High-light Pro','Low-light Pro',
                                                                   'Crocosphaera','Viruses','Eukaryota',
                                                                   'Synechococcus','Cyanothece','Euryarchaeota','--'))
full_c_cycle_frame<-data.frame(rbind(metabolite_scaled,
                                     lipid_scaled,
                                     photo_genes_scaled,
                                     tca_genes_scaled,
                                     sugar_transporters_scaled,
                                     as.numeric(optics_data_matched)),type=rep(c('Metabolite','TAG','Photosynthesis','TCA Cycle','Sugar Transporters','POC'),
                                                                               c(nrow(metabolite_scaled),
                                                                                 nrow(lipid_scaled),
                                                                                 nrow(photo_genes_scaled),
                                                                                 nrow(tca_genes_scaled),nrow(sugar_transporters_scaled),1))) %>%
  rownames_to_column(var='Specific') %>%
  gather(key='Time',value='Intensity',contains('X')) %>%
  mutate(Time=gsub('X','',Time)) %>%
  mutate(number_time=factor(Time,levels=unique(Time))) %>%
  mutate(broad_type=ifelse(type %in% c('Metabolite','TAG'),'Molecule',ifelse(type=='POC','POC','Transcript'))) %>%
  mutate(broad_type=factor(broad_type,levels=c('POC','Molecule','Transcript'))) %>%
  mutate(type=factor(type,levels=c('POC','TAG','Metabolite','Photosynthesis','TCA Cycle','Sugar Transporters'))) %>%
  mutate(special_type=ifelse(type!='Sugar Transporters',as.character(type),ifelse((gsub('_.*$','',Specific)%in%heterotroph_names),'Sugar_Heterotroph','Sugar_Autotroph'))) %>%
  mutate(special_type=ifelse(type!='TCA Cycle',special_type,ifelse((gsub('_.*$','',Specific)%in%heterotroph_names),'TCA_Heterotroph','TCA_Autotroph'))) %>%
  filter(special_type!='Sugar_Autotroph')
ggplot()+
  geom_smooth(data=full_c_cycle_frame,aes(x=number_time,y=Intensity,group=special_type,col=type,
                                          linetype=(special_type %in% c('TCA_Heterotroph',
                                                                        'Sugar_Heterotroph'))),
              size=1.5)+
  facet_wrap(~broad_type,ncol=1,scales='free_y',strip.position='right')+
  theme_bw()+
  scale_x_discrete(breaks=c('0600_727',
                            '1800_727',
                            '0600_728',
                            '1800_728',
                            '0600_729'),labels=c('7/27 6AM',
                                                 '7/27 6PM',
                                                 '7/28 6AM',
                                                 '7/28 6PM',
                                                 '7/29 6AM'),
                   expand=c(0,0))+
  theme(strip.background=element_blank(),
        text=element_text(size=18),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,vjust=0),
        axis.title.y=element_blank(),
        legend.position='bottom')+
  scale_color_brewer(name=NULL,guide=guide_legend(override.aes=list(fill='white'),
                                                  nrow=3),palette='Dark2')+
  guides(linetype='none')+
  annotate('rect',xmax=c('0600_727','0600_728','0600_729'),
           xmin=c('0200_727','1800_727','1800_728'),
           ymin=-Inf,ymax=Inf,alpha=0.25,fill='darkgrey')+
  geom_segment(data=data.frame(x='1000_728',xend='1400_728',y=1.25,yend=1.15,broad_type='Transcript'),
               aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(type='closed',length=unit(0.1,'inches')))+
  geom_segment(data=data.frame(x='0200_729',xend='2200_728',y=1.25,yend=0.35,broad_type='Transcript'),
               aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(type='closed',length=unit(0.1,'inches')))+
  geom_segment(data=data.frame(x='0200_729',xend='0600_729',y=1.25,yend=0.75,broad_type='Transcript'),
               aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(type='closed',length=unit(0.1,'inches')))+
  geom_segment(data=data.frame(x='0600_728',xend='0200_728',y=1.25,yend=0.8,broad_type='Transcript'),
               aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(type='closed',length=unit(0.1,'inches')))+
  geom_text_repel(data=data.frame(x='1000_728',xend='1400_728',y=1.25,yend=1.2,broad_type='Transcript',label='Bacterial\nHeterotrophs'),
                  aes(x=x,y=y,label=label),size=3,
                  force=0,
                  nudge_x=-0.55,
                  nudge_y=0.1)+
  geom_text_repel(data=data.frame(x='0200_729',xend='1400_728',y=1.25,yend=1.2,broad_type='Transcript',label='Photoautotrophs'),
                  aes(x=x,y=y,label=label),size=3,
                  force=0,
                  nudge_x=0,
                  nudge_y=0.2)
ggsave(device='pdf',filename='f5a.pdf')

