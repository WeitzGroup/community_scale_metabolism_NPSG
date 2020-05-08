## The goal of this script will be to produce a summatory figure showing the ideas behind the
## N partitioning and contributions to daily C cycle stuff
## Figuring out what data we want

## Time series?

## POC Data
optics_data<-read.csv('~/Dropbox (GaTech)/Weitz Group Team Folder/Manuscripts/Member papers/dmuratore/2019mahalo_growth_met/growth_met_source/intermediate_data/optics_data.csv')
optics_data_matched<-scale(optics_data[c(1,5,9,
                                   13,18,21,
                                   25,29,33,
                                   37,41,45,
                                   49,53,57),c('POC')])
## hNB/AA data

## Metabolite (Fixed C) + TAG time series data? Units or no units?
## Let's try the smoothing thing so we go for scaled data
metabolite_scaled<-lined_up_full_mat[c('UDP-glucose',
                                       'UDP-glucosamine',
                                       'Trehalose',
                                       'Choline',
                                       'Glucosylglycerol',
                                       'Ribose 5 phosphate',
                                       'DHPS','Gluconic Acid',
                                       'Sucrose','Vanillic Acid','Chitobiose',
                                       'trans Retinal','DMSPE','HomarineE'),]
lipid_scaled<-lined_up_full_mat[grep('TAG',rownames(lined_up_full_mat)),]

## Getting photosynthesis and TCA cycle genes
photo_gene_kos<-as.character(unique(filter(graphing_frame_rotated,path=='Photosynthesis')$kos))
photo_genes_scaled<-lined_up_full_mat[(gsub('^.*_','',rownames(lined_up_full_mat)) %in% photo_gene_kos),]
tca_gene_kos<-as.character(unique(filter(graphing_frame_rotated,path=='Citrate cycle (TCA cycle)')$kos))
tca_genes_scaled<-lined_up_full_mat[(gsub('^.*_','',rownames(lined_up_full_mat)) %in% tca_gene_kos),]
heterotroph_names<-setdiff(unique(gsub('_.*$','',clade_diels$X)),c('High-light Pro','Low-light Pro',
                                                            'Crocosphaera','Viruses','Eukaryota',
                                                            'Synechococcus','Cyanothece','Euryarchaeota','--'))
full_c_cycle_frame<-data.frame(rbind(metabolite_scaled,
                          lipid_scaled,
                          photo_genes_scaled,
                          tca_genes_scaled,
                         as.numeric(optics_data_matched)),type=rep(c('Metabolite','TAG','Photosynthesis','TCA Cycle','POC'),
                                                     c(nrow(metabolite_scaled),
                                                       nrow(lipid_scaled),
                                                       nrow(photo_genes_scaled),
                                                     nrow(tca_genes_scaled),1))) %>%
  rownames_to_column(var='Specific') %>%
  gather(key='Time',value='Intensity',contains('X')) %>%
  mutate(Time=gsub('X','',Time)) %>%
  mutate(number_time=factor(Time,levels=unique(Time))) %>%
  mutate(broad_type=ifelse(type %in% c('Metabolite','TAG'),'Molecule',ifelse(type=='POC','POC','Transcript'))) %>%
  mutate(broad_type=factor(broad_type,levels=c('POC','Molecule','Transcript'))) %>%
  mutate(type=factor(type,levels=c('POC','TAG','Metabolite','Photosynthesis','TCA Cycle'))) %>%
  mutate(special_type=ifelse(type!='TCA Cycle',as.character(type),ifelse((gsub('_.*$','',Specific)%in%heterotroph_names),'Heterotroph','Autotroph')))
ggplot()+
  geom_smooth(data=full_c_cycle_frame,aes(x=number_time,y=Intensity,group=special_type,col=type,
                                          linetype=(special_type=='Heterotroph')),
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
  geom_segment(data=data.frame(x='0200_729',xend='2200_728',y=1.25,yend=0.5,broad_type='Transcript'),
               aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(type='closed',length=unit(0.1,'inches')))+
  geom_segment(data=data.frame(x='0200_729',xend='0600_729',y=1.25,yend=0.75,broad_type='Transcript'),
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
ggsave(device='pdf',filename='../figures/f5_test.pdf')
## Gene time series data chunked up by pathway association? What pathways do we use

## Ribosome
## The genes we pulled for the NMDS subset
## Photosynthesis

nbs_aas<-read.csv("zscore_timeseries_n_related_529.csv",row.names=1)
just_nbaa<-nbs_aas[c(59,60),]
colnames(just_nbaa)<-gsub('X','',colnames(just_nbaa))
glutamine<-lined_up_full_mat['Glutamine',]
granular_annotes<-read.csv('specific_gene_information_517.csv')
gs_ones<-lined_up_full_mat[which(gsub('^.*_','',rownames(lined_up_full_mat)) %in% granular_annotes$kos[grep('gln',granular_annotes$gene_function)]),]
gogat_ones<-lined_up_full_mat[which(gsub('^.*_','',rownames(lined_up_full_mat)) %in% granular_annotes$kos[grep('glt',granular_annotes$gene_function)]),]
uptake_ones<-lined_up_full_mat[which(gsub('^.*_','',rownames(lined_up_full_mat)) %in% granular_annotes$kos[grep('Transport',granular_annotes$path)]),]
ribosome_ones<-lined_up_full_mat[which(gsub('^.*_','',rownames(lined_up_full_mat)) %in% graphing_frame_rotated$kos[grep('Ribosome',graphing_frame_rotated$path)]),]
tRNA_biosynth<-lined_up_full_mat[which(gsub('^.*_','',rownames(lined_up_full_mat)) %in% graphing_frame_rotated$kos[grep('Aminoacyl-tRNA biosynthesis',graphing_frame_rotated$path)]),]
peptide_aa_transporters<-lined_up_full_mat[which(gsub('^.*_','',rownames(lined_up_full_mat)) %in% new_transporters$kos[which(new_transporters$General_category=='amino_acids')]),]
new.categories <- read_csv("kos_pathway_geneFunction_plus.csv") %>% select(-X1) %>%
  filter(kos!="K04488",kos!="K02588")
nif_ones<-lined_up_full_mat[which(gsub('^.*_','',rownames(lined_up_full_mat)) %in% new.categories$kos[grep('niF',new.categories$Gene_Function_2)]),]
general_categories<-c(rep('Molecule',3),rep('Assimilation\nand\nTranslation',sum(nrow(gs_ones),nrow(gogat_ones),nrow(ribosome_ones),nrow(tRNA_biosynth))),
                      rep('Uptake',nrow(uptake_ones)+nrow(nif_ones)),rep('Uptake',nrow(peptide_aa_transporters)))
more_specific_categories<-c('Nitrogenous Bases','Amino Acids','Glutamine',
                            rep('Glutamine Synthetase',nrow(gs_ones)),
                            rep('Glutamate Synthase',nrow(gogat_ones)),
                            rep('Ribosomal',nrow(ribosome_ones)),
                            rep('tRNA Synthesis',nrow(tRNA_biosynth)),
                            as.character(granular_annotes$gene_function[match(rownames(uptake_ones),granular_annotes$full_ids)]),
                            rep('N2 Fixation',nrow(nif_ones)),
                            rep('AA+Peptide Transporters',nrow(peptide_aa_transporters)))
more_specific_categories[which(more_specific_categories=='amt_ammonia_transporter')]<-'Ammonia'
more_specific_categories[which(more_specific_categories=='NRT_NO34_transport')]<-'Nitrate'
more_specific_categories[grep('nitT_',more_specific_categories)]<-'Nitrite'
more_specific_categories[grep('niT_',more_specific_categories)]<-'Nitrite'
more_specific_categories[grep('urtA',more_specific_categories)]<-'Urea'

nitrogen_everything<-cbind(rbind(just_nbaa,glutamine,gs_ones,gogat_ones,ribosome_ones,tRNA_biosynth,uptake_ones,nif_ones,peptide_aa_transporters),
                           general_categories,more_specific_categories) %>%
  rownames_to_column(var='specific_id') %>%
  mutate(specific_id=ifelse(specific_id=='mean_hydrolyzable_nbs_ugL','NB',
                            ifelse(specific_id=='mean_hydrolyzable_aas_ugL','AA',specific_id))) %>%
  mutate(taxon=gsub('_.*$','',specific_id)) %>%
  mutate(tax_cat=ifelse(taxon %in% heterotroph_names, 'Heterotroph',
                        ifelse(taxon %in% c('High-light Pro',
                                            'Low-light Pro',
                                            'Crocosphaera',
                                            'Synechococcus',
                                            'Cyanothece'), 'Cyanobacteria','Eukaryote'))) %>%
  gather(key='Time',value='Intensity',contains('00')) %>%
  mutate(number_time=factor(Time,levels=unique(Time))) %>%
  mutate(finest_grain=paste(more_specific_categories,tax_cat)) %>%
  mutate(general_categories=factor(general_categories,levels=c('Molecule','Assimilation\nand\nTranslation',
                                                               'Uptake'))) %>%
mutate(more_specific_categories=factor(more_specific_categories,levels=c('Amino Acids',
                                                                         'Nitrogenous Bases',
                                                                         'Glutamine',
                                                                         'Ribosomal',
                                                                         'tRNA Synthesis',
                                                                         'Glutamate Synthase',
                                                                         'Glutamine Synthetase',
                                                                         'N2 Fixation',
                                                                         'Ammonia','Nitrate','Nitrite',
                                                                         'Urea','AA+Peptide Transporters')))



summary_n_frame<-nitrogen_everything %>%
  filter(more_specific_categories!='tRNA Synthesis') %>%
  group_by(general_categories,more_specific_categories,tax_cat,finest_grain,number_time) %>%
  summarize(mean_int=mean(Intensity))

ggplot()+
  #geom_smooth(data=nitrogen_everything,aes(x=number_time,y=Intensity,col=more_specific_categories,
  #                                       group=more_specific_categories),
  #            size=1.5,method='loess')+
  #geom_line(data=summary_n_frame,aes(x=number_time,y=mean_int,linetype=tax_cat,group=finest_grain,
  #                                   col=more_specific_categories),
  #          size=1.5)+
  geom_point(data=summary_n_frame,
             aes(x=number_time,y=mean_int,col=more_specific_categories),size=2)+
  geom_line(data=summary_n_frame,
            aes(x=number_time,y=mean_int,col=more_specific_categories,group=finest_grain))+
  facet_grid(tax_cat~general_categories,scales='free_y')+
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
  scale_color_discrete(name='Data Type',guide=guide_legend(override.aes=list(fill='white'),
                                                         ncol=4))+
  #scale_shape_discrete(guide=guide_legend(nrow=3))+
  #guides(linetype='none')+
  scale_linetype(guide=guide_legend(nrow=3),name='Broad Taxonomy',palette=)+
  annotate('rect',xmax=c('0600_727','0600_728','0600_729'),
           xmin=c('0200_727','1800_727','1800_728'),
           ymin=-Inf,ymax=Inf,alpha=0.25,fill='darkgrey')
ggsave(device='pdf',filename='../figures/figure_5b_try.pdf')


more_filtered_n<-nitrogen_everything %>%
  filter(more_specific_categories %in% c('Ammonia','Glutamate Synthase','Glutamine Synthetase',
                                         'Amino Acids','Nitrogenous Bases','Ribosomal','AA+Peptide Transporters')) %>%
  mutate(updated_cat=ifelse(more_specific_categories %in% c('Glutamate Synthase','Glutamine Synthetase'),
                            'GS-GOGAT',as.character(more_specific_categories))) %>%
  mutate(updated_cat=factor(updated_cat,levels=c('Amino Acids','Nitrogenous Bases','Ribosomal','GS-GOGAT','Ammonia','AA+Peptide Transporters'))) %>%
  mutate(broader_cat=ifelse(updated_cat %in% c('Amino Acids','Nitrogenous Bases'),'Molecule',
                            as.character(updated_cat))) %>%
  mutate(tax_cat=ifelse(broader_cat=='Molecule',as.character(updated_cat),as.character(tax_cat))) %>%
  mutate(new_finest_grain=paste(updated_cat,tax_cat)) %>%
  group_by(general_categories,broader_cat,updated_cat,tax_cat,new_finest_grain,number_time) %>%
  summarize(mean_int=mean(Intensity)) %>%
  mutate(facet_cat=ifelse(broader_cat=='Molecule',
                          'Molecule',
                          ifelse(broader_cat=='Ribosomal','Ribosomal',
                          ifelse(broader_cat=='GS-GOGAT',
                                 'N Assimilation','N Uptake')))) %>%
  mutate(facet_cat=factor(facet_cat,levels=c('Molecule','Ribosomal',
                                             'N Assimilation',
                                             'N Uptake'))) %>%
  ungroup() %>%
  mutate(tax_cat=factor(tax_cat,levels=c('Amino Acids','Nitrogenous Bases','Cyanobacteria','Eukaryote','Heterotroph'))) %>%
  mutate(tax_cat2=ifelse(tax_cat %in% c('Amino Acids','Nitrogenous Bases'),'Eukaryote',as.character(tax_cat)),
         updated_cat2=ifelse(broader_cat=='Molecule','Molecule',as.character(updated_cat)),
         updated_cat2=factor(updated_cat2,c('Molecule','Ribosomal','GS-GOGAT','Ammonia','AA+Peptide Transporters')))


ggplot()+
  geom_line(data=filter(more_filtered_n,updated_cat2!='Molecule'),aes(x=number_time,
                                     y=mean_int,
                                     group=new_finest_grain,
                                     color=updated_cat),
            size=1.5)+
  facet_grid(facet_cat~tax_cat2,
             scale='free_y',drop=TRUE)+
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
        text=element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,vjust=0),
        axis.title.y=element_blank(),
        legend.position='left',
        strip.text=element_text(size=14))+
  scale_linetype(guide=guide_legend(nrow=3),name='Broad Taxonomy',palette=)+
  annotate('rect',xmax=c('0600_727','0600_728','0600_729'),
           xmin=c('0200_727','1800_727','1800_728'),
           ymin=-Inf,ymax=Inf,alpha=0.25,fill='darkgrey')+
  scale_color_brewer(name=NULL,guide=guide_legend(override.aes=list(fill='white'),
                                                         nrow=6),palette='Dark2',
                     labels=c('Ribosomal','GS-GOGAT','Ammonia Transporters','Amino Acid Transporters'))#+
  #ggrepel::geom_label_repel(data=data.frame(updated_cat2='Molecule',tax_cat2='Eukaryote',x='0600_727',y=c(0.57531460,
  #                                                                                                        0.06036468),label=c('Amino Acids','Nitrogenous Bases')),
  #          aes(x=x,y=y,label=label),segment.color='darkred',nudge_y=0.25,nudge_x=0.25)


ggplot()+
  geom_line(data=filter(more_filtered_n,updated_cat2!='Molecule') %>%
              filter(!(updated_cat=='AA+Peptide Transporters' & tax_cat2 %in% c('Eukaryote','Cyanobacteria'))),aes(x=number_time,
                                                                      y=mean_int,
                                                                      group=new_finest_grain,
                                                                      color=updated_cat),
            size=1.5)+
  facet_grid(facet_cat~tax_cat2,
             scale='free_y',drop=TRUE)+
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
        text=element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,vjust=0),
        axis.title.y=element_blank(),
        legend.position='left',
        strip.text=element_text(size=14))+
  scale_linetype(guide=guide_legend(nrow=3),name='Broad Taxonomy',palette=)+
  annotate('rect',xmax=c('0600_727','0600_728','0600_729'),
           xmin=c('0200_727','1800_727','1800_728'),
           ymin=-Inf,ymax=Inf,alpha=0.25,fill='darkgrey')+
  scale_color_brewer(name=NULL,guide=guide_legend(override.aes=list(fill='white'),
                                                  nrow=6),palette='Dark2',
                     labels=c('Ribosomal','GS-GOGAT','Ammonia Transporters','Amino Acid Transporters'))

f5b_grob<-ggplotGrob(fig_5bplot)
idx <- which(f5b_grob$layout$name %in% c("panel-1-1", "panel-3-3"));
for (i in idx) f5b_grob$grobs[[i]] <- grid::nullGrob();
grid::grid.newpage()
grid::grid.draw(f5b_grob)
ggsave(device='pdf',filename='f5b_newtry_transcripts.pdf')
moew<-ggplot()+
  geom_line(data=filter(more_filtered_n,updated_cat2=='Molecule'),aes(x=number_time,
                                                                      y=mean_int,
                                                                      group=new_finest_grain,
                                                                      color=updated_cat),
            size=1.5)+
  facet_wrap(~updated_cat2,
             scale='free_y',drop=TRUE)+
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
        strip.text=element_blank(),
        legend.position='bottom')+
  scale_linetype(guide=guide_legend(nrow=3),name='Broad Taxonomy')+
  annotate('rect',xmax=c('0600_727','0600_728','0600_729'),
           xmin=c('0200_727','1800_727','1800_728'),
           ymin=-Inf,ymax=Inf,alpha=0.25,fill='darkgrey')+
  scale_color_manual(name=NULL,guide=guide_legend(override.aes=list(fill='white'),
                                                  nrow=3),values=c('black','magenta'))
meow2<-ggplotGrob(moew)
ggsave(meow2,device='pdf',filename='f5b_newtry_chemicals.pdf')
