## Command to generate figure 3 (data come from script_3b_transcript_ordination.R)

#### NMDS
labeling_function<-function(string,string2){
  current_tax_labels<-unique(as.character(graphing_frame_rotated$new_tax))
  fixed_tax_labels<-c('Metabolite',
                      'Chloroplast Membrane',
                      'TAG',
                      'Quinone',
                      'Cell Membrane',
                      'Carotenoid',
                      'Chloropigment',
                      'Dinophyta (E)',
                      'Other',
                      'Bacillariophyta (E)',
                      'Bicosoecida (E)',
                      'Cercozoa (E)',
                      'Chlorophyta (E)',
                      'Cryptophyta (E)',
                      'Haptophyta (E)',
                      'Ochrophyta (E)',
                      'Rhodophyta (E)',
                      'Sarcomastigophora (E)',
                      'High-Light Prochlorococcus (BA)',
                      'SAR92 (BH)',
                      'Crocosphaera (BA)',
                      'SAR11 (BH)')
  swap_table<-cbind(current_tax_labels,fixed_tax_labels)
  output_index<-swap_table[which(swap_table[,1]==as.character(string)),2]
  if(output_index=='Other'){
    if(string2=='eukaryote'){
      output_index<-'Other (E)'
    }
    if(string2=='prokaryote_non_photoauto'){
      output_index<-'Other (BH)'
    }
    if(string2=='prokaryote_photoauto'){
      output_index<-'Other (BA)'
    }
  }
  return(output_index)
}
new_labels<-c()
for(i in 1:nrow(graphing_frame_rotated)){
  new_labels[i]<-labeling_function(graphing_frame_rotated$new_tax[i],graphing_frame_rotated$tax_group[i])
}
plot_graph_frame<-data.frame(graphing_frame_rotated,plotting_labels=new_labels) %>%
  filter(plotting_labels!='Other')
clock_figure_r<-ggplot(plot_graph_frame,
                       aes(x=x_rot,y=y_rot,col=factor(time_rank)))+
  geom_point()+
  coord_fixed()+
  scale_colour_manual(name='Peak Expression',
                      labels=c('2200','0000','0200','0400','0600','0800','1000','1200','1400',
                               '1600','1800','2000'),
                      values=c('blue4','navy','royalblue3','lightblue1',
                               'khaki1','gold','darkorange','red','salmon',
                               'maroon','magenta4','dodgerblue4'))+
  facet_wrap(~plotting_labels,labeller = label_wrap_gen(width=6),
             ncol=5)+
  theme(panel.background=element_blank(),
        legend.text=element_text(size=16),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        strip.background = element_blank(),
        text=element_text(size=16),
        legend.key=element_rect(fill='white'),
        legend.key.size = unit(2, "lines"),
        legend.position='bottom',
        strip.text.y=element_text(angle=180),
        strip.text=element_text(size=11))+
  guides(colour = guide_legend(override.aes = list(size=4),position='bottom'))
ggsave(plot=clock_figure_r,filename='../figures/nmds_figure_try_fixed_labels.png',height=11,width=8,units='in',
       device='png',dpi=400)

###### Pulling out sensitive vs agnostic top 10 # taxa kos
sen_top10<-c('K02276','K04078','K03686','K02706','K04079','K00600','K00123','K04487','K03568','K00134')
ag_top10<-c('K00239','K03798','K04077','K04043','K02275','K02274','K00240','K02358','K03403','K02703')
sen_frame<-graphing_frame_rotated %>%
  filter(kos %in% sen_top10,
         tax_group %in% c('eukaryotic','prokaryotic_non_photoauto','prokaryotic_photoauto'))
ag_frame<-graphing_frame_rotated %>%
  filter(kos %in% ag_top10,
         tax_group %in% c('eukaryotic','prokaryotic_non_photoauto','prokaryotic_photoauto'))
ggplot(sen_frame,aes(x=x_rot,y=y_rot,col=tax_group))+
  geom_point(size=3)+
  facet_wrap(~kos)+
  ylab('NMDS2')+
  xlab('NMDS1')+
  theme(panel.background=element_blank(),
        text=element_text(size=18))+
  ggtitle('Asynchronous')
ggplot(ag_frame,aes(x=x_rot,y=y_rot,col=tax_group))+
  geom_point(size=3)+
  facet_wrap(~kos)+
  ylab('NMDS2')+
  xlab('NMDS1')+
  theme(panel.background=element_blank(),
        text=element_text(size=18))+
  ggtitle('Not asynchronous')
