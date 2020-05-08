### Doing some transporter nonsense
transporter_subset<-filter(graphing_frame_rotated,path=='Transport')
transporter_prevalence<-transporter_subset %>% count(kos) %>% arrange(desc(n))
het_trans<-filter(transporter_subset,tax_group=='prokaryotic_non_photoauto')
het_trans_summary<-het_trans %>%
   group_by(kos,cluster) %>%
   summarize(ntax=n()) %>%
  ungroup() %>%
  spread('cluster','ntax',fill=0) %>%
  arrange(desc(Night),desc(Dusk))
View(het_trans_summary)
ggplot(transporter_subset)+
  geom_point(aes(x=x_rot,y=y_rot,col=tax_group))+
  facet_wrap(~kos)
write.csv(transporter_subset,'intermediate_data/transporter_unannotated.csv')
new_transporters<-read.csv('../intermediate_data/transporter_annotated.csv')[,-1]
ggplot(filter(new_transporters,!(tax_group %in% c('','archaea'))) %>% filter(General_category %in% c('amino_acids','inorganic','organic')))+
  geom_point(aes(x=x_rot,y=y_rot,col=Specific_category))+
  facet_grid(tax_group~General_category)+
  coord_fixed()+
  theme_bw()+
  ggforce::geom_circle(aes(x0=0.35,y0=0.25,r=3.5))+
  theme(text=element_text(size=16),
        strip.background=element_blank(),
        axis.title=element_blank())

ggplot(filter(new_transporters,!(tax_group %in% c('','archaea'))))+
  geom_point(aes(x=x_rot,y=y_rot,col=General_category),size=4)+
  facet_grid(~tax_group,labeller = as_labeller(c('eukaryotic'='Eukaryote',
                                                                 'prokaryotic_non_photoauto'='Heterotroph',
                                                                 'prokaryotic_photoauto'='Cyanobacteria',
                                                                 'amino_acids'='Amino Acids',
                                                                 'inorganic'='Inorganics',
                                                                 'metal'='Metals',
                                                                 'organic'='Organic Molecules',
                                                                 'other'='Other Junk',
                                                                 'sugar'='Sugars')))+
  coord_fixed()+
  theme_bw()+
  ggforce::geom_circle(aes(x0=0.35,y0=0.25,r=3.5))+
  #ggrepel::geom_label_repel(aes(x=x_rot,y=y_rot,label=Specific_category))+
  theme(text=element_text(size=16),
        strip.background=element_blank(),
        axis.title=element_blank())


ggplot(filter(new_transporters,!(tax_group %in% c('','archaea')) & General_category=='amino_acids'))+
  geom_point(aes(x=x_rot,y=y_rot),size=2)+
  facet_grid(~tax_group,labeller = as_labeller(c('eukaryotic'='Eukaryote',
                                                                 'prokaryotic_non_photoauto'='Heterotroph',
                                                                 'prokaryotic_photoauto'='Cyanobacteria',
                                                                 'amino_acids'='Amino Acids',
                                                                 'inorganic'='Inorganics',
                                                                 'metal'='Metals',
                                                                 'organic'='Organic Molecules',
                                                                 'other'='Other Junk',
                                                                 'sugar'='Sugars')))+
  coord_fixed()+
  theme_bw()+
  ggforce::geom_circle(aes(x0=0.35,y0=0.25,r=3.5))+
  ggrepel::geom_label_repel(aes(x=x_rot,y=y_rot,label=Specific_category))+
  theme(text=element_text(size=16),
        strip.background=element_blank(),
        axis.title=element_blank())

View(filter(graphing_frame_rotated,cluster=='Morning') %>% 
       filter(taxa %in% c('Bacillariophyta','Dinophyta')) %>%
       arrange(path))


## Working now on B12 stuff

meow<-c('K00798','K19221',
        'K02232',
        'K02225','K02227',
        'K02231',
        'K00768','K02226','K22316','K02233')
graphing_frame_rotated[which(graphing_frame_rotated$kos %in% meow),]
