## Making new figure 5b component

diel_amino_acids<-c('Tyrosine',
                    'Serine',
                    'Methionine',
                    'Taurine','Valine',
                    'Alanine','Glutamine',
                    'Asparagine','Aspartic acid',
                    'Proline','Tryptophan','Phenylalanine',
                    'Betaine','Isoleucine-Leucine')
diel_aa_transporters<-c('K01999','K02000','K02002',
                        'K02030','K02049','K02050',
                        'K02051','K11954','K15551','K09969','K11928')
updated_titles<-c('Tyr','Ser','Met','Taurine','Val','Ala','Gln','Asn',
                  'Asp','Pro','Trp','Phe','Betaine','Ile/Leu',
                  'Val/Ile/Leu','Betaine/Pro','Betaine/Pro',
                  'Polar AAs',
                  'Taurine','Taurine','Taurine',
                  'Neutral AAs','Taurine','General AAs','Proline'
                  )
title_frame<-data.frame(informal=c(diel_amino_acids,diel_aa_transporters),formal=updated_titles,
                        stringsAsFactors = FALSE)
set.seed(1509)
amino_acid_data<-graphing_frame_rotated %>%
  filter(taxa %in% diel_amino_acids | kos %in% diel_aa_transporters) %>%
  mutate(kos=ifelse(kos=='Metabolite',as.character(taxa),as.character(kos))) %>%
  full_join(title_frame,by=c('kos'='informal')) %>%
  mutate(tax_group=factor(tax_group,levels=c('eukaryotic','prokaryotic_photoauto','prokaryotic_non_photoauto','Chemical'))) %>%
  mutate(theta_point=atan2(y_rot,x_rot),
         label_y=4.25*sin(theta_point),
         label_x=4.25*cos(theta_point),
         label_y_jit=jitter(label_y,amount=1.5),
         label_x_jit=jitter(label_x,amount=1.5),
         h=ifelse(x_rot>0.35,1,0),
         xnudge=ifelse(x_rot>0.35,4-x_rot,-4-x_rot)) %>%
  filter(!(tax_group %in% c('eukaryotic','prokaryotic_photoauto'))) %>%
  mutate(plot_tax=ifelse(taxa %in% diel_amino_acids, 'Metabolite',as.character(taxa)),
         plot_tax2=factor(plot_tax,labels=c('Chloroflexi','Metabolite','Alphaproteobacteria',
                                            'Gammaproteobacteria','Roseobacter','SAR11-like','SAR116'))) %>%
  mutate(plot_tax2=factor(plot_tax2,levels=c('Alphaproteobacteria','SAR11-like','Roseobacter','Gammaproteobacteria','SAR116','Chloroflexi','Metabolite')))

transporters_color_label<-ggplot(amino_acid_data)+
  geom_point(aes(x=x_rot,y=y_rot,col=plot_tax2),size=5)+
  coord_fixed()+
  theme_bw()+
  scale_color_brewer(palette='Paired',name='',guide=guide_legend(ncol=2))+
  scale_fill_brewer(palette='Paired',guide='none')+
  ggforce::geom_circle(aes(x0=0.35,y0=0.25,r=3.5))+
  ggrepel::geom_label_repel(mapping=aes(x=x_rot,y=y_rot,label=formal,
                                fill=plot_tax2),
                            alpha=0.75,
                            force=25,
                            size=4,
                            point.padding=1)+
  facet_grid(~tax_group,labeller = as_labeller(c('eukaryotic'='Eukaryote',
                                               'prokaryotic_non_photoauto'='Heterotroph',
                                                 'prokaryotic_photoauto'='Cyanobacteria',
                                                 'Chemical'='Amino Acid')),
             switch='y')+
  #geom_label(aes(x=label_x_jit,y=label_y_jit,label=formal,fill=plot_tax2),
  #                          alpha=0.75)+
  #geom_segment(aes(x=label_x_jit*0.925,y=label_y_jit*0.925,xend=x_rot,yend=y_rot))+
  theme(text=element_text(size=18),
        strip.background=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position='bottom',
        strip.text=element_blank())
transporters_bw_label<-ggplot(filter(amino_acid_data,!(tax_group %in% c('eukaryotic','prokaryotic_photoauto'))) %>%
                                mutate(plot_tax=ifelse(taxa %in% diel_amino_acids, 'Metabolite',as.character(taxa)),
                                       plot_tax2=factor(plot_tax,labels=c('Chloroflexi','Metabolite','Alphaproteobacteria',
                                                                          'Gammaproteobacteria','Roseobacter','SAR11-like','SAR116'))) %>%
                                mutate(plot_tax2=factor(plot_tax2,levels=c('Alphaproteobacteria','SAR11-like','Roseobacter','Gammaproteobacteria','SAR116','Chloroflexi','Metabolite'))))+
  geom_point(aes(x=x_rot,y=y_rot,col=plot_tax2),size=3.5)+
  facet_grid(tax_group~.,labeller = as_labeller(c('eukaryotic'='Eukaryote',
                                                  'prokaryotic_non_photoauto'='Heterotroph',
                                                  'prokaryotic_photoauto'='Cyanobacteria',
                                                  'Chemical'='Amino Acid')),
             switch='y')+
  coord_fixed()+
  theme_bw()+
  scale_color_brewer(palette='Paired',name='',guide=guide_legend(ncol=2))+
  ggforce::geom_circle(aes(x0=0.35,y0=0.25,r=3.5))+
  ggrepel::geom_label_repel(aes(x=x_rot,y=y_rot,label=formal),force=3,size=4.5)+
  theme(text=element_text(size=18),
        strip.background=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position='bottom',
        strip.text=element_blank())

ggsave(plot=transporters_color_label,filename='../figures/transporters_color4.pdf',device='pdf')
ggsave(plot=transporters_bw_label,filename='../figures/transpoerters_bw.pdf',device='pdf')
## NVM but we rolling

