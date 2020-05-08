### Doing some summary analysis wrt heterotroph metabolism using organic substrates
playing_frame<-graphing_frame_rotated # in case we want to mess with some stuff
heterotroph_subset<-filter(graphing_frame_rotated,tax_group=='prokaryotic_non_photoauto')
## Total 819 bacterial heterotroph diel transcripts
## Let's see the pathway distribution
path_distribution<-table(heterotroph_subset$path)
## Things that we can work with
# Ala/Asp/Glu metabolism - 4 transcripts
# Cys/Met metabolism - 8 transcripts
# Fatty acid biosynthesis - 5 transcripts
# Fatty acid degredation - 9 transcripts
# Glycolysis stuff - 19 transcripts
# Nicotinate/nicotinamide metabolism 5 transcripts
# Gly/Ser/Thr metabolism - 19 transcripts
# Phe/Tyr/Trp metabolism - 6 transcripts
# N metabolism - 20 transcripts
# Purine metabolism - 14 transcripts
# Pyrimidine metabolism - 8 transcripts
# 2-component systems 10 transcripts
# Transporters 71 transcripts

## Of the amino acid metabolism ones, Gly/Ser/Thr metabolism appears to be really well represented
## so let's look at what those are
gly_ser_thr<-filter(heterotroph_subset,path=='Glycine, serine and threonine metabolism') %>%
  arrange(kos,taxa)
## Kos that appear a lot are K00281, K00315, K00643
## Those are glycine dehydrogenase, dimethylglycine dehydrogenase, 5-aminolevulinate synthase
## Furthermore, we see diel expression of gcvT (aminomethyltransferase) K00605 in the same
## heterotrophs which according to KEGG catalyzes the rxn 	
## [Protein]-S8-aminomethyldihydrolipoyllysine + Tetrahydrofolate <=> 
## Dihydrolipoylprotein + 5,10-Methylenetetrahydrofolate + Ammonia
## Maybe not bad starting point
## Furthermore, we have precursor to DMGly Betaine in the diel metabolites :D 
## However comma, we also have DMGly *and* Sarcosine in the metabolite data and it was found to not be diel D:
## boo. 


## Figuring out pathway
## Choline -> Betaine
cho_betainealdehyde<-c('K60499',
           'K17755',
           'K00108',
           'K11440')
betaine_aldehyde_betaine<-c('K00130',
                            'K14085',
                            'K17755')
betaine_proline_uptake<-c('K02000','K02002')
betaine_dmg<-c('K00544')
dmg_sarcosine<-c('K00309','K00315')
sarcosine_gly<-c('K00301','K00302','K00303','K00304','K00305','K00306','K00314')
aminomethyltransferase<-c('K00605')
gly_lipoylprotein<-c('K00281','K00282','K00283')
metabolites<-c('BetaineE','Proline','Choline')
pathway_list<-list(cho_betainealdehyde=cho_betainealdehyde,
                   betaine_aldehyde_betaine=betaine_aldehyde_betaine,
                   betaine_proline_uptake=betaine_proline_uptake,
                   betaine_dmg=betaine_dmg,
                   dmg_sarcosine=dmg_sarcosine,
                   sarcosine_gly=sarcosine_gly,
                   aminomethyltransferase=aminomethyltransferase,
                   gly_lipoylprotein=gly_lipoylprotein)
full_list<-lapply(pathway_list,function(x) graphing_frame_rotated[which(graphing_frame_rotated$kos %in% x),])
n_hits<-do.call(rbind,lapply(full_list,nrow))
path_names<-rep(names(pathway_list),n_hits)
betaine_frame<-data.frame(do.call(rbind,full_list),path_names=path_names)
betaine_mets<-data.frame(graphing_frame_rotated[which(graphing_frame_rotated$full_ids %in% metabolites),],path_names=metabolites)
betaine_frame<-rbind(betaine_frame,betaine_mets)
ggplot(betaine_frame)+
  ggforce::geom_circle(aes(x0=0,y0=0,r=3))+
  geom_point(data=betaine_frame,aes(x=x,y=y,color=taxa),size=3)+
  coord_fixed()+
  #ggrepel::geom_label_repel(data=betaine_frame,aes(x=x,y=y,label=path_names))+
  facet_wrap(~path_names)+
  theme_bw()

betaine_het_only<-filter(betaine_frame,tax_group %in% c('prokaryotic_non_photoauto','Chemical')) %>%
  filter(!(full_ids %in% c('Choline','Proline'))) %>%
  mutate(path_names=factor(path_names,levels=c('Choline',
                                               'betaine_proline_uptake',
                                               'dmg_sarcosine',
                                               'sarcosine_gly',
                                               'gly_lipoylprotein',
                                               'aminomethyltransferase')))
ggplot(betaine_het_only)+
  ggforce::geom_circle(aes(x0=0,y0=0,r=3))+
  geom_point(aes(x=x_rot,y=y_rot,color=path_names),size=3)+
  coord_fixed()+
  #ggrepel::geom_label_repel(data=betaine_frame,aes(x=x,y=y,label=path_names))+
  facet_wrap(~taxa)+
  theme_bw()+
  scale_color_discrete(name='Enzyme',labels=c('Aminomethyltransferase',
                                              'Betaine Uptake',
                                              'DMG -> Sarcosine',
                                              'Glycine -> Lipoylprotein',
                                              'Sarcosine -> Glycine','Betaine')[c(6,2,3,5,4,1)])+
  theme(strip.background=element_blank(),
        axis.title=element_blank())
