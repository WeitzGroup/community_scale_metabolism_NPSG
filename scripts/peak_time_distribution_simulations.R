### Updated analysis for identifying KOs that peak at the same time
### across taxa with diel expression of that KO

## Setting up necessary functions
## Calculating the circular distance between two peak times
calculate_med_dist<-function(pair){
  difference<-abs(pair[1]-pair[2])
  return(min(difference,6-difference))
}
## Finding the mean pairwise circular distance between all peaktimes in a set
generate_mean_stat<-function(vec){
  all_pairwise<-t(combn(vec,2))
  all_differences<-apply(all_pairwise,1,calculate_med_dist)
  return(mean(all_differences))
}
## Using the number of diel transcripts for one KO, simulate a distribution of
## average pairwise circular distance between peaktimes of KOs drawn at random from all diel KOs
## We do this because the peaktimes are not uniformly distributed over the 6 peak times
generate_null_ensemble<-function(vec,pop,niter=1e4){
  number_samples<-length(vec)
  bootstrap_ensemble<-t(replicate(niter,sample(x=pop,size=number_samples,replace=TRUE)))
  null_dn<-apply(bootstrap_ensemble,1,generate_mean_stat)
  return(null_dn)
}
## Use the null ensemble to simulate a p-value
simulate_p<-function(vec,nd){
  sample_statistic<-generate_mean_stat(vec)
  simulated_p<-length(which(nd<=sample_statistic))/length(nd)
  return(simulated_p)
}
## A function to do all of these things and keep you updated on what's happening
full_wrapper<-function(vec,pop,niter=1e4){
  print('Reading Data')
  vec<-vec
  pop<-pop
  print('Monte Carlo Simulating Null Distribution (the slow part)...')
  nd<-generate_null_ensemble(vec,pop,niter=niter)
  print('Null distribution generated')
  print('Calculating empirical p-value')
  pval<-simulate_p(vec,nd)
  print('Test complete, gathering outputs')
  output_list<-list(pval=pval,md=generate_mean_stat(vec),null=nd)
  return(output_list)
}

## Implementing on our data

pop_times<-floor(graphing_frame_rotated$time_rank)

ko_counts<-graphing_frame_rotated %>%
  group_by(kos) %>%
  count() %>%
  arrange(desc(n))

## Only test on KOs which had at least 4 taxa with diel transcription of that KO
keeping_kos<-ko_counts[ko_counts$n>=4,1]

output_storage<-list()
for(i in 1:length(keeping_kos$kos)){
  ko_of_i<-keeping_kos$kos[i]
  sub_vec<-floor(graphing_frame_rotated$time_rank[which(graphing_frame_rotated$kos==ko_of_i)])
  output_storage[[i]]<-full_wrapper(sub_vec,pop_times)
  print(i)
}
## Do 10% ABH
## Get p-values
output_ps<-do.call(rbind,lapply(output_storage,function(x) x$pval))
output_diffs<-do.call(rbind,lapply(output_storage, function(x) x$md))
q<-0.1
i<-1:length(output_ps)
m<-length(output_ps)
ordered_ps<-output_ps[order(output_ps,decreasing=FALSE)]
stop_reject<-min(which(ordered_ps>=i*q/m))
decisions<-rep(c('reject','fail'),c(stop_reject,length(ordered_ps)-stop_reject))
plot(1:m,ordered_ps)
lines(1:m,i*q/m,col='red',lwd=2)
abline(v=stop_reject)
pval_frame<-data.frame(pval=output_ps,
                       diffs=output_diffs,
                       ko=keeping_kos$kos,
                       n_tax=ko_counts$n[ko_counts$n>=4]) %>%
  arrange(pval) %>%
  mutate(rejects=decisions)
## writing to output
write.csv(pval_frame,quote=FALSE,file='./outputs_and_summaries/peak_time_difference_no_annote.csv')
## KEGG Annotations were added manually from KEGG.jp web server

ggplot(pval_frame)+
  geom_point(aes(x=diffs,y=pval,col=factor(n_tax),shape=rejects))+
  geom_line(aes(x=diffs,y=pval,group=factor(n_tax),color=factor(n_tax)))

ggplot(pval_frame)+
  geom_boxplot(aes(x=factor(n_tax),y=diffs*4))+
  geom_jitter(data=filter(pval_frame,!(ko %in% c('K04752','K00266','K01915','K00284','K00265',
                                                 'K00239','K03798','K04043','K02706','K02703','K02274'))),
              aes(x=factor(n_tax),y=diffs*4,col=rejects))+
  scale_color_brewer(palette='Set2',name='Orthologue Synchronous?',labels=c('No','Yes'),
                    guide=guide_legend(override.aes = list(shape=16,size=3)))+
  geom_point(data=filter(pval_frame,ko %in% c('K04752','K00266','K01915','K00284','K00265',
                                              'K00239','K03798','K04043','K02706','K02703','K02274')),
             aes(x=factor(n_tax),y=diffs*4,fill=rejects),size=4,shape=24,color='black')+
  theme_bw()+
  xlab('# Taxa with Diel Transcript of Orthologue')+
  ylab('Mean difference in peak time [hrs]')+
geom_text_repel(
  data          = filter(pval_frame,ko %in% c('K04752','K00266','K01915','K00284','K00265',
                                              'K02274','K00239','K03798','K04043','K02706','K02703')),
  aes(x=factor(n_tax),
      y=diffs*4,
      label=str_wrap(c('sdhA (Succinate dehydrogenase)',
                       'ftsH (Cell division protease)',
                       'psbA (Photosystem II)',
                       'coxA (Cytochrome C oxidase)',
                       'dnaK (Molecular chaperone)',
                       'psbD (Photosystem II)',
                       'glnA (GS)',
                       'gltS (GOGAT)',
                       'gltD (GOGAT small chain)',
                       'gltB (GOGAT large chain)',
                       'glnK (N regulatory protein PII)'
      ),
      width=8)),
  nudge_y       = c(-4*filter(pval_frame,ko %in% c('K00239','K03798','K02703'))$diffs,
                    (10 - filter(pval_frame,ko %in% c('K04752','K00266','K01915','K00284','K00265',
                                                      'K02274','K02706','K04043'))$diffs)),
  segment.size  = 0.5,
  segment.color = "darkred",
  direction     = "x",
  force=10
)+
  ylim(c(-0.1,10))+
  theme(legend.position='bottom',
        text=element_text(size=16))+
  scale_fill_manual(values=c('#a6cee3','darkorange'))+
  guides(fill=FALSE)

ggsave(filename='figures/co_peaking_kos_lowlabels.pdf',device='pdf')


### Pointing out new things
dist_plot<-ggplot(pval_frame)+
  geom_boxplot(aes(x=factor(n_tax),y=diffs*4))+
  geom_jitter(data=filter(pval_frame,!(ko %in% c('K04752','K00266','K01915','K00284','K00265',
                                                 'K00239','K03798','K04043','K02706','K02703','K02274'))),
              aes(x=factor(n_tax),y=diffs*4,col=rejects))+
  scale_color_brewer(palette='Set2',name='Orthologue Synchronous?',labels=c('No','Yes'),
                     guide=guide_legend(override.aes = list(shape=16,size=3)))+
  geom_point(data=filter(pval_frame,ko %in% c('K04752','K00266','K01915','K00284','K00265',
                                              'K00239','K03798','K04043','K02706','K02703','K02274')),
             aes(x=factor(n_tax),y=diffs*4,fill=rejects),size=4,shape=24,color='black')+
  theme_bw()+
  xlab('# Taxa with Diel Transcript of Orthologue')+
  ylab('Mean difference in peak time [hrs]')+
  geom_text_repel(
    data          = filter(pval_frame,ko %in% c('K04752','K00266','K01915','K00284','K00265',
                                                'K02274','K00239','K03798','K04043','K02706','K02703')),
    aes(x=factor(n_tax),
        y=diffs*4,
        label=str_wrap(c('sdhA\n(Succinate dehydrogenase)',
                         'ftsH\n(Cell division protease)',
                         'psbA\n(Photosystem II)',
                         'coxA\n(Cytochrome C oxidase)',
                         'dnaK\n(Molecular chaperone)',
                         'psbD\n(Photosystem II)',
                         'glnA\n(GS)',
                         'gltS\n(GOGAT)',
                         'gltD\n(GOGAT small chain)',
                         'gltB\n(GOGAT large chain)',
                         'glnK\n(N regulatory protein PII)'
        ),
        width=8)),
    nudge_y       = c(-4*filter(pval_frame,ko %in% c('K00239','K03798','K02703'))$diffs,
                      (10 - filter(pval_frame,ko %in% c('K04752','K00266','K01915','K00284','K00265',
                                                        'K02274','K02706','K04043'))$diffs)),
    segment.size  = 0.5,
    segment.color = "darkred",
    direction     = "x",
    force=10
  )+
  ylim(c(-0.1,10))+
  theme(legend.position='bottom',
        text=element_text(size=16))+
  scale_fill_manual(values=c('#66c2a5','darkorange'))+
  guides(fill=FALSE)

ggsave(dist_plot,filename='../figures/f5a_04_17_20.pdf',device='pdf')

ggplot(pval_frame)+
  geom_boxplot(aes(x=factor(n_tax),y=diffs*4))+
  geom_jitter(data=filter(pval_frame,!(ko %in% c('K04752','K00266','K01915','K00284','K00265'))),
              aes(x=factor(n_tax),y=diffs*4,col=rejects))+
  scale_color_brewer(palette='Set2',name='Orthologue Synchronous?',labels=c('No','Yes'),
                     guide=guide_legend(override.aes = list(shape=16,size=3)))+
  geom_point(data=filter(pval_frame,ko %in% c('K04752','K00266','K01915','K00284','K00265')),
             aes(x=factor(n_tax),y=diffs*4,fill=rejects),size=4,shape=24,color='black')+
  theme_bw()+
  xlab('# Taxa with Diel Transcript of Orthologue')+
  ylab('Mean difference in peak time [hrs]')+
  geom_text_repel(
    data          = filter(pval_frame,ko %in% c('K04752','K00266','K01915','K00284','K00265')),
    aes(x=factor(n_tax),
        y=diffs*4,
        label=str_wrap(c('glnA (GS)',
                         'gltS (GOGAT)',
                         'gltD (GOGAT small chain)',
                         'gltB (GOGAT large chain)',
                         'glnK (N regulatory protein PII)'
        ),
        width=8)),
    nudge_y       = c(
                      (10 - filter(pval_frame,ko %in% c('K04752','K00266','K01915','K00284','K00265'))$diffs)),
    segment.size  = 0.5,
    segment.color = "darkred",
    direction     = "x",
    force=10
  )+
  ylim(c(-0.1,10))+
  theme(legend.position='bottom',
        text=element_text(size=16))+
  scale_fill_manual(values=c('#a6cee3','darkorange'))+
  guides(fill=FALSE)
ggsave(filename='figures/co_peaking_kos.pdf',device='pdf')

difference_with_context_data<-graphing_frame_rotated %>%
  filter(kos %in% pval_frame$ko) %>%
  filter(tax_group %in% c('eukaryotic','prokaryotic_non_photoauto','prokaryotic_photoauto')) %>%
  mutate(new_tr=floor(time_rank),
         new_group=factor(as.numeric(tax_group),levels=c(3,5,4),
                          labels=c('Eukaryote','Bacteria Autotroph','Bacteria Heterotroph')),
         new_finetax=factor(new_tax,levels=levels(new_tax)[c(1,2,3,4,6,7,8,10,12,15,
                                                             5,9,13,14,11)])) %>%
  arrange(kos,desc(new_tr)) %>%
  full_join(pval_frame,by=c('kos'='ko')) %>%
  arrange(desc(n_tax),desc(diffs))


## New strategy
## Just show all 4+ taxa ones like 50 at a time or something organized by average peak time
## Estimating best fit circle using (http://www.cs.bsu.edu/homepages/kjones/kjones/circles.pdf)
## Closed form estimators
get_best_fit_circle<-function(data){
  x<-data[,1]
  y<-data[,2]
  n<-nrow(data)
  a_hat_elements<-rep(NA,n-2)
  b_hat_elements<-rep(NA,n-2)
  for(i in 1:(n-2)){
    ## Estimating x coordinate of circle center
    omega_ijk<-(x[i]*(y[i+2]-y[i+1]))+(x[i+1]*(y[i]-y[i+2]))+(x[i+2]*(y[i+1]-y[i]))
    omega_bar_ijk<-((x[i]^2)*(y[i+2]-y[i+1]))+((x[i+1]^2)*(y[i]-y[i+2]))+((x[i+2]^2)*(y[i+1]-y[i]))
    ybar_ijk<-(y[i+1]-y[i])*(y[i+2]-y[i+1])*(y[i]-y[i+2])
    a_hat_elements[i]<-(omega_bar_ijk-ybar_ijk)/(omega_ijk)
    ## Estimating y coordinate of circle center
    z_ijk<-(y[i]*(x[i+2]-x[i+1]))+(y[i+1]*(x[i]-x[i+2]))+(y[i+2]*(x[i+1]-x[i]))
    z_bar_ijk<-((y[i]^2)*(x[i+2]-x[i+1]))+((y[i+1]^2)*(x[i]-x[i+2]))+((y[i+2]^2)*(x[i+1]-x[i]))
    xbar_ijk<-(x[i+1]-x[i])*(x[i+2]-x[i+1])*(x[i]-x[i+2])
    b_hat_elements[i]<-(z_bar_ijk-xbar_ijk)/(z_ijk)
  }
  a_hat<-sum(a_hat_elements)/(2*choose(n,3))
  b_hat<-sum(b_hat_elements)/(2*choose(n,2))
  r_hat<-sum((sqrt((x-a_hat)^2+(y-b_hat)^2)))/n
  return(c(a_hat,b_hat,r_hat))
}
pars<-get_best_fit_circle(graphing_frame_rotated[,c('x_rot','y_rot')])

cleaned_plot_frame<-difference_with_context_data %>%
  filter(tax_group %in% c('eukaryotic',
                          'prokaryotic_non_photoauto',
                          'prokaryotic_photoauto')) %>%
  filter(kos!='--') %>%
  filter(n_tax>=4)
small_diffs<-filter(cleaned_plot_frame,
       diffs %in% sort(unique(diffs))[1:5])
small_diffs2<-filter(cleaned_plot_frame,
                     kos %in% c('K02274','K00239','K03798','K04043','K02706','K02703'))
taxa_coverage<-unique(small_diffs$new_finetax)
shape_vec<-c(rep(16,sum(taxa_coverage %in% c('Bacillariophyta',
                                             'Bicosoecida','Haptophyta','Ochrophyta',
                                             'Rhodophyta','Chlorophyta','Cryptophyta',
                                             'Cercozoa','Sarcomastigophora','Dinophyta'))),
             rep(17,sum(taxa_coverage %in% c('Crocosphaera','High-light Pro'))),
             rep(15,sum(taxa_coverage %in% c('SAR11','SAR92'))),
             rep(16,sum(taxa_coverage=='other')))

small_ko_names<-c('K02935'='rplL (1.31 hr)',
                  'K01682'='acnB (0.889 hr)',
                  'K00958'='met3 (1.00 hr)',
                  'K02183'='CALM (0.00 hr)',
                  'K01899'='LSC1 (0.00 hr)',
                  'K07374'='TUBA (0.00 hr)',
                  'K11251'='H2A (0.00 hr)',
                  'K00031'='IDH1 (1.14 hr)',
                  'K02115'='ATPF1G (1.14 hr)',
                  'K00855'='prkB (0.00 hr)',
                  'K00325'='pntB (0.00 hr)')
ggplot(small_diffs)+
  geom_point(aes(x=x_rot,y=y_rot,col=new_finetax,shape=new_group),size=3)+
  ggforce::geom_circle(data=data.frame(x0=pars[1]+1,
                                       y0=pars[2],
                                       r=pars[3]),aes(x0=x0,y0=y0,r=r))+
  coord_fixed()+
  facet_wrap(~kos,ncol=6,labeller=as_labeller(small_ko_names))+
  scale_shape_discrete(name=NA,
                       guide=guide_legend(title.hjust=1.75,
                                          override.aes=list(size=3),
                                          order=1))+
  scale_color_discrete(name='Taxonomy',
                       guide=guide_legend(ncol=2,title.hjust=0.5,
                                          override.aes=list(size=3,
                                                            shape=shape_vec),
                                          order=2))+
  theme(panel.background=element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14,hjust=0.8),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        strip.background = element_blank(),
        text=element_text(size=16),
        legend.key=element_rect(fill='white',size=0.5),
        legend.key.size = unit(1, "lines"),
        legend.position='bottom',
        legend.direction='vertical',
        panel.border=element_rect(fill=NA,color='black'),
        plot.title=element_text(hjust=0.5),
        strip.text=element_text(size=11))+
  ggtitle('Minimum Average Peak Time Difference')

ggsave(filename='figures/min_diff_circles.pdf',device='pdf')
  
#### Using updated transcripts
big_diffs2<-filter(cleaned_plot_frame,kos %in% c('K04752','K00266','K01915','K00284','K00265'))
full_diffs<-rbind(big_diffs2,small_diffs2) %>%
  mutate(kos=factor(kos,levels=c('K04752',
                                 'K00266',
                                 'K01915',
                                 'K00284','K00265','K02274','K00239',
                                 'K03798','K04043','K02706','K02703')),
         new_tax2=ifelse(as.character(new_tax)=='High-light Pro','HL Prochlorococcus',as.character(new_tax)))
updated_names<-c('K04752'='glnK',
                 'K00266'='gltD',
                 'K01915'='glnA',
                 'K00284'='gltS',
                 'K00265'='gltB',
                 'K02274'='coxA',
                 'K00239'='sdhA',
                 'K03798'='ftsH',
                 'K04043'='dnaK',
                 'K02706'='psbD',
                 'K02703'='psbA')
tax_cov<-unique(full_diffs$new_tax2)
shape_vec_update<-c(rep(16,sum(tax_cov %in% c('Bacillariophyta',
                                                     'Bicosoecida','Haptophyta','Ochrophyta',
                                                     'Rhodophyta','Chlorophyta','Cryptophyta',
                                                     'Sarcomastigophora','Dinophyta'))),
                 rep(17,sum(tax_cov %in% c('Crocosphaera','HL Prochlorococcus'))),
                 rep(15,sum(tax_cov %in% c('SAR11','SAR92'))),
                 rep(16,sum(tax_cov=='other')))
top_plot<-ggplot(full_diffs)+
  geom_point(aes(x=x_rot,y=y_rot,col=new_finetax,shape=new_group),size=3)+
  ggforce::geom_circle(data=data.frame(x0=pars[1]+1,
                                       y0=pars[2],
                                       r=pars[3]),aes(x0=x0,y0=y0,r=r))+
  coord_fixed()+
  facet_wrap(~kos,
             labeller = as_labeller(updated_names),ncol=5)+
  scale_shape_discrete(name=NA,
                       guide=guide_legend(title.hjust=1.75,
                                          override.aes=list(size=3),
                                          order=1))+
  scale_color_discrete(name='Taxonomy',
                       guide=guide_legend(ncol=2,title.hjust=0.5,
                                          override.aes=list(size=3,
                                                            shape=shape_vec_update),
                                          order=2))+
  theme(panel.background=element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14,hjust=0.8),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        strip.background = element_blank(),
        text=element_text(size=16),
        legend.key=element_rect(fill='white',size=0.5),
        legend.key.size = unit(1, "lines"),
        legend.position='bottom',
        legend.direction='vertical',
        panel.border=element_rect(fill=NA,color='black'),
        plot.title=element_text(hjust=0.5))
ggsave(plot=top_plot,filename='../figures/updated_f4b_0406.pdf',device='pdf')
top_plot2<-ggplot(mutate(full_diffs,kos=factor(kos,levels=levels(kos)[c(6:11,1:5)])))+
  geom_point(aes(x=x_rot,y=y_rot,col=new_finetax,shape=new_group),size=3)+
  ggforce::geom_circle(data=data.frame(x0=pars[1]+1,
                                       y0=pars[2],
                                       r=pars[3]),aes(x0=x0,y0=y0,r=r))+
  coord_fixed()+
  facet_wrap(~kos,
             labeller = as_labeller(updated_names),ncol=6)+
  scale_shape_discrete(name=NA,
                       guide=guide_legend(title.hjust=1.75,
                                          override.aes=list(size=3),
                                          order=1))+
  scale_color_discrete(name='Taxonomy',
                       guide=guide_legend(ncol=2,title.hjust=0.5,
                                          override.aes=list(size=3,
                                                            shape=shape_vec_update),
                                          order=2))+
  theme(panel.background=element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14,hjust=0.8),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        strip.background = element_blank(),
        text=element_text(size=16),
        legend.key=element_rect(fill='white',size=0.5),
        legend.key.size = unit(1, "lines"),
        legend.position='bottom',
        legend.direction='vertical',
        panel.border=element_rect(fill=NA,color='black'),
        plot.title=element_text(hjust=0.5))

ggsave(plot=top_plot2,filename='../figures/updated_f4b_0417.pdf',device='pdf')

big_diffs<-filter(cleaned_plot_frame,
                    diffs %in% sort(unique(diffs),decreasing=TRUE)[1:5])
taxa_coverage_big<-unique(big_diffs$new_finetax)

shape_vec_big<-c(rep(16,sum(taxa_coverage_big %in% c('Bacillariophyta',
                                             'Bicosoecida','Haptophyta','Ochrophyta',
                                             'Rhodophyta','Chlorophyta','Cryptophyta',
                                             'Cercozoa','Sarcomastigophora','Dinophyta'))),
             rep(17,sum(taxa_coverage_big %in% c('Crocosphaera','High-light Pro'))),
             rep(15,sum(taxa_coverage_big %in% c('SAR11','SAR92'))),
             rep(16,sum(taxa_coverage_big=='other')))

big_diff_ko_names<-c('K00266'='gltD (6.14 hr)',
                     'K00265'='gltB (6.48 hr)',
                     'K00382'='pdhD (6.48 hr)',
                     'K03797'='ctpA (6.48 hr)',
                     'K00384'='trxB (6.67 hr)',
                     'K02112'='atpD (6.67 hr)',
                     'K00281'='gcvP (6.40 hr)',
                     'K02109'='atpF (6.40 hr)',
                     'K02992'='rpsG (6.40 hr)',
                     'K00284'='gltS (5.87 hr)')

ggplot(big_diffs)+
  geom_point(aes(x=x_rot,y=y_rot,col=new_finetax,shape=new_group),size=3)+
  ggforce::geom_circle(data=data.frame(x0=pars[1]+1,
                                       y0=pars[2],
                                       r=pars[3]),aes(x0=x0,y0=y0,r=r))+
  coord_fixed()+
  facet_wrap(~kos,ncol=5,labeller=as_labeller(big_diff_ko_names))+
  scale_shape_discrete(name=NA,
                       guide=guide_legend(title.hjust=1.75,
                                          override.aes=list(size=3),
                                          order=1))+
  scale_color_discrete(name='Taxonomy',
                       guide=guide_legend(ncol=2,title.hjust=0.5,
                                          override.aes=list(size=3,
                                                            shape=shape_vec_big),
                                          order=2))+
  theme(panel.background=element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14,hjust=0.8),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        strip.background = element_blank(),
        text=element_text(size=16),
        legend.key=element_rect(fill='white',size=0.5),
        legend.key.size = unit(1, "lines"),
        legend.position='bottom',
        legend.direction='vertical',
        panel.border=element_rect(fill=NA,color='black'),
        plot.title=element_text(hjust=0.5))+
  ggtitle('Maximum Average Peak Time Difference')
ggsave(filename='figures/max_diff_circle.pdf',device='pdf')

for_angie<-pval_frame %>%
  arrange(desc(n_tax),desc(diffs))
write.csv(for_angie,'peak_time_difference.csv')
