## Notes for figure scripts
library(ggplot2)
library(gridExtra)
library(grid)
library(gplots)
library(fields)
pie_fig<-ggplot(clust_plot_frame,aes(x=cluster,y=numsig,fill=tax))+geom_col()+
  #scale_fill_manual(name="Signal Type",
  #                  labels=c("Cilliate","Crocosphaera","Diatom",
  #                           "Dinoflagellate","Haptophyte","High-Light Pro",
  #                           "Lipid","Metabolite","Ochrophyte",
  #                           "Other Transcript","SAR11","SAR92"),
  #                  values=c('lightpink','darkorange','springgreen2','darkgreen',
  #                           'darkcyan','blue4','yellow','green','lightcoral','red','black','gray'))+
  theme_bw()+xlab("Cluster ID")+ylab('# Signals')+ggtitle('Cluster Signal Distribution')+
  coord_flip()+
  scale_x_discrete(breaks=c(1,3,4,2),labels=c('Morning','Afternoon','Dusk','Night'),
                   limits=rev(levels(clust_plot_frame$cluster)))

################


samplers<-som4$unit.classif
set.seed(98014)
clust_ex_mat<-matrix(data=NA,ncol=15)
colnames(clust_ex_mat)<-colnames(lined_up_full_mat)
for(i in 1:4){
  clust_sample<-sample(which(samplers==i),100,replace=FALSE)
  clust_timeseries<-lined_up_full_mat[clust_sample,]
  clust_ex_mat<-rbind(clust_ex_mat,clust_timeseries)
}
full_code_mat<-rbind(clust_ex_mat[-1,],as.matrix(som4$codes[[1]]))
code_frame<-reshape(full_code_mat,direction='long',varying=list(1:15)) %>%
  mutate(time_label=factor(x=time,
                           levels=as.character(1:15),
                           labels=colnames(full_code_mat)),
         clust_label=factor(x=rep(c(rep(1:4,each=200),1,2,3,4),15),
                            levels=c(1,3,4,2),
                            labels=c('Morning','Afternoon','Dusk','Night')),
         archetype=rep(c(rep('no',200*4),rep('yes',4)),15))
ggplot()+
  geom_rect(data=data.frame(xmin=c(1,5,11),xmax=c(2,8,14),ymin=rep(-Inf,3),ymax=rep(Inf,3)),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.5,fill='grey')+
  geom_line(data=code_frame,aes(x=time,y=`0200_727`,group=id),alpha=0.2)+
  geom_line(data=filter(code_frame,archetype=='yes'),aes(x=time,y=`0200_727`,group=id,col=clust_label),lwd=1.5)+
  facet_wrap(~clust_label,ncol=1,strip.position='left')+
  scale_x_continuous(name=NULL,
                   breaks=seq(2,15,by=3),
                   labels=c('6AM','6PM','6AM','6PM','6AM'),expand=c(0,0))+
  scale_color_manual(labels=c('Morning','Afternoon','Dusk','Night'),
                     values=rev(c("lightslategrey", "darkslategray4", "gold", "lightgoldenrod1")))+
  scale_y_continuous(limits=c(-2,2))+
  theme_bw()+
  theme(strip.background=element_blank(),
        text=element_text(size=16),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=60,vjust=0.5),
        panel.border=element_blank(),
        strip.placement='outside',
        panel.spacing.y=unit(2,'lines'),
        panel.grid=element_blank(),
        legend.position='none')

no_color_line<-ggplot()+
  geom_rect(data=data.frame(xmin=c(1,5,11),xmax=c(2,8,14),ymin=rep(-Inf,3),ymax=rep(Inf,3)),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.5,fill='grey')+
  geom_line(data=code_frame,aes(x=time,y=`0200_727`,group=id),alpha=0.2)+
  geom_line(data=filter(code_frame,archetype=='yes'),aes(x=time,y=`0200_727`,group=id),col='beige',lwd=1.5)+
  facet_wrap(~clust_label,ncol=1,strip.position='left')+
  scale_x_continuous(name=NULL,
                     breaks=seq(2,15,by=3),
                     labels=c('6AM','6PM','6AM','6PM','6AM'),expand=c(0,0))+
  scale_color_manual(labels=c('Morning','Afternoon','Dusk','Night'),
                     values=rev(c("lightslategrey", "darkslategray4", "gold", "lightgoldenrod1")))+
  scale_y_continuous(limits=c(-2,2))+
  theme_bw()+
  theme(strip.background=element_blank(),
        text=element_text(size=16),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=60,vjust=0.5),
        panel.border=element_blank(),
        strip.placement='outside',
        panel.spacing.y=unit(1,'lines'),
        panel.grid=element_blank(),
        legend.position='none')+
  coord_fixed(ratio=1.7)
ggsave(filename='../figures/som_archetypes.pdf',device='pdf',width=6,units='in')
ggsave(no_color_line,filename='../figures/som_no_color.pdf',device='pdf',width=6,units='in')

##### Options for lipid data
lip_intro<-ggplot(lipid_plot_frame)+
  geom_rect(data=data.frame(ymin=rep(-Inf,3),
                            ymax=rep(Inf,3),
                            xmin=c(1,seq(5,15,by=6)),
                            xmax=c(2,seq(8,15,by=6))),
            aes(ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),
            fill='grey',
            alpha=0.5)+
  geom_line(aes(x=time,y=scaled_conc,group=signal),alpha=0.2)+
  geom_line(aes(x=time,y=avg_ts,group=type),lwd=1.5,col='darkred')+
  facet_wrap(~type,ncol=1,strip.position='left')+
  scale_x_discrete(name=NULL,
                   breaks=levels(lipid_plot_frame$time)[seq(2,15,by=3)],
                   labels=c('6AM','6PM','6AM','6PM','6AM'),expand=c(0,0))+
  #scale_y_continuous(name=NULL,breaks=c(-1,0,1),labels=c(-1,0,1),expand=c(0,0))+
  scale_y_continuous(name=NULL,expand=c(0,0),breaks=c(-2,0,2),labels=c(-2,0,2))+
  theme(axis.text.x=element_text(angle=75,vjust=0.7),
        axis.ticks.x=element_blank(),
        strip.background=element_blank(),
        strip.placement='outside',
        strip.text.y=element_text(angle=180),
        panel.spacing.y=unit(1.5,'lines'),
        text=element_text(size=16),
        panel.background=element_blank())

ggsave(filename='../figures/lipid_timeseries.png',
      device='png',
      width=6,units='in')

ggplot(lipid_plot_frame)+
  geom_line(aes(x=time,y=sum_conc,group=type),lwd=1.5)+
  geom_line(aes(x=time,y=July.27..2015,group=signal),alpha=0.2)+
  facet_wrap(~type,ncol=1,strip.position='left',scale='free_y')+
  scale_x_discrete(name=NULL,
                   breaks=levels(lipid_plot_frame$time)[seq(2,15,by=6)],
                   labels=rep('6AM',3),expand=c(0,0))+
  geom_rect(data=data.frame(ymin=rep(-Inf,3),
                            ymax=rep(Inf,3),
                            xmin=c(1,seq(5,15,by=6)),
                            xmax=c(2,seq(8,15,by=6))),
            aes(ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),
            fill='darkgrey',
            alpha=0.3)+
  theme(axis.text.x=element_text(angle=75,vjust=0.7),
        axis.ticks.x=element_blank(),
        strip.background=element_blank(),
        strip.placement='outside',
        strip.text.y=element_text(angle=240),
        panel.spacing.y=unit(2,'lines'),
        axis.title.y=element_blank())

## We need to get an 'other' category in there
metab_intro<-ggplot(mean.df)+
  geom_rect(data=data.frame(ymin=rep(-Inf,3),
                            ymax=rep(Inf,3),
                            xmin=c(1,seq(5,15,by=6)),
                            xmax=c(2,seq(8,15,by=6))),
            aes(ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),
            fill='grey',
            alpha=0.5)+
  geom_line(aes(x=SampID_new,y=scaled_conc,group=Compound.Name),alpha=0.2)+
  geom_line(aes(x=SampID_new,y=avg_conc,group=new_compound),lwd=1.5,col='darkred')+
  facet_wrap(~new_compound,ncol=2,nrow=6,strip.position='left')+
  scale_x_continuous(name=NULL,
                   breaks=seq(2,15,by=3),
                   labels=c('6AM','6PM','6AM','6PM','6AM'),expand=c(0,0))+
  scale_y_continuous(name=NULL,breaks=c(-2,0,2),labels=c(-2,0,2),expand=c(0,0))+
  theme(axis.text.x=element_text(angle=75,vjust=0.7),
        axis.ticks.x=element_blank(),
        strip.background=element_blank(),
        strip.placement='outside',
        strip.text.y=element_text(angle=180),
        panel.spacing.y=unit(1.5,'lines'),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        text=element_text(size=16))
ggsave(filename='../figures/metabolite_timeseries.png',
       device='png',
       width=12,units='in')

## Need to do same thing for transcript counts
small_intro<-ggplot(clade_units_plotting)+
  geom_rect(data=data.frame(ymin=rep(-Inf,3),
                            ymax=rep(Inf,3),
                            xmin=c(1,seq(5,15,by=6)),
                            xmax=c(2,seq(8,15,by=6))),
            aes(ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),
            fill='grey',
            alpha=0.5)+
  geom_line(aes(x=time,y=scaled_conc,group=id),alpha=0.2)+
  geom_line(aes(x=time,y=avg_conc,group=new_taxon),lwd=1.5,col='darkred')+
  facet_wrap(~new_taxon,strip.position='left',ncol=1)+
  theme(axis.text.x=element_text(angle=75,vjust=0.7),
        axis.ticks.x=element_blank(),
        strip.background=element_blank(),
        strip.placement='outside',
        strip.text.y=element_text(angle=180),
        panel.spacing.y=unit(1.5,'lines'),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        text=element_text(size=16))+
  scale_x_continuous(name=NULL,
                     breaks=seq(2,15,by=3),
                     labels=c('6AM','6PM','6AM','6PM','6AM'),expand=c(0,0))+
  scale_y_continuous(name=NULL,breaks=c(-2,0,2),labels=c(-2,0,2),expand=c(0,0))
ggsave(filename='../figures/small_transcripts_timeseries.png',
       device='png',
       width=6,units='in')

## Need to do same thing for euk transcript counts
## Need to do same thing for transcript counts
big_intro<-ggplot(big_units_plotting)+
  geom_rect(data=data.frame(ymin=rep(-Inf,3),
                            ymax=rep(Inf,3),
                            xmin=c(1,seq(5,15,by=6)),
                            xmax=c(2,seq(8,15,by=6))),
            aes(ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax),
            fill='grey',
            alpha=0.5)+
  geom_line(aes(x=time,y=scaled_conc,group=id),alpha=0.1)+
  geom_line(aes(x=time,y=avg_conc,group=new_taxon),lwd=2,col='darkred')+
  facet_wrap(~new_taxon,strip.position='left',nrow=6)+
  theme(axis.text.x=element_text(angle=75,vjust=0.7),
        axis.ticks.x=element_blank(),
        strip.background=element_blank(),
        strip.placement='outside',
        strip.text.y=element_text(angle=180),
        panel.spacing.y=unit(1.5,'lines'),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        text=element_text(size=16))+
  scale_x_continuous(name=NULL,
                     breaks=seq(2,15,by=3),
                     labels=c('6AM','6PM','6AM','6PM','6AM'),expand=c(0,0))+
  scale_y_continuous(name=NULL,breaks=c(-2,0,2),labels=c(-2,0,2),expand=c(0,0))

ggsave(filename='../figures/big_transcripts_timeseries.png',
       device='png',
       width=12,units='in')




## Making updated ODI at good resolution we hope
pdf('../figures/distance_odi.pdf')
image(rearranged_dmat,col=colorRampPalette(c('Blue','White','Gold'))(n=100),
axes=FALSE,
main='Distance Matrix',
cex.main=2)
abline(h=length(nights)/length(som4$unit.classif),lwd=2,lty=5)
abline(h=(length(nights)+length(dusks))/length(som4$unit.classif),lwd=2,lty=5)
abline(h=(length(nights)+length(dusks)+length(afternoons))/length(som4$unit.classif),lwd=2,lty=5)
abline(v=length(nights)/length(som4$unit.classif),lwd=2,lty=5)
abline(v=(length(nights)+length(dusks))/length(som4$unit.classif),lwd=2,lty=5)
abline(v=(length(nights)+length(dusks)+length(afternoons))/length(som4$unit.classif),lwd=2,lty=5)
text(x=c(length(nights)/(2*length(som4$unit.classif)),
         (length(nights)+(0.5*length(dusks)))/length(som4$unit.classif),
         (length(nights)+length(dusks)+(0.5*length(afternoons)))/length(som4$unit.classif),
         (length(som4$unit.classif)-(0.5*length(mornings)))/length(som4$unit.classif)),
     y=c(length(nights)/(2*length(som4$unit.classif)),
         (length(nights)+(0.5*length(dusks)))/length(som4$unit.classif),
         (length(nights)+length(dusks)+(0.5*length(afternoons)))/length(som4$unit.classif),
         (length(som4$unit.classif)-(0.5*length(mornings)))/length(som4$unit.classif)),
     labels=c('Night','Dusk','Afternoon','Morning'),cex=1.5
)
image.plot(rearranged_dmat[seq(1,nrow(rearranged_dmat),by=100),
                      seq(1,nrow(rearranged_dmat),by=100)],col=colorRampPalette(c('Blue','White','Gold'))(n=100),
      axes=FALSE,
      horizontal=TRUE,
      main='Distance Matrix',
      cex=10,legend.only=TRUE,
      axis.args=list(cex.axis=2))
dev.off()




#### NMDS
labeling_function<-function(string,string2){
  current_tax_labels<-unique(as.character(graphing_frame_rotated$new_tax))
  fixed_tax_labels<-c('Metabolite (M)',
                      'Chloroplast Membrane (L)',
                      'TAG (L)',
                      'Quinone (L)',
                      'Cell Membrane (L)',
                      'Carotenoid (L)',
                      'Chloropigment (L)',
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
  filter(plotting_labels!='Other') %>%
  mutate(plotting_labels=factor(plotting_labels,levels(plotting_labels)[c(14,
                                                                          3,7,8,17,22,
                                                                          4,13,9,19,20,
                                                                          1,2,5,6,10,11,12,15,18,21,16)]))
clock_figure_r_transcript<-ggplot(filter(plot_graph_frame,big_class=='Transcript'),
                       aes(x=x_rot,y=y_rot,col=factor(floor(time_rank))))+
  ggforce::geom_circle(aes(x0=0.35,y0=0.25,r=3.5),color='black',inherit.aes=FALSE)+
  geom_point(size=2)+
  coord_fixed()+
  scale_colour_manual(name='Peak Expression',
                     labels=c('2200','0000','0200','0400','0600','0800','1000','1200','1400',
                              '1600','1800','2000')[c(seq(3,12,by=2),1)],
                     values=c('blue4','navy','royalblue3','lightblue1',
                              'yellow3','gold','darkorange','red','salmon',
                              'maroon','magenta4','dodgerblue4')[c(seq(3,12,by=2),1)],
                     guide=guide_legend(nrow=3))+
  facet_wrap(~plotting_labels,labeller = label_wrap_gen(width=6),
             ncol=5)+
  theme(panel.background=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14,hjust=0.8),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        strip.background = element_blank(),
        text=element_text(size=16),
        legend.key=element_rect(fill='white',size=0.25),
        legend.key.size = unit(1, "lines"),
        legend.position=c(0.91,0.16),
        #legend.position=c(0.7,0.12) with chem,
        strip.text.y=element_text(angle=180),
        strip.text=element_text(size=10),
        panel.border=element_rect(color='black',fill=NA))+
  guides(colour = guide_legend(override.aes = list(size=4),position='bottom',nrow=3))
ggsave(plot=clock_figure_r_transcript,filename='../figures/nmds_figure_try_6col_nochem_darkyellow.pdf',height=11,width=8,units='in',
       device='pdf')
clock_figure_r_chem<-ggplot(filter(plot_graph_frame,big_class!='Transcript') %>%
                              mutate(new_label=ifelse(plotting_labels %in% c('Chloropigment (L)','Carotenoid (L)'),'Pigment (L)',as.character(plotting_labels))),
                                  aes(x=x_rot,y=y_rot,col=factor(floor(time_rank))))+
  ggforce::geom_circle(aes(x0=0.35,y0=0.25,r=3.5),color='black',inherit.aes=FALSE)+
  geom_point(size=2)+
  coord_fixed()+
  scale_colour_manual(name='Peak Expression',
                      labels=c('2200','0000','0200','0400','0600','0800','1000','1200','1400',
                               '1600','1800','2000')[c(seq(3,12,by=2),1)],
                      values=c('blue4','navy','royalblue3','lightblue1',
                               'yellow3','gold','darkorange','red','salmon',
                               'maroon','magenta4','dodgerblue4')[c(seq(3,12,by=2),1)],
                      guide=guide_legend(nrow=3))+
  facet_wrap(~new_label,labeller = label_wrap_gen(width=6),
             ncol=5)+
  theme(panel.background=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14,hjust=0.8),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        strip.background = element_blank(),
        text=element_text(size=16),
        legend.key=element_rect(fill='white',size=0.25),
        legend.key.size = unit(1, "lines"),
        legend.position=c(0.55,0.2),
        #legend.position=c(0.7,0.12) with chem,
        strip.text.y=element_text(angle=180),
        strip.text=element_text(size=10),
        panel.border=element_rect(color='black',fill=NA))+
  guides(colour = guide_legend(override.aes = list(size=4),position='bottom',nrow=3))
ggsave(plot=clock_figure_r_chem,filename='../figures/nmds_figure_try_6col_justchem_darkeryellow_42020.pdf',width=8,units='in',
       device='pdf')
## Figuring out top 50 genes by ''abundance''
## Strategy 1: Sum read counts over all KOs for each species
summed_counts_big<-big_units_plotting %>%
  group_by(new_taxon,id) %>%
  summarize(total_count=sum(V4)) %>%
  group_by(new_taxon) %>%
  mutate(gene_rank=rank(total_count),
         ko=gsub('^.*_','',id))
summed_counts_small<-clade_units_plotting %>%
  group_by(new_taxon,id) %>%
  summarize(total_count=sum(S14C001)) %>%
  group_by(new_taxon) %>%
  mutate(gene_rank=rank(total_count),
         ko=gsub('^.*_','',id))
full_transcript_count_frame<-rbind(summed_counts_big,summed_counts_small) %>%
  group_by(new_taxon) %>%
  mutate(rel_level=total_count/sum(total_count))
full_count_summary<- full_transcript_count_frame %>%
  group_by(new_taxon) %>%
  mutate(normed_rank=gene_rank/sum(gene_rank)) %>%
  group_by(ko) %>%
  summarize(avg_rank=sum(normed_rank),
            number_taxa=n(),
            mean_rel=mean(rel_level)) %>%
  filter(number_taxa>=4) %>%
  arrange(desc(mean_rel))
kos_to_mess_with<-full_count_summary$ko[2:51]
overall_summary<-full_transcript_count_frame %>%
  group_by(new_taxon) %>%
  summarize(avg_count=mean(total_count),
            sd_count=sd(total_count),
            med_count=median(total_count))
by_gene_nmds_frame<-graphing_frame_rotated %>%
  filter(kos %in% kos_to_mess_with,
         tax_group %in% c('prokaryotic_non_photoauto',
                          'prokaryotic_photoauto',
                          'eukaryotic'))
ggplot(by_gene_nmds_frame)+
  facet_wrap(~kos,nrow=10)+
  geom_point(aes(x=x_rot,y=y_rot,col=new_tax,shape=tax_group))+
  scale_shape_discrete(name='group',
                       guide=guide_legend(nrow=3))+
  scale_color_discrete(name='',guide=guide_legend(nrow=6))+
  theme(strip.background=element_blank(),
        legend.position='bottom',
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  coord_fixed()

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
pars<-get_best_fit_circle(plot_graph_frame[,c('x_rot','y_rot')])
## We have a best fit circle... ish, roughly centered at 0 with radius 3.5, we'll see if it looks
## acceptable
tax_sen_kos<-read.csv('../outputs_and_summaries/tax_sen_ko_info.csv')
tax_sen_ko_coords<-graphing_frame_rotated %>%
  filter(kos %in% tax_sen_kos$tax_sen_kos) %>%
  #filter(tax_group %in% c('eukaryotic','prokaryotic_non_photauto','prokaryotic_photoauto')) %>%
  mutate(new_tr=floor(time_rank),
         new_group=factor(as.numeric(tax_group),levels=c(3,5,4),
                          labels=c('Eukaryote','Bacteria Autotroph','Bacteria Heterotroph')),
         new_finetax=factor(new_tax,levels=levels(new_tax)[c(1,2,3,4,6,7,8,10,12,15,
                                                             5,9,13,14,11)])) %>%
  arrange(kos,desc(new_tr))
tax_sen_ko_modes<-tax_sen_ko_coords %>%
  group_by(kos) %>%
  summarize(mode_time=pracma::Mode(new_tr),
            number_genes=n())

## Adding in gene names here we go


tax_ag_kos<-read.csv('../outputs_and_summaries/tax_ag_ko_info.csv')
tax_ag_ko_coords<-graphing_frame_rotated %>%
  filter(kos %in% tax_ag_kos$tax_ag_kos) %>%
  filter(tax_group %in% c('eukaryotic','prokaryotic_non_photoauto','prokaryotic_photoauto')) %>%
  mutate(new_tr=floor(time_rank),
         new_group=factor(as.numeric(tax_group),levels=c(3,5,4),
                          labels=c('Eukaryote','Bacteria Autotroph','Bacteria Heterotroph')),
         new_finetax=factor(new_tax,levels=levels(new_tax)[c(1,2,3,4,6,7,8,10,12,15,
                                                         5,9,13,14,11)])) %>%
  arrange(kos,desc(new_tr))
## Id mode peak times(? we'll see)
tax_ag_ko_modes<-tax_ag_ko_coords %>%
  group_by(kos) %>%
  summarize(mode_time=pracma::Mode(new_tr),
            number_genes=n())
make_facet_plot<-function(mt,df,modedf){
ggplot()+
  geom_point(data=filter(df,kos %in% subset(modedf,mode_time %in% mt)$kos) %>%
               filter(kos %in% subset(modedf,number_genes>=4)$kos),
             aes(x=x_rot,y=y_rot,col=factor(new_tr),shape=tax_group),
             size=1.75)+
  geom_circle(data=data.frame(x0=pars[1],y0=pars[2],r=pars[3]),aes(x0=x0,y0=y0,r=r))+
  facet_wrap(~kos)+
  scale_colour_manual(name='Peak Expression',
                      labels=c('2200','0000','0200','0400','0600','0800','1000','1200','1400',
                               '1600','1800','2000')[seq(1,12,by=2)],
                      values=c('blue4','navy','royalblue3','lightblue1',
                               'khaki1','gold','darkorange','red','salmon',
                               'maroon','magenta4','dodgerblue4')[seq(1,12,by=2)],
                      guide=guide_legend(nrow=3))+
    scale_shape_discrete(name='Broad Taxonomy',
                         labels=c('Eukaryote','Bacteria\nPhotoautotroph',
                                  'Bacteria\nHeterotroph'))+
    coord_fixed()+
    theme(panel.background=element_blank(),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14,hjust=0.8),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          axis.title=element_blank(),
          strip.background = element_blank(),
          text=element_text(size=16),
          legend.key=element_rect(fill='white',size=0.5),
          legend.key.size = unit(1, "lines"))
}

make_facet_plot_colortax<-function(mt,df,modedf,n_tax){
  shapeprep<-df %>%
    filter(tax_group %in% c('eukaryotic',
                            'prokaryotic_non_photoauto',
                            'prokaryotic_photoauto')) %>%
    filter(kos %in% subset(modedf,number_genes>=n_tax)$kos)
  taxa_coverage<-unique(shapeprep$new_finetax)
  shape_vec<-c(rep(16,sum(taxa_coverage %in% c('Bacillariophyta',
                                                  'Bicosoecida','Haptophyta','Ochrophyta',
                                                  'Rhodophyta','Chlorophyta','Cryptophyta',
                                                  'Cercozoa','Sarcomastigophora','Dinophyta'))),
               rep(17,sum(taxa_coverage %in% c('Crocosphaera','High-light Pro'))),
               rep(15,sum(taxa_coverage %in% c('SAR11','SAR92'))),
               rep(16,sum(taxa_coverage=='other')))
  ggplot()+
    geom_point(data=filter(df,kos %in% subset(modedf,mode_time %in% mt)$kos) %>%
                 filter(kos %in% subset(modedf,number_genes>=n_tax)$kos) %>%
                 filter(tax_group %in% c('eukaryotic','prokaryotic_non_photoauto','prokaryotic_photoauto')),
               aes(x=x_rot,
                   y=y_rot,
                   shape=new_group,
                   col=new_finetax),
               size=2.3)+
    geom_circle(data=data.frame(x0=pars[1],y0=pars[2],r=pars[3]),aes(x0=x0,y0=y0,r=r))+
    facet_wrap(~kos)+
    scale_shape_discrete(name='Broad Taxonomy',
                         guide=guide_legend(title.hjust=1.75,
                                            override.aes=list(size=3),
                                            order=1))+
    scale_color_discrete(name='Taxonomy',
                         guide=guide_legend(ncol=2,title.hjust=0.5,
                                            #override.aes=list(size=3,
                                            #                  shape=c(rep(16,9),
                                            #                          rep(17,2),
                                            #                          rep(15,2),
                                            #                          16)),
                                            #override.aes=list(size=3,
                                            #                  shape=c(rep(16,10),
                                            #                          rep(17,2),
                                            #                          16)),
                                            override.aes=list(size=3,
                                                              shape=shape_vec),
                                            order=2))+
    coord_fixed()+
    theme(panel.background=element_blank(),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14,hjust=0.8),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          axis.title=element_blank(),
          strip.background = element_blank(),
          text=element_text(size=16),
          legend.key=element_rect(fill='white',size=0.5),
          legend.key.size = unit(1, "lines"))
}
make_facet_plot(6,tax_ag_ko_coords,tax_ag_ko_modes)
make_facet_plot(5,tax_ag_ko_coords,tax_ag_ko_modes)
make_facet_plot(4,tax_ag_ko_coords,tax_ag_ko_modes)
make_facet_plot(3,tax_ag_ko_coords,tax_ag_ko_modes)
make_facet_plot(2,tax_ag_ko_coords,tax_ag_ko_modes)
make_facet_plot(1,tax_ag_ko_coords,tax_ag_ko_modes)
make_facet_plot(1:6,tax_sen_ko_coords,tax_sen_ko_modes)
make_facet_plot_colortax(6,tax_ag_ko_coords,tax_ag_ko_modes)
make_facet_plot_colortax(2,tax_ag_ko_coords,tax_ag_ko_modes)
plota<-make_facet_plot_colortax(1:6,tax_sen_ko_coords,tax_sen_ko_modes,7)+
  ggtitle('Asynchronous (p<0.05) 7+ Taxa')
ggsave(plota,filename='asynchronous_shapelegend.pdf',device='pdf')
plotb<-make_facet_plot_colortax(1:6,tax_ag_ko_coords,tax_ag_ko_modes,8)+
  ggtitle('Not Asynchronous (p>0.05) 8+ Taxa')
ggsave(plotb,filename='non_asynchronous_shapelegend.pdf',device='pdf')



### Doing a nitrogen focus
n_focus<-graphing_frame_rotated %>%
  filter(kos %in% c('K03320',
         'K00284',
         'K00265',
         'K00266',
         'K01915',
         'K04571',
         'K04572',
         'K02575',
         'K02598',
         'K02049',
         'K02050',
         'K02051',
         'K10534',
         'K11959',
         'K02587')) %>%
  filter(tax_group %in% c('eukaryotic','prokaryotic_non_photoauto','prokaryotic_photoauto')) %>%
  mutate(new_tr=floor(time_rank),
         new_group=factor(as.numeric(tax_group),levels=c(3,5,4),
                          labels=c('Eukaryote','Bacteria Autotroph','Bacteria Heterotroph')),
         new_finetax=factor(new_tax,levels=levels(new_tax)[c(1,2,3,4,6,7,8,10,12,15,
                                                             5,9,13,14,11)]),
         function_annote=factor(kos,levels=c('K03320',
                                             'K00284',
                                             'K00265',
                                             'K00266',
                                             'K01915',
                                             'K04571',
                                             'K04572',
                                             'K02575',
                                             'K02598',
                                             'K02049',
                                             'K02050',
                                             'K02051',
                                             'K10534',
                                             'K11959',
                                             'K02587'),
                                labels=str_wrap(c('NH3 Uptake',
                                         'gltS Glutamate Synthase',
                                         'gltB Glutamate Synthase',
                                         'gltD Glutamate Synthase',
                                         'glnA Glutamine Synthetase',
                                         'glnB Glutamine Synthetase Regulator PII',
                                         'glnK Glutamine Synthetase Regulator PII',
                                         'NOx Uptake',
                                         'NO3 Uptake',
                                         'NO2 Uptake',
                                         'NO2 Uptake',
                                         'NO2 Uptake',
                                         'NO3 Reductase',
                                         'Urea Uptake',
                                         'N2 Fix'),width=12)),
         function_broad=c(rep('GS',21),
                          rep('GOGAT',7),
                          rep('Uptake',7),
                          rep('Fixation'),
                          rep('Uptake',16))) %>%
  arrange(kos,desc(new_tr))


ggplot()+
  geom_point(data=n_focus %>%
               filter(tax_group %in% c('eukaryotic','prokaryotic_non_photoauto','prokaryotic_photoauto')),
             aes(x=x_rot,
                 y=y_rot,
                 shape=new_group,
                 col=new_finetax),
             size=2.3)+
  geom_circle(data=data.frame(x0=pars[1],y0=pars[2],r=pars[3]),aes(x0=x0+1,y0=y0,r=r))+
  facet_wrap(~function_annote,nrow=4)+
  scale_shape_discrete(name='Broad Taxonomy',
                       guide=guide_legend(title.hjust=1.75,
                                          override.aes=list(size=3),
                                          order=1))+
  scale_color_discrete(name='Taxonomy',
                       guide=guide_legend(ncol=2,title.hjust=0.5,
                                          override.aes=list(size=3,
                                                            shape=c(rep(16,9),
                                                                    rep(17,2),
                                                                    rep(15,2),
                                                                    16))
                                          #override.aes=list(size=3,
                                          #                  shape=c(rep(16,10),
                                          #                          rep(17,2),
                                          #                          16)),
                                          #override.aes=list(size=3,
                                          #                  shape=shape_vec),
                                          #order=2))+
                       ))+
  coord_fixed()+
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
        panel.border=element_rect(fill=NA,color='black'))
ggsave(filename='~/Desktop/n_metabolism_subplots_fixed_labels.pdf',device='pdf',width=12,height=12)


### Map
new_samp_loc<-fread('../raw_data/km1513_nav') %>%
  .[seq(2,nrow(.),by=25),] %>%
  mutate(correct_day=between(V2,208,210),
         correct_start=ifelse(V2==208,
                              ifelse(V3>=12,TRUE,FALSE),FALSE),
         correct_end=ifelse(V2==210,
                            ifelse(V3<22,TRUE,FALSE),FALSE),
         in_diel=ifelse(correct_day==TRUE,
                        ifelse(correct_start==TRUE,
                               TRUE,
                               ifelse(correct_end==TRUE,
                                      TRUE,
                                      ifelse(V2==209,
                                             TRUE,
                                             FALSE))),
                        FALSE))
A<-read.csv('~/Dropbox (GaTech)/Weitz Group Team Folder/Projects/mahalo_gt_internal/SEAFLOW_ribalet/SeaFlow-phytoplanktonAbundance-June2016.csv') %>%
  mutate(dates=as.Date(time_GMT),
         hour=as.numeric(gsub(':.*$','',gsub('^.* ','',time_GMT))),
         correct_day=between(dates,as.Date('2015-07-27'),as.Date('2015-07-29')),
         correct_start=ifelse(dates==as.Date('2015-07-27'),
                              ifelse(hour>=12,TRUE,FALSE),FALSE),
         correct_end=ifelse(dates==as.Date('2015-07-29'),
                            ifelse(hour<22,TRUE,FALSE),FALSE),
         in_diel=ifelse(correct_day==TRUE,ifelse(correct_start==TRUE,TRUE,ifelse(correct_end==TRUE,TRUE,ifelse(dates==as.Date('2015-07-28'),TRUE,FALSE))),FALSE))
wmap<-getMap(resolution='high')
wmap_d<-map_data(wmap)
small_world<-wmap_d %>%
  filter(long>=-162,lat>=19) %>%
  filter(long<=-150,lat<=25)
zoom_newfig<-ggplot()+
  geom_polygon(data=small_world,aes(x=long,y=lat,group=group),fill='lightgoldenrod1',col='darkseagreen1')+
  coord_fixed(1.3)+
  geom_path(data=new_samp_loc,
            aes(x=V9,y=V8),
            col=c(ifelse(new_samp_loc$in_diel==TRUE,'darkblue','darkorange')),
            size=1.5)+
  theme(axis.line=element_blank(),
        panel.grid.major=element_line(color='black',size=0.1),
        panel.background=element_rect(fill='azure'),
        text=element_text(size=18),
        panel.border=element_rect(color='darkgrey',fill=NA),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  geom_point(aes(x=-158,y=22.75),size=2,col='red')+
  geom_label(aes(y=22.75,x=-158.75,label='Station ALOHA'))+
  geom_segment(aes(x=-161.5,xend=-161.5+(100/102),y=19.5,yend=19.5),
               col='darkgrey')+
  geom_text(aes(x=-161,y=19.3,label='100km'),
            col='darkgrey')+
  #geom_point(data=samp_loc,aes(x=Longitude..degrees_east.,
  #                             y=Latitude..degrees_north.),
  #           size=5,shape=17)+
  #geom_point(data=samp_loc,aes(x=Longitude..degrees_east.,
  #                             y=Latitude..degrees_north.),
  #           size=5,shape=2,col='black')+
  geom_text(data=samp_loc,aes(x=Longitude..degrees_east.[1]+0.0725,y=Latitude..degrees_north.[1],
                              label='7/27 0200'),
            size=6)+
  geom_point(data=samp_loc,aes(x=Longitude..degrees_east.[15],y=Latitude..degrees_north.[15]),
             size=5)+
  geom_point(data=samp_loc,aes(x=Longitude..degrees_east.[1],y=Latitude..degrees_north.[1]),
             size=5,shape=1,stroke=1.5)+
  geom_text(data=samp_loc,aes(x=Longitude..degrees_east.[15],y=Latitude..degrees_north.[15]-0.025,
                              label='7/29 1000'),
            size=6)+
  geom_path(aes(x=c(-156.25,-156.25,-156.95,-156.95,-156.25),
                y=c(24.15,24.65,24.65,24.15,24.15)),col='navy',size=2)+
  scale_x_continuous(breaks=seq(-156.45,-156.95,by=-0.2),labels=c('156.45','156.65','156.85'),
                     limits=c(-156.95,-156.25),expand=c(0,0))+
  scale_y_continuous(limits=c(24.15,24.65),expand=c(0,0))+
  scale_color_manual(values=c('navy','navy'),guide=FALSE)
ggsave(zoom_newfig,filename='~/Desktop/fixed_map_zoom_newnav.pdf',device='pdf')

## Making full map:
full_map<-ggplot()+
  geom_polygon(data=small_world,aes(x=long,y=lat,group=group),fill='cornsilk',col='darkseagreen1')+
  coord_fixed(1.3)+
  #geom_path(data=A,aes(x=longitude_degW,
  #                     y=latitude_degN,
  #                     col=in_diel),
  #          col=c(ifelse(A$in_diel==TRUE,'darkblue','darkorange')),size=1.5)+
  geom_path(data=new_samp_loc,
            aes(x=V9,y=V8),
            col=c(ifelse(new_samp_loc$in_diel==TRUE,'darkblue','darkorange')),size=1.5)+
  theme(axis.line=element_blank(),
        panel.grid.major=element_line(color='black',size=0.1),
        panel.background=element_rect(fill='azure'),
        text=element_text(size=18),
        panel.border=element_rect(color='darkgrey',fill=NA),
        axis.ticks=element_blank())+
  geom_point(aes(x=-158,y=22.75),size=5,col='black')+
  geom_text(aes(y=22.75,x=-159.15,label='Station ALOHA'),size=5)+
  geom_segment(aes(x=-161.5,xend=-161.5+(100/102),y=19.5,yend=19.5),
               col='darkgrey')+
  geom_text(aes(x=-161,y=19.3,label='100km'),
            col='darkgrey')+
  #geom_point(data=samp_loc,aes(x=Longitude..degrees_east.,
  #                             y=Latitude..degrees_north.,
  #                             col=light_status),
  #           size=3,shape=17)+
  geom_path(aes(x=c(-156.25,-156.25,-156.95,-156.95,-156.25),
                y=c(24.15,24.65,24.65,24.15,24.15)),col='navy')+
  xlab('\u00B0W')+ylab('\u00B0N')+
  scale_x_continuous(breaks=seq(-162,-156,by=2),labels=c('162','160','158','156'))+
  scale_color_manual(values=c('navy','navy'),guide=FALSE)
ggsave(full_map,filename='~/Desktop/fullmap_newnav.pdf',device='pdf')