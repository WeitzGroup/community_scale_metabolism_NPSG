
## Generating Figure 1 line chart componenets
## Loading libraries
library(ggplot2)
## All data generated in script_1b_units_dataprep.R
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

lipid_nolabels<-ggplot(lipid_plot_frame)+
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
  scale_y_continuous(name=NULL,expand=c(0,0),breaks=c(-2,0,2),labels=c(-2,0,2))+
  theme(axis.text.x=element_text(angle=75,vjust=0.7),
        axis.ticks.x=element_blank(),
        strip.background=element_blank(),
        strip.placement='outside',
        strip.text.y=element_blank(),
        panel.spacing.y=unit(1.5,'lines'),
        text=element_text(size=16),
        panel.background=element_blank())
ggsave('figures/lipid_nolabel_timeseries.png',device='png',plot=lipid_nolabels,width=6,units='in')

## Metabolomics data 
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

metab_nolabel<-ggplot(mean.df)+
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
        strip.text.y=element_blank(),
        panel.spacing.y=unit(1.5,'lines'),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        text=element_text(size=16))
ggsave('figures/metabolite_nolabel_timeseries.png',device='png',plot=metab_nolabel,width=12,units='in')

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

## Need to do same thing for transcript counts
small_nolabel<-ggplot(clade_units_plotting)+
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
        strip.text.y=element_blank(),
        panel.spacing.y=unit(1.5,'lines'),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        text=element_text(size=16))+
  scale_x_continuous(name=NULL,
                     breaks=seq(2,15,by=3),
                     labels=c('6AM','6PM','6AM','6PM','6AM'),expand=c(0,0))+
  scale_y_continuous(name=NULL,breaks=c(-2,0,2),labels=c(-2,0,2),expand=c(0,0))
ggsave(plot=small_nolabel,
       filename='figures/small_transcripts_nolabel.png',
       device='png',
       width=6,
       units='in')

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

big_nolabel<-ggplot(big_units_plotting)+
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
        strip.text.y=element_blank(),
        panel.spacing.y=unit(1.5,'lines'),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        text=element_text(size=16))+
  scale_x_continuous(name=NULL,
                     breaks=seq(2,15,by=3),
                     labels=c('6AM','6PM','6AM','6PM','6AM'),expand=c(0,0))+
  scale_y_continuous(name=NULL,breaks=c(-2,0,2),labels=c(-2,0,2),expand=c(0,0))
ggsave('figures/big_transcript_nolabel_timeseries.png',device='png',plot=big_nolabel,width=12,units='in')


### New figure concept for right panel
## The idea is to sample 25 metabs; 25 lips; 10 transcripts per taxon (do them by lowest p-val??)
## Scale so max-min=1
## Rearrange in 'peak time order'
## Somehow plot for cascade
## Add POC on bottom
## Add some kind of color code for what type of thing it is
