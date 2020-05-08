library(ggplot2)
require(gridExtra)
library(grid)
library(fields)
## Commands to generate components of figure 2, data come from script_2_som_clustering.R
## Making updated ODI at good resolution we hope

## Making codebook figure for conceptual demonstration
set.seed(89781)
big_som<-som(as.matrix(lined_up_full_mat),grid=somgrid(6,4,'hexagonal'))
big_som_hc<-hcut((object.distances(big_som,'codes')),k=4)
png('../figures/codebook_figure.png',
    units='in',
    width=6,
    height=6,
    res=400)
plot(big_som,type='codes',shape='straight',main='')
add.cluster.boundaries(big_som,big_som_hc$cluster)
dev.off()

png('../figures/distance_odi.png',width=6,height=6,units='in',res=275)
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
#text(x=(length(nights)+length(dusks)+(0.5*length(afternoons)))/length(som4$unit.classif),
#y=(length(nights)+length(dusks)+(0.5*length(afternoons)))/length(som4$unit.classif),
#labels='Night',
#cex=1.2)
image.plot(rearranged_dmat[seq(1,nrow(rearranged_dmat),by=100),
                           seq(1,nrow(rearranged_dmat),by=100)],col=colorRampPalette(c('Blue','White','Gold'))(n=100),
           axes=FALSE,
           horizontal=TRUE,
           main='Distance Matrix',
           cex=10,legend.only=TRUE,
           axis.args=list(cex.axis=2))
dev.off()

png('../figures/distance_odi_som5.png',
    units='in',
    width=6,
    height=6,
    res=400)
image(dmat_5,col=colorRampPalette(c('Blue','White','Gold'))(n=100),
      axes=FALSE,
      main='Distance Matrix',
      cex.main=2)
abline(h=length(clust_5)/length(som5$unit.classif),lwd=2,lty=5)
abline(h=(length(clust_5)+length(clust_4))/length(som5$unit.classif),lwd=2,lty=5)
abline(h=(length(clust_5)+length(clust_4)+length(clust_3))/length(som5$unit.classif),lwd=2,lty=5)
abline(h=(length(clust_5)+length(clust_4)+length(clust_3)+length(clust_2))/length(som5$unit.classif),lwd=2,lty=5)
abline(v=length(clust_5)/length(som5$unit.classif),lwd=2,lty=5)
abline(v=(length(clust_5)+length(clust_4))/length(som5$unit.classif),lwd=2,lty=5)
abline(v=(length(clust_5)+length(clust_4)+length(clust_3))/length(som5$unit.classif),lwd=2,lty=5)
abline(v=(length(clust_5)+length(clust_4)+length(clust_3)+length(clust_2))/length(som5$unit.classif),lwd=2,lty=5)
image.plot(rearranged_dmat[seq(1,nrow(rearranged_dmat),by=100),
                           seq(1,nrow(rearranged_dmat),by=100)],col=colorRampPalette(c('Blue','White','Gold'))(n=100),
           axes=FALSE,
           horizontal=TRUE,
           main='Distance Matrix',
           cex=10,legend.only=TRUE,
           axis.args=list(cex.axis=2))
dev.off()


samplers<-som4$unit.classif
set.seed(98014)
clust_ex_mat<-matrix(data=NA,ncol=15)
colnames(clust_ex_mat)<-colnames(lined_up_full_mat)
for(i in 1:4){
  clust_sample<-sample(which(samplers==i),200,replace=FALSE)
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
  theme_bw()+
  theme(strip.background=element_blank(),
        text=element_text(size=16),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=60,vjust=0.5),
        panel.border=element_blank(),
        strip.placement='outside',
        panel.spacing.y=unit(2,'lines'),
        panel.grid=element_blank(),
        legend.position='none',
        strip.text=element_blank())
ggsave(filename='../figures/som_archetypes_nolabels.png',device='png',width=6,units='in')


## Bottom panel
dir<-'../intermediate_data/'
dta <- read.table(paste0(dir,'fullcluster.txt'),
                  sep="\t", header=TRUE)
dta2 <- read.table(file=paste0(dir,'fullcluster_count.txt'),
                   sep="\t", header=TRUE)
dta2$analyte<-gsub('lipid','Lipid',dta2$analyte)
dta2$analyte<-gsub('metabolite','Metabolite',dta2$analyte)
dta$analyte<-factor(dta$analyte,levels=dta2$analyte[order(dta2$logcount,decreasing=TRUE)])
dta2$analyte<-factor(dta2$analyte,levels=dta2$analyte[order(dta2$logcount,decreasing=TRUE)])
cbPalette <- c("lightslategrey", "darkslategray4", "gold", "lightgoldenrod1")
p1<-ggplot(subset(dta,analyte!='Other prokaryote'), aes(x=analyte, y=proportion, fill=time)) +
  geom_bar(stat="identity", colour="black", size=.3)+
  scale_fill_manual(values=cbPalette,
                    labels=c('Night','Dusk','Afternoon','Morning'),
                    name='Cluster')+
  theme_classic() +
  theme(text = element_text(size=20),
        axis.title.y=element_text(hjust=0.66,size=10),
        axis.title=element_blank())+
  scale_y_continuous(expand = c(0,0))+ 
  facet_grid(.~taxa,scales="free",space = "free_x",
             switch='x',
             labeller=as_labeller(c('euk'='>5 Micron',
                                    'pro'='<5 Micron',
                                    'xlipid'='Molecule',
                                    'ymetab'=''),multi_line=TRUE))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position='bottom',
        strip.placement='outside',
        strip.background.x=element_blank(),
        strip.text=element_text(angle=90))+
  theme(strip.text=element_blank(),
        axis.text.x=element_blank())+
  guides(fill=FALSE)

p1_labels<-p1 +
  theme(axis.text.x=element_text(angle=90,hjust=0.95,size=14),
        axis.text.y=element_text(size=12,angle=60,vjust=0.5))+
  guides(fill=guide_legend())
ggsave(p1_labels,filename='../figures/fingerprint_with_labels.pdf',device='pdf')

cbPalette2 <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
p2<-ggplot(subset(dta2,analyte!='Other prokaryote'),aes(x=analyte, y=logcount, fill=taxa))+
  geom_bar(stat="identity", colour="black", size=.3)+
  scale_fill_manual(values=cbPalette2,name='Data Type',
                    labels=c('Euk','Pro','Lipid','Metabolite'))+
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),strip.text.x=element_blank())+
  facet_grid(.~taxa,scales="free",space = "free_x")+
  theme(text = element_text(size=12))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  labs(y=expression(Log[10]~Count))+
  theme(axis.text=element_blank())


p2<-ggplot(dta2,aes(x=analyte, y=logcount))+
  geom_bar(stat="identity", colour="black", size=.3)+
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),strip.text.x=element_blank())+
  facet_grid(.~taxa,scales="free",space = "free_x")+
  theme(text = element_text(size=15))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  labs(y=expression(Log[10]~Count))+
  theme(#axis.text=element_blank(),
        axis.title.y=element_text(size=14))

plot_ob<-rbind(ggplotGrob(p2),ggplotGrob(p1_labels),size='first')
actual_plots<-unique(plot_ob$layout$t[grep('panel',plot_ob$layout$name)])
plot_ob$heights[actual_plots[1]]<-unit(1,'null')
plot_ob$heights[actual_plots[2]]<-unit(4,'null')

grid.newpage()
grid.draw(plot_ob)
ggsave(plot_ob,filename='../figures/new_fingerprint_plot_withlabels.pdf',
       device='pdf',units='in',width=12,height=12)
