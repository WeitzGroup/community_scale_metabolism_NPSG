### This script is designed to generate data products dealing with data in raw units as opposed to scaled detrended data for cluster
## analyses etc. Here we also prepare additional proteomics data that was generated to complement our analysis

## PART ONE: GETTING OBSERVATIONAL CONTEXTUAL DATA (Optics-derived POC, Chl, associated PAR)
## We need to mention that we downloaded this from Station ALOHA portal and provide citations
raw_spreadsheet<-read.csv('../raw_data/underway_optics.csv')
## Okay say all measurements are 1 hr apart using NA rows as spacers
no_units<-raw_spreadsheet[-1,]
measures<-rep(1,nrow(no_units))
gaps<-which(no_units$Hour=='NaN')
measures[gaps]<-0
plot_times<-seq(0,nrow(no_units)-1)
plot_frame<-cbind(no_units,plot_times)
keeper_keys<-seq(51,108)
plot_sub<-plot_frame[keeper_keys,]
write.csv(plot_sub,'../intermediate_data/optics_data.csv',row.names=FALSE)

## Now getting units-level data for lipids
lipids_to_use<-read.csv('../intermediate_data/diel_lipid_units.csv',row.names=1)
# Subsetting for the measurements that are fully overlapping
lipid_to_use_match<-lipids_to_use[,1:15]
lipid_types<-c('TAG','UQ','PQ','PG','PE',
               'PC','DGTS','SQDG','MGDG',
               'DGDG','DGCC')
type_annotes<-c('Storage','Quinone',
                'Quinone','Cell Membrane',
                'Cell Membrane','Cell Membrane',
                'Cell Membrane','Chloroplast Membrane',
                'Chloroplast Membrane','Chloroplast Membrane',
                'Chloroplast Membrane')
lipid_annotes<-sapply(lipid_types,
                      function(x) grep(x,rownames(lipid_to_use_match)))
mol_types<-rep(NA,nrow(lipid_to_use_match))
for(i in 1:length(lipid_annotes)){
  mol_types[lipid_annotes[[i]]]<-type_annotes[i]
}
rownames(lipid_to_use_match)[which(is.na(mol_types)==TRUE)]
mol_types[which(is.na(mol_types)==TRUE)]<-c('Carotenoid',
                                            'Carotenoid',
                                            'Carotenoid',
                                            'Carotenoid',
                                            'Chloropigment',
                                            'Chloropigment',
                                            'Chloropigment',
                                            'Chloropigment',
                                            'Carotenoid',
                                            'Carotenoid',
                                            'Chloropigment')

lipid_plot_frame<-data.frame(signal=rownames(lipid_to_use_match),
                             type=mol_types,
                             lipid_to_use_match)
lipid_plot_frame<-reshape(direction='long',
                          lipid_plot_frame,
                          varying=list(3:17)) %>%
  group_by(signal) %>%
  mutate(scaled_conc=scale(July.27..2015)) %>%
  group_by(type,time) %>%
  mutate(avg_ts=mean(scaled_conc),sum_conc=sum(July.27..2015))
lipid_plot_frame$time<-factor(lipid_plot_frame$time,
                              labels=c('07_27_0200',
                                       '07_27_0600',
                                       '07_27_1000',
                                       '07_27_1400',
                                       '07_27_1800',
                                       '07_27_2200',
                                       '07_28_0200',
                                       '07_28_0600',
                                       '07_28_1000',
                                       '07_28_1400',
                                       '07_28_1800',
                                       '07_28_2200',
                                       '07_29_0200',
                                       '07_29_0600',
                                       '07_29_1000'))

## Doing the same for metabolites
plot.Conc <- read.csv("../raw_data/metabolite_concentration_estimates.csv",row.names=1,as.is=T) %>%
  filter(day.time>500,
         day.time<8201) %>%
  arrange(desc(pM.Val.mean)) %>%
  mutate(Compound.Name = ifelse(grepl("QE",Compound.Name),str_replace(Compound.Name,"_QE",""),Compound.Name),
         realDate = as.POSIXct(Date.Time,format = "%Y-%d-%b-%H%M"))
Diel.Compounds <- read.csv("../outputs_and_summaries/met_cluster_details.csv", row.names = 1, as.is = T)
Diel.Compounds <- Diel.Compounds %>%
  mutate(Compound.Name = ifelse(grepl("eE",signal),str_replace(signal,"eE","e"),signal)) %>%
  select(home_cluster, Compound.Name)
plot.Conc.Diel <- left_join(plot.Conc, Diel.Compounds) %>%
  filter(!is.na(home_cluster)) %>%
  mutate(new_compound=ifelse(Compound.Name %in% c('Serine',
                                                'Arachidonic Acid',
                                                'Homarine',
                                                'DHPS',
                                                'Glutamine',
                                                'Aspartic acid',
                                                'Betaine',
                                                'Alanine',
                                                'Adenosine',
                                                'Glucosylglycerol'),Compound.Name,'Other Compounds'))

## Mean
mean.df <- plot.Conc.Diel %>%
  select(Compound.Name, pM.Val.mean, realDate, SampID,new_compound)%>%
  filter(SampID %in% 6:20) %>%
  mutate(SampID_new=SampID-5,
         Value  = pM.Val.mean/1000) %>%
  group_by(Compound.Name) %>%
  mutate(scaled_conc=scale(Value)) %>%
  group_by(new_compound,SampID) %>%
  mutate(avg_conc=mean(scaled_conc))


## Doing the same for transcripts <5um fraction
clade_diels<-read.csv('../intermediate_data/clade_diel_data.csv')
clade_diels_matching<-clade_diels[,(5:25)[2:16]]
rownames(clade_diels_matching)<-clade_diels[,1]
clade_taxa<-gsub('_.*$','',clade_diels[,1])
clade_frame<-data.frame(taxon=clade_taxa,id=clade_diels[,1],clade_diels_matching)

## Now we need a df that has these components
clade_units_prelim<-reshape(direction='long',
                            clade_frame,varying=list(3:17))
clade_units_plotting<-clade_units_prelim %>%
  mutate(new_taxon=ifelse(taxon %in% c('High-light Pro',
                                       'Crocosphaera',
                                       'SAR11',
                                       'SAR92'),as.character(taxon),'Other <5micron')) %>%
  group_by(id) %>%
  mutate(scaled_conc=scale(S14C001)) %>%
  group_by(new_taxon,time) %>%
  mutate(avg_conc=mean(scaled_conc))


## Now for >5um transcripts
big_transcripts<-read.csv('../intermediate_data/diel_big_transcripts.csv')
big_taxa<-gsub('_.*$','',big_transcripts$X)
big_transcript_frame<-data.frame(taxon=big_taxa,id=big_transcripts$X,big_transcripts[,c(2:16)])
big_units_prelim<-reshape(direction='long',
                            big_transcript_frame,varying=list(3:17))
big_units_plotting<-big_units_prelim %>%
  mutate(new_taxon=ifelse(taxon %in% c('Haptophyta',
                                       'Bacillariophyta',
                                       'Ochrophyta',
                                       'Chlorophyta',
                                       'Dinophyta',
                                       'Cryptophyta',
                                       'Rhodophyta',
                                       'Bicosoecida',
                                       'Cercozoa',
                                       'Sarcomastigophora'),as.character(taxon),'Other >5micron')) %>%
  group_by(id) %>%
  mutate(scaled_conc=scale(V4)) %>%
  group_by(new_taxon,time) %>%
  mutate(avg_conc=mean(scaled_conc))
  