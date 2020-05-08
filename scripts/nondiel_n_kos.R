## Not sure where this should go/what we should do with this
kos_to_check<-read.csv('~/Downloads/kos_pathway_geneFunction_plus.csv')

## Looking up non-diel data
raw_pro_transcripts<-read.delim('../raw_data/transcripts_table_ko_byClade.txt')
raw_euk_transcripts<-read.csv('../intermediate_data/big_transcript_signals.csv')

## Looking up diel data
diel_trans<-read.csv('../intermediate_data/full_diel_transcripts.csv')

## Getting all measured transcripts
all_trans<-unique(c(paste(raw_pro_transcripts$clade_curated,raw_pro_transcripts$ko,sep='_'),
             as.character(raw_euk_transcripts$X)))
diel_ids<-as.character(diel_trans$X)

nondiel_measured<-setdiff(all_trans,diel_ids)

nondiel_ko_hits<-sapply(as.character(kos_to_check$kos),function(x) grep(x,nondiel_measured))
nondiel_ko_ids<-do.call(c,lapply(nondiel_ko_hits,function(x) as.character(nondiel_measured[x])))
nondiel_ko_taxa<-gsub('_.*$','',nondiel_ko_ids)
nondiel_ko_numbers<-gsub('^.*_','',nondiel_ko_ids)
nondiel_pros<-sapply(nondiel_ko_taxa,function(x) x%in%raw_pro_transcripts$clade_curated)
euk_taxa<-as.character(gsub('_.*$','',raw_euk_transcripts$X))
nondiel_euks<-sapply(nondiel_ko_taxa,function(x) x%in%euk_taxa)
tax_group<-rep(NA,length(nondiel_ko_ids))
tax_group[nondiel_euks]<-'Euk'
tax_group[nondiel_pros]<-'Pro'

new_tax_group<-rep(NA,length(nondiel_ko_ids))
new_tax_group[nondiel_euks]<-'Euk'
new_tax_group[nondiel_ko_taxa %in% c('Cyanothece',
                                     'Rivulariaceae',
                                     'Synechococcus',
                                     'Crocosphaera',
                                     'High-light Pro',
                                     'Low-light Pro')]<-'B_auto'
new_tax_group[is.na(new_tax_group)]<-'B_het'

output<-data.frame(nondiel_ko_ids,nondiel_ko_taxa,nondiel_ko_numbers,tax_group,new_tax_group)
write.csv(output,'measured_nondiel_nitrogen_kos_with_metab.csv',quote=FALSE)
