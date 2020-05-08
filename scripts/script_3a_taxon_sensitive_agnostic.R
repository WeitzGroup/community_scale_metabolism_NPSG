## if data are taxonomically agnostic, the covariance matrix of the time series
## should be compound symmetric (ie if all time series are drawn from some 1 gene m-dimensional d'n)
## with constant \Sigma. 
## So we can test that using a permutational test described by https://digitalcommons.wayne.edu/cgi/viewcontent.cgi?referer=https://www.google.com/&httpsredir=1&article=1323&context=jmasm
## which estimates the covariance matrix, compares it to the observed, and then uses permutations
## to determine the significance of the deviation
## We estimate the covariance matrix by taking
##mean(sample variance)((1-mean(correlation))I_p-r1_p1_p')

library(permute)
## Reading in appropriate information from som clustering script
setwd('intermediate_data/')
lined_up_full_mat<-read.csv('../intermediate_data/lined_up_full_mat.csv',row.names=1)
transcipt_signals<-read.csv('../intermediate_data/full_diel_transcripts.csv',row.names=1)
ko_names<-gsub('^.*_','',rownames(transcript_signals))


## Now we examine the KOs and look at taxonomic agnostic v sensitive KOs
## Find all KOs represented in our dataset
all_kos<-unique(ko_names)
## Figure out how often each one appears
ko_counts<-sapply(all_kos,function(x) length(grep(x,rownames(lined_up_full_mat))))
## Order them by # appearances
ko_ordered<-ko_counts[order(ko_counts,decreasing = TRUE)[2:length(ko_counts)]]
## Figure out which appear at least 4 times 
more_than_3<-which(ko_ordered>=4)

## Writing a function to perform this calculation
calc_stat<-function(m){
  cov_obs<-cov(t(m))
  samp_v<-mean(diag(cov_obs))
  cor_obs<-cov2cor(cov_obs)
  mean_cor<-mean(cor_obs[upper.tri(cor_obs)])
  i_p<-diag(nrow(cov_obs))
  one_p<-rep(1,nrow(i_p))
  sigma_hat<-samp_v*((1-mean_cor)*i_p+mean_cor*one_p%*%t(one_p))
  diff_mat<-cov_obs-sigma_hat
  d_stat<-t(rep(1,nrow(m)*(nrow(m)-1)/2))%*%diff_mat[upper.tri(diff_mat)]
  return(as.numeric(d_stat))}
## Writing a function to perform permutations to get d'n of test statistic
## First we need a function to turn a vector into matrix of appropriate dimensions
vec_2_mat<-function(v){
  v<-matrix(v,nrow=15)
  return(t(v))
}
get_dn<-function(m){
  shuff_plot<-Plots(strata=gl(nrow(m),15))
  cont<-how(plots=shuff_plot,within=Within(type="series"))
  shuffles<-replicate(5000,shuffle(nrow(m)*15,cont))
  shuffled_mat<-apply(shuffles,2,function(x) m[x])
  value<-c()
  ## this could be WAY optimized but it does the job as it stands
  for(i in 1:ncol(shuffled_mat)){
    mat<-vec_2_mat(shuffled_mat[,i])
    value[i]<-calc_stat(mat)
  }
  return(value)
}

## Perform this test on all KOs with at least 4 taxa represented
p_outputs<-c()
dn_list<-list()
stat_list<-c()
k<-1
## Again can be optimized but WARNING: EXTREMELY SLOW
set.seed(4629)
for(i in 1:length(more_than_3)){
  using_matrix<-as.matrix(lined_up_full_mat[grep(names(ko_ordered)[more_than_3[i]],rownames(lined_up_full_mat)),])
  dn<-get_dn(using_matrix)
  stat_val<-calc_stat(using_matrix)
  if(stat_val<0){
    p_outputs[k]<-length(which(dn<=stat_val))/length(dn)
  } else p_outputs[k]<-length(which(dn>=stat_val))/length(dn)
  dn_list[[k]]<-dn
  stat_list[k]<-stat_val
  k<-k+1
  ## Printing progress if you're impatient like me
  print(k)
}

## Summarizing outputs
## If we just say 0.05 significance level, then we can identify agnostic and
## not KOs
tax_ag_kos_n<-ko_ordered[more_than_3[which(p_outputs>=0.05)]]
tax_sen_kos_n<-ko_ordered[more_than_3[which(p_outputs<0.05)]]
tax_ag_kos<-names(tax_ag_kos_n)
tax_sen_kos<-names(tax_sen_kos_n)
## finding which taxa those are
tax_ag_tax<-sapply(tax_ag_kos,function(x) paste(rownames(transcript_signals)[which(ko_names %in% x)],collapse=','))
tax_sen_tax<-sapply(tax_sen_kos,function(x) paste(rownames(transcript_signals)[which(ko_names %in% x)],collapse=','))
## I manually searched these on KEGG's website for annotations
narrow_diels_one<-names(ko_ordered[which(ko_ordered==1)])
narrow_diels_one_n<-rep(1,length(narrow_diels_one))
narrow_diels_twothree<-names(ko_ordered[which(ko_ordered %in% c(2,3))])
narrow_diels_twothree_n<-ko_ordered[which(ko_ordered %in% c(2,3))]
## Finding which taxa those are
one_tax<-sapply(narrow_diels_one,function(x) paste(rownames(transcript_signals)[which(ko_names %in% x)],collapse=','))
twothree_tax<-sapply(narrow_diels_twothree,function(x) paste(rownames(transcript_signals)[which(ko_names %in% x)],collapse=','))
## Writing these files to output
tax_ag_out<-cbind(tax_ag_kos,tax_ag_kos_n,tax_ag_tax)
tax_sen_out<-cbind(tax_sen_kos,tax_sen_kos_n,tax_sen_tax)
narrow_out<-cbind(c(narrow_diels_one,narrow_diels_twothree),
                  c(narrow_diels_one_n,narrow_diels_twothree_n),
                  c(one_tax,twothree_tax))
write.csv(tax_ag_out,'../outputs_and_summaries/tax_ag_ko_info.csv')
write.csv(tax_sen_out,'../outputs_and_summaries/tax_sen_ko_info.csv')
write.csv(narrow_out,'../outputs_and_summaries/tax_narrow_out.csv')
