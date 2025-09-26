
#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

#install.packages("ggplot2")

#library
library(TwoSampleMR)
library(ggplot2)
library(foreach)

setwd("//step1")

iddf=read.table("id.txt",header =T,sep = "\t")


bioid=as.vector(iddf$id)

result=data.frame()

foreach(i=bioid, .errorhandling = "pass") %do%{
expo_rt<- read.table("expo_rt_out.txt",header = T,sep = "\t")
#expo_rt$samplesize.exposure=40234

outc_rt <- read_outcome_data(
  snps = expo_rt$SNP,
  filename = paste0("bmsnp/",i,".txt.gz"),
  sep = "\t",
  snp_col = "rsids2",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value")

harm_rt <- harmonise_data(
  exposure_dat =  expo_rt, 
  outcome_dat = outc_rt,action=1)
harm_rt<-subset(harm_rt,mr_keep==TRUE)
harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2)) /
  (2 * (harm_rt$beta.exposure^2) +
     2 * harm_rt$samplesize.exposure * harm_rt$se.exposure^2)
harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
harm_rt$meanf<- mean( harm_rt$f)
harm_rt<-harm_rt[harm_rt$f>10,]

harm_rt_steiger<-steiger_filtering(harm_rt)
harm_rt_steiger<-subset(harm_rt_steiger,steiger_dir==TRUE)

mr_result<- mr(harm_rt)
result_or=generate_odds_ratios(mr_result) 
result_or$FDR=p.adjust(result_or$pval,method = "fdr",n=338)

  filename=paste0("r_result/",i)
  dir.create(filename) 
  write.table(harm_rt, file =paste0(filename,"/harmonise.txt"),row.names = F,sep = "\t",quote = F)
  write.table(result_or[,5:ncol(result_or)],file =paste0(filename,"/OR.txt"),row.names = F,sep = "\t",quote = F)
  pleiotropy=mr_pleiotropy_test(harm_rt)
  write.table(pleiotropy,file = paste0(filename,"/pleiotropy.txt"),sep = "\t",quote = F)
  heterogeneity=mr_heterogeneity(harm_rt)
  write.table(heterogeneity,file = paste0(filename,"/heterogeneity.txt"),sep = "\t",quote = F)
  
  p1 <- mr_scatter_plot(mr_result, harm_rt)
  ggsave(p1[[1]], file=paste0(filename,"/scatter.pdf"), width=8, height=8)
  

  singlesnp_res<- mr_singlesnp(harm_rt)
  singlesnpOR=generate_odds_ratios(singlesnp_res)
  write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.txt"),row.names = F,sep = "\t",quote = F)
  p2 <- mr_forest_plot(singlesnp_res)
  ggsave(p2[[1]], file=paste0(filename,"/forest.pdf"), width=8, height=8)
  sen_res<- mr_leaveoneout(harm_rt)
  p3 <- mr_leaveoneout_plot(sen_res)
  ggsave(p3[[1]], file=paste0(filename,"/sensitivity-analysis.pdf"), width=8, height=8)
  res_single <- mr_singlesnp(harm_rt)
  p4 <- mr_funnel_plot(singlesnp_res)
  ggsave(p4[[1]], file=paste0(filename,"/funnelplot.pdf"), width=8, height=8)
  presso=run_mr_presso(harm_rt,NbDistribution = 1000)
  capture.output(presso,file = paste0(filename,"/presso.txt"))
}
}
write.table(result,"r_result.txt",sep = "\t",quote = F,row.names = F)






