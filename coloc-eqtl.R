
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
}

BiocManager::install("org.Hs.eg.db")


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

mythreshold <- 1000000 #1MB 1000kb #500kb #100kb

library(data.table)
library(dplyr)
library(coloc)
library(foreach)
library(org.Hs.eg.db)

setwd("E:\\AA")

genelist=read.table("id_sig.txt",header = T,sep = "\t")

genelist1=as.vector(genelist$id)
result=data.frame()


foreach(i=genelist1, .errorhandling = "pass") %do%{
df1=fread(file = paste0("eqtlclump/",i,".txt"))

head(df1)
topsnp1<-df1%>%arrange(p)
topsnp=topsnp1$SNP[1]
df1=df1[df1$eaf>0.05,] 
df1=df1[df1$p<5e-8,] 


df1=data.frame(SNP=df1$SNP,chrom=df1$chr,
               pos=df1$pos,A1=df1$A1,A2=df1$A2,beta=df1$beta,se=df1$se,MAF=df1$eaf,N=df1$n)

df2=fread("finngen_R10_E4_HYPOPIT.gz",header = T)
head(df2)
df2=data.frame(SNP=df2$rsids,chrom= df2$'#chrom',
               pos=df2$pos,A1=df2$alt, A2=df2$ref,beta=df2$beta,se=df2$sebeta,p=df2$pval) 


dfall=merge(df1,df2,by="SNP")
dfall = dfall[!duplicated(dfall$SNP),]




dfall = dfall %>% filter((A1.x==A1.y&A2.x==A2.y)|(A1.x==A2.y&A2.x==A1.y)) 
dfall = dfall %>% mutate(beta.y = ifelse(A1.x==A1.y,beta.y,-beta.y))


dfall$VAR1 = dfall$se.x^2
dfall$VAR2 =dfall$se.y^2
dfall = dfall[dfall$VAR1!=0 & dfall$VAR2!=0 ,]


cdf1 =data.frame(beta=dfall$beta.x,varbeta=dfall$VAR1,snp=dfall$SNP,MAF=dfall$MAF,N=dfall$N)
cdf2= data.frame(beta=dfall$beta.y,varbeta=dfall$VAR2,snp=dfall$SNP)

cdf1=na.omit(cdf1)
cdf2=na.omit(cdf2)

cdf1 = as.list(cdf1)
cdf2 = as.list(cdf2)



cdf1$type = "quant"#quant
cdf2$type = "cc"

res = coloc.abf(cdf1,cdf2,p1=1e-4,p2=1e-4,p12=1e-5)
color_result1=res$summary
color_result2=res$results
write.table(color_result2,file = paste0(i,"_color_result2.txt"),row.names = F,quote = F,sep = "\t")

colorsnpdf <- color_result2 %>% dplyr::arrange(desc(SNP.PP.H4))
colorsnp=colorsnpdf[1,1]

result=rbind(result,cbind(id=i,nsnps=color_result1[[1]],
                          topsnp=topsnp,
                          colorsnp=colorsnp,
                          PP.H0.abf=color_result1[[2]],
                          PP.H1.abf=color_result1[[3]],
                          PP.H2.abf=color_result1[[4]],
                          PP.H3.abf=color_result1[[5]],
                          PP.H4.abf=color_result1[[6]]))

}

write.table(result,"cloc_res_eqtl.txt",sep = "\t",row.names = F,quote = F)

