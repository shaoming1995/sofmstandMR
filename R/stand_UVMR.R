#' @title 标准单变量孟德尔随机化分析
#' @param expgwas 输入暴露的GWAS摘要数据
#' @param outgwas 输入结局的GWAS摘要数据
#' @param clump_p1 输入工具变量的选择P值,默认5e-08
#' @param clump_r2 输入工具变量的选择的r2,默认0.001
#' @param clump_kb 输入工具变量的选择的距离,默认10000
#' @param pop 输入工具变量的选择的人群,默认EUR
#' @param outfile 输入分析结果的文件夹
#' @param steiger 是否进行反向过滤，默认是TURE
#' @param Fvalue 是否计算F值，默认是TURE
#' @param pt 是否进行绘图，默认是TURE
#' @export

stand_UVMR<-function(expgwas,outgwas,clump_p1=5e-08,clump_r2=0.001,clump_kb=10000,pop="EUR",outfile="MR结果",steiger=T,Fvalue=T,pt=T){
  A_temp <- c()#
  B_temp <- c()#
  C_temp <- c()#
  D_temp <- c()#
  dir.create(outfile)
  EXP<-expgwas[,c("SNP",
                  "effect_allele.exposure",
                  "other_allele.exposure",
                  "eaf.exposure",
                  "beta.exposure",
                  "se.exposure",
                  "pval.exposure",
                  "id.exposure",
                  "exposure",
                  "samplesize.exposure"
  )]
  OUT<-outgwas[,c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome",
                  "beta.outcome","se.outcome","pval.outcome","id.outcome","outcome",
                  "samplesize.outcome")]
  expiv<-subset(EXP,pval.exposure<clump_p1)
  expiv<- clump_data(expiv,clump_kb = clump_kb,clump_r2 = clump_r2,clump_p1 = 1,clump_p2 = 1,pop = pop)
  if(Fvalue==T){
    expiv$R2<-expiv$beta.exposure*expiv$beta.exposure*2*(expiv$eaf.exposure)*(1-expiv$eaf.exposure)
    expiv$Fvalue<-(expiv$samplesize.exposure-2)*expiv$R2/(1-expiv$R2)
    expiv<-subset(expiv,Fvalue>10)}
  #在结局GWAS summary中寻找与暴露对应的SNPs
  if(dim(expiv)[[1]]!=0){
    total1<-merge(OUT,expiv,by.x="SNP",by.y="SNP",all = F)
    #去除与结局有gwas显著性的SNPs以及可能重复的SNP
    total1<-subset(total1,pval.outcome>5e-08)
    total1<-total1[!duplicated(total1$SNP),]
    if(dim(total1)[[1]]!=0){
      #分别取出暴露与结局的数据
      EXP1<-total1[,c("SNP","effect_allele.exposure","other_allele.exposure", "eaf.exposure",
                      "beta.exposure","se.exposure", "pval.exposure","id.exposure","exposure",
                      "samplesize.exposure")]
      OUT1<-total1[,c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome",
                      "beta.outcome","se.outcome","pval.outcome","id.outcome","outcome",
                      "samplesize.outcome")]
      #去除回文
      dat1<-harmonise_data(exposure_dat=EXP1,outcome_dat=OUT1,action=2)
      #Steiger过滤
      if(steiger==T){
        dat1<-steiger_filtering(dat1)
        dat1<-subset(dat1,steiger_dir==TRUE)}
      res <- mr(dat1)
      mr_OR<-generate_odds_ratios(res)
      mr_OR$or<-round(mr_OR$or,3)
      mr_OR$or_lci95<-round(mr_OR$or_lci95,3)
      mr_OR$or_uci95 <- round(mr_OR$or_uci95,3)
      mr_OR$OR_CI <- paste0(mr_OR$or,"(",mr_OR$or_lci95,"-",mr_OR$or_uci95,")")
      het <- mr_heterogeneity(dat1)
      ple <- mr_pleiotropy_test(dat1)
      data_h_TableS1 <- dat1
      data_h_TableS1$R2<-data_h_TableS1$beta.exposure*data_h_TableS1$beta.exposure*2*(data_h_TableS1$eaf.exposure)*(1-data_h_TableS1$eaf.exposure)
      data_h_TableS1$Fvalue<-(data_h_TableS1$samplesize.exposure-2)*data_h_TableS1$R2/(1-data_h_TableS1$R2)
      path1<-paste0(getwd(),"/",outfile,"/MR.csv")
      path2<-paste0(getwd(),"/",outfile,"/het.csv")
      path3<-paste0(getwd(),"/",outfile,"/ple.csv")
      path4<-paste0(getwd(),"/",outfile,"/IV.csv")
      write.csv(mr_OR, path1, row.names = F)
      write.csv(het, path2, row.names = F)
      write.csv(ple, path3, row.names = F)
      write.csv(data_h_TableS1, path4, row.names = F)
      if(pt==T){
        #散点图
        path5<-paste0(getwd(),"/",outfile,"/散点图.pdf")
        pdf(path5, width = 10, height = 10)
        p1 <- mr_scatter_plot(res[1:5,], dat1)
        print(p1[[1]])
        dev.off()
        #敏感性分析Leave-one-out plot
        path6<-paste0(getwd(),"/",outfile,"/留一法.pdf")
        pdf(path6, width = 10, height = 10)
        mr_outcome_loo <- mr_leaveoneout(dat1)
        p3 <- mr_leaveoneout_plot(mr_outcome_loo)
        print(p3[[1]])
        dev.off()
        #森林图
        path7<-paste0(getwd(),"/",outfile,"/森林图.pdf")
        pdf(path7, width = 10, height = 10)
        mr_outcome_single <- mr_singlesnp(dat1)
        p2 <- mr_forest_plot(mr_outcome_single)
        print(p2[[1]])
        dev.off()
        #漏斗图
        path8<-paste0(getwd(),"/",outfile,"/漏斗图.pdf")
        pdf(path8,width = 10,height = 10)
        mr_outcome_single <- mr_singlesnp(dat1)
        p4 <- mr_funnel_plot(mr_outcome_single)
        print(p4[[1]])
        dev.off()}
    }else{cat("当前阈值可能严格，未找到工具变量")}
  }else{cat("当前阈值可能严格，未找到工具变量")}
  cat("此R包由作者邵明个人编制供MR爱好者免费开放试用一周，请关注抖音号793742981（医小研）")
}
