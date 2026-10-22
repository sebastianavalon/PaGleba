suppressPackageStartupMessages({library(vegan);library(permute);library(jsonlite)})
a<-commandArgs(TRUE);arg<-function(k)a[match(k,a)+1];jobs<-fromJSON(arg('--jobs'),simplifyVector=FALSE);OUT<-arg('--out');results<-list();dispersion<-list();alpha_all<-list()
wr<-function(x,p)write.table(x,p,sep='\t',row.names=FALSE,quote=FALSE,na='NA')
for(j in jobs){
 if(isTRUE(j$status_only)){results[[j$id]]<-data.frame(Task=j$task,Block=j$block,N=j$N,Status=j$inferential_status,Residual_df=NA_integer_,R2=NA_real_,F=NA_real_,P=NA_real_,Permutations=0,Reason=j$note);next}
 folder<-j$folder;md<-read.delim(j$metadata,check.names=FALSE,stringsAsFactors=FALSE);D0<-read.delim(file.path(folder,'native_distance.tsv'),check.names=FALSE);rownames(D0)<-D0$Sample;D<-as.dist(as.matrix(D0[,-1,drop=FALSE]));md$Group<-factor(md[[j$group]])
 covs<-unlist(j$covariates);if(length(covs))for(v in covs)md[[v]]<-factor(md[[v]])
 covs<-covs[vapply(covs,function(v)nlevels(md[[v]])>1,logical(1))];set.seed(j$seed)
 if(all(is.finite(D)) && nrow(md)>=3 && max(D)>0){
  pc<-tryCatch(wcmdscale(D,k=min(2,nrow(md)-1),eig=TRUE,add='lingoes'),error=function(e)NULL)
  if(!is.null(pc)){
   pts<-as.matrix(pc$points);if(ncol(pts)<2)pts<-cbind(pts,0);ev<-pc$eig;den<-sum(ev[ev>0]);wr(data.frame(Sample=md$sample_id,Group=md$Group,PCoA1=pts[,1],PCoA2=pts[,2]),file.path(folder,'PCoA.tsv'));wr(data.frame(Axis=seq_along(ev),Eigenvalue=ev,Positive_fraction=ifelse(ev>0,ev/den,0)),file.path(folder,'PCoA_eigenvalues.tsv'))
  }
 }
 base<-data.frame(Task=j$task,Block=j$block,N=nrow(md),Status=j$inferential_status,Residual_df=NA_integer_,R2=NA_real_,F=NA_real_,P=NA_real_,Permutations=0,Reason='')
 if(j$inferential_status!='eligible_for_inference'){base$Reason<-'Unreplicated group, missing/empty planned sample, single group or undefined distance';results[[j$id]]<-base;next}
 C<-model.matrix(if(length(covs))reformulate(c(covs,'Group'))else~Group,md);df<-nrow(md)-qr(C)$rank;base$Residual_df<-df
 if(df<2 || qr(C)$rank<ncol(C)){base$Status<-'descriptive_only';base$Reason<-'Saturated, confounded or rank-deficient design';results[[j$id]]<-base;next}
 ctrl<-how(nperm=j$permutations)
 if(!is.null(j$strata)&&nzchar(j$strata))setBlocks(ctrl)<-factor(md[[j$strata]])
 P<-tryCatch(shuffleSet(nrow(md),j$permutations,control=ctrl),error=function(e)matrix(integer(),0,nrow(md)))
 if(nrow(P)<1){base$Status<-'descriptive_only';base$Reason<-'No exchangeable permutations';results[[j$id]]<-base;next}
 # Distance variation within nuisance groups is permuted according to the supplied strata, never invented pairing.
 form<-if(length(covs))reformulate(c(covs,'Group'),response='D')else D~Group
 fit<-tryCatch(adonis2(form,data=md,permutations=P,by='margin'),error=function(e)e)
 if(inherits(fit,'error')){base$Status<-'error';base$Reason<-conditionMessage(fit);results[[j$id]]<-base;next}
 wr(data.frame(Term=rownames(fit),fit,check.names=FALSE),file.path(folder,'PERMANOVA_terms.tsv'));g<-fit['Group',];base$R2<-g$R2;base$F<-g$F;base$P<-g$`Pr(>F)`;base$Permutations<-nrow(P);base$Status<-'tested_exploratory';base$Reason<-'Exact/sampled permutation resolution limits inference for small groups';results[[j$id]]<-base
 bd<-tryCatch(betadisper(D,md$Group,type='median',bias.adjust=TRUE,add='lingoes'),error=function(e)NULL)
 if(!is.null(bd)){pd<-permutest(bd,permutations=P)$tab;dispersion[[j$id]]<-data.frame(Task=j$task,Block=j$block,P=pd[1,'Pr(>F)']);wr(data.frame(Sample=md$sample_id,Group=md$Group,Distance_to_centroid=bd$distances),file.path(folder,'PERMDISP_distances.tsv'))}
 # Freedman-Lane feature screening, excluding locally absent/constant features before CLR.
 tab<-read.delim(file.path(folder,'matrix.tsv'),check.names=FALSE);Y0<-as.matrix(tab[,-1,drop=FALSE]);keep<-colSums(Y0)>0 & apply(Y0,2,sd)>1e-12;Y0<-Y0[,keep,drop=FALSE]
 if(ncol(Y0)>0){
  Y<-if(j$kind=='taxonomy'){z<-log(Y0+1e-6);z-rowMeans(z)}else sqrt(Y0)
  C0<-model.matrix(if(length(covs))reformulate(covs)else~1,md);C1<-model.matrix(if(length(covs))reformulate(c(covs,'Group'))else~Group,md);q0<-qr(C0);q1<-qr(C1);d1<-q1$rank-q0$rank;d2<-nrow(md)-q1$rank
  statistic<-function(Z){s0<-colSums(qr.resid(q0,Z)^2);s1<-colSums(qr.resid(q1,Z)^2);pmax(0,(s0-s1)/d1)/pmax(s1/d2,1e-15)}
  obs<-statistic(Y);fitted0<-qr.fitted(q0,Y);res0<-qr.resid(q0,Y);hits<-rep(0,ncol(Y))
  for(k in seq_len(nrow(P)))hits<-hits+(statistic(fitted0+res0[P[k,],,drop=FALSE])>=obs-1e-10)
  pvals<-(hits+1)/(nrow(P)+1);means<-sapply(levels(md$Group),function(g)colMeans(Y0[md$Group==g,,drop=FALSE]));if(is.null(dim(means)))means<-matrix(means,nrow=ncol(Y0));maxgroup<-levels(md$Group)[max.col(means,ties.method='first')]
  wr(data.frame(Feature=colnames(Y),F=obs,P=pvals,Q_BH=p.adjust(pvals,'BH'),Greatest_mean_group=maxgroup,Method='Freedman_Lane_conditional_permutation'),file.path(folder,'differential_features.tsv'))
 }
 if(j$kind=='taxonomy')for(alpha_file in c('alpha_observed_counts.tsv','alpha_rarefied.tsv')){
  aa<-read.delim(file.path(folder,alpha_file),check.names=FALSE,na.strings='NA');aa<-aa[match(md$sample_id,aa$Sample),,drop=FALSE];metrics<-if(alpha_file=='alpha_observed_counts.tsv')c('Observed_Feature_Rows','Shannon','Simpson')else c('Observed_Mean','Shannon_Mean','Simpson_Mean')
  for(metric in metrics){
   yy<-aa[[metric]];ar<-data.frame(Task=j$task,Block=j$block,Alpha_Source=alpha_file,Metric=metric,N=length(yy),F=NA_real_,P=NA_real_,Status='not_estimable_missing_or_saturated')
   if(all(is.finite(yy))){
    am<-md;am$logReadDepth<-if('Reads'%in%names(aa))log1p(aa$Reads)else 0;ac<-c(covs,if(sd(am$logReadDepth)>1e-9)'logReadDepth'else character());A0<-model.matrix(if(length(ac))reformulate(ac)else~1,am);A1<-model.matrix(if(length(ac))reformulate(c(ac,'Group'))else~Group,am);q0a<-qr(A0);q1a<-qr(A1);df1<-q1a$rank-q0a$rank;df2<-nrow(am)-q1a$rank
    if(df1>0&&df2>=2&&q1a$rank==ncol(A1)&&sd(yy)>1e-12){
     stat<-function(v){s0<-sum(qr.resid(q0a,v)^2);s1<-sum(qr.resid(q1a,v)^2);max(0,(s0-s1)/df1)/max(s1/df2,1e-15)};ff<-stat(yy);f0<-qr.fitted(q0a,yy);r0<-qr.resid(q0a,yy);vals<-apply(P,1,function(ii)stat(f0+r0[ii]));ar$F<-ff;ar$P<-(1+sum(vals>=ff-1e-10))/(1+nrow(P));ar$Status<-'tested_exploratory'
    }
   }
   alpha_all[[paste(j$id,alpha_file,metric)]]<-ar
  }
 }
 # Every pair of actual groups is a prespecified secondary family; skip missing/empty/confounded subsets.
 pairrows<-list()
 for(pair in combn(levels(md$Group),2,simplify=FALSE)){
  ix<-which(md$Group%in%pair);dm<-droplevels(md[ix,,drop=FALSE]);ds<-as.dist(as.matrix(D)[ix,ix]);cv<-covs[vapply(covs,function(v)nlevels(dm[[v]])>1,logical(1))];CC<-model.matrix(if(length(cv))reformulate(c(cv,'Group'))else~Group,dm)
  row<-data.frame(Task=j$task,Block=j$block,Group_A=pair[1],Group_B=pair[2],N_A=sum(dm$Group==pair[1]),N_B=sum(dm$Group==pair[2]),R2=NA_real_,F=NA_real_,P=NA_real_,Status='descriptive_only')
  if(min(table(dm$Group))>=2 && nrow(dm)-qr(CC)$rank>=2 && qr(CC)$rank==ncol(CC)){
   ct<-how(nperm=j$permutations);if(!is.null(j$strata)&&nzchar(j$strata))setBlocks(ct)<-factor(dm[[j$strata]]);PP<-tryCatch(shuffleSet(nrow(dm),j$permutations,control=ct),error=function(e)NULL)
   if(!is.null(PP)&&nrow(PP)>0){ff<-if(length(cv))reformulate(c(cv,'Group'),response='ds')else ds~Group;z<-tryCatch(adonis2(ff,data=dm,permutations=PP,by='margin'),error=function(e)NULL);if(!is.null(z)){row$R2<-z['Group','R2'];row$F<-z['Group','F'];row$P<-z['Group','Pr(>F)'];row$Status<-'tested_exploratory'}}
  }
  pairrows[[paste(pair,collapse='_')]]<-row
 }
 if(length(pairrows))wr(do.call(rbind,pairrows),file.path(folder,'pairwise_PERMANOVA.tsv'))
}
r<-do.call(rbind,results);r$Q_BH<-p.adjust(r$P,'BH');wr(r,file.path(OUT,'PERMANOVA_summary.tsv'))
if(length(dispersion)){d<-do.call(rbind,dispersion);d$Q_BH<-p.adjust(d$P,'BH');wr(d,file.path(OUT,'PERMDISP_summary.tsv'))}
pp<-list.files(OUT,'pairwise_PERMANOVA.tsv$',recursive=TRUE,full.names=TRUE)
if(length(pp)){x<-do.call(rbind,lapply(pp,read.delim,check.names=FALSE));x$Q_BH_all_pairs<-p.adjust(x$P,'BH');wr(x,file.path(OUT,'pairwise_PERMANOVA_all_tasks.tsv'))}
if(length(alpha_all)){aa<-do.call(rbind,alpha_all);aa$Q_BH<-p.adjust(aa$P,'BH');wr(aa,file.path(OUT,'alpha_group_tests.tsv'))}
sink(file.path(OUT,'R_sessionInfo.txt'));print(sessionInfo());sink()
if(any(r$Status=='error'))stop('One or more estimable task analyses failed; inspect summary and R log')
