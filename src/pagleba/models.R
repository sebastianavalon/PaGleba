suppressPackageStartupMessages({library(jsonlite);library(igraph);library(vegan);library(permute)})
a<-commandArgs(TRUE);arg<-function(k)a[match(k,a)+1];cfg<-fromJSON(arg('--settings'),simplifyVector=FALSE);jobs<-fromJSON(cfg$jobs,simplifyVector=FALSE);OUT<-cfg$out;dir.create(OUT,recursive=TRUE,showWarnings=FALSE);allstatus<-list();concordance_all<-list();wr<-function(x,p)write.table(x,p,sep='\t',row.names=FALSE,quote=FALSE,na='NA')
if(!is.null(cfg$r_library)) .libPaths(c(cfg$r_library,.libPaths()))
for(task in unique(vapply(jobs,function(j)j$task,''))){
 J<-Filter(function(j)j$task==task&&!isTRUE(j$status_only),jobs);wanted<-if(!is.null(cfg$blocks))unlist(cfg$blocks)else vapply(J,function(j)j$block,'');J<-Filter(function(j)j$block%in%wanted,J);folder<-file.path(OUT,task);dir.create(folder,showWarnings=FALSE)
 if(length(J)<2){allstatus[[task]]<-data.frame(Task=task,N=0,Network='not_configured_fewer_than_two_blocks',PLS_PM='not_configured_fewer_than_two_blocks');next}
 ids<-Reduce(intersect,lapply(J,function(j)read.delim(j$metadata,check.names=FALSE)$sample_id));md<-read.delim(J[[1]]$metadata,check.names=FALSE);md<-md[match(ids,md$sample_id),,drop=FALSE];n<-length(ids);covs<-unlist(J[[1]]$covariates);for(v in covs)md[[v]]<-factor(md[[v]]);covs<-covs[vapply(covs,function(v)nlevels(md[[v]])>1,logical(1))];C<-model.matrix(if(length(covs))reformulate(covs)else~1,md);df<-n-qr(C)$rank-1
 if(any(vapply(J,function(j)length(j$missing_samples)+length(j$empty_samples)>0,logical(1)))){allstatus[[task]]<-data.frame(Task=task,N=n,Network='not_inferred_missing_planned_samples',PLS_PM='not_fitted_missing_planned_samples');next}
 Xs<-list();nd<-list();ind<-list();loadings<-list();pcvar<-list();node_i<-0
 for(j in J){
  tab<-read.delim(file.path(j$folder,'matrix.tsv'),check.names=FALSE);Y<-as.matrix(tab[match(ids,tab$Sample),-1,drop=FALSE]);X<-if(j$kind=='taxonomy'){z<-log(Y+1e-6);z-rowMeans(z)}else sqrt(Y);Xs[[j$block]]<-X
  pc<-tryCatch(prcomp(X),error=function(e)NULL)
  if(!is.null(pc)&&pc$sdev[1]>1e-9){ind[[j$block]]<-pc$x[,1];loadings[[j$block]]<-data.frame(Block=j$block,Feature=rownames(pc$rotation),PC1=pc$rotation[,1]);pcvar[[j$block]]<-data.frame(Block=j$block,PC1_fraction=pc$sdev[1]^2/sum(pc$sdev^2))}
  keep<-which(apply(Y,2,sd)>1e-10)
  if(j$kind=='taxonomy'){keep<-keep[colMeans(Y[,keep,drop=FALSE]>0)>=.2 & colMeans(Y[,keep,drop=FALSE])>=.0001];keep<-head(keep[order(colMeans(Y[,keep,drop=FALSE]),decreasing=TRUE)],20)}
  Z<-if(j$kind=='taxonomy')X[,keep,drop=FALSE]else log1p(100*Y[,keep,drop=FALSE]);nodeids<-paste0('N',node_i+seq_len(ncol(Z)));node_i<-node_i+ncol(Z);labels<-read.delim(file.path(j$folder,'feature_labels.tsv'),check.names=FALSE);nd[[j$block]]<-data.frame(ID=nodeids,Block=j$block,Feature=colnames(Z),Label=labels$Label[match(colnames(Z),labels$Feature)],Kind=j$kind);colnames(Z)<-nodeids;Xs[[paste0(j$block,'__network')]]<-Z
 }
 if(length(J)>=2&&n>=6&&n-qr(C)$rank>=3){
  set.seed(cfg$seed);ctrl<-how(nperm=cfg$permutations);st<-J[[1]]$strata;if(!is.null(st)&&nzchar(st))setBlocks(ctrl)<-factor(md[[st]]);PERM<-tryCatch(shuffleSet(n,cfg$permutations,control=ctrl),error=function(e)NULL)
  if(!is.null(PERM)&&nrow(PERM)>0)for(pair in combn(vapply(J,function(j)j$block,''),2,simplify=FALSE)){
   A<-Xs[[pair[1]]];B0<-Xs[[pair[2]]];RA<-qr.resid(qr(C),A);RB<-qr.resid(qr(C),B0)
   z<-tryCatch({raw<-mantel(dist(A),dist(B0),method='spearman',permutations=PERM);adj<-mantel(dist(RA),dist(RB),method='spearman',permutations=PERM);pa<-prcomp(RA)$x[,1:min(2,ncol(RA)),drop=FALSE];pb<-prcomp(RB)$x[,1:min(2,ncol(RB)),drop=FALSE];if(ncol(pa)!=ncol(pb)){k<-min(ncol(pa),ncol(pb));pa<-pa[,1:k,drop=FALSE];pb<-pb[,1:k,drop=FALSE]};pro<-protest(pa,pb,permutations=PERM);data.frame(Task=task,Block_A=pair[1],Block_B=pair[2],Mantel_raw_r=as.numeric(raw$statistic),Mantel_raw_P=raw$signif,Mantel_conditional_r=as.numeric(adj$statistic),Mantel_conditional_P=adj$signif,Procrustes_conditional_r=pro$t0,Procrustes_conditional_P=pro$signif)},error=function(e)NULL)
   if(!is.null(z))concordance_all[[paste(task,pair[1],pair[2])]]<-z
  }
 }
 status<-data.frame(Task=task,N=n,Network='not_inferred_small_n_or_residual_df',PLS_PM='not_estimable');wr(md,file.path(folder,'samples.tsv'))
 if(n>=10&&df>=8){
  nodes<-do.call(rbind,nd);Z<-do.call(cbind,Xs[paste0(vapply(J,function(j)j$block,''),'__network')]);R<-qr.resid(qr(C),apply(Z,2,rank));keep<-apply(R,2,sd)>1e-8;R<-R[,keep,drop=FALSE];nodes<-nodes[keep,,drop=FALSE]
  if(ncol(R)>=2){
   rho<-cor(R);ix<-which(upper.tri(rho),arr.ind=TRUE);rr<-rho[ix];p<-2*pt(abs(rr)*sqrt(df/pmax(1-rr^2,1e-15)),df,lower.tail=FALSE);e<-data.frame(Source=nodes$ID[ix[,1]],Target=nodes$ID[ix[,2]],BlockA=nodes$Block[ix[,1]],BlockB=nodes$Block[ix[,2]],Rho=rr,P=p,Q_BH=p.adjust(p,'BH'));e$Selected<-abs(e$Rho)>=cfg$network_rho&e$Q_BH<=cfg$network_q
   selected<-which(e$Selected)
   if(length(selected)&&cfg$network_bootstrap>0){
    M<-Z[,keep,drop=FALSE];edgeix<-ix[selected,,drop=FALSE];stratum<-if(!is.null(cfg$bootstrap_strata)&&cfg$bootstrap_strata%in%names(md))factor(md[[cfg$bootstrap_strata]])else factor(rep('all',n));parts<-split(seq_len(n),stratum);set.seed(cfg$seed)
    seeds<-sample.int(1e8,cfg$network_bootstrap);one<-function(seed){set.seed(seed);ii<-unlist(lapply(parts,function(v)sample(v,length(v),replace=TRUE)));tryCatch(cor(qr.resid(qr(C[ii,,drop=FALSE]),apply(M[ii,,drop=FALSE],2,rank)))[edgeix],error=function(e)rep(NA_real_,length(selected)))}
    BB<-if(any(lengths(parts)>1))do.call(rbind,parallel::mclapply(seeds,one,mc.cores=cfg$threads))else matrix(NA_real_,0,length(selected));e$Bootstrap_sign_agreement<-NA_real_;e$Bootstrap_CI_low<-NA_real_;e$Bootstrap_CI_high<-NA_real_;e$Bootstrap_valid_draws<-NA_integer_
    for(k in seq_along(selected)){v<-BB[,k];v<-v[is.finite(v)];j<-selected[k];e$Bootstrap_valid_draws[j]<-length(v);if(length(v)>=.8*cfg$network_bootstrap){e$Bootstrap_sign_agreement[j]<-mean(sign(v)==sign(e$Rho[j]));e$Bootstrap_CI_low[j]<-quantile(v,.025);e$Bootstrap_CI_high[j]<-quantile(v,.975)}}
   }
   wr(e,file.path(folder,'network_all_pairs.tsv'));wr(nodes,file.path(folder,'network_nodes.tsv'));es<-e[e$Selected,,drop=FALSE];wr(es,file.path(folder,'network_selected_edges.tsv'));g<-graph_from_data_frame(es[,1:2],FALSE,vertices=data.frame(name=nodes$ID,nodes[,-1]));E(g)$rho<-es$Rho;E(g)$q_BH<-es$Q_BH;g<-induced_subgraph(g,V(g)[degree(g)>0]);set.seed(cfg$seed)
   if(vcount(g)){lay<-layout_with_fr(g,weights=abs(E(g)$rho),niter=3000);V(g)$x<-lay[,1];V(g)$y<-lay[,2];V(g)$degree<-degree(g);V(g)$betweenness<-betweenness(g,normalized=TRUE,weights=NA);V(g)$module<-membership(cluster_louvain(g,weights=abs(E(g)$rho)));wr(data.frame(ID=V(g)$name,Block=V(g)$Block,Feature=V(g)$Feature,Label=V(g)$Label,X=lay[,1],Y=lay[,2],Degree=V(g)$degree,Betweenness=V(g)$betweenness,Module=V(g)$module),file.path(folder,'network_layout.tsv'))}
   write_graph(g,file.path(folder,'network.graphml'),'graphml');status$Network<-paste('exploratory',nrow(es),'selected edges')
  }else status$Network<-'not_inferred_fewer_than_two_variable_nodes'
 }
 if(length(ind)>=3&&n>=6){
  requireNamespace('plspm',quietly=TRUE)||stop('Install CRAN plspm for explicitly configured PLS-PM models')
  X<-do.call(cbind,ind);blocks<-colnames(X);edges<-cfg$paths
  if(is.null(edges)||!length(edges)){status$PLS_PM<-'not_configured_explicit_paths_required';allstatus[[task]]<-status;next}
  order<-unlist(cfg$block_order);if(is.null(order)||length(order)==0)order<-blocks
  if(!all(order%in%blocks)){status$PLS_PM<-'not_fitted_missing_requested_blocks';allstatus[[task]]<-status;next}
  X<-X[,order,drop=FALSE];blocks<-order;path<-matrix(0,length(blocks),length(blocks),dimnames=list(blocks,blocks))
  for(e in edges){stopifnot(e$from%in%blocks,e$to%in%blocks,match(e$from,blocks)<match(e$to,blocks));path[e$to,e$from]<-1};ix<-which(path>0,arr.ind=TRUE)
  fitfun<-function(z)suppressWarnings(plspm::plspm(as.data.frame(z),path,as.list(seq_along(blocks)),modes=rep('B',length(blocks)),scheme='path',maxiter=1000,tol=1e-7))
  results<-list();wr(data.frame(Sample=ids,X),file.path(folder,'PLS_local_PC1_indicators.tsv'));wr(do.call(rbind,loadings),file.path(folder,'PLS_local_PC1_loadings.tsv'));wr(do.call(rbind,pcvar),file.path(folder,'PLS_axis_variance.tsv'))
  for(mode in c('raw','adjusted')){
   A<-if(mode=='raw')X else qr.resid(qr(C),X);effective<-n-(if(mode=='raw')1 else qr(C)$rank)-max(rowSums(path));pf<-file.path(folder,paste0('PLS_',mode));dir.create(pf,showWarnings=FALSE)
   if(effective<2||any(apply(A,2,sd)<1e-9))next
   fit<-tryCatch(fitfun(A),error=function(e)NULL);if(is.null(fit))next
   orient<-function(f,anchor){sgn<-vapply(seq_along(blocks),function(i)if(cor(f$scores[,i],anchor[,i])<0)-1 else 1,1);f$path_coefs*outer(sgn,sgn)};coef<-orient(fit,A);B<-matrix(numeric(),0,nrow(ix));strata<-if(!is.null(cfg$bootstrap_strata)&&cfg$bootstrap_strata%in%names(md))factor(md[[cfg$bootstrap_strata]])else factor(rep('all',n));stratum_indices<-split(seq_len(n),strata)
   eligible<-n>=10&&effective>=8&&any(lengths(stratum_indices)>1)
   if(eligible){set.seed(cfg$seed);seeds<-sample.int(1e8,cfg$bootstrap);bootone<-function(seed){set.seed(seed);ii<-unlist(lapply(stratum_indices,function(v)sample(v,length(v),replace=TRUE)));tryCatch(orient(fitfun(A[ii,,drop=FALSE]),A[ii,,drop=FALSE])[ix],error=function(e)rep(NA_real_,nrow(ix)))};B<-do.call(rbind,parallel::mclapply(seeds,bootone,mc.cores=cfg$threads));B<-B[complete.cases(B),,drop=FALSE]}
   r<-data.frame(Path=paste(blocks[ix[,2]],blocks[ix[,1]],sep=' -> '),Estimate=coef[ix],CI_low=NA_real_,CI_high=NA_real_,Simultaneous_low=NA_real_,Simultaneous_high=NA_real_,Successful_bootstraps=nrow(B),Status='descriptive_point_estimates_only')
   if(nrow(B)>=.8*cfg$bootstrap&&nrow(B)>1){r$CI_low<-apply(B,2,quantile,.025);r$CI_high<-apply(B,2,quantile,.975);se<-apply(B,2,sd);ma<-apply(abs(sweep(sweep(B,2,r$Estimate,'-'),2,pmax(se,1e-12),'/')),1,max);crit<-quantile(ma,.95);r$Simultaneous_low<-r$Estimate-crit*se;r$Simultaneous_high<-r$Estimate+crit*se;r$Status<-'exploratory_conditional_bootstrap'}
   wr(r,file.path(pf,'paths.tsv'));wr(data.frame(Block=rownames(fit$inner_summary),fit$inner_summary),file.path(pf,'inner_summary.tsv'));saveRDS(list(fit=fit,bootstrap=B,path=path),file.path(pf,'model.rds'));wr(data.frame(Sample=ids,fit$scores),file.path(pf,'latent_scores.tsv'));results[[mode]]<-r$Status[1]
  }
  status$PLS_PM<-paste(names(results),unlist(results),collapse='; ')
 }
 allstatus[[task]]<-status
}
if(length(concordance_all)){cc<-do.call(rbind,concordance_all);cc$Mantel_raw_Q<-p.adjust(cc$Mantel_raw_P,'BH');cc$Mantel_conditional_Q<-p.adjust(cc$Mantel_conditional_P,'BH');cc$Procrustes_conditional_Q<-p.adjust(cc$Procrustes_conditional_P,'BH');wr(cc,file.path(OUT,'cross_block_concordance.tsv'))}
wr(do.call(rbind,allstatus),file.path(OUT,'MODEL_STATUS.tsv'));sink(file.path(OUT,'R_sessionInfo.txt'));print(sessionInfo());sink()
