library(RColorBrewer)
###########################################
#' Prepare dataset for HapICE calculation. 
#'
#' \code{HapICE.prepare_dataset} Read in read-region content matrix, 
#' reference information files, and output required dataset for calculation
#'
#' @param mat data.frame, original read-region content matrix.
#' @param mod_region data.frame, region position file in bed format. 
#' @param mod_info data.frame, position transfer file. 
#' @param int_pos character, intrested positions in genome position. 
#' 
#' @return A list contain: all_mat, group_info, group_name and pre_define. 
#'
#' @examples
#' \dontrun{
#'  load('demo/HapICE_demo_20201204.RData')
#'  res <- HapICE.prepare_dataset(mat=demo_mat,mod_region=mod_region,
#'         mod_info=mod_info,int_pos=int_pos)
#'  all_mat <- res$all_mat
#'  group_info <- res$group_info
#'  group_name <- res$group_name
#'  pre_define <- res$pre_define
#' }
#' @export
HapICE.prepare_dataset <- function(mat=NULL,mod_region=NULL,mod_info=NULL,int_pos=NULL){
  #if(mod_info$genome_pos[1] < mod_info$genome_pos[nrow(mod_info)]){strand='+'}else{strand='-'ult/$sample_name.res.pdf  $result/$sample_name.res.txt
  int_pos <- as.numeric(unlist(strsplit(int_pos,',')))
  int_region <- sprintf('TID_%s',int_pos)
  r1 <- mod_region; d1 <- mat; f1 <- mod_info;
  w1 <- which(int_region %in% r1$V4)
  int_pos <- int_pos[w1]; int_region <- int_region[w1]
  if(length(w1)==0) stop('no matched interesting position !')
  r1$colname <- r1$V4
  r1$colname <- gsub('-','.',r1$colname)
  rownames(r1) <- r1$colname
  r1 <- r1[which(r1$colname %in% colnames(d1)),];
  d2 <- d1[,r1$colname]
  w1 <- which(int_region %in% r1$V4);
  int_pos <- int_pos[w1]; int_region <- int_region[w1]
  ##
  int_region_len <- r1[int_region,3]-r1[int_region,2]
  ref_nt <- do.call(rbind,lapply(1:length(int_region),function(x){
    xp <- int_pos[x]
    xl <- int_region_len[x]
    x1 <- f1[which(f1$genome_pos %in% xp:(xp+xl-1)),]
    #if(strand == "+") {x1 <- f1[which(f1$genome_pos %in% (xp-xl+1):xp),]}
    #if(strand == "-") {x1 <- f1[which(f1$genome_pos %in% xp:(xp+xl-1)),]}
    x2 <- c(int_region[x],paste(x1[,4],collapse=''),paste(x1[,5],collapse=''),x1[1,3])
    names(x2) <- c('region',colnames(f1)[c(4,5,3)])
    x2
  }))
  ref_nt <- as.data.frame(ref_nt,stringsAsFactors=FALSE)
  all_mat <- d2[,int_region]
  group_info <- ref_nt[,2:(ncol(ref_nt)-1)]
  group_name <- ref_nt[,ncol(ref_nt)]
  pre_define <- colorRampPalette(brewer.pal(8,'Spectral'))(length(group_info))
  names(pre_define) <- names(group_info)
  rownames(group_info) <- colnames(all_mat)
  names(group_name) <- colnames(all_mat)
  return(list(all_mat=all_mat,group_info=group_info,group_name=group_name,pre_define=pre_define))
}
###########################################
## funcitons for HapICE.plot_readsComponent(); HapICE.plot_adjacentCombination()
get.class.color <- function(x,use_color=NULL,pre_define=NULL) {
  if(is.null(pre_define)==FALSE & is.null(names(pre_define))==TRUE){
    message('No class name for the color vector, please check and re-try !');return(FALSE);
  }
  if(is.null(use_color)==TRUE){
    use_color <- brewer.pal(9, 'Set1')
  }
  if (base::length(base::intersect(x, names(pre_define))) == 0) {
    w1 <- base::length(base::unique(x))
    if(w1 < length(use_color)){
      cc2 <- use_color[1:w1]
    }else{
      cc2 <- grDevices::colorRampPalette(use_color)(base::length(base::unique(x)))
    }
    names(cc2) <- base::unique(x)
    cc2 <- cc2[x]
  } else{
    x1 <- base::unique(x)
    x2 <- base::setdiff(x1, names(pre_define))
    cc1 <- NULL
    w1 <- base::length(x2)
    if (w1 > 0) {
      if(w1 < length(use_color)){
        cc1 <- use_color[1:w1]
      }else{
        cc1 <- grDevices::colorRampPalette(use_color)(w1)
      }
      names(cc1) <- x2
    }
    cc2 <- c(pre_define, cc1)
    cc2 <- cc2[x]
  }
  return(cc2)
}
plot_new <- function(xlim=c(),ylim=c()){
  plot(1,xlim=xlim,ylim=ylim,
       xaxt='n',yaxt='n',xlab='',ylab='',bty='n',col='white')
}
curve_sankey <- function (x0, x1, y0, y1,nsteps = 1000){
  xx <- seq(-pi/2, pi/2, length.out = nsteps)
  yy <- y0 + (y1 - y0) * (sin(xx) + 1)/2
  xx <- seq(x0, x1, length.out = nsteps)
  list(x=xx,y=yy)
}
sankey_polygon <- function(Fx0, Fx1, Fy0, Fy1,Sx0, Sx1, Sy0, Sy1,
                           nsteps = 100,col='red',border=NA,...){
  curve_pos1 <- curve_sankey(Fx0,Fx1,Fy0,Fy1)
  curve_pos2 <- curve_sankey(Sx0,Sx1,Sy0,Sy1)
  polygon(x=c(curve_pos1$x,rev(curve_pos2$x)),
          y=c(curve_pos1$y,rev(curve_pos2$y)),
          col=col,border=border,xpd=TRUE,...)
}

###########################################
#' Plot read positional-components
#'
#' \code{HapICE.plot_readsComponent} plot reads component for each position.
#'
#' @param all_mat data.frame, read-region content matrix with interested positions.
#' @param remove_char character, character for none matching sequence content, default is ".". 
#' @param group_info data.frame, sequence context at each position for parental and pseudogene. 
#' @param group_name character, genome position for each interested positions. 
#'
#' @examples
#' \dontrun{
#' HapICE.plot_readsComponent(all_mat,group_info=group_info,
#'                           remove_char='.',
#'                           group_name=group_name)
#' }
#' @export
HapICE.plot_readsComponent <- function(all_mat,group_info=NULL,
                        remove_char='.',
                        group_name=NULL){
  n <- ncol(all_mat); 
  all_r <- paste(rep(remove_char,n),collapse = '_')
  tmp1 <- apply(all_mat,1,function(x)paste(x,collapse = '_'))
  tmp2 <- sort(table(tmp1),decreasing = T)
  tmp2 <- tmp2[which(names(tmp2)!=all_r)]
  m <- length(group_info[[1]])
  par(mar=c(2,2,3,2))
  plot_new(xlim=c(0,m),ylim=c(-0.5,1.25))
  text(x=1:m,y=0,group_info[[1]],xpd=TRUE,cex=0.9)
  text(x=1:m,y=1,group_info[[2]],xpd=TRUE,cex=0.9)
  text(x=0,y=0,names(group_info)[1],pos=2,xpd=TRUE,cex=0.7)
  text(x=0,y=1,names(group_info)[2],pos=2,xpd=TRUE,cex=0.7)
  text(x=1:m,y=1.2,group_name,xpd=TRUE,adj=0,srt=90,cex=0.5)
  cc <- rev(colorRampPalette(c('black','red'))(length(tmp2)))
  dyy <- 0.8/length(tmp2)
  for(j in 1:length(tmp2)){
    x1 <- unlist(strsplit(names(tmp2)[j],'_'))
    w0 <- which(x1!='.')
    if(length(w0)>1){
      for(i in 1:(length(w0)-1)){
        w1 <- which(group_info[w0[i],] == x1[w0[i]])
        w2 <- which(group_info[w0[i+1],] == x1[w0[i+1]])
        if(length(w1)>0 & length(w2)>0){
          segments(x0=w0[i],x1=w0[i+1],
                   y0=w1-1-dyy*j,y1=w2-1-dyy*j,
                   col=adjustcolor(cc[j],0.5),
                   lwd=3*tmp2[j]/max(tmp2)+1,xpd=TRUE,
                   lty=ifelse(w0[i+1]-w0[i]==1,1,2))
          points(w0[i],w1-1-dyy*j,xpd=TRUE,cex=0.5)
          points(w0[i+1],w2-1-dyy*j,xpd=TRUE,cex=0.5)
        }
      }
    }
  }
}

###########################################
#' Plot adjacent combination of genotypes
#'
#' \code{HapICE.plot_adjacentCombination} plot reads component for each position.
#'
#' @param all_mat data.frame, read-region content matrix with interested positions.
#' @param top_n, numeric, number of top sequence context at each position. Default is 3.
#' @param each_width, numeric, bar width for each position. Default is 0.25. 
#' @param remove_char character, character for none matching sequence content, default is ".". 
#' @param count_link_thre numeric, threshold for the percentage of each intersected sequence context to the total number of valid sequence.   
#' Default is 0.1.
#' @param part_thre numeric, threshold for the percentage of each intersected sequence context to the total number of start-side valid sequence. 
#' Default is 0.1.
#' @param group_info data.frame, sequence context at each position for parental and pseudogene. 
#' @param only_use_group logical, if true, sequence not match group_info will not be plotted. Default is FALSE.
#' @param group_name character, genome position for each interested positions. 
#' @param strategy character, choose from 'percentage' and 'count'. 
#' If "percentage", the height for each position is normalized to the total number of valid sequence. 
#' If "count", the height for each position is the true read number for each sequence context. 
#' Default is 'percentage'. 
#'
#' @examples
#' \dontrun{
#' HapICE.plot_adjacentCombination(all_mat,remove_char='.',
#'                                group_info=group_info,
#'                                group_name=group_name,only_use_group=FALSE,
#'                                count_link_thre=0,part_thre=0,
#'                                top_n=2)
#' }
#' @export
HapICE.plot_adjacentCombination <- function(all_mat,top_n=3,each_width=0.25,
                       remove_char='.',
                       count_link_thre=0.1,part_thre=0.1,
                       group_info=NULL,
                       only_use_group=TRUE,
                       group_name=NULL,strategy='percentage'){
  # strategy = 'count'
  n <- ncol(all_mat)
  if(is.null(group_info)==FALSE){
    par(mar=c(5,3,ifelse(is.null(group_name)==TRUE,2,5),5))
  }else{
    par(mar=c(5,3,ifelse(is.null(group_name)==TRUE,2,5),2))
  }
  plot_new(xlim=c(1.1,n+0.25),ylim=c(0,1))
  hap_start <- list(); hap_end <- list(); all_top <- list();
  all_valid <- list(); all_topc <- list()
  if(strategy=='count'){
    w1 <- apply(all_mat,2,function(x)length(which(x!=remove_char)))
    max_w1 <- max(w1); scale_w1 <- w1/max_w1;
  }else{
    scale_w1 <- rep(1,n)
  }
  for(i in 1:n){
    w1 <- which(all_mat[,i]!=remove_char)
    all_valid[[i]] <- w1
    tmp1 <- sort(table(all_mat[w1,i]),decreasing = T)
    tmp1 <- tmp1/sum(tmp1)
    tmp2 <- as.numeric(tmp1); names(tmp2) <- names(tmp1);
    tmp1 <- tmp2
    if(is.null(group_info)==FALSE){
      w1 <- unlist(lapply(group_info,function(x)x[[i]]))
      if(only_use_group==TRUE){
        tmp2 <- tmp1[c(w1)]; names(tmp2) <- w1; tmp1 <- tmp2;
        tmp1c <- names(group_info)
        w2 <- which(tmp1>0);tmp1<-tmp1[w2];tmp1c<-tmp1c[w2] ## 2022-06-21: remove type with count 0
      }else{
        tmp2 <- tmp1[c(w1,setdiff(names(tmp1),w1))]
        names(tmp2) <- c(w1,setdiff(names(tmp1),w1))
        tmp1 <- tmp2
        tmp1c <- c(names(group_info),setdiff(names(tmp1),w1))
        w2 <- which(tmp1>0);tmp1<-tmp1[w2];tmp1c<-tmp1c[w2] ## 2022-06-21: remove type with count 0
      }
    }
    tmp1[which(is.na(tmp1)==T)] <- 0;
    if(length(tmp1)>top_n){
      tmp1 <- c(tmp1[1:top_n],other=sum(tmp1)-sum(tmp1[1:top_n]))
      top_k <- top_n+1
    }else{
      top_k <- length(tmp1)
    }
    tmp1 <- tmp1*scale_w1[i]
    tmp2 <- cumsum(rev(tmp1))
    tmp31 <- c(0,tmp2[1:(top_k-1)]); tmp32 <- c(tmp2[1:top_k])
    hap_start[[i]] <- rev(tmp31);
    hap_end[[i]]   <- rev(tmp32);
    all_top[[i]] <- names(tmp1);
    cc <- get.class.color(tmp1c,use_color=brewer.pal(8,'Pastel2'),
                                  pre_define = pre_define)
    names(cc) <- names(tmp1)
    all_topc[[i]] <- cc;
    for(j in 1:top_k){
      rect(xleft=i,xright = i+each_width,
           ybottom = tmp31[j],ytop=tmp32[j],col=cc[names(tmp32[j])],xpd=TRUE)
    }
    legend(x=i,y=-0.05,names(tmp1),cc,
           cex=0.5,xpd=T,border = NA,bty='n',yjust = 1)
  }
  if(is.null(group_name)==FALSE){
    text(1:n+each_width/2,1.05,group_name,srt=90,adj = 0,xpd=TRUE,cex=0.7)
  }
  mat1 <- as.matrix(all_mat)
  for(i in 1:(n-1)){
    start_left <- rep(0,length(all_top[[i]]))
    end_right <- rep(0,length(all_top[[i+1]]))
    for(j1 in 1:length(all_top[[i]])){
      for(j2 in 1:length(all_top[[i+1]])){
        inter_count <- length(which(mat1[,i]==all_top[[i]][j1] & 
                                      mat1[,i+1]==all_top[[i+1]][j2]))
        count_link <- inter_count/length(intersect(all_valid[[i]],all_valid[[i+1]]))
        per1_part <- inter_count/length(which(mat1[,i]==all_top[[i]][j1]))
        per1 <- inter_count/length(all_valid[[i]])*scale_w1[i]
        per2 <- inter_count/length(all_valid[[i+1]])*scale_w1[i+1]
 	#print(c(i,j1,j2,count_link,per1_part,per1,per2))
		if(is.nan(count_link) == T | is.nan(per1_part) == T){next}
        if(count_link > count_link_thre & per1_part>part_thre){
          Fy0 <- hap_end[[i]][j1]-start_left[j1]
          Fy1 <- hap_end[[i+1]][j2]-end_right[j2]
          Sy0 <- Fy0-per1
          Sy1 <- Fy1-per2
          sankey_polygon(i+each_width,i+1,Fy0,Fy1,i+each_width,i+1,Sy0,Sy1,
                         col = adjustcolor(all_topc[[i]][all_top[[i]][j1]],0.3),
                         border = 1,lwd=0.4)
          start_left[j1] <- start_left[j1]+per1
          end_right[j2] <- end_right[j2]+per2
        }
      }
    }
  }
  if(is.null(group_info)==FALSE){
    legend(x=n+0.5,y=0.5,names(group_info),
           fill=cc[1:length(names(group_info))],xpd=T,cex=0.75,
           border = NA,bty='n')
  }
  if(strategy == 'count'){
    ss <- round(seq(1,max_w1,length.out = 11))
    axis(side=2,at=c(0:10)/10,labels = ss,
         las=2,cex.axis=0.7,cex=0.7)
  }else{
    axis(side=2,at=c(0:10)/10,labels = sprintf('%s%s',c(0:10)*10,'%'),
         las=2,cex.axis=0.7,cex=0.7)
  }
}

###########################################
#' Summarize inference results. 
#'
#' \code{HapICE.summ_readsComponent} plot reads component for each position.
#'
#' @param all_mat data.frame, read-region content matrix with interested positions.
#' @param group_info data.frame, sequence context at each position for parental and pseudogene. 
#' @param group_name character, genome position for each interested positions. 
#'
#' @examples
#' \dontrun{
#'res_readsComponent <- HapICE.summ_readsComponent(all_mat,group_info,group_name)
#' @export
HapICE.summ_readsComponent <- function(all_mat,group_info,group_name){
  d4 <- sapply(1:ncol(all_mat),function(i){
    x1 <- all_mat[,i]
    w1 <- which(x1==group_info[i,1])
    w2 <- which(x1==group_info[i,2])
    x1[w1] <- sprintf("%s:%s",colnames(group_info)[1],x1[w1])
    x1[w2] <- sprintf("%s:%s",colnames(group_info)[2],x1[w2])
    x1
  })
  tmp1 <- apply(d4,1,function(x)paste(x,collapse = '_'))
  tmp2 <- sort(table(tmp1),decreasing = T)
  #tmp2 <- tmp2[which(tmp2>=nrow(d3)*0.001)]
  all_res <- list()
  for(j in 1:length(tmp2)){
    x1 <- unlist(strsplit(names(tmp2)[j],'_'))
    w0 <- which(x1!='.')
    if(length(w0)>1){
      tmp11 <- apply(d4[,w0],1,function(x)paste(x,collapse = '_'))
      r_count <- length(which(tmp11==paste(x1[w0],collapse = '_')))
      all_res[[j]] <- c(x1,r_count)
    }
  }
  all_res <- do.call(rbind,all_res)
  all_res <- all_res[order(as.numeric(all_res[,ncol(all_res)]),decreasing = T),]
  all_res <- rbind(cbind(t(group_info),'.'),all_res)
  all_res <- data.frame(all_res,stringsAsFactors = F)
  colnames(all_res)[-ncol(all_res)] <- group_name
  colnames(all_res)[ncol(all_res)] <- 'Support_Read_Number'
  rownames(all_res) <- c(colnames(group_info),
                         sprintf('Top%s',1:(nrow(all_res)-ncol(group_info)))) 
  all_res
}

###########################################
## funcitons for HapICE.infer_Haplotype()
get_valid <- function(all_mat,x,remove_character='.'){
  all_mat[which(all_mat[,x]!=remove_character),]
}
get_from <- function(group_info,each_t,x){
  x1 <- unlist(lapply(group_info,function(xx)xx[each_t]))
  x2 <- names(group_info)[which(x1==x)]
  if(length(which(x1==x))==0) x2 <- 'other'
  x2
}
compare_twoSeq <- function(x1,x2,min_overlap=2,remove_character='.'){ # x1-ref, x2-obs
  x1 <- unlist(strsplit(x1,'_'))
  x2 <- unlist(strsplit(x2,'_'))
  u1 <- which(x1!=remove_character); u2 <- which(x2!=remove_character)
  x1 <- sprintf('P%s:%s',1:length(x1),x1)[u1]
  x2 <- sprintf('P%s:%s',1:length(x2),x2)[u2]
  w0 <- intersect(x1,x2)
  w12 <- setdiff(x1,x2)
  w21 <- setdiff(x2,x1)
  if(length(w21)==0 & length(w0)>=min_overlap){
    return(1)
  }else{
    return(0)
  }
}
compare_twoSeq_count <- function(x1,x2,min_overlap=2,remove_character='.'){ # x1-ref, x2-obs
  x1 <- unlist(strsplit(x1,'_'))
  x2 <- unlist(strsplit(x2,'_'))
  u1 <- which(x1!=remove_character); u2 <- which(x2!=remove_character)
  x1 <- sprintf('P%s:%s',1:length(x1),x1)[u1]
  x2 <- sprintf('P%s:%s',1:length(x2),x2)[u2]
  w0 <- intersect(x1,x2)
  w12 <- setdiff(x1,x2)
  w21 <- setdiff(x2,x1)
  if(length(w21)==0 & length(w0)>=min_overlap){
    return(length(w0))
  }else{
    return(0)
  }
}



###########################################
#' infer for haplotype
#'
#' \code{HapICE.infer_Haplotype} perform inference for haplotypes. 
#'
#' @param all_mat data.frame, read-region content matrix with interested positions.
#' @param group_info data.frame, sequence context at each position for parental and pseudogene. 
#' @param top_each, numeric, number of top haplotypes to leave at each step. Default is 3.
#' @param use_site, character, if not NULL, only infer for the input sites. Default is NULL. 
#' @param min_overlap, numeric, minimun matching_score for a read assigned to a haplotype. Default is 2.  
#' @param check_thre, numeric, threshold for the remaining possibility when only consider top_each haplotypes at each step. 
#' If the remaining possibility is higher than check_thre, we will include more haplotypes until the remaining is possibility less than check_thre. 
#' Default is 0.05.
#' @param perm_k numeric, maximum number of permutation, default is 5000. 
#' @param perm_strategy character, choose from 'sample' and 'frequency'. 
#' If "sample", the read’s haplotype belonging is random-sampled with matching_score*haplotype frequency as probability weights.
#' If "frequency", read’s haplotype belonging is distributed according to the weight value.
#' Default is 'sample'. 
#' @param remove_char character, character for none matching sequence content, default is ".". 
#' @param min_frequency numeric, minimun frequency for the final haplotype to output. Default is 0.01.
#' @param adjust_per logical, whether to normalize output haplotype frequency to 1.Default is TRUE. 
#'
#' @examples
#' \dontrun{
#' res_trace <- HapICE.infer_Haplotype(all_mat=all_mat,
#'                                     group_info=group_info,
#'                                     top_each=3)
#' }
#' @export
HapICE.infer_Haplotype <- function(all_mat,group_info,use_site=NULL,
                         min_overlap=2,check_thre=0.05,top_each=3,
                         adjust_per=TRUE,perm_k=5000,perm_strategy='frequency',
                         remove_char='.',min_diff=1e-30,min_frequency=0.01){
  ##
  ori_group_info <- group_info;
  w1 <- apply(all_mat,1,function(x)length(setdiff(unique(x),remove_char)))
  all_mat <- all_mat[which(w1>0),]
  if(is.null(use_site)==FALSE){
    all_mat <- all_mat[,use_site]
    group_info <- lapply(group_info,function(x)x[use_site])
  }
  ##### each_t: for each Site
  ## S1, S2, ... Smax_t
  max_t <- length(group_info[[1]]) ## each_t
  result_t <- list()
  check_prob <- c()
  ## site=1
  ## get possible haplotypes
  each_t <- 1
  use_mat <- get_valid(all_mat,each_t)
  tmp1 <- table(use_mat[,each_t])/nrow(use_mat)
  ## use top
  o1 <- order(tmp1,decreasing = T)[1:min(top_each,length(tmp1))]
  check_prob[[each_t]] <- sum(tmp1)-sum(tmp1[o1])
  tmp1 <- tmp1[o1]
  ## normalize
  tmp1 <- tmp1/sum(tmp1)
  res <- tmp1
  result_t[[each_t]] <- list(seq=lapply(names(res),c),
                             trace=lapply(names(res),function(x)get_from(group_info,each_t,x)),
                             record=lapply(1:length(res),function(x)x),
                             prob=as.numeric(res),
                             check_prob=check_prob[[each_t]])
  ##
  em_part <- function(read2hap,min_diff=min_diff,perm_k=perm_k,perm_strategy=perm_strategy){
    new_prob <- colSums(read2hap)/nrow(read2hap); 
    new_prob <- new_prob/sum(new_prob)
    all_prob <- list()
    for(k in 1:perm_k){
      read2hap_new <- t(apply(read2hap,1,function(x){
        if(perm_strategy=='sample'){
          x1 <- sample(1:length(x),size=1,prob=x*new_prob)
          x[x1] <- 1; x[setdiff(1:length(x),x1)] <- 0; return(x)
        }
        if(perm_strategy=='frequency'){
          x <- x*new_prob; x <- x/sum(x);
        }
      }))
      new_prob1 <- colSums(read2hap_new)/nrow(read2hap_new)
      s1 <- sum(abs(new_prob1-new_prob))
      if(sum(abs(new_prob1-new_prob))<min_diff){
        break
      }else{
        new_prob <- new_prob1
      }
      all_prob[[k]] <- new_prob
    }
    return(list(new_prob,read2hap_new,all_prob))
  }
  ##
  for(each_t in 2:max_t){
    use_mat <- get_valid(all_mat,each_t)
    use_res <- result_t[[each_t-1]]
    tmp1 <- use_res$seq ## previous position
    tmp2 <- sort(table(use_mat[,each_t]),decreasing = T) ## real observed seq in the position
    tmp2 <- names(tmp2)
    res <- list(seq=list(),trace=list(),record=list(),prob=c(),check_prob=c(),support_read = c())
    ## get all possible haplotype --> get read to hap matrix
    read2hap <- matrix(0,nrow(all_mat),ncol=length(tmp1)*length(tmp2)); 
    rownames(read2hap) <- rownames(all_mat)
    hap_count <- 0
    for(i in 1:length(tmp1)){
      for(j in 1:length(tmp2)){
        hap_count <- hap_count+1
        ref <- c(tmp1[[i]],tmp2[j])
        obs <- apply(all_mat[,1:each_t,drop=F],1,function(x){
          compare_twoSeq_count(ref,x,min_overlap=min(min_overlap,each_t))
        })
        obs1 <- obs[which(obs>0)]
		if(length(obs1)==0){ ## for conditions with no read support 1st and 2nd, start from 2nd
        	obs <- apply(all_mat[,1:each_t,drop=F],1,function(x){
          		compare_twoSeq_count(ref,x,min_overlap=1)
       		})
        	obs1 <- obs[which(obs>0)]
		}
        read2hap[names(obs1),hap_count] <- obs1
        res$seq[[hap_count]] <- ref
        res$trace[[hap_count]] <- c(result_t[[each_t-1]]$trace[[i]],
                                    get_from(group_info,each_t,tmp2[j]))
        res$record[[hap_count]] <- result_t[[each_t-1]]$record[[i]]
      }
    }
    ## filter reads
    w1 <- which(rowSums(read2hap)>0)
    read2hap <- read2hap[w1,,drop=F]
    ## filter hap
    w1 <- which(colSums(read2hap)>0)
    read2hap <- read2hap[,w1,drop=F]
    res$seq <- res$seq[w1]; res$trace <- res$trace[w1]; 
    res$record <- res$record[w1]
    ## do em to get res$prob
    res1 <- em_part(read2hap,min_diff=min_diff,perm_k=perm_k,perm_strategy=perm_strategy)
    prob <- res1[[1]]; read2hap_new <- res1[[2]]
    res$prob <- prob; res$support_read <- round(colSums(read2hap_new))
    w1 <- order(res$prob,decreasing = T)
    res$seq <- res$seq[w1]; res$trace <- res$trace[w1]; 
    res$prob <- res$prob[w1]; res$record <- res$record[w1];
    res$support_read <- res$support_read[w1]
    if(length(res$record)>0){
    	for(i in 1:length(res$record)){
      	  res$record[[i]] <- c(res$record[[i]],i)
    	}
	}
    ## filter top
    ori_o <- order(res$prob,decreasing = T)
    prob_o <- 1-cumsum(res$prob[ori_o])
    if(length(which(prob_o<=check_thre))>0){
      w1 <- min(which(prob_o<=check_thre))
    }else{
      w1 <- length(prob_o)
    }
    o1 <- ori_o[1:max(w1,min(top_each,length(res$prob)))]
    check_prob <- sum(res$prob)-sum(res$prob[o1])
    res$seq <- res$seq[o1];res$prob <- res$prob[o1]
    res$trace <- res$trace[o1];res$record <- res$record[o1]; 
    res$support_read <- res$support_read[o1]
    res$check_prob <- check_prob
    result_t[[each_t]] <- res
  }
  ## for max_t
  w1 <- which(result_t[[max_t]]$prob >= min_frequency)
  result_t[[max_t]]$seq <- result_t[[max_t]]$seq[w1];
  result_t[[max_t]]$prob <- result_t[[max_t]]$prob[w1]
  result_t[[max_t]]$trace <- result_t[[max_t]]$trace[w1];
  result_t[[max_t]]$record <- result_t[[max_t]]$record[w1]; 
  result_t[[max_t]]$support_read <- result_t[[max_t]]$support_read[w1]
  result_t[[max_t]]$check_prob <- 1-sum(result_t[[max_t]]$prob)
  if(adjust_per==TRUE){
    result_t[[max_t]]$prob <- result_t[[max_t]]$prob/sum(result_t[[max_t]]$prob)
  }
  #if(is.null(res$record[[1]])==TRUE){ ## if no reads support 1st and 2nd position,remove 1st position
  #	message('No reads support the 1st and 2nd position, remove the 1st position and re-run!')
  #  if(is.null(use_site)==TRUE) use_site <- colnames(all_mat)[2:ncol(all_mat)]
  #  if(is.null(use_site)==FALSE) use_site <- setdiff(use_site,colnames(all_mat)[1])
  #	result_t <- HapICE.infer_Haplotype(all_mat=all_mat[,2:ncol(all_mat),drop=F],
  #                       group_info=ori_group_info[2:nrow(group_info),,drop=F],
  #						 use_site=use_site,
  #                       min_overlap=min_overlap,check_thre=check_thre,top_each=top_each,
  #                       adjust_per=adjust_per,perm_k=perm_k,perm_strategy=perm_strategy,
  #                       remove_char=remove_char,min_diff=min_diff,min_frequency=min_frequency)
  #  return(result_t)
  #}
  ## output
  names(result_t) <- rownames(ori_group_info)
  message(sprintf('Output for position: %s',paste(names(result_t),collapse=',')))
  return(result_t)
}


###########################################
#' Plot for inference results
#'
#' \code{HapICE.plot_inferHaplotype} plot the inference results. 
#'
#' @param result_t list, output object from HapICE.infer_Haplotype(). 
#' @param prob_thre, numeric, threshold for the haplotype frequency to show. Default is 0.01.
#' @param group_name character, genome position for each interested positions. 
#' @param consider_site_per logical, whether to display each site according to the possibility at each step. 
#' Default is FALSE.
#'
#' @examples
#' \dontrun{
#'       HapICE.plot_inferHaplotype(res_trace,prob_thre=0.01,
#'                           group_name=group_name)
#' }
#' @export
HapICE.plot_inferHaplotype <- function(result_t,prob_thre=0.01,group_name=NULL,
                              consider_site_per=FALSE){
  par(mar=c(2,2,5,8))
  group_name <- group_name[names(result_t)]
  if(consider_site_per==TRUE){
    tmp1 <- lapply(result_t, function(x){
      w1 <- which(x$prob>=prob_thre)
      list(seq=x$seq[w1],trace=x$trace[w1],prob=x$prob[w1],record=x$record[w1])
    })
    use_result_t <- result_t
    use_result_t[max(length(result_t))] <- tmp1[max(length(result_t))]
    max_x <- length(use_result_t)
    max_y <- max(unlist(lapply(use_result_t,function(x)length(x$prob))))
    all_m <- unique(unlist(lapply(use_result_t,function(x)unlist(x$trace))))
    cc <- get.class.color(all_m,use_color=brewer.pal(8,'Pastel2'),
                          pre_define = pre_define)
    # points
    plot_new(xlim = c(1,max_x),ylim=c(1,max_y))
    for(i in 2:max_x){
      each_r <- use_result_t[[i]]
      each_t <- each_r$trace
      each_d <- each_r$record
      for(j in 1:length(each_t)){
        for(j1 in 1:(i-1)){
          use_t <- each_t[[j]][j1:(j1+1)]
          segments(x0=j1,x1=j1+1,
                   y0=1+max_y-each_d[[j]][j1],
                   y1=1+max_y-each_d[[j]][j1+1],
                   lwd=3+5*each_r$prob[j],
                   col=adjustcolor('dark grey',0.5))
        }
      }
    }
    ##
    for(i in 1:max_x){
      each_r <- use_result_t[[i]]
      each_b <- each_r$prob
      each_p <- unlist(lapply(each_r$trace,function(x)x[i]))
      each_n <- unlist(lapply(each_r$seq,function(x)x[i]))
      yy <- 1+max_y-(1:length(each_p))
      points(x=rep(i,length.out=length(yy)),y=yy,
             col=cc[each_p],pch=16,cex=each_b+2,xpd=TRUE)
      text(x=i,y=yy,each_n,xpd=TRUE,adj=0.5)
    }
    legend(x=max_x+max_x/5,y=max_y/2,all_m,
           fill=cc[all_m],xpd=T,cex=0.75,
           border = NA,bty='n',xjust=0,yjust=0.5)
    text(x=1:max_x,y=max_y+max_y/10,group_name,srt=90,adj=0,cex=0.7,xpd=TRUE)
    text(x=i,y=yy,sprintf('%s%s',round(each_b*10000)/100,'%'),xpd=TRUE,pos=4)
  }else{
    max_x <- length(result_t)
    final_t <- result_t[[max_x]]
    w1 <- which(final_t$prob>=prob_thre)
    final_t_mat <- do.call(rbind,final_t$seq)[w1,]
    final_s_mat <- do.call(rbind,final_t$trace)[w1,]
    max_x <- ncol(final_t_mat)
    max_y <- nrow(final_t_mat)
    all_m <- unique(unlist(final_t$trace[w1]))
    cc <- get.class.color(all_m,use_color=brewer.pal(8,'Pastel2'),
                          pre_define = pre_define)
    # points
    plot_new(xlim = c(1,max_x),ylim=c(1,max_y))
    for(i in 1:max_x){
      each_p <- final_s_mat[,i]
      each_n <- final_t_mat[,i]
      yy <- 1+max_y-c(1:length(each_p))
      points(x=rep(i,length.out=length(yy)),y=yy,
             col=cc[each_p],pch=16,cex=2,xpd=TRUE)
      text(x=i,y=yy,each_n,xpd=TRUE,adj=0.5)
    }   
    each_b <- final_t$prob[w1]
    for(i in 1:max_y){
      segments(x0=1,x1=max_x,
             y0=i,y1=i,
             lwd=3+5*each_b[i],
             col=adjustcolor('dark grey',0.5),xpd=T)
    }
    legend(x=max_x+max_x/5,y=max_y/2,all_m,
           fill=cc[all_m],xpd=T,cex=0.75,
           border = NA,bty='n',xjust=0,yjust=0.5)
    text(x=1:max_x,y=max_y+max_y/10,group_name,srt=90,adj=0,cex=0.7,xpd=TRUE)
    text(x=max_x,y=max_y:1,
         sprintf('%s%s',round(each_b*10000)/100,'%'),xpd=TRUE,pos=4)
  }
  ##
}

###########################################
#' Calculate allele frequency for each site. 
#'
#' \code{HapICE.cal_siteProb} calculate allele frequency for each site.
#'
#' @param final_t list, one list object from HapICE.infer_Haplotype(). 
#'
#' @examples
#' \dontrun{
#'   res_finalHaplotype <- res_trace[[max(length(res_trace))]]
#'   site_prob <- HapICE.cal_siteProb(res_finalHaplotype)
#' }
#' @export
HapICE.cal_siteProb <- function(final_t){
  n <- length(final_t$seq[[1]])
  prob <- final_t$prob
  s1 <- do.call(rbind,final_t$seq)
  t1 <- do.call(rbind,final_t$trace)
  all_r <- list()
  for(i in 1:n){
    r2 <- aggregate(prob,list(t1[,i]),sum)
    rownames(r2) <- r2$Group.1
    all_r[[i]] <- r2
  }
  all_n <- unique(unlist(final_t$trace))
  all_r <- lapply(all_r,function(x){
    x1 <- x$x; names(x1) <- x$Group.1; x1[all_n]
  })
  all_r <- do.call(cbind,all_r)
  rownames(all_r) <- all_n
  all_r <- t(t(all_r)/colSums(all_r,na.rm=T))
  all_r
}

###########################################
#' Output statistics to text files. 
#'
#' \code{HapICE.output_result} output result statistics to text files. 
#' 
#' @param res_readsComponent list, object from HapICE.summ_readsComponent(). 
#' @param res_finalHaplotype list, one list object from HapICE.infer_Haplotype(). 
#' @param group_name character, genome position for each interested positions. 
#' @param output_file character, output file name. 
#'
#' @examples
#' \dontrun{
#' res_readsComponent <- HapICE.summ_readsComponent(all_mat,group_info,group_name)
#' res_finalHaplotype <- res_trace[[max(length(res_trace))]]
#' HapICE.output_result(res_readsComponent=res_readsComponent,
#'                     res_finalHaplotype=res_finalHaplotype,
#'                     group_name=group_name,
#'                     output_file=output_file)
#' }
#' @export
HapICE.output_result <- function(res_readsComponent,
                                 res_finalHaplotype,group_name,
                                 output_file){
  ## output reads component results
  write.table(res_readsComponent,file=gsub(".txt$","_readsComponent.txt",output_file),
              row.names = T,col.names = T,quote = F,sep='\t')
  ## output site probability
  final_res_trace <- do.call(rbind,lapply(1:length(res_finalHaplotype$seq),function(i){
    cbind(rbind(res_finalHaplotype$seq[[i]],res_finalHaplotype$trace[[i]]),
          signif(res_finalHaplotype$prob[[i]],3),res_finalHaplotype$support_read[i])
  }))
  colnames(final_res_trace)[ncol(final_res_trace)] <- 'Support_Read'
  colnames(final_res_trace)[ncol(final_res_trace)-1] <- 'Support_Percentage'
  colnames(final_res_trace)[1:(ncol(final_res_trace)-2)] <- group_name
  write.table(final_res_trace,file=gsub(".txt$","_haplotype.txt",output_file),
              row.names = T,col.names = T,quote = F,sep='\t')
  ## output site probability
  site_prob <- HapICE.cal_siteProb(res_finalHaplotype)
  colnames(site_prob) <- colnames(final_res_trace)[1:ncol(site_prob)]
  write.table(site_prob,file=gsub(".txt$","_siteFreq.txt",output_file),
              row.names = T,col.names = T,quote = F,sep='\t')
  return(TRUE)
}
