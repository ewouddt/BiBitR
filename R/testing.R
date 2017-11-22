
# Including autochoice based on clust_sel of Javier
# Notes: 
# - Doesn't really seem to be working
# - row coverage is not necessarily (strict) monotome -> problem?


# - verhouding van hoe fel iets niet meer stijgt ten opzicht van stijging (percentage) -> 1 is cte lijn
#   - vermindering van stijging/stijging => percentage => m.a.w. als dit gelijk is aan 1 dan zijn we op een cte lijn uitgekomen
# - Question: indicatie dat stijging gestopt is, maar de grootte van de stijging wordt niet echt meer in achting genomen

# AUTOGUESS will try and find starts of "plateau's". 
# Add in documentation that it can however still be interesting to check these points -1 (to see what changed/got added)
# Mention caution cause a large increase could also be because of a small pattern (with noise) that has not that much meaning.
# Close inspection is advised (first selected might not always be the best one!!)
# Add in documentation that an extra column is added to the output


# Autotype=1: find where the increasing stops ->  stopping of increasement / current increasment = measure of constant we have gotten
# Autotype=2: find where there was a big jump that stopped

# Add ylim?

ClusterRowCoverage_new <- function(result,matrix,maxCluster=20,rangeCluster=NULL,
                               noise=0.1,noise_select = 0,
                               plots=c(1:3),
                               autoguess=0,
                               verbose=TRUE,
                               plot.type="device",filename="RowCoverage"){
  
  if(!all(plots%in%c(1:3))){stop("plots should be part of c(1,2,3)")}
  if(length(plot.type)!=1){stop("plot.type should be of length 1",call.=FALSE)}
  if(!(plot.type %in% c("device","file","other"))){stop("plot.type should be 'device', 'file' or 'other'",call.=FALSE)}
  FIRSTPLOT <- TRUE
  
  ## PARAMETER CHECKS ##
  if(class(result)!="BiBitWorkflow"){stop("result needs to be of class 'BiBitWorkflow'")}  
  if(class(matrix)!="matrix"){stop("matrix parameter should contain a matrix object",call.=FALSE)}
  if(!identical(as.numeric(as.vector(matrix)),as.numeric(as.logical(matrix)))){stop("matrix is not a binary matrix!",call.=FALSE)}
  if(is.null(rownames(matrix))){rownames(matrix) <- paste0("Row",c(1:nrow(matrix)))}
  if(is.null(colnames(matrix))){colnames(matrix) <- paste0("Col",c(1:ncol(matrix)))}
  biclust_correctdim(result=result,matrix=matrix)
  
  
  
  
  if(is.null(rangeCluster)){
    rangeCluster <- 2:maxCluster
  }
  nrow <- length(rangeCluster)
  if(any(rangeCluster<2)){stop("The minimum number of clusters is 2.")}
  
  cov_df <- data.frame(clusters=rep(NA,nrow),NumberRows=rep(NA,nrow),RowPerc=rep(NA,nrow),FinalBC=rep(NA,nrow))
  
  if(verbose){
    cat("Merging clusters and growing rows:\n")
    pb <- txtProgressBar(min=0,max=sum(1:nrow),initial=0,style=3)
  }
  
  
  for(i in 1:nrow){
    # cat(i,"\n")
    # temp <- capture.output({
    #   out_temp <- BiBitWorkflow(matrix=matrix,
    #                             cut_type="number",cut_pm=i+1,
    #                             noise=noise,noise_select = noise_select,
    #                             plots=c(),
    #                             BCresult = result$info$BiclustInitial,
    #                             simmatresult = result$info$BiclustSimInitial,
    #                             treeresult = result$info$Tree)
    #   
    # })
    # cov_df[i,] <- c(i+1,out_temp$info$Coverage$RowCoverage,out_temp$Biclust@Number)
    
    temp <- capture.output({
      
      result2 <- workflow_mergeBC(result=result$info$BiclustInitial,tree=result$info$Tree,JI=NULL,number=rangeCluster[i])
      
      BC.Merge <- result2@info$BC.Merge
      MergedColPatterns <- lapply(as.list(seq_len(result2@Number)),FUN=function(pattern){
        return(which(result2@NumberxCol[pattern,]))
      })
      names(MergedColPatterns) <- paste0("Pat",seq_len(result2@Number))
      
      NoisexNumber <- apply(result2@NumberxCol,MARGIN=1,FUN=function(pattern){
        return(apply(matrix[,pattern],MARGIN=1,FUN=function(row){return(sum(row==0))}))
      })
      data_noisescree <- lapply(seq_len(ncol(NoisexNumber)),FUN=function(pattern){
        tab <- table(NoisexNumber[,pattern])
        return(data.frame(Noise=as.numeric(names(tab)),Total=as.numeric(tab)))
      })
      noise_threshold <- workflow_noisethreshold(noise,noise_select,data_noisescree)
      result2@info$Noise.Threshold <- noise_threshold
      
      result3 <-  workflow_UpdateBiclust_RowNoise(result=result2,matrix=matrix, noise=result2@info$Noise.Threshold,removeBC=TRUE,NoisexNumber=NoisexNumber)
    })
    
    BC_rownames <- unlist(apply(result3@RowxNumber,MARGIN=2,FUN=function(x){return(rownames(matrix)[x])}))
    nBC_rownames <- length(unique(BC_rownames))
    coverage <- c(NumberRows=round(nBC_rownames,0),RowPerc=round((nBC_rownames/nrow(matrix))*100,2))
    
    
    cov_df[i,] <- c(rangeCluster[i],coverage,result3@Number)
    
    if(verbose){
      setTxtProgressBar(pb,sum(1:i))
    }
    
  }
  
  if(verbose){
    close(pb)
  }
  
  
  ## AUTO CHOICE ##
  # if(autoguess>0){
  #   if(verbose){cat("\nAutoguessing...")}
  #   
    # # gaps <- -diff(cov_df$NumberRows[nrow(cov_df):1])
    # # gaps <- diff(cov_df$NumberRows) + 0.01
    # gaps <- -diff(cov_df$NumberRows[nrow(cov_df):1])+0.01
    # 
    # diffcominggap_gap <- diff(gaps)/gaps[-length(gaps)]*100
    # # diffcominggap_gap <- sapply(diffcominggap_gap,FUN=function(x){
    # #   if(is.na(x)){return(NA)}
    # #   if(x==Inf | x==-Inf){return(0)} else {return(x)}
    # # })
    # 
    # selection <- sort.list(diffcominggap_gap)[1:autoguess]
    # # selection <- cov_df$clusters[nrow(cov_df):1][]
    # # selection <- sort.list(abs(diffcominggap_gap),decreasing=TRUE)[1:autoguess]
    # # selection <- sort.list(-diffcominggap_gap)[1:autoguess]
    # 
    # 
    # cov_df$autoguess <- 0
    # cov_df$autoguess[nrow(cov_df):1][selection] <- 1:autoguess
    # # cov_df$autoguess[selection] <- 1:autoguess
  #   
  #   
  #   
  #   # minus <- c(2,2,1,10,8,1,3,1,3,1)
  #   # sm1 <- rep(60,length(minus)) - cumsum(minus)
  #   # sm1 <- sm1[length(sm1):1]
  #   # gaps <- -diff(sm1)
  #   # col <- rep("black",length(sm1))
  #   # col[selection] <- "red"
  #   # plot(sm1,col=col,pch=19)
  #   
  #   if(verbose){cat("Done")}
  # }
  
  #################
  if(autoguess>0){
    # if(verbose){cat("\nAutoguessing...")}
    
    # Auto Guess - Try 2
    gaps <- diff(cov_df$NumberRows)+0.01 # +0.01 in order to avoid /0  OR do not this and delete Inf -Inf later
    

    
    # fixing decreasing? # 
    gaps[gaps<0] <- 0.01 # or 0
    # Notes on it:
    # - Could work, but any increase after the drop would be ignored and not chosen as a good point
    # - this may however be what we want, because why choose a later point of the same height if we can do it with less clusters?
    
    
    
    # Prevent the small increase and then cte always winning 
    # Can not put all to NA, cause that would also exclude large increases, then cte
    # BUT what to do with decreasing? -> handle seperately -> probably okay to just see it as cte's always 
    # g0 <- which(gaps==0) # or 0.01
    # Either percentage based on max jump OR a minimum (e.g. 10) 
    # gaps[g0[which(gaps[g0-1]<0.1*max(gaps))]] <- NA
    # gaps[g0[which(gaps[g0-1]<10)]] <- NA
    
    # OR can we come up with a percentage score that adds more weight to the jump than the fact it levels out (goes cte)
    
    
    # In case of not 0.01 for 0:
    # temp[temp==Inf | temp==-Inf] <- NA
    
    
    # autotype=1
    # selection <- order(diff(gaps))+1  # This difference will be very negative if the `increasing/climbing` is stopping
    
    # autotype=2
    selection2 <- order(diff(gaps)/gaps[-length(gaps)])+1  
    # Why do we divide? In order to avoid selecting points 'being in the increase'? (No)
    # Dividing gives preference to points after which the growth is supersmall (this is why having constant is a problem!)
  
    
    # rbind(diff(gaps),gaps[-length(gaps)])
    
    cov_df$autoguess <- 0
    cov_df$autoguess[selection2[1:autoguess]] <- 1:autoguess
    # if(verbose){cat("Done")}
  }
  ###
  
  col <- rep("black",nrow(cov_df))
  if(autoguess>0){col[cov_df$autoguess>0] <- "red"}
  
  if(1 %in% plots){
    if(plot.type=="device"){
      dev.new()
    }else if(plot.type=="file" & FIRSTPLOT){
      pdf(paste0(filename,".pdf"))
      FIRSTPLOT <- FALSE
    }
    
    plot(cov_df$clusters,cov_df$RowPerc,col=col,main="Clusters vs Row Coverage Perc.",xlab="n clusters",ylab="Row CovPerc",pch=19)
    points(cov_df$clusters,cov_df$RowPerc,type="l")
    text(cov_df$clusters,cov_df$RowPerc,as.character(round(cov_df$RowPerc,2)),pos=3)
  }
  
  if(2 %in% plots){
    if(plot.type=="device"){
      dev.new()
    }else if(plot.type=="file" & FIRSTPLOT){
      pdf(paste0(filename,".pdf"))
      FIRSTPLOT <- FALSE
    }
    plot(cov_df$clusters,cov_df$NumberRows,col=col,main="Clusters vs Total Number Rows",xlab="n clusters",ylab="Total Number Rows",pch=19)
    points(cov_df$clusters,cov_df$NumberRows,type="l")
    text(cov_df$clusters,cov_df$NumberRows,as.character(cov_df$NumberRows),pos=3)
  }
  
  if(3 %in% plots){
    if(plot.type=="device"){
      dev.new()
    }else if(plot.type=="file" & FIRSTPLOT){
      pdf(paste0(filename,".pdf"))
      FIRSTPLOT <- FALSE
    }
    plot(cov_df$clusters,cov_df$FinalBC,main="Clusters vs Final Number BC",xlab="n clusters",ylab="Final BC",pch=19)
    points(cov_df$clusters,cov_df$FinalBC,type="l")
    text(cov_df$clusters,cov_df$FinalBC,as.character(cov_df$FinalBC),pos=3)
    
  }
  if(plot.type=="file" & length(plots)>0){
    dev.off()
  }
  
  return(cov_df)
}
