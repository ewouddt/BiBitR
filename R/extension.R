
.EnvBIBIT <- new.env()


.GetEnvBIBIT <- function(x){
  if(!exists(x,envir=.EnvBIBIT,inherits=FALSE)){
    return(NULL)
  }
  else{
    return(get(x=x,envir=.EnvBIBIT,inherits=FALSE))
  }
  
}

.AssignEnvBIBIT <- function(x,value){
  assign(x=x,value=value,envir=.EnvBIBIT)
}




check_candidates <- function(data,included,candidates,noise){
  
  if(length(candidates)>0){
    
    noise <- ifelse(((noise<1)&(noise>0)),ceiling(noise*(length(included)+1)),noise)
    
    noise_in_rows <- ncol(data[,included,drop=FALSE])-rowSums(data[,included,drop=FALSE])
    
    rows_noise_allowed <- which((noise - noise_in_rows)>0)
    
    compatible_candidates <- sapply(candidates,FUN=function(x){
      
      zero_rows <- which(data[,x]==0)
      
      if(all(zero_rows %in% rows_noise_allowed)){
        return(TRUE)
      }else{
        return(FALSE)
      }
    })
    return(candidates[compatible_candidates])
  }else{
    return(candidates)
  }
  

}


BC_column_extension_recursive <- function(result,data,noise,extend_limitcol=1,extend_mincol=1){
  
  if(extend_limitcol<=0){stop("Unallowed extend_limitcol")}
  if(extend_mincol<1){stop("Unallowed extend_mincol")}
  extend_mincol <- as.integer(extend_mincol)-1
  
  time_extend <- round(proc.time()['elapsed']/60,2)
  
  BC.extended_list <- rowxnumber_list <-  numberxcol_list <- vector("list",result@Number)
  
  
  for(i.BC in 1:result@Number){
    
    if(extend_limitcol<1){
      extend_limitcol2 <- ceiling(extend_limitcol*sum(result@RowxNumber[,i.BC]))
    }else{
      extend_limitcol2 <- extend_limitcol
    }

    .AssignEnvBIBIT(x="extensions",value=list())
    
    included_temp <- which(result@NumberxCol[i.BC,])
    candidates_temp <- which(!result@NumberxCol[i.BC,])
    
    # Before going into recursion which checks candidates, already delete candidates which have a full-zero column when extending this BC
    
    candidates_temp <- candidates_temp[(which(colSums(data[result@RowxNumber[,i.BC],candidates_temp])>=extend_limitcol2))]
    
    
    temp <- extension_recursive(data=data[result@RowxNumber[,i.BC],],included=included_temp,candidates=candidates_temp,noise=noise,startlength = length(included_temp)+extend_mincol)
    
    # Have we found extensions?
    extensions <- .GetEnvBIBIT("extensions")
    
    # assign(paste0("TEST",i.BC),extensions,envir = .GlobalEnv)
    
    
    if(length(extensions)>0){
      number_columns <- unlist(lapply(extensions,FUN=function(x){x$number_columns}))
      selected_ext <- which(number_columns==max(number_columns))
      
      numberxcol <- matrix(FALSE,nrow=length(selected_ext)+1,ncol=ncol(data))
      
      numberxcol[1,] <- result@NumberxCol[i.BC,]
      for(i.ext in 1:length(selected_ext)){numberxcol[i.ext+1,extensions[[selected_ext[i.ext]]]$included] <- TRUE}
      
      rowxnumber <- matrix(rep(result@RowxNumber[,i.BC],length(selected_ext)+1),ncol=length(selected_ext)+1)
      
      names_temp <- c(paste0("BC",i.BC),paste0("BC",i.BC,"_Ext",1:length(selected_ext)))
      rownames(numberxcol) <- names_temp
      colnames(rowxnumber) <- names_temp
      
      rowxnumber_list[[i.BC]] <- rowxnumber
      numberxcol_list[[i.BC]] <- numberxcol
      
      BC.extended_list[[i.BC]] <- data.frame(BC_Original=i.BC,Number_Extended=length(selected_ext),Number_AddedColumns=(sum(numberxcol[2,])-sum(result@NumberxCol[i.BC,])))
      
    }else{
      rowxnumber_list[[i.BC]] <- result@RowxNumber[,i.BC,drop=FALSE]
      colnames(rowxnumber_list[[i.BC]]) <- paste0("BC",i.BC)
      numberxcol_list[[i.BC]] <- result@NumberxCol[i.BC,,drop=FALSE]
      rownames(numberxcol_list[[i.BC]]) <- paste0("BC",i.BC)
      BC.extended_list[[i.BC]] <- NULL
    }
    
    
  }
  
  #make a object with which BC's were extended, how many resulted BC's (when equal adding length), how many EXTRA columns
  BC.extended <- do.call(rbind,BC.extended_list)
  
  # do a do call rbind/cbind on lists to make final result
  RowxNumber <- do.call(cbind,rowxnumber_list)
  NumberxCol <- do.call(rbind,numberxcol_list)
  
  
  time_extend <- round(proc.time()['elapsed']/60-time_extend,2)
  
  info_temp <- result@info
  info_temp$Time_Min$extend <- time_extend
  info_temp$Time_Min$full <- info_temp$Time_Min$full + time_extend
  info_temp$BC.Extended <- BC.extended
  
  OUT <- new("Biclust",Parameters=result@Parameters,RowxNumber=RowxNumber,NumberxCol=NumberxCol,Number=ncol(RowxNumber),info=info_temp)
  
  return(OUT)

}

extension_recursive <- function(data,included,candidates,noise,startlength){
  
  candidates <- check_candidates(data=data,included=included,candidates,noise)
    
  
  if(length(candidates)>0){
    for(i.candidates in candidates){
      
      included_new <- included
      included_new[length(included_new)+1] <- i.candidates
      candidates_new <- candidates[!candidates==i.candidates]
      
      temp <- extension_recursive(data=data,included=included_new,candidates=candidates_new,noise=noise,startlength=startlength)
      
    }
    
  }else if(length(included)>startlength){
    
    
    extensions_list <- .GetEnvBIBIT("extensions")
    extensions_list[[length(extensions_list)+1]] <- list(
      number_columns=length(included),
      included=included
    )
    .AssignEnvBIBIT(x="extensions",value=extensions_list)
    
  }
  
}





BC_column_extension <- function(result,data,noise,extend_mincol=1,extend_limitcol=1){
  
  time_extend <- round(proc.time()['elapsed']/60,2)
  
  if(extend_limitcol<=0){stop("Unallowed extend_limitcol")}
  if(extend_mincol<1){stop("Unallowed extend_mincol")}
  extend_mincol <- as.integer(extend_mincol)
  
  BC.extended_list <- rowxnumber_list <-  numberxcol_list <- vector("list",result@Number)
  
  
  for(i.BC in 1:result@Number){
    included_columns <- result@NumberxCol[i.BC,]
    
    colsums_temp <- colSums(data[result@RowxNumber[,i.BC],])
    column_candidates <- order(colsums_temp,decreasing=TRUE)
    
    if(extend_limitcol<1){
      extend_limitcol2 <- ceiling(extend_limitcol*sum(result@RowxNumber[,i.BC]))
    }else{
      extend_limitcol2 <- extend_limitcol
    }
    
    column_candidates <- column_candidates[which(colsums_temp[column_candidates]>=extend_limitcol2)]
    

    GO <- TRUE
    i.candidate <- 1
    
    
    while(GO & (i.candidate<=length(column_candidates))){
      
      if(!included_columns[column_candidates[i.candidate]]){
        
        included_columns_temp <- included_columns
        included_columns_temp[column_candidates[i.candidate]] <- TRUE
        
        zeros_allowed <- ifelse(((noise<1)&(noise>0)),ceiling(noise*sum(included_columns_temp)),noise)
        
        zeros_in_rows <- apply(data[result@RowxNumber[,i.BC],included_columns_temp],MARGIN=1,FUN=function(x){sum(x==0)})
        
        if(all(zeros_in_rows<=zeros_allowed)){
          
          included_columns <- included_columns_temp
          i.candidate <- i.candidate+1
          
          
        }else{
          GO <- FALSE
        }
      }else{
        i.candidate <- i.candidate+1
      }
      
    }
    
    n_addedcolumns <- sum(included_columns)-sum(result@NumberxCol[i.BC,])
    
    if(n_addedcolumns>=extend_mincol){
      
      rowxnumber <- matrix(rep(result@RowxNumber[,i.BC],2),ncol=2)
      numberxcol <- matrix(FALSE,nrow=2,ncol=ncol(data))
      numberxcol[1,] <- result@NumberxCol[i.BC,]
      numberxcol[2,included_columns] <- TRUE
      names_temp <- c(paste0("BC",i.BC),paste0("BC",i.BC,"_Ext1"))
      colnames(rowxnumber) <- names_temp
      rownames(numberxcol) <- names_temp
      
      rowxnumber_list[[i.BC]] <- rowxnumber
      numberxcol_list[[i.BC]] <- numberxcol
      
      BC.extended_list[[i.BC]] <- data.frame(BC_Original=i.BC,Number_Extended=1,Number_AddedColumns=n_addedcolumns)
      
      
    }else{
      
      rowxnumber_list[[i.BC]] <- result@RowxNumber[,i.BC,drop=FALSE]
      colnames(rowxnumber_list[[i.BC]]) <- paste0("BC",i.BC)
      numberxcol_list[[i.BC]] <- result@NumberxCol[i.BC,,drop=FALSE]
      rownames(numberxcol_list[[i.BC]]) <- paste0("BC",i.BC)
      BC.extended_list[[i.BC]] <- NULL
    }
    
    

  }
  BC.extended <- do.call(rbind,BC.extended_list)
  
  RowxNumber <- do.call(cbind,rowxnumber_list)
  NumberxCol <- do.call(rbind,numberxcol_list)
  
  
  time_extend <- round(proc.time()['elapsed']/60-time_extend,2)
  
  info_temp <- result@info
  info_temp$Time_Min$extend <- time_extend
  info_temp$Time_Min$full <- info_temp$Time_Min$full + time_extend
  info_temp$BC.Extended <- BC.extended
  
  OUT <- new("Biclust",Parameters=result@Parameters,RowxNumber=RowxNumber,NumberxCol=NumberxCol,Number=ncol(RowxNumber),info=info_temp)
  
  return(OUT)
}






# library(BiBitR)
# set.seed(1)
# data <- matrix(sample(c(0,1),100*100,replace=TRUE,prob=c(0.9,0.1)),nrow=100,ncol=100)
# data[1:10,1:10] <- 1 # BC1
# data[11:20,11:20] <- 1 # BC2
# data[21:30,21:30] <- 1 # BC3
# data <- data[sample(1:nrow(data),nrow(data)),sample(1:ncol(data),ncol(data))]
# result <- bibit(data,minr=5,minc=5)
# result