bibit3_alt <- function(matrix=NULL,minr=1,minc=2,noise=0,pattern_matrix=NULL,subpattern=TRUE,pattern_combinations=FALSE,arff_row_col=NULL,
                   extend_columns="none",extend_mincol=1,extend_limitcol=1,extend_noise=noise,extend_contained=FALSE){
  
  pm <- match.call()
  minr <- minr + 2
  
  ###
  # Legacy compatibility for GUI
  if(is.logical(extend_columns)){
    extend_columns <- ifelse(extend_columns,"naive","none")
  }
  
  if(noise<0){stop("noise parameter can not be negative",call.=FALSE)}
  if(noise>=1){noise <- as.integer(noise)}
  
  ## Extend parameters
  if(extend_noise<0){stop("extend_noise parameter can not be negative",call.=FALSE)}
  if(extend_noise>=1){extend_noise <- as.integer(extend_noise)}
  if(extend_noise<noise){stop("extend_noise can't be lower than noise",call.=FALSE)}
  if(length(extend_columns)!=1){stop("extend_columns needs 1 input",call.=FALSE)}
  if(!(extend_columns)%in%c("none","naive","recursive")){stop("extend_columns should be \"none\", \"naive\" or \"recursive\"",call.=FALSE)}
  if(extend_limitcol<=0){stop("extend_limitcol should be larger than 0",call.=FALSE)}
  if(extend_mincol<1){stop("extend_mincol should be larger than or equal to 1",call.=FALSE)}
  
  ###
  
  
  if(is.null(arff_row_col)){
    
    # Check if matrix is binary (DISCRETIZED NOT YET IMPLEMENTED!)
    if(class(matrix)!="matrix"){stop("matrix parameter should contain a matrix object",call.=FALSE)}
    if(!identical(as.numeric(as.vector(matrix)),as.numeric(as.logical(matrix)))){stop("matrix is not a binary matrix!",call.=FALSE)}
    
    if(is.null(rownames(matrix))){rownames(matrix) <- paste0("Row",c(1:nrow(matrix)))}
    if(is.null(colnames(matrix))){colnames(matrix) <- paste0("Col",c(1:ncol(matrix)))}
    
    # Check if rownames & colnames contain ; or ,  -> should be deleted and give warnings it was deleted
    rowdot <- grepl(",",rownames(matrix))
    if(sum(rowdot)>0){
      rownames(matrix) <- gsub(",","",rownames(matrix))
      warning(paste0("Row names ",paste0(which(rowdot),collapse = ",")," contained a ',' which was deleted."),call.=FALSE)
    }
    rowsc <- grepl(";",rownames(matrix))
    if(sum(rowsc)>0){
      rownames(matrix) <- gsub(";","",rownames(matrix))
      warning(paste0("Row names ",paste0(which(rowsc),collapse = ",")," contained a ';' which was deleted."),call.=FALSE)
    }
    coldot <- grepl(",",colnames(matrix))
    if(sum(coldot)>0){
      colnames(matrix) <- gsub(",","",colnames(matrix))
      warning(paste0("Column names ",paste0(which(coldot),collapse = ",")," contained a ',' which was deleted."),call.=FALSE)
    }
    colsc <- grepl(";",colnames(matrix))
    if(sum(colsc)>0){
      colnames(matrix) <- gsub(";","",colnames(matrix))
      warning(paste0("Column names ",paste0(which(colsc),collapse = ",")," contained a ';' which was deleted."),call.=FALSE)
    }
    
    # No duplicate row names allowed!
    if(sum(table(rownames(matrix))>1)){stop("No duplicate row names allowed!")}
    
    # Check pattern matrix
    if(is.null(pattern_matrix)){stop("pattern_matrix needs to be provided",call.=FALSE)}
    if(class(pattern_matrix)!="matrix"){stop("pattern_matrix parameter should contain a matrix object",call.=FALSE)}
    if(!identical(as.numeric(as.vector(pattern_matrix)),as.numeric(as.logical(pattern_matrix)))){stop("pattern_matrix is not a binary matrix!",call.=FALSE)}
    if(is.null(rownames(pattern_matrix))){rownames(pattern_matrix) <- paste0("Pattern",1:nrow(pattern_matrix))}
    if(ncol(pattern_matrix)!=ncol(matrix)){stop("matrix and pattern_matrix have a different number of columns",call.=FALSE)}
    
    
    # If combinations required, add to pattern!
    if(pattern_combinations & nrow(pattern_matrix)>1){
      cat("Computing pattern combinations...")
      comb_temp <- combn(1:nrow(pattern_matrix),2)
      comb_matrix <- matrix(NA,nrow=ncol(comb_temp),ncol=ncol(pattern_matrix),dimnames=list(paste0("comb",1:ncol(comb_temp))))
      
      for(i.comb in 1:ncol(comb_temp)){
        comb_matrix[i.comb,] <- ((pattern_matrix[comb_temp[1,i.comb],]+pattern_matrix[comb_temp[2,i.comb],])==2)+0
        rownames(comb_matrix)[i.comb] <- paste0(rownames(pattern_matrix)[comb_temp[1,i.comb]],"_",rownames(pattern_matrix)[comb_temp[2,i.comb]])
      }
      
      pattern_matrix <- rbind(pattern_matrix,comb_matrix)
      
      cat("DONE\n\n")
    }
    
    # Delete zero-rows
    zero_rows <- which(rowSums(pattern_matrix)==0)
    if(length(zero_rows)>0){
      pattern_matrix <- pattern_matrix[-zero_rows,,drop=FALSE]
    }
    
    
    nPatterns <- nrow(pattern_matrix)
    if(nPatterns==0){stop("No viable patterns in pattern_matrix, all zero values.")}
    
  }else{
    time_arff <- 0
    
    if(length(arff_row_col)!=3){stop("arff_row_col should contain 3 elements",call.=FALSE)}
    bibitdata_path <- arff_row_col[1]
    bibitrows_path <- arff_row_col[2]
    bibitcols_path <- arff_row_col[3]
    
    pattern_matrix <- matrix(NA,nrow=1,ncol=1,dimnames=list("arff_Pattern"))
    nPatterns <- 1
  }
  
  #############################################
  ## PREPARE BASIC ARFF FILE & READ IN LINES ##
  #############################################
  
  cat("Transform matrix into arff format...")
  
  bibitbasic_path <- tempfile("bibitbasic",fileext=".arff")
  write.arff(t(matrix),file=bibitbasic_path)
  basic_file <- file(bibitbasic_path)
  basic_lines <- readLines(basic_file)
  close(basic_file)
  
  number_white <- nrow(matrix)+2
  
  cat("DONE\n\n")
  
  ######################################
  ## START FOR LOOP OVER ALL PATTERNS ##
  ######################################
  FINAL_RESULT <- vector("list",nPatterns)
  names(FINAL_RESULT) <- rownames(pattern_matrix)
  
  
  for(i.pattern in 1:nPatterns){
    
    if(i.pattern>1){cat("\n=============================================================================\n\n")}
    cat(toupper(rownames(pattern_matrix)[i.pattern]),"\n\n")
    
    if(is.null(arff_row_col)){
      time_arff <- round(proc.time()['elapsed']/60,2)
      
      # Add patterns to matrix
      matrix_with_pattern <- rbind(matrix(rep(pattern_matrix[i.pattern,],2),nrow=2,byrow=TRUE,dimnames = list(paste0(rownames(pattern_matrix)[i.pattern],"_Art",c(1,2)))),matrix)
      
      # Transform data into arff format
      cat("Changing arff file...",rownames(pattern_matrix)[i.pattern],"...")
      
      bibitdata_path <- tempfile("bibitdata",fileext=".arff")
      bibitrows_path <- tempfile("bibitrows",fileext=".csv")
      bibitcols_path <- tempfile("bibitcols",fileext=".csv")
      
      
      new_lines_meta <- basic_lines[1:number_white]
      new_lines_data <- basic_lines[(number_white+1):length(basic_lines)]
      
      pattern <- apply(cbind(matrix_with_pattern[1,],matrix_with_pattern[2,]),MARGIN=1,FUN=paste0,collapse=",")
      new_rownames <- rownames(matrix_with_pattern)[c(1,2)]
      
      meta1 <- new_lines_meta[1]
      new_lines_meta <- new_lines_meta[-1]
      new_lines_meta <- c(meta1,paste0("@attribute ",new_rownames," numeric"),new_lines_meta)
      new_lines_data <- apply(cbind(pattern,new_lines_data),MARGIN=1,FUN=paste0,collapse=",")
      
      new_file <- file(bibitdata_path)
      writeLines(c(new_lines_meta,new_lines_data),new_file)
      close(new_file)
       
      write.table(matrix(rownames(matrix_with_pattern),ncol=1),quote=FALSE,row.names=FALSE,col.names=FALSE,file=bibitrows_path)
      write.table(matrix(colnames(matrix_with_pattern),ncol=1),quote=FALSE,row.names=FALSE,col.names=FALSE,file=bibitcols_path)
      
      cat("DONE\n")
      cat("\n")
      
      time_arff <- round(proc.time()['elapsed']/60-time_arff,2)
      
    }else{
      matrix_with_pattern <- NULL
      if(extend_columns!="none"){
        
        matrix_with_pattern <- read.arff(bibitdata_path)
        rownames.data <- as.character(read.table(bibitrows_path,header=FALSE)[,1])
        colnames.data <- as.character(read.table(bibitcols_path,header=FALSE)[,1])
        if(length(rownames.data)!=nrow(matrix_with_pattern)){
          matrix_with_pattern <- t(matrix_with_pattern)
        }
        rownames(matrix_with_pattern) <- rownames.data
        colnames(matrix_with_pattern) <- colnames.data
        
      }
    }
    
    # Apply BiBit Algorithm
    
    cat("Initiate BiBit for",rownames(pattern_matrix)[i.pattern],"...\n")
    cat("\n")
    
    bibitoutput_path <- tempfile("bibitoutput",fileext = "")
    
    
    time_bibit <- proc.time()['elapsed']/60
    
    
    javaloc <- paste0(find.package("BiBitR")[1],"/java/BiBit3.jar")
    # javaloc <- paste0(getwd(),"/inst/java/BiBit3.jar")
    
    subpat <- ifelse(subpattern,1,0)
    
    # BiBit.jar location needs to be standardized for package location! # .libPaths()
    command <- paste("java -jar -Xmx1000M",paste0("\"",javaloc,"\""),paste0("\"",bibitdata_path,"\""),"1",minr,minc,paste0("\"",bibitoutput_path,"\""),paste0("\"",bibitrows_path,"\""),paste0("\"",bibitcols_path,"\""),1,paste0(" ",noise),paste0(" ",subpat))
    # cat(command,"\n")
    system(command)
    
    time_bibit <- round(proc.time()['elapsed']/60-time_bibit,2)
    
    
    cat("\n")
    cat("Transforming into biclust output...")
    
    time_biclust <- round(proc.time()['elapsed']/60,2)
    result <- bibit2biclust(data=matrix_with_pattern,resultpath=paste0(bibitoutput_path,"_1.txt"),arff_row_col = arff_row_col)
    cat("DONE\n")
    time_biclust <- round(proc.time()['elapsed']/60-time_biclust,2)
    
    
    # Small prep
    if(!is.null(arff_row_col)){
      rownames.data <- as.character(read.table(arff_row_col[2],header=FALSE)[,1])
      colnames.data <- as.character(read.table(arff_row_col[3],header=FALSE)[,1])
      nrow.data <- length(rownames.data)
      ncol.data <- length(colnames.data)
    }else{
      nrow.data <- nrow(matrix_with_pattern)
      ncol.data <- ncol(matrix_with_pattern)
    }
    
    
    
    # Look for and label the Biclusters (Full Pattern (zero or not)/Sub Pattern)
    
    
    if(!is.null(result)){
      result2 <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                     RowxNumber=result$RowxNumber,
                     NumberxCol=result$NumberxCol,
                     Number=result$Number,
                     info=list(Time_Min=list(arff=time_arff,bibit=time_bibit,biclust=time_biclust,full=time_arff+time_bibit+time_biclust)))
      
      
      
      FullPattern <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                         RowxNumber=result2@RowxNumber[,1,drop=FALSE],
                         NumberxCol=result2@NumberxCol[1,,drop=FALSE],
                         Number=1,
                         info=list())
      
      
      if(subpattern & result2@Number>1){
        SubPattern <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                          RowxNumber=result2@RowxNumber[,2:result2@Number,drop=FALSE],
                          NumberxCol=result2@NumberxCol[2:result2@Number,,drop=FALSE],
                          Number=result2@Number-1,
                          info=list())
      }else{
        SubPattern <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                          RowxNumber=matrix(FALSE,nrow=nrow.data,ncol=1),
                          NumberxCol=matrix(FALSE,nrow=1,ncol=ncol.data),
                          Number=0,
                          info=list())
      }
      
      
      if(extend_columns!="none"){
        
        # Reduce matrix and result only for Extended part (Artificial rows may not influence extension procedure)
        ########################
        ########################
        result2_temp <- result2
        
        result2_temp@RowxNumber <- result2_temp@RowxNumber[-c(1,2),,drop=FALSE]
        
        if(result2_temp@Number>0){
          deleteBC_index <- which(colSums(result2_temp@RowxNumber)==0)
          
          if(length(deleteBC_index)>0){
            if(length(deleteBC_index)==result2_temp@Number){
              result2_temp <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                                  RowxNumber=matrix(FALSE,nrow=nrow.data,ncol=1),
                                  NumberxCol=matrix(FALSE,nrow=1,ncol=ncol.data),
                                  Number=0,
                                  info=list())
            }else{
              result2_temp@RowxNumber <- result2_temp@RowxNumber[,-deleteBC_index]
              result2_temp@NumberxCol <- result2_temp@NumberxCol[-deleteBC_index,]
              result2_temp@Number <- ncol(result2_temp@RowxNumber)
            }
          }
        }
        ########################
        ########################
        
        # Use extension_procedure, delete original BC's, check if there were extensions...
        # check for BC.Extender, if NULL, then make similar object below, otherwise delete parts
        
        Extended <- extension_procedure(result2=result2_temp,data=matrix_with_pattern[-c(1,2),],extend_noise=extend_noise,extend_mincol=extend_mincol,extend_limitcol=extend_limitcol,extend_columns=extend_columns,extend_contained=extend_contained)
        
        if(!is.null(Extended@info$BC.Extended)){
          
          original_index <- which(!grepl("_Ext",colnames(Extended@RowxNumber)))
          Extended@Number <- Extended@Number - length(original_index)
          Extended@RowxNumber <- Extended@RowxNumber[,-original_index,drop=FALSE]
          Extended@NumberxCol <- Extended@NumberxCol[-original_index,,drop=FALSE]
          
          time_extend <- Extended@info$Time_Min$extend
          Extended@info$Time_Min <- NULL
          Number_Extended <- Extended@Number
          
        }else{
          Extended <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                          RowxNumber=matrix(FALSE,nrow=nrow.data,ncol=1),
                          NumberxCol=matrix(FALSE,nrow=1,ncol=ncol.data),
                          Number=0,
                          info=list())  
          time_extend <- 0
          Number_Extended <- 0
        }
        
        
        # TO DO: take time extend from result + TO DO: add and check parameters + add documentaiton
        
      }else{
        Extended <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                        RowxNumber=matrix(FALSE,nrow=nrow.data,ncol=1),
                        NumberxCol=matrix(FALSE,nrow=1,ncol=ncol.data),
                        Number=0,
                        info=list())  
        time_extend <- 0
        Number_Extended <- 0
      }
      
      time_final <- list(arff=result2@info$Time_Min$arff,bibit= result2@info$Time_Min$bibit,biclust=result2@info$Time_Min$biclust,extend=time_extend,full=result2@info$Time_Min$full+time_extend)
      
      FINAL_RESULT[[i.pattern]] <- list(Number=result2@Number,Number_Extended=Number_Extended,FullPattern=FullPattern,SubPattern=SubPattern,Extended=Extended,info=list(Time_Min=time_final))
      
      
    }else{
      
      
      result2 <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                     RowxNumber=matrix(FALSE,nrow=nrow.data,ncol=1),
                     NumberxCol=matrix(FALSE,nrow=1,ncol=ncol.data),
                     Number=0,
                     info=list())
      
      FINAL_RESULT[[i.pattern]] <- list(Number=0,Number_Extended=0,FullPattern=result2,SubPattern=result2,Extended=result2,info=list(Time_Min=list(arff=time_arff,bibit=time_bibit,biclust=time_biclust,extend=0,full=time_arff+time_bibit+time_biclust)))
      
    }
    
  }
  
  # DELETE ARTIFICIAL ROWS FROM BC RESULTS , if no other rows remain, go to empty result
  
  for(i.list in 1:length(FINAL_RESULT)){
    for(j.list in c("FullPattern","SubPattern")){
      
      FINAL_RESULT[[i.list]][[j.list]]@RowxNumber <- FINAL_RESULT[[i.list]][[j.list]]@RowxNumber[-c(1,2),,drop=FALSE]
      
      if(FINAL_RESULT[[i.list]][[j.list]]@Number>0){
        deleteBC_index <- which(colSums(FINAL_RESULT[[i.list]][[j.list]]@RowxNumber)==0)
        
        if(length(deleteBC_index)>0){
          if(length(deleteBC_index)==FINAL_RESULT[[i.list]][[j.list]]@Number){
            FINAL_RESULT[[i.list]][[j.list]] <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                                                    RowxNumber=matrix(FALSE,nrow=nrow.data,ncol=1),
                                                    NumberxCol=matrix(FALSE,nrow=1,ncol=ncol.data),
                                                    Number=0,
                                                    info=list())
          }else{
            FINAL_RESULT[[i.list]][[j.list]]@RowxNumber <- FINAL_RESULT[[i.list]][[j.list]]@RowxNumber[,-deleteBC_index]
            FINAL_RESULT[[i.list]][[j.list]]@NumberxCol <- FINAL_RESULT[[i.list]][[j.list]]@NumberxCol[-deleteBC_index,]
            FINAL_RESULT[[i.list]][[j.list]]@Number <- ncol(FINAL_RESULT[[i.list]][[j.list]]@RowxNumber)
          }
        }
      }
    }
  }
  
  # END RESULT
  FINAL_RESULT$pattern_matrix <- pattern_matrix
  class(FINAL_RESULT) <- "bibit3"
  return(FINAL_RESULT)
}

