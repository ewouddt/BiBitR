
## IMPORTS ##

#' @importFrom foreign write.arff 
#' @import biclust
NULL


#' @title The BiBit Algorithm
#' 
#' @description A R-wrapper which directly calls the original Java code for the BiBit algorithm (\url{http://eps.upo.es/bigs/BiBit.html} and transforms it to the output format of \code{Biclust}.
#' 
#' @details This function uses the original Java code directly (with the intended input and output). Because the Java code was not refactored, the \code{rJava} package could not be used.
#' The \code{bibit} function does the following:
#' \enumerate{
#' \item Convert R matrix to a \code{.arff} output file.
#' \item Use the \code{.arff} file as input for the Java code which is called by \code{system()}.
#' \item The outputted \code{.txt} file from the Java BiBit algorithm is read in and transformed to a \code{Biclust} object.
#' }
#' Because of this, there is a chance of \emph{overhead} when applying the algorithm on large datasets. Make sure your machine has enough RAM available when applying to big data.
#' 
#' @author Ewoud De Troyer
#' 
#' @references Domingo S. Rodriguez-Baena, Antonia J. Perez-Pulido and Jesus S. Aguilar-Ruiz (2011), "A biclustering algorithm for extracting bit-patterns from binary datasets", \emph{Bioinformatics}
#' 
#' @export
#' @param matrix The binary input matrix.
#' @param minr The minimum number of rows of the Biclusters.
#' @param minc The minimum number of columns of the Biclusters.
#' 
#' @return A Biclust S4 Class object.
#' 
#' @examples 
#' \dontrun{
#' data <- matrix(sample(c(0,1),100*100,replace=TRUE,prob=c(0.9,0.1)),nrow=100,ncol=100)
#' data[1:10,1:10] <- 1 # BC1
#' data[11:20,11:20] <- 1 # BC2
#' data[21:30,21:30] <- 1 # BC3
#' data <- data[sample(1:nrow(data),nrow(data)),sample(1:ncol(data),ncol(data))]
#' result <- bibit(data,minr=5,minc=5)
#' result
#' }
bibit <- function(matrix,minr=2,minc=2){
  
  pm <- match.call()
  time_arff <- round(proc.time()['elapsed']/60,2)
    
  # Check if matrix is binary (DISCRETIZED NOT YET IMPLEMENTED!)
  if(!identical(as.vector(matrix),as.numeric(as.logical(matrix)))){stop("matrix is not a binary matrix!",call.=FALSE)}
  
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
  
  
  # Transform data into arff format
  cat("Transform matrix into arff format...")

  bibitdata_path <- tempfile("bibitdata",fileext=".arff")
  bibitrows_path <- tempfile("bibitrows",fileext=".csv")
  bibitcols_path <- tempfile("bibitcols",fileext=".csv")
  
  write.arff(t(matrix),file=bibitdata_path)
  write.table(matrix(rownames(matrix),ncol=1),quote=FALSE,row.names=FALSE,col.names=FALSE,file=bibitrows_path)
  write.table(matrix(colnames(matrix),ncol=1),quote=FALSE,row.names=FALSE,col.names=FALSE,file=bibitcols_path)
  
  bibitoutput_path <- tempfile("bibitoutput",fileext = "")
  cat("DONE\n")
  cat("\n")
  
  time_arff <- round(proc.time()['elapsed']/60-time_arff,2)
  time_bibit <- proc.time()['elapsed']/60
  
  javaloc <- paste0(.libPaths(),"/BiBitR/java/BiBit.jar")
  # javaloc <- gsub("/","\\\\",javaloc)
  
  # BiBit.jar location needs to be standardized for package location! # .libPaths()
  # command <- paste("java -jar -Xmx1000M",javaloc,bibitdata_path,"1",minr,minc,bibitoutput_path,bibitrows_path,bibitcols_path,1)
  command <- paste("java -jar -Xmx1000M",paste0("\"",javaloc,"\""),paste0("\"",bibitdata_path,"\""),"1",minr,minc,paste0("\"",bibitoutput_path,"\""),paste0("\"",bibitrows_path,"\""),paste0("\"",bibitcols_path,"\""),1)
  
  system(command)
  
  time_bibit <- round(proc.time()['elapsed']/60-time_bibit,2)
  
  cat("\n")
  cat("Transforming into biclust output...")
  
  time_biclust <- round(proc.time()['elapsed']/60,2)
  result <- bibit2biclust(data=matrix,resultpath=paste0(bibitoutput_path,"_1.txt"))
  cat("DONE\n")
  time_biclust <- round(proc.time()['elapsed']/60-time_biclust,2)
  
  result$info$time_minutes <- list(arff=time_arff,bibit=time_bibit,biclust=time_biclust,full=time_arff+time_bibit+time_biclust)
  
  result2 <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                 RowxNumber=result$RowxNumber,
                 NumberxCol=result$NumberxCol,
                 Number=result$Number,
                 info=list(Time_Min=result$info$time_minutes))
                 
  return(result2)
}



bibit2biclust <- function(data,resultpath){
  result <- read.table(resultpath,header=TRUE,sep=";")
  
  if(dim(result)[1]>0){
    
    result$Rows <- as.character(result$Rows)
    result$Columns <- as.character(result$Columns)
    
    Number <- nrow(result)
    
    rowlist <- strsplit(result$Rows,",")
    # for(i in 1:length(rowlist)){
    #   rowlist[[i]] <- rowlist[[i]][1:result$NumOfRows[i]]
    # }
    
    collist <- strsplit(result$Columns,", ")
    # for(i in 1:length(collist)){
    #   collist[[i]] <- collist[[i]][1:result$NumOfColumns[i]]
    # }
    
    # Let's add a quick to avoid problems...
    if(!identical(result$NumOfRows,unlist(lapply(rowlist,FUN=length)))){warning("Issue reading row names...")}
    if(!identical(result$NumOfColumns,unlist(lapply(collist,FUN=length)))){warning("Issue reading column names...")}
    
    
    rowlist_index <- lapply(rowlist,FUN=function(x){rownames(data) %in%  x})
    collist_index <- lapply(collist,FUN=function(x){colnames(data) %in%  x})
    
    RowxNumber <- matrix(unlist(rowlist_index),byrow=FALSE,nrow=nrow(data),ncol=Number)
    NumberxCol <- matrix(unlist(collist_index),byrow=TRUE,nrow=Number,ncol=ncol(data))
    
    # again quick BC dimension check 
    if(!identical(result$NumOfRows,as.integer(colSums(RowxNumber)))){warning("Issue row BC dimension")}
    if(!identical(result$NumOfColumns,as.integer(rowSums(NumberxCol)))){warning("Issue column BC dimension")}
    
    
    
    # Temporart list output, needs to be changed to biclust object
    return(list(Parameters=list(),Number=Number,RowxNumber=RowxNumber,NumberxCol=NumberxCol,info=list()))
    
  }else{
    return(NULL)
  }
}

#' Finding Maximum Size Biclusters
#' 
#' Simple function which scans a \code{Biclust} result and returns which biclusters have maximum row, column or size (row*column).
#' 
#' @author Ewoud De Troyer
#' 
#' @export
#' @param bicresult A \code{Biclust} result. (e.g. The return object from \code{bibit})
#' 
#' @return A list containing:
#' \itemize{
#' \item \code{$row}: A matrix containing in the columns the Biclusters which had maximum rows, and in the rows the Row Dimension, Column Dimension and Size.
#' \item \code{$column}: A matrix containing in the columns the Biclusters which had maximum columns, and in the rows the Row Dimension, Column Dimension and Size.
#' \item \code{$size}: A matrix containing in the columns the Biclusters which had maximum size, and in the rows the Row Dimension, Column Dimension and Size.
#' }
#' 
#' @examples 
#' \dontrun{
#' data <- matrix(sample(c(0,1),100*100,replace=TRUE,prob=c(0.9,0.1)),nrow=100,ncol=100)
#' data[1:10,1:10] <- 1 # BC1
#' data[11:20,11:20] <- 1 # BC2
#' data[21:30,21:30] <- 1 # BC3
#' data <- data[sample(1:nrow(data),nrow(data)),sample(1:ncol(data),ncol(data))]
#' result <- bibit(data,minr=2,minc=2)
#' 
#' MaxBC(result)
#' 
#' }
MaxBC <- function(bicresult){
  if(class(bicresult)!="Biclust"){stop("bicresult needs to be of class 'Biclust'")}
  
  colsum <- colSums(bicresult@RowxNumber)
  rowsum <- rowSums(bicresult@NumberxCol)
  sizesum <- rowsum*colsum
  
  ind.colmax <- which(max(colsum)==colsum)
  ind.rowmax <- which(max(rowsum)==rowsum)
  ind.sizemax <- which(max(sizesum)==sizesum)
  

  row <- rbind(RowDim=rowsum[ind.rowmax],ColDim=colsum[ind.rowmax],SizeDim=sizesum[ind.rowmax])
  colnames(row) <- paste0("BC",ind.rowmax)
  
  column <- rbind(RowDim=rowsum[ind.colmax],ColDim=colsum[ind.colmax],SizeDim=sizesum[ind.colmax])
  colnames(column) <- paste0("BC",ind.colmax)
  
  size <- rbind(RowDim=rowsum[ind.sizemax],ColDim=colsum[ind.sizemax],SizeDim=sizesum[ind.sizemax])
  colnames(size) <- paste0("BC",ind.sizemax)
  
  return(list(row=row,column=column,size=size))
}


