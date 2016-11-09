
## IMPORTS ##

#' @importFrom foreign write.arff 
#' @import biclust
NULL


#' @title The BiBit Algorithm
#' 
#' @description A R-wrapper which directly calls the original Java code for the BiBit algorithm (\url{http://eps.upo.es/bigs/BiBit.html}) and transforms it to the output format of the \code{Biclust} R package.
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
#' @param arff_row_col If you want to circumvent the internal R function to convert the matrix to \code{.arff} format, provide the pathname of this file here. Additionally, two \code{.csv} files should be provided containing 1 column of row and column names. These two files should not contain a header or quotes around the names, simply 1 column with the names.\cr 
#' (\emph{Example}: \code{arff_row_col=c("...\\\\data\\\\matrix.arff","...\\\\data\\\\rownames.csv","...\\\\data\\\\colnames.csv")})\cr
#' \emph{Note:} These files can be generated with the \code{make_arff_row_col} function.
#' @param output_path If as output, the original txt output of the Java code is desired, provide the outputh path here (without extension). In this case the \code{bibit} function will skip the transformation to a Biclust class object and simply return \code{NULL}.\cr 
#' (\emph{Example}: \code{output_path="...\\\\out\\\\bibitresult"})
#' \cr
#' (\emph{Description Output}: The following information about every bicluster generated will be printed in the output file: number of rows, number of columns, name of rows and name of columns.
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
bibit <- function(matrix=NULL,minr=2,minc=2,arff_row_col=NULL,output_path=NULL){
  
  pm <- match.call()
  
  
  if(is.null(arff_row_col)){
    time_arff <- round(proc.time()['elapsed']/60,2)
    
    # Check if matrix is binary (DISCRETIZED NOT YET IMPLEMENTED!)
    if(class(matrix)!="matrix"){stop("matrix parameter should contain a matrix object",call.=FALSE)}
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
    
    cat("DONE\n")
    cat("\n")
    
    time_arff <- round(proc.time()['elapsed']/60-time_arff,2)
    
    
  }else{
    time_arff <- 0
    
    if(length(arff_row_col)!=3){stop("arff_row_col should contain 3 elements",call.=FALSE)}
    bibitdata_path <- arff_row_col[1]
    bibitrows_path <- arff_row_col[2]
    bibitcols_path <- arff_row_col[3]
    
  }

  if(is.null(output_path)){
    bibitoutput_path <- tempfile("bibitoutput",fileext = "")
  }else{
    bibitoutput_path <- output_path
  }

  
  time_bibit <- proc.time()['elapsed']/60
  
  javaloc <- paste0(.libPaths(),"/BiBitR/java/BiBit.jar")
  # javaloc <- gsub("/","\\\\",javaloc)
  
  # BiBit.jar location needs to be standardized for package location! # .libPaths()
  # command <- paste("java -jar -Xmx1000M",javaloc,bibitdata_path,"1",minr,minc,bibitoutput_path,bibitrows_path,bibitcols_path,1)
  command <- paste("java -jar -Xmx1000M",paste0("\"",javaloc,"\""),paste0("\"",bibitdata_path,"\""),"1",minr,minc,paste0("\"",bibitoutput_path,"\""),paste0("\"",bibitrows_path,"\""),paste0("\"",bibitcols_path,"\""),1)
  
  system(command)
  
  time_bibit <- round(proc.time()['elapsed']/60-time_bibit,2)
  
  
  if(is.null(output_path)){
    cat("\n")
    cat("Transforming into biclust output...")
    
    time_biclust <- round(proc.time()['elapsed']/60,2)
    result <- bibit2biclust(data=matrix,resultpath=paste0(bibitoutput_path,"_1.txt"),arff_row_col = arff_row_col)
    cat("DONE\n")
    time_biclust <- round(proc.time()['elapsed']/60-time_biclust,2)
    

    if(!is.null(result)){
      result2 <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                     RowxNumber=result$RowxNumber,
                     NumberxCol=result$NumberxCol,
                     Number=result$Number,
                     info=list(Time_Min=list(arff=time_arff,bibit=time_bibit,biclust=time_biclust,full=time_arff+time_bibit+time_biclust)))
      
    }else{
      
      if(!is.null(arff_row_col)){
        rownames.data <- as.character(read.table(arff_row_col[2],header=FALSE)[,1])
        colnames.data <- as.character(read.table(arff_row_col[3],header=FALSE)[,1])
        nrow.data <- length(rownames.data)
        ncol.data <- length(colnames.data)
      }else{
        nrow.data <- nrow(matrix)
        ncol.data <- ncol(matrix)
      }
      
      
      result2 <- new("Biclust",Parameters=list(Call=pm,Method="BiBit"),
                     RowxNumber=matrix(FALSE,nrow=nrow.data,ncol=1),
                     NumberxCol=matrix(FALSE,nrow=1,ncol=ncol.data),
                     Number=0,
                     info=list(Time_Min=list(arff=time_arff,bibit=time_bibit,biclust=time_biclust,full=time_arff+time_bibit+time_biclust)))
    }

    return(result2)
    
  }else{
    return(NULL)
  }

}



bibit2biclust <- function(data,resultpath,arff_row_col){
  result <- read.table(resultpath,header=TRUE,sep=";")
  
  if(is.null(arff_row_col)){
    rownames.data <- rownames(data)
    colnames.data <- colnames(data)
    nrow.data <- nrow(data)
    ncol.data <- ncol(data)
  }else{
    rownames.data <- as.character(read.table(arff_row_col[2],header=FALSE)[,1])
    colnames.data <- as.character(read.table(arff_row_col[3],header=FALSE)[,1])
    nrow.data <- length(rownames.data)
    ncol.data <- length(colnames.data)
  }
  
  
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
    
    
    rowlist_index <- lapply(rowlist,FUN=function(x){rownames.data %in%  x})
    collist_index <- lapply(collist,FUN=function(x){colnames.data %in%  x})
    
    RowxNumber <- matrix(unlist(rowlist_index),byrow=FALSE,nrow=nrow.data,ncol=Number)
    NumberxCol <- matrix(unlist(collist_index),byrow=TRUE,nrow=Number,ncol=ncol.data)
    
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
#' @param top The number of top row/col/size dimension which are searched for. (e.g. default \code{top=1} gives only the maximum)
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
MaxBC <- function(bicresult,top=1){
  if(class(bicresult)!="Biclust"){stop("bicresult needs to be of class 'Biclust'")}
  
  rowsum <- colSums(bicresult@RowxNumber)
  colsum <- rowSums(bicresult@NumberxCol)
  sizesum <- rowsum*colsum
  
  top.col <- sort(unique(colsum),decreasing=TRUE)[1:top]
  top.row <- sort(unique(rowsum),decreasing=TRUE)[1:top]
  top.size <- sort(unique(sizesum),decreasing=TRUE)[1:top]
  
  
  for(i in 1:top){
   
    ind.colmax <- which(top.col[i]==colsum)
    ind.rowmax <- which(top.row[i]==rowsum)
    ind.sizemax <- which(top.size[i]==sizesum)
    
    if(i==1){
      
      row <- rbind(RowDim=rowsum[ind.rowmax],ColDim=colsum[ind.rowmax],SizeDim=sizesum[ind.rowmax])
      colnames(row) <- paste0("BC",ind.rowmax)
      
      column <- rbind(RowDim=rowsum[ind.colmax],ColDim=colsum[ind.colmax],SizeDim=sizesum[ind.colmax])
      colnames(column) <- paste0("BC",ind.colmax)
      
      size <- rbind(RowDim=rowsum[ind.sizemax],ColDim=colsum[ind.sizemax],SizeDim=sizesum[ind.sizemax])
      colnames(size) <- paste0("BC",ind.sizemax)
      
    }else{
      row.temp <- rbind(RowDim=rowsum[ind.rowmax],ColDim=colsum[ind.rowmax],SizeDim=sizesum[ind.rowmax])
      colnames(row.temp) <- paste0("BC",ind.rowmax)
      
      column.temp <- rbind(RowDim=rowsum[ind.colmax],ColDim=colsum[ind.colmax],SizeDim=sizesum[ind.colmax])
      colnames(column.temp) <- paste0("BC",ind.colmax)
      
      size.temp <- rbind(RowDim=rowsum[ind.sizemax],ColDim=colsum[ind.sizemax],SizeDim=sizesum[ind.sizemax])
      colnames(size.temp) <- paste0("BC",ind.sizemax)
    
      row <- cbind(row,row.temp)
      column <- cbind(column,column.temp)
      size <- cbind(size,size.temp)
    }
    
  }
  

  
  return(list(row=row,column=column,size=size))
}






#' Transform R matrix object to BiBit input files.
#' 
#' Transform the R matrix object to 1 \code{.arff} for the data and 2 \code{.csv} files for the row and column names. These are the 3 files required for the original BiBit Java algorithm
#' The path of these 3 files can then be used in the \code{arff_row_col} parameter of the \code{bibit} function.
#' 
#' @author Ewoud De Troyer
#' 
#' @export
#' @param matrix The binary input matrix.
#' @param name Basename for the 3 input files.
#' @param path Directory path where to write the 3 input files to.
#' 
#' @return 3 input files for BiBit:
#' \itemize{
#' \item 1 \code{.arff} file containing the data.
#' \item 1 \code{.csv} file for the row names. The file contains 1 column of names without quotation.
#' \item 1 \code{.csv} file for the column names. The file contains 1 column of names without quotation.
#' }
#' 
#' @examples 
#' \dontrun{
#' data <- matrix(sample(c(0,1),100*100,replace=TRUE,prob=c(0.9,0.1)),nrow=100,ncol=100)
#' data[1:10,1:10] <- 1 # BC1
#' data[11:20,11:20] <- 1 # BC2
#' data[21:30,21:30] <- 1 # BC3
#' data <- data[sample(1:nrow(data),nrow(data)),sample(1:ncol(data),ncol(data))]
#' 
#' make_arff_row_col(matrix=data,name="data",path="")
#' }
make_arff_row_col <- function(matrix,name="data",path=""){
  if(class(matrix)!="matrix"){stop("matrix parameter should contain a matrix object",call.=FALSE)}
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
  
  write.arff(t(matrix),file=paste0(getwd(),"/",path,"/",name,"_arff.arff"))
  write.table(matrix(rownames(matrix),ncol=1),quote=FALSE,row.names=FALSE,col.names=FALSE,file=paste0(getwd(),"/",path,"/",name,"_rownames.csv"))
  write.table(matrix(colnames(matrix),ncol=1),quote=FALSE,row.names=FALSE,col.names=FALSE,file=paste0(getwd(),"/",path,"/",name,"_colnames.csv"))
  
}

