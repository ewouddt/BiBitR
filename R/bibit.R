
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






#' @title The BiBit Algorithm with Noise Allowance
#' 
#' @description Same function as \code{\link{bibit}} with an additional new noise parameter which allows 0's in the discovered biclusters (See Details for more info).
#' 
#' @details \code{bibit2} follows the same steps as described in the Details section of \code{\link{bibit}}.\cr
#' Following the general steps of the BiBit algorithm, the allowance for noise in the biclusters is inserted in the original algorithm as such:
#' \enumerate{
#' \item Binary data is encoded in bit words.
#' \item Take a pair of rows as your starting point.
#' \item Find the maximal overlap of 1's between these two rows and save this as a pattern/motif. You now have a bicluster of 2 rows and N columns in which N is the number of 1's in the motif.
#' \item Check all remaining rows if they match this motif, \emph{however} allow a specific amount of 0's in this matching as defined by the \code{noise} parameter. Those rows that match completely or those within the allowed noise range are added to bicluster.
#' \item Go back to \emph{Step 2} and repeat for all possible row pairs.
#' }
#' \emph{Note:} Biclusters are only saved if they satisfy the \code{minr} and \code{minc} parameter settings and if the bicluster is not already contained completely within another bicluster.\cr
#' \cr
#' What you will end up with are biclusters not only consisting out of 1's, but biclusters in which 2 rows (the starting pair) are all 1's and in which the other rows could contain 0's (= noise).\cr
#' \cr
#' \emph{Note:} Because of the extra checks involved in the noise allowance, using noise might increase the computation time a little bit.
#' 
#' @author Ewoud De Troyer
#' 
#' @references Domingo S. Rodriguez-Baena, Antonia J. Perez-Pulido and Jesus S. Aguilar-Ruiz (2011), "A biclustering algorithm for extracting bit-patterns from binary datasets", \emph{Bioinformatics}
#' 
#' @export
#' @param matrix The binary input matrix.
#' @param minr The minimum number of rows of the Biclusters.
#' @param minc The minimum number of columns of the Biclusters.
#' @param noise Noise parameter which determines the amount of zero's allowed in the bicluster (i.e. in the extra added rows to the starting row pair).
#' \itemize{
#' \item \code{noise=0}: No noise allowed. This gives the same result as using the \code{\link{bibit}} function.
#' \item \code{0<noise<1}: The \code{noise} parameter will be a noise percentage. The number of allowed 0's in a (extra) row in the bicluster will depend on the column size of the bicluster. 
#' More specifically \code{zeros_allowed = ceiling(noise * columnsize)}. For example for \code{noise=0.10} and a bicluster column size of \code{5}, the number of allowed 0's would be \code{1}.
#' \item \code{noise>=1}: The \code{noise} parameter will be the number of allowed 0's in a (extra) row in the bicluster independent from the column size of the bicluster. In this noise option, the noise parameter should be an integer.
#' }
#' 
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
#' 
#' result1 <- bibit2(data,minr=5,minc=5,noise=0.2)
#' result1
#' result2 <- bibit2(data,minr=5,minc=5,noise=3)
#' result2
#' }
bibit2 <- function(matrix=NULL,minr=2,minc=2,noise=0,arff_row_col=NULL,output_path=NULL){
  
  pm <- match.call()
  
  
  if(noise<0){stop("noise parameter can not be negative",call.=FALSE)}
  if(noise>=1){noise <- as.integer(noise)}
  
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
  
  javaloc <- paste0(.libPaths(),"/BiBit2R/java/BiBit2.jar")
  # javaloc <- gsub("/","\\\\",javaloc)
  
  # BiBit.jar location needs to be standardized for package location! # .libPaths()
  command <- paste("java -jar -Xmx1000M",paste0("\"",javaloc,"\""),paste0("\"",bibitdata_path,"\""),"1",minr,minc,paste0("\"",bibitoutput_path,"\""),paste0("\"",bibitrows_path,"\""),paste0("\"",bibitcols_path,"\""),1,paste0(" ",noise))
  # cat(command,"\n")
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





