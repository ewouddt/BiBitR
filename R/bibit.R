## IMPORTS ##

#' @importFrom foreign write.arff read.arff
#' @import biclust
#' @importFrom methods new
#' @importFrom utils read.table write.table combn
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
#' \emph{Note:} These files can be generated with the \code{\link{make_arff_row_col}} function.
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
#' MaxBC(result)
#' }
bibit <- function(matrix=NULL,minr=2,minc=2,arff_row_col=NULL,output_path=NULL){
  
  pm <- match.call()
  
  
  if(is.null(arff_row_col)){
    time_arff <- round(proc.time()['elapsed']/60,2)
    
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
  
  javaloc <- paste0(.libPaths()[1],"/BiBitR/java/BiBit.jar")
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
#' \emph{Note:} These files can be generated with the \code{\link{make_arff_row_col}} function.
#' @param output_path If as output, the original txt output of the Java code is desired, provide the outputh path here (without extension). In this case the \code{bibit} function will skip the transformation to a Biclust class object and simply return \code{NULL}.\cr 
#' (\emph{Example}: \code{output_path="...\\\\out\\\\bibitresult"})
#' \cr
#' (\emph{Description Output}: The following information about every bicluster generated will be printed in the output file: number of rows, number of columns, name of rows and name of columns.
#' @param extend_columns (EXPERIMENTAL!) Boolean value which applies a column extension procedure to the result of the BiBit algorithm. Columns will be sequentially added, keeping the noise beneath the allowed level. The procedure is the same as in \code{\link{bibit3}}, but now no artificial rows have to be ignored in the noise levels. 
#' \cr Note: The \code{@info} slot will also contain a \code{BC.Extended} value which contains the indices of which Biclusters's columns were extended.
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
#' MaxBC(result1,top=1)
#' 
#' result2 <- bibit2(data,minr=5,minc=5,noise=3)
#' result2
#' MaxBC(result2,top=2)
#' }
bibit2 <- function(matrix=NULL,minr=2,minc=2,noise=0,arff_row_col=NULL,output_path=NULL,extend_columns=FALSE){
  
  pm <- match.call()
  
  
  if(noise<0){stop("noise parameter can not be negative",call.=FALSE)}
  if(noise>=1){noise <- as.integer(noise)}
  
  if(is.null(arff_row_col)){
    time_arff <- round(proc.time()['elapsed']/60,2)
    
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
  
  javaloc <- paste0(.libPaths()[1],"/BiBitR/java/BiBit2.jar")
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
      
      if(extend_columns){
        result2 <- BC_column_extension(result=result2,data=matrix,noise=noise)
      }
      
      
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




#' @title The BiBit Algorithm with Noise Allowance guided by Provided Patterns.
#' 
#' @description Same function as \code{\link{bibit2}} but only aims to discover biclusters containing the (sub) pattern of provided patterns or their combinations.
#' @details The goal of the \code{\link{bibit3}} function is to provide one or multiple patterns in order to only find those biclusters exhibiting those patterns.
#' Multiple patterns can be given in matrix format, \code{pattern_matrix}, and their pairwise combinations can automatically be added to this matrix by setting \code{pattern_combinations=TRUE}.
#' All discovered biclusters are still subject to the provided \code{noise} level.
#' 
#' Three types of Biclusters can be discovered:
#' \describe{
#' \item{\emph{Full Pattern: }}{Bicluster which overlaps completely (within allowed noise levels) with the provided pattern. The column size of this bicluster is always equal to the number of 1's in the pattern.}
#' \item{\emph{Sub Pattern: }}{Biclusters which overlap with a part of the provided pattern within allowed noise levels. Will only be given if \code{subpattern=TRUE} (default). Setting this option to \code{FALSE} decreases computation time.}
#' \item{\emph{Extended: }}{Using the resulting biclusters from the full and sub patterns, other columns will be attempted to be added to the biclusters while keeping the noise as low as possible (the number of rows in the BC stays constant). Naturally the articially added pattern rows will not be taken into account with the noise levels as they are 0 in each other column.
#' \cr The question which is attempted to be answered here is \emph{`Do the rows, which overlap partly or fully with the given pattern, have other similarities outside the given pattern?`}}
#' } 
#' 
#' \emph{How?}
#' \cr The BiBit algorithm is applied to a data matrix that contains 2 identical artificial rows at the top which contain the given pattern. 
#' The default algorithm is then slightly altered to only start from this articial row pair (=Full Pattern) or from 1 artificial row and 1 other row (=Sub Pattern).
#' 
#' \emph{Note 1 - Large Data:}
#' \cr The \code{arff_row_col} can still be provided in case of large data matrices, but the \code{.arff} file should already contain the pattern of interest in the first two rows. Consequently not more than 1 pattern at a time can be investigated with a single call of \code{bibit3}.
#' 
#' \emph{Note 2 - Viewing Results:}
#' \cr A \code{print} and \code{summary} method has been implemented for the output object of \code{bibit3}. It gives an overview of the amount of discovered biclusters and their dimensions
#' \cr Additionally, the \code{\link{bibit3_patternBC}} function can extract a Bicluster and add the artificial pattern rows to investigate the results.
#' 
#' @author Ewoud De Troyer
#' 
#' @references Domingo S. Rodriguez-Baena, Antonia J. Perez-Pulido and Jesus S. Aguilar-Ruiz (2011), "A biclustering algorithm for extracting bit-patterns from binary datasets", \emph{Bioinformatics}
#' 
#' @export
#' @param matrix The binary input matrix.
#' @param minr The minimum number of rows of the Biclusters. (Note that in contrast to \code{\link{bibit}} and \code{\link{bibit2}}, this can be be set to 1 since we are looking for additional rows to the provided pattern.)
#' @param minc The minimum number of columns of the Biclusters.
#' @param noise Noise parameter which determines the amount of zero's allowed in the bicluster (i.e. in the extra added rows to the starting row pair).
#' \itemize{
#' \item \code{noise=0}: No noise allowed. This gives the same result as using the \code{\link{bibit}} function.
#' \item \code{0<noise<1}: The \code{noise} parameter will be a noise percentage. The number of allowed 0's in a (extra) row in the bicluster will depend on the column size of the bicluster. 
#' More specifically \code{zeros_allowed = ceiling(noise * columnsize)}. For example for \code{noise=0.10} and a bicluster column size of \code{5}, the number of allowed 0's would be \code{1}.
#' \item \code{noise>=1}: The \code{noise} parameter will be the number of allowed 0's in a (extra) row in the bicluster independent from the column size of the bicluster. In this noise option, the noise parameter should be an integer.
#' }
#' @param pattern_matrix Matrix (Number of Patterns x Number of Data Columns) containing the patterns of interest.
#' @param subpattern Boolean value if sub patterns are of interest as well (default=TRUE).
#' @param extend_columns Boolean value if columns of Biclusters should also be extended for additional results (default=TRUE). See Details Section for more info.
#' @param pattern_combinations Boolean value if the pairwise combinations of patterns (the intersecting 1's) should also used as starting points (default=FALSE).
#' @param arff_row_col Same argument as in \code{\link{bibit}} and \code{\link{bibit2}}. However you can only provide 1 pattern by using this option. For \code{bibit3} to work, the pattern has to be added 2 times on top of the matrix (= identical first 2 rows).
#' @return A S3 list object, \code{"bibit3"} in which each element (apart from the last one) corresponds with a provided pattern or combination thereof. \cr
#' Each element is a list containing:
#' \describe{
#' \item{\code{Number}: }{Number of Initially found BC's by applying BiBit with the provided pattern.} 
#' \item{\code{Number_Extended}: }{Number of additional discovered BC's by extending the columns.}
#' \item{\code{FullPattern}: }{Biclust S4 Class Object containing the Bicluster with the Full Pattern.}
#' \item{\code{SubPattern}: }{Biclust S4 Class Object containing the Biclusters showing parts of the pattern.}
#' \item{\code{Extended}: }{Biclust S4 Class Object containing the additional Biclusters after extending the biclusters (column wise) of the full and sub patterns}
#' \item{\code{info}: }{Contains \code{Time_Min} element which includes the elapsed time of parts and the full analysis.}
#' }
#' The last element in the list is a matrix containing all the investigated patterns.
#' 
#' @examples 
#' \dontrun{ 
#' set.seed(1)
#' data <- matrix(sample(c(0,1),100*100,replace=TRUE,prob=c(0.9,0.1)),nrow=100,ncol=100)
#' data[1:10,1:10] <- 1 # BC1
#' data[11:20,11:20] <- 1 # BC2
#' data[21:30,21:30] <- 1 # BC3
#' colsel <- sample(1:ncol(data),ncol(data))
#' data <- data[sample(1:nrow(data),nrow(data)),colsel]
#' 
#' pattern_matrix <- matrix(0,nrow=3,ncol=100)
#' pattern_matrix[1,1:7] <- 1
#' pattern_matrix[2,11:15] <- 1
#' pattern_matrix[3,13:20] <- 1
#' 
#' pattern_matrix <- pattern_matrix[,colsel]
#' 
#' 
#' out <- bibit3(matrix=data,minr=2,minc=2,noise=0.1,pattern_matrix=pattern_matrix,
#'               subpattern=TRUE,extend_columns=TRUE,pattern_combinations=TRUE)
#' out  # OR print(out) OR summary(out)
#' 
#' 
#' bibit3_patternBC(result=out,matrix=data,pattern=c(1),type=c("full","sub","ext"),BC=c(1,2))
#' }
bibit3 <- function(matrix=NULL,minr=1,minc=2,noise=0,pattern_matrix=NULL,subpattern=TRUE,extend_columns=TRUE,pattern_combinations=FALSE,arff_row_col=NULL){
  
  pm <- match.call()
  minr <- minr + 2
  
  ###
  if(noise<0){stop("noise parameter can not be negative",call.=FALSE)}
  if(noise>=1){noise <- as.integer(noise)}
  
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
      cat("Transform matrix into arff format...",rownames(pattern_matrix)[i.pattern],"...")

      bibitdata_path <- tempfile("bibitdata",fileext=".arff")
      bibitrows_path <- tempfile("bibitrows",fileext=".csv")
      bibitcols_path <- tempfile("bibitcols",fileext=".csv")

      write.arff(t(matrix_with_pattern),file=bibitdata_path)
      write.table(matrix(rownames(matrix_with_pattern),ncol=1),quote=FALSE,row.names=FALSE,col.names=FALSE,file=bibitrows_path)
      write.table(matrix(colnames(matrix_with_pattern),ncol=1),quote=FALSE,row.names=FALSE,col.names=FALSE,file=bibitcols_path)

      cat("DONE\n")
      cat("\n")

      time_arff <- round(proc.time()['elapsed']/60-time_arff,2)
      
    }else{
      matrix_with_pattern <- NULL
      if(extend_columns){
        matrix_with_pattern <- t(foreign::read.arff(bibitdata_path))
      }
    }
    
    # Apply BiBit Algorithm
    
    cat("Initiate BiBit for",rownames(pattern_matrix)[i.pattern],"...\n")
    cat("\n")
    
    bibitoutput_path <- tempfile("bibitoutput",fileext = "")
    
    
    time_bibit <- proc.time()['elapsed']/60
    
    javaloc <- paste0(.libPaths()[1],"/BiBitR/java/BiBit3.jar")
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
      
      
      if(extend_columns){
        time_extend <- round(proc.time()['elapsed']/60,2)
        Extended <- BC_column_extension_pattern(result=result2,data=matrix_with_pattern,noise=noise)
        
        
        time_extend <- round(proc.time()['elapsed']/60-time_extend,2)
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
    for(j.list in c("FullPattern","SubPattern","Extended")){
      
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





BC_column_extension_pattern <- function(result,data,noise){
  
  
  BC_extended <- rep(FALSE,result@Number)
  
  for(i.BC in 1:result@Number){
    included_columns <- result@NumberxCol[i.BC,]
    
    column_candidates <- order(colSums(data[result@RowxNumber[,i.BC],]),decreasing=TRUE)
    
    GO <- TRUE
    i.candidate <- 1
    
    
    while(GO & (i.candidate<=length(column_candidates))){

      if(!included_columns[column_candidates[i.candidate]]){
        
        included_columns_temp <- included_columns
        included_columns_temp[column_candidates[i.candidate]] <- TRUE
        
        zeros_allowed <- ifelse(((noise<1)&(noise>0)),ceiling(noise*sum(included_columns_temp)),noise)
        
        zeros_in_rows_withoutpattern <- apply(data[result@RowxNumber[,i.BC],included_columns_temp],MARGIN=1,FUN=function(x){sum(x==0)})[-c(1,2)]
        
        if(all(zeros_in_rows_withoutpattern<=zeros_allowed)){
          
          included_columns <- included_columns_temp
          i.candidate <- i.candidate+1
          
          
        }else{
          GO <- FALSE
        }
      }else{
        i.candidate <- i.candidate+1
      }
      
    }
    
    if(sum(included_columns)>sum(result@NumberxCol[i.BC,])){BC_extended[i.BC] <- TRUE}
    result@NumberxCol[i.BC,] <- included_columns
    
  }
  
  if(sum(BC_extended)>=1){
    result@RowxNumber <- result@RowxNumber[,BC_extended,drop=FALSE]
    result@NumberxCol <- result@NumberxCol[BC_extended,,drop=FALSE]
    result@Number <- sum(BC_extended)
    result@info <- list()
  }else{
    nrow.data <- nrow(result@RowxNumber)
    ncol.data <- ncol(result@NumberxCol)
    
    result <- new("Biclust",Parameters=result@Parameters,
                  RowxNumber=matrix(FALSE,nrow=nrow.data,ncol=1),
                  NumberxCol=matrix(FALSE,nrow=1,ncol=ncol.data),
                  Number=0,
                  info=list()) 
  }

  
  return(result)
}



BC_column_extension <- function(result,data,noise){
  
  time_extend <- round(proc.time()['elapsed']/60,2)
  
  
  BC_extended <- rep(FALSE,result@Number)
  
  for(i.BC in 1:result@Number){
    included_columns <- result@NumberxCol[i.BC,]
    
    column_candidates <- order(colSums(data[result@RowxNumber[,i.BC],]),decreasing=TRUE)
    
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
    
    if(sum(included_columns)>sum(result@NumberxCol[i.BC,])){BC_extended[i.BC] <- TRUE}
    result@NumberxCol[i.BC,] <- included_columns
    
  }
  time_extend <- round(proc.time()['elapsed']/60-time_extend,2)
  
  result@info$Time_Min$extend <- time_extend
  result@info$Time_Min$full <- result@info$Time_Min$full + time_extend
  result@info$BC.Extended <- which(BC_extended)
  
  # info=list(Time_Min=list(arff=time_arff,bibit=time_bibit,biclust=time_biclust,full=time_arff+time_bibit+time_biclust)))

  
  return(result)
}


