#' Finding Maximum Size Biclusters
#' 
#' Simple function which scans a \code{Biclust} result and returns which biclusters have maximum row, column or size (row*column).
#' 
#' @author Ewoud De Troyer
#' 
#' @export
#' @param bicresult A \code{Biclust} result. (e.g. The return object from \code{bibit} or \code{bibit2})
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
#' 
#' result <- bibit(data,minr=5,minc=5,
#'                 arff_row_col=c("data_arff.arff","data_rownames.csv","data_colnames.csv"))
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




#' Finding BC's with specific rows.
#' 
#' Simple function which scans a \code{Biclust} result and returns which biclusters contain all rows given in the \code{rows} parameter.
#' 
#' @author Ewoud De Troyer
#' 
#' @export
#' @param bicresult A \code{Biclust} result. (e.g. The return object from \code{bibit} or \code{bibit2})
#' @param rows A vector containing containing the row numbers which should be in the bicluster.
#' 
#' @return A matrix containing the biclusters in the columns and the row, column and size dimensions on the rows.
#' 
#' @examples 
#' \dontrun{
#' data <- matrix(sample(c(0,1),100*100,replace=TRUE,prob=c(0.9,0.1)),nrow=100,ncol=100)
#' data[1:10,1:10] <- 1 # BC1
#' data[11:20,11:20] <- 1 # BC2
#' data[21:30,21:30] <- 1 # BC3
#' result <- bibit(data,minr=2,minc=2)
#' 
#' rows_in_BC(result,rows=c(21,22,23))
#' 
#' }
rows_in_BC <- function(bicresult,rows){
  
  if(class(bicresult)!="Biclust"){stop("bicresult is not a Biclust class object",call.=FALSE)}
  
  BC.boolean <- sapply(1:bicresult@Number,FUN=function(x){
    return(all(rows%in%which(bicresult@RowxNumber[,x])))
  })
  
  BC.sel <- which(BC.boolean)
  
  rowdim <- colSums(bicresult@RowxNumber[,BC.sel,drop=FALSE])
  coldim <- rowSums(bicresult@NumberxCol[BC.sel,,drop=FALSE])
  sizedim <- rowdim*coldim
  
  out <- matrix(c(rowdim,coldim,sizedim),byrow=TRUE,nrow=3,ncol=length(BC.sel),dimnames=list(c("RowDim","ColDim","SizeDim"),paste0("BC",BC.sel)))
  out <- out[,order(sizedim,decreasing=TRUE)]
  
  return(out)
}






#' @title Finding BC's with specific rows which only 1's in the BC.
#' 
#' @description Simple function which scans a \code{Biclust} result and returns which biclusters contain all rows given in the \code{rows} parameter,
#' but only if these rows only contain 1's in the bicluster. This can be particularly helpful after having added articial row-pairs with a pattern of interest. 
#' With this function you can retrieve the biclusters that grew from these pairs from all the discovered biclusters.
#' 
#' @author Ewoud De Troyer
#' 
#' @export
#' @param matrix The binary input matrix.
#' @param bicresult A \code{Biclust} result. (e.g. The return object from \code{bibit} or \code{bibit2})
#' @param rows A vector containing containing the row numbers which should be in the bicluster.
#' 
#' @return A matrix containing the biclusters in the columns and the row, column and size dimensions on the rows.
#' 
#' @examples 
#' \dontrun{
#' data <- matrix(sample(c(0,1),100*100,replace=TRUE,prob=c(0.9,0.1)),nrow=100,ncol=100)
#' data[1:10,1:10] <- 1 # BC1
#' data[11:20,11:20] <- 1 # BC2
#' data[21:30,21:30] <- 1 # BC3
#' 
#' extra_rows <- rep(0,100)
#' extra_rows[11:25] <- 1
#' 
#' data <- rbind(data,rbind(extra_rows,extra_rows)) 
#' rownames(data) <- NULL
#' 
#' result <- bibit2(data,minr=2,minc=2,noise=0.2)
#' 
#' rows_full1_in_BC(matrix=data,bicresult=result,rows=c(101,102))
#' 
#' }
rows_full1_in_BC <- function(matrix,bicresult,rows){
  
  if(class(bicresult)!="Biclust"){stop("bicresult is not a Biclust class object",call.=FALSE)}
  if(class(matrix)!="matrix"){stop("matrix parameter should contain a matrix object",call.=FALSE)}
  
  
  BC.boolean <- sapply(1:bicresult@Number,FUN=function(x){
    
    if(all(rows%in%which(bicresult@RowxNumber[,x]))){
      
      submat <- matrix[rows,bicresult@NumberxCol[x,]]
      return(all(submat==1))
      
    }else{
      return(FALSE)
    }
  })
  
  BC.sel <- which(BC.boolean)
  
  rowdim <- colSums(bicresult@RowxNumber[,BC.sel,drop=FALSE])
  coldim <- rowSums(bicresult@NumberxCol[BC.sel,,drop=FALSE])
  sizedim <- rowdim*coldim
  
  out <- matrix(c(rowdim,coldim,sizedim),byrow=TRUE,nrow=3,ncol=length(BC.sel),dimnames=list(c("RowDim","ColDim","SizeDim"),paste0("BC",BC.sel)))
  out <- out[,order(sizedim,decreasing=TRUE)]
  
  return(out)
}


