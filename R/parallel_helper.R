
simmat_par_double <- function(x,worker_rows,BCresult1,BCresult2,type){
  
  worker_rows <- worker_rows[[x]]
  
  simmat <- matrix(0,nrow=length(worker_rows),ncol=BCresult2@Number,dimnames=list(paste0("Result1_BC",worker_rows),paste0("Result2_BC",1:BCresult2@Number)))
  
  if(type=="both"){
    for(i in 1:nrow(simmat)){
      for(j in 1:ncol(simmat)){
        row_contain_temp <- sum(which(BCresult1@RowxNumber[,worker_rows[i]])%in%which(BCresult2@RowxNumber[,j]))
        col_contain_temp <- sum(which(BCresult1@NumberxCol[worker_rows[i],])%in%which(BCresult2@NumberxCol[j,]))
        m1 <- sum(BCresult1@RowxNumber[,worker_rows[i]])*sum(BCresult1@NumberxCol[worker_rows[i],])
        m2 <- sum(BCresult2@RowxNumber[,j])*sum(BCresult2@NumberxCol[j,])
        m12 <- m1+m2-row_contain_temp*col_contain_temp
        simmat[i,j] <- (m1+m2-(m12))/(m12)
      }
    }
  }else if(type=="row"){
    for(i in 1:nrow(simmat)){
      for(j in 1:ncol(simmat)){
        x1 <- BCresult1@RowxNumber[,worker_rows[i]]
        x2 <- BCresult2@RowxNumber[,j]
        m1 <- sum(x1)
        m2 <- sum(x2)
        m12 <- sum(as.logical(x1+x2))
        simmat[i,j] <- (m1+m2-m12)/m12
      }
    }
  }else if(type=="col"){
    for(i in 1:nrow(simmat)){
      for(j in 1:ncol(simmat)){
        x1 <- BCresult1@NumberxCol[worker_rows[i],]
        x2 <- BCresult2@NumberxCol[j,]
        m1 <- sum(x1)
        m2 <- sum(x2)
        m12 <- sum(as.logical(x1+x2))
        simmat[i,j] <- (m1+m2-m12)/m12
      }
    }
  }
  
  return(simmat)
}




simmat_par_single <- function(x,worker_rows,BCresult1,type){
  
  workflow_jaccard_bc <- function(result,BC1,BC2,type="both"){
    if(type=="both"){
      row_contain_temp <- sum(which(result@RowxNumber[,BC1])%in%which(result@RowxNumber[,BC2]))
      col_contain_temp <- sum(which(result@NumberxCol[BC1,])%in%which(result@NumberxCol[BC2,]))
      m1 <- sum(result@RowxNumber[,BC1])*sum(result@NumberxCol[BC1,]) 
      m2 <- sum(result@RowxNumber[,BC2])*sum(result@NumberxCol[BC2,])
      m12 <- m1+m2-row_contain_temp*col_contain_temp
      JI <- (m1+m2-(m12))/(m12)
      return(JI)
    }else{
      if(type=="row"){
        x1 <- result@RowxNumber[,BC1]
        x2 <- result@RowxNumber[,BC2]
      }
      if(type=="col"){
        x1 <- result@NumberxCol[BC1,]
        x2 <- result@NumberxCol[BC2,]
      }
      m1 <- sum(x1)
      m2 <- sum(x2)
      m12 <- sum(as.logical(x1+x2))
      JI <- (m1+m2-m12)/m12
      return(JI)
    }
  }
  
  
  final_job <- length(worker_rows)
  worker_rows <- worker_rows[[x]]
  simmat <- matrix(0,nrow=length(worker_rows),ncol=BCresult1@Number,dimnames=list(paste0("Result1_BC",worker_rows),paste0("Result1_BC",1:BCresult1@Number)))
  
  if(x==final_job){
    for(i in 1:(nrow(simmat)-1)){
      for(j in (worker_rows[[i]]+1):(ncol(simmat))){
        simmat[i,j] <- workflow_jaccard_bc(BCresult1,worker_rows[i],j,type=type)
      }
    }
  }else{
    for(i in 1:nrow(simmat)){
      for(j in (worker_rows[[i]]+1):(ncol(simmat))){
        simmat[i,j] <- workflow_jaccard_bc(BCresult1,worker_rows[i],j,type=type)
      }
    }
  }

  return(simmat)
}



simmat_par_single_bibitworkflow <- function(x,worker_rows,BCresult1,type){
  
  workflow_jaccard_bc <- function(result,BC1,BC2,type="both"){
    if(type=="both"){
      row_contain_temp <- sum(which(result@RowxNumber[,BC1])%in%which(result@RowxNumber[,BC2]))
      col_contain_temp <- sum(which(result@NumberxCol[BC1,])%in%which(result@NumberxCol[BC2,]))
      m1 <- sum(result@RowxNumber[,BC1])*sum(result@NumberxCol[BC1,]) 
      m2 <- sum(result@RowxNumber[,BC2])*sum(result@NumberxCol[BC2,])
      m12 <- m1+m2-row_contain_temp*col_contain_temp
      JI <- (m1+m2-(m12))/(m12)
      return(JI)
    }else{
      if(type=="row"){
        x1 <- result@RowxNumber[,BC1]
        x2 <- result@RowxNumber[,BC2]
      }
      if(type=="col"){
        x1 <- result@NumberxCol[BC1,]
        x2 <- result@NumberxCol[BC2,]
      }
      m1 <- sum(x1)
      m2 <- sum(x2)
      m12 <- sum(as.logical(x1+x2))
      JI <- (m1+m2-m12)/m12
      return(JI)
    }
  }
  
  
  final_job <- length(worker_rows)
  worker_rows <- worker_rows[[x]]
  simmat <- matrix(0,nrow=length(worker_rows),ncol=BCresult1@Number)
  
  if(x==final_job){
    for(i in 1:(nrow(simmat)-1)){
      for(j in (worker_rows[[i]]+1):(ncol(simmat))){
        simmat[i,j] <- workflow_jaccard_bc(BCresult1,worker_rows[i],j,type=type)
      }
    }
  }else{
    for(i in 1:nrow(simmat)){
      for(j in (worker_rows[[i]]+1):(ncol(simmat))){
        simmat[i,j] <- workflow_jaccard_bc(BCresult1,worker_rows[i],j,type=type)
      }
    }
  }
  
  return(simmat)
}


