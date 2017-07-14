
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


simmat_par_single <- function(){}