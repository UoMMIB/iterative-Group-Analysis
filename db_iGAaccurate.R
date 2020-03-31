#ACCURATE VERSION
#function for the "double boundary iterative Group Analysis" (iGA).
# The main idea is to build a wrapper function that, using various
# parameters calls this function or the other.


#inputs:
### metric: (N x 1) vector containing the metric used to evaluate the
###         significance of the variables, not necessarly ordered
### group.membership: (N x M) binary matrix containing the group membership
###                   if the (i,j) element is equal to 1 it means that the ith
###                   variable is a member of the jth group
### var.names: an optional input containing the names of the variables
### groups: (M x 1) vector containing the names of the groups



"db.iGA_acc" <- function(metric, group.membership, groups, var.names = NULL, decreasing=FALSE){
    library(Rmpfr) 
    # cheching inputs
    if(length(metric) != dim(group.membership)[1]){
        cat('\n wrong inputs')
        stop()
    }
    if (!is.null(var.names) & (length(metric) != length(var.names))){
        cat('\n wrong var.names')
        stop()
    }
    
    # counting number of groups
    N.groups <- dim(group.membership)[2]
    N <- length(metric)
    
    # counting the number of elements for each group
    group.elements <- rep(NA, N.groups)
    for (g in 1:N.groups){
        group.elements[g] <- length(which(group.membership[,g] == 1))
    }
    
    # ordering metrics, group.mebership and var.names
    ord <- order(metric, decreasing= decreasing)
    metric <- metric[ord]
    group.membership <- group.membership[ord,]
    if (!is.null(var.names)) var.names <- var.names[ord]
    rm(ord)
    
    PC.list <- list() #this list will contain a vector per each group,
    #such vector contains all the PC-values
    min.PCs <- rep(NA, N.groups) # this vector will contain the min PC-value for
    # per each group
    min.PCs.pos <- matrix(NA, N.groups, 2)
    colnames(min.PCs.pos) <- c("start", "end")
    rownames(min.PCs.pos) <- groups
    for (g in 1:N.groups){
        x <- group.elements[g]
        PC.mat <- matrix(NA,x,N)
        position <- 1:N
        rows<-position[which(group.membership[,g] == 1)]
        for(start in 1:N){
            position <- 1:N
            tmp <- group.membership[,g]
            TMP <-tmp
            for(z in 1:x){
            next.group.member <- which(tmp == 1)[1]
            t <- position[next.group.member]
            bef.start <- 0
            if(start>1) bef.start <- length(which(TMP[1:(start-1)]== 1))
            Z <- z-bef.start
            if (t<start){
                p <- NA
            }else{
                 #p <- hyper.geom(Z,N,t-start+1,x)
                 p <- iga_acc(Z, N, t-start+1, x)
            }            
            if (next.group.member < length(tmp)){
                tmp <- tmp[(next.group.member+1):length(tmp)]
                position <- position[(next.group.member+1):length(position)]
            }
            PC.mat[z,start]<- as.numeric(p)
          }
        }
        row.names(PC.mat) <- rows
    PC.list[[g]] <- PC.mat
    min.PCs[g] <- min(PC.mat,na.rm=TRUE)
    pos <- which(PC.mat==min(PC.mat, na.rm = TRUE), arr.ind = TRUE) 
    if(is.matrix(pos)){
     min.PCs.pos[g,1] <- as.numeric(pos[1,2])
     min.PCs.pos[g,2] <- as.numeric(rownames(PC.mat)[pos[1,1]])
     }else{
     min.PCs.pos[g,1] <- as.numeric(pos[2])
     min.PCs.pos[g,2] <- as.numeric(rownames(PC.mat)[pos[1]])    
     }
    }
    #preparing output
    var.selected <- rep(NA,N.groups)
    position <- 1:N
    if (!is.null(var.names)) var.selected.names <- list()
    for(g in 1:N.groups){
        ind<-which(group.membership[min.PCs.pos[g,1]:min.PCs.pos[g,2],g] == 1)
        ind<-position[min.PCs.pos[g,1]:min.PCs.pos[g,2]][ind]
        var.selected[g] <- length(ind)
        if (!is.null(var.names)) var.selected.names[[g]] <- var.names[ind] 
    }
    
    summary <- matrix(NA, N.groups, 5)
    rownames(summary) <- groups
    summary[,1] <- min.PCs
    summary[,2] <- min.PCs.pos[,1]
    summary[,3] <- min.PCs.pos[,2]
    summary[,4] <- var.selected
    summary[,5] <- group.elements
    colnames(summary) <- c("minPC","start", "end", "N.var.selected","N.var.group")
    

    
    if (is.null(var.names)){
        out <- list(PC.list = PC.list, minPCs = min.PCs, minPCs.pos = min.PCs.pos,
                    summary = summary)
    }else{
        out <- list(PC.list = PC.list, minPCs = min.PCs, minPCs.pos = min.PCs.pos,
                    summary = summary, var.sel.list = var.selected.names)    
    }
    out
    
    
    
}

"hyper.geom_acc" <- function(z,n,t,x){
    if((x-z)>(n-t)){
        p<-0
    }else{
        p <- (chooseZ(t,z)*chooseZ(n-t, x-z))/chooseZ(n, x)
    }
 p
}

"iga_acc" <- function(z,n,t,x){
    if (z==0){
        p<-1
    } else {
        p <- iga_acc(z-1, n,t,x) - hyper.geom_acc(z-1,n,t,x)
    }
 p
}




#"hyper.geom" <- function(z,n,t,x){
#    if (z<=0 | t==n | t==0){
#        p <- 1
#    } else{
#        p <- 0
#        for (i in 1:z-1){
#            p <- p + ((choose(t,i) * choose(n-t, x-i)) / choose(n,x))
#        }
#        p <- 1-p
#    }
#    p 
#}
