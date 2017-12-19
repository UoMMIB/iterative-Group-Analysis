#ACCURATE VERSION

#function for the "iterative Group Analysis" (iGA).
# for the "double boundary iterative group analysis" (db-iGA) there is another
# function. The main idea is to build a wrapper function that, using various
# parameters calls this function or the other.

#inputs:
### metric: (N x 1) vector containing the metric used to evaluate the
###         significance of the variables, not necessarly ordered
### group.membership: (N x M) binary matrix containing the group membership
###                   if the (i,j) element is equal to 1 it means that the ith
###                   variable is a member of the jth group
### var.names: an optional input containing the names of the variables
### groups: (M x 1) vector containing the names of the groups
# outputs:
### 


"iGA_acc" <- function(metric, group.membership, groups, var.names = NULL, decreasing=FALSE){
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
        group.elements[g] <- length(which(group.membership[,g]==1))
    }
    
    # ordering metrics, group.mebership and var.names
    ord <- order(metric, decreasing= decreasing)
    metric <- metric[ord]
    group.membership <- group.membership[ord,]
    group.membership <- t(t(group.membership))
    if (!is.null(var.names)) var.names <- var.names[ord]
    rm(ord)
    
    PC.list <- list() #this list will contain a vector per each group,
                      #such vector contains all the PC-values
    min.PCs <- rep(NA, N.groups) # this vector will contain the min PC-value for
                                 # per each group
    min.PCs.pos <- rep(NA,N.groups) # this vector will contain number of group
                                   # elements used for the min PC-value
    
    # for cycle that considers each group
    for (g in 1:N.groups){
        x <- group.elements[g]  # number of elements in the gth group
        position <- 1:N         # vector from which we extract t
        PC.vec <- rep(NA, x)
        tmp <- group.membership[,g]
        store.t<-rep(NA, x)
        for (z in 1:x){
            next.group.member <- which(tmp==1)[1]
            t <- position[next.group.member]
            p <- iga_acc(z,N,t,x)
            PC.vec[z] <- as.numeric(p)
            if (next.group.member < length(tmp)){
                tmp <- tmp[(next.group.member+1):length(tmp)]
                position <- position[(next.group.member+1):length(position)]
            }
            store.t[z] <-t
        }
        PC.list[[g]] <- PC.vec
        min.PCs[g] <- min(PC.vec)
        min.PCs.pos[g] <- store.t[which(PC.vec==min(PC.vec))[1]]        
    }
    
#preparing the output
var.selected <- rep(NA,N.groups)
if (!is.null(var.names)) var.selected.names <- list()
#####
for(g in 1:N.groups){
    ind<-which(group.membership[1:min.PCs.pos[g],g] == 1)
    var.selected[g] <- length(ind)
    if (!is.null(var.names)) var.selected.names[[g]] <- var.names[ind] 
}
summary <- matrix(NA, N.groups, 4)
rownames(summary) <- groups
summary[,1] <- min.PCs
summary[,2] <- min.PCs.pos
summary[,3] <- var.selected
summary[,4] <- group.elements
colnames(summary) <- c("minPC","list.position", "N.var.selected","N.var.group")
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
