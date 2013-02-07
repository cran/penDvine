Dvine <- function(data,K=8,lambda=100,order.Dvine=TRUE,pen=1,base="Bernstein",m=2,cores=NULL,q=2) {
  if(require(doParallel)) doMC <- TRUE else doMC <- FALSE
    if(doMC) {
      library(doParallel)
      if(is.null(cores)) registerDoParallel() else registerDoParallel(cores=cores)
      mcoptions <- list(preschedule=FALSE)
    }

  #if(base=="Bernstein") q <- NULL
  if(!is.matrix(data)) data <- as.matrix(data)
  U <- data
S <- seq(1:(dim(U)[2]))
N <-  dim(U)[1]
k <-  1
Dvine <- list()
SSi <-  list()

#order options
#pairwise AIC-order in the first level
help.env <- new.env()
assign("doMC",doMC,help.env)
assign("S",S,help.env)
assign("U",U,help.env)
assign("lambda",lambda,help.env)
assign("K",K,help.env)
assign("base",base,help.env)
assign("m",m,help.env)

if(order.Dvine) {
  library(TSP)
  assign("order.stat","cAIC",help.env)
  order.Dvine(help.env)
  U <- U[,get("order",help.env)]
}
else assign("order",seq(1,length(S)),help.env)

for ( i in 1:length(S))
  {
    vine.knot <-  list(j1=S[i],j2=NULL,D=NULL,v=NULL, U = U[,S[i]] )
    SSi <- append(SSi, list(vine.knot))
  }
Dvine <- append(Dvine, list(SSi))

level <-  1
log.like <- AIC <- cAIC <- 0

SS <-  Dvine[[level]]
nSS <- length(SS)
while (nSS>1)
  {
    SSi <-  list()
    for (i in 2:nSS)
      {
        S1 <- c(SS[[i-1]]$j1, SS[[i-1]]$j2,SS[[i-1]]$D)
        S2 <- c(SS[[i]]$j1, SS[[i]]$j2,SS[[i]]$D)
        index1 <- rep(TRUE,length(S1))
        index2 <- rep(TRUE,length(S2))
        S3 <- c()
        for (j in 1:length(S1))
          {
            indexi <- S1[j]==S2
            if (sum(indexi)>0)
              {
                S3 <- c(S3, S1[j])
                index1[j] <- FALSE
                index2[indexi] <- FALSE                
              }
          }
        if(!is.null(S3)) S3 <- sort(S3)
        S1 <- S1[index1]
        S2 <- S2[index2]
        vine.knot <-  list(j1=S1,j2=S2,D=S3,v=NULL, U = NULL, log.like = NULL, AIC = NULL)
        SSi <- append(SSi, list(vine.knot))
      }
    Dvine <- append(Dvine, list(SSi))
    level <- level + 1
    SS <-  Dvine[[level]]
    nSS <- length(SS)
  }

level <-  2
Tree.l <- Dvine[[level]]

if(order.Dvine) {
  ord <- get("order",help.env)
  pairs <- get("pairs",help.env)
  if(!doMC) {
    Tree.l.temp <- list()
    for(i in 1:length(Tree.l)) {
      he <- ord[c(i,i+1)]
      #line <- c()
      #for(mm in 1:dim(pairs)[1])  {
      #  if(all(he%in%pairs[mm,])) line <- mm
      #  if(all(he[c(2,1)]%in%pairs[mm,])) line <- mm    
      #}
      index.j1 <- Tree.l[[i]]$j1
      index.j2 <- Tree.l[[i]]$j2 
      #model.l <- get("cal.order",help.env)[[line]]
      model.l <- paircopula(U[,c(he)],K=K,lambda=lambda,pen=pen,base=base,m=m,q=q)
      Tree.l.temp[[i]] <- list(j1=index.j1,j2=index.j2,D=NULL,v=get("ck.val",model.l),U=U[,c(index.j1,index.j2)],log.like = get("log.like",model.l),cop=model.l , AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l))
    }
  }
  if(doMC) {
    Tree.l.temp <- foreach(i=1:length(Tree.l),.combine=list,.multicombine=TRUE,.options.multicore=mcoptions) %dopar%   {
      he <- ord[c(i,i+1)]
      #line <- c()
      #for(mm in 1:dim(pairs)[1])  {
      #  if(all(he%in%pairs[mm,])) line <- mm
      #  if(all(he[c(2,1)]%in%pairs[mm,])) line <- mm    
      #}
      index.j1 <- Tree.l[[i]]$j1
      index.j2 <- Tree.l[[i]]$j2 
      #model.l <- get("cal.order",help.env)[[line]]
      model.l <- paircopula(U[,c(he)],K=K,lambda=lambda,pen=pen,base=base,m=m,q=q)
      list(j1=index.j1,j2=index.j2,D=NULL,v=get("ck.val",model.l),U=U[,c(index.j1,index.j2)],log.like = get("log.like",model.l),cop=model.l , AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l))
    }
  }
}

if(!order.Dvine) {
  if(!doMC) {
    Tree.l.temp <- foreach(i=1:length(Tree.l),.combine=list,.multicombine=TRUE) %do%   {
      index.j1 <- Tree.l[[i]]$j1
      index.j2 <- Tree.l[[i]]$j2
      model.l <- paircopula(U[,c(index.j1,index.j2)],K=K,lambda=lambda,pen=pen,base=base,m=m,q=q)
      list(j1=index.j1,j2=index.j2,D=NULL,v=get("ck.val",model.l),U=U[,c(index.j1,index.j2)],log.like = get("log.like",model.l),AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l))
    }
  }
 if(doMC) {
    Tree.l.temp <- foreach(i=1:length(Tree.l),.combine=list,.multicombine=TRUE,.options.multicore=mcoptions) %dopar%   {
      index.j1 <- Tree.l[[i]]$j1
      index.j2 <- Tree.l[[i]]$j2
  
      model.l <- paircopula(U[,c(index.j1,index.j2)],K=K,lambda=lambda,pen=pen,base=base,m=m,q=q)
      list(j1=index.j1,j2=index.j2,D=NULL,v=get("ck.val",model.l),U=U[,c(index.j1,index.j2)],log.like = get("log.like",model.l),AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l))
    }
  }
}
  if(base=="Bernstein") ddb <- K+1
  if(base=="B-spline" & q==2) ddb <- K+q-1
  if(base=="B-spline" & q==1) ddb <- K+q-2
 
  Index.basis.D <- matrix(NA,ddb^2,2)
  Index.basis.D[,1] <- rep(seq(1,ddb),ddb)
  Index.basis.D[,2] <- sort(Index.basis.D[,1])

  if(level<length(S)) {
    for(i in 1:length(Tree.l.temp)) {
      log.like <- log.like+Tree.l.temp[[i]]$log.like
      AIC <- AIC+Tree.l.temp[[i]]$AIC
      cAIC <- cAIC+Tree.l.temp[[i]]$cAIC
    }
  }
  else {
    log.like <- log.like+Tree.l.temp$log.like
    AIC <- AIC+Tree.l.temp$AIC
    cAIC <- cAIC+Tree.l.temp[[i]]$cAIC
    Tree.l.temp <- list(Tree.l.temp)
  }

Dvine[[level]] <- Tree.l.temp
if(doMC) {
  for(level in 3:length(S)) {
    Tree.l <- Dvine[[level]] # current tree in vine
    Tree.l1 <- Dvine[[level-1]] # previous tree in vine, one level up
    t.length <- length(Tree.l)
    Tree.l.temp <- foreach(i=1:t.length,.combine=list,.multicombine=TRUE,.options.multicore=mcoptions) %dopar%  {
      UU <-  c()
      index <- list(c(Tree.l[[i]]$j1,Tree.l[[i]]$D),c(Tree.l[[i]]$j2,Tree.l[[i]]$D)) # list of involved indices in current knot, i.e {j1,D} and {j2,D}
      for(j in 1:2) # Loop over both index sets, i.e. {j1,D} and {j2,D}
        {
          indexi <- sort(index[[j]])
          index.ancestor <- c()
          for (ml in 1:length(Tree.l1))
            {
              index.ancestor <- c(index.ancestor, all(indexi==sort(c(Tree.l1[[ml]]$j1,Tree.l1[[ml]]$j2,Tree.l1[[ml]]$D))))
            }
          ancestor.knot <- Tree.l1[index.ancestor][[1]] # Ancestor knot of j-th Element in Index set, i.e. knot with indices {j1,D} for j=1 and {j2,D} for j=2, respectively.
          if(j==1) UU <-cbind(UU,cond.cop(ancestor.knot$U,coef=ancestor.knot$v,K=K,diff="u2",Index.basis.D,base=base,q=q))
          if(j==2) UU <-cbind(UU,cond.cop(ancestor.knot$U,coef=ancestor.knot$v,K=K,diff="u1",Index.basis.D,base=base,q=q))
        }
      model.l <- paircopula(UU,K=K,lambda=lambda,pen=pen,base=base,m=m,q=q)
      list(j1=Tree.l[[i]]$j1,j2=Tree.l[[i]]$j2,D=Tree.l[[i]]$D,U=UU,v=get("ck.val",model.l),log.like = get("log.like",model.l),AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l))
    }
    
    if(level<length(S)) {
      for(i in 1:length(Tree.l.temp)) {
        log.like <- log.like+Tree.l.temp[[i]]$log.like
        AIC <- AIC+Tree.l.temp[[i]]$AIC
        cAIC <- cAIC+Tree.l.temp[[i]]$cAIC
      }
    }
    else {
      log.like <- log.like+Tree.l.temp$log.like
      AIC <- AIC+Tree.l.temp$AIC
      cAIC <- cAIC+Tree.l.temp$cAIC
      Tree.l.temp <- list(Tree.l.temp)
    }
    Dvine[[level]] <- Tree.l.temp
  }
}  

if(!doMC) {
  for(level in 3:length(S)) {
    Tree.l <- Dvine[[level]] # current tree in vine
    Tree.l1 <- Dvine[[level-1]] # previous tree in vine, one level up
    t.length <- length(Tree.l)
    Tree.l.temp <- list()
    for(i in 1:t.length) {
      UU <-  c()
      index <- list(c(Tree.l[[i]]$j1,Tree.l[[i]]$D),c(Tree.l[[i]]$j2,Tree.l[[i]]$D)) # list of involved indices in current knot, i.e {j1,D} and {j2,D}
      for(j in 1:2) # Loop over both index sets, i.e. {j1,D} and {j2,D}
        {
          indexi <- index[[j]]
          index.ancestor <- c()
          for (ml in 1:length(Tree.l1))
            {
              index.ancestor <- c(index.ancestor, all( sort(indexi)==sort(c(Tree.l1[[ml]]$j1,Tree.l1[[ml]]$j2,Tree.l1[[ml]]$D))))
            }
          ancestor.knot <- Tree.l1[index.ancestor][[1]] # Ancestor knot of j-th Element in Index set, i.e. knot with indices {j1,D}
                   # for j=1 and {j2,D} for j=2, respectively. 
          if(j==1) UU <-cbind(UU,cond.cop(ancestor.knot$U,coef=ancestor.knot$v,K=K,diff="u2",Index.basis.D,base=base,q=q))
          if(j==2) UU <-cbind(UU,cond.cop(ancestor.knot$U,coef=ancestor.knot$v,K=K,diff="u1",Index.basis.D,base=base,q=q))
        }
      model.l <- paircopula(UU,K=K,lambda=lambda,pen=pen,base=base,m=m,q=q)
      Tree.l.temp[[i]] <- list(j1=Tree.l[[i]]$j1,j2=Tree.l[[i]]$j2,D=Tree.l[[i]]$D,U=UU,v=get("ck.val",model.l),log.like = get("log.like",model.l),AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l))
    }
    if(level<length(S)) {
      for(i in 1:length(Tree.l.temp)) {
        log.like <- log.like+Tree.l.temp[[i]]$log.like
        AIC <- AIC+Tree.l.temp[[i]]$AIC
        cAIC <- cAIC+Tree.l.temp[[i]]$cAIC
      }
    }
    else {
      log.like <- log.like+Tree.l.temp[[i]]$log.like
      AIC <- AIC+Tree.l.temp[[i]]$AIC
      cAIC <- cAIC+Tree.l.temp[[i]]$cAIC
    }
    Dvine[[level]] <- Tree.l.temp
  }
}
  class(Dvine) <- "penDvine"
return(list(Dvine=Dvine,log.like=log.like,AIC=AIC,cAIC=cAIC,K=K,order=get("order",help.env),S=S,N=N,base=base))
}
