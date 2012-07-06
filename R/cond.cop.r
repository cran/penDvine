cond.cop <- function(data,coef,K,diff="u2",Index.basis.D,base,q=NULL) {
  p <- 2
  penden.env <- new.env()
  assign("K",K,penden.env)
  int.bernstein.help(penden.env)
  ddb <- K+1
  assign("u1",data[,1],penden.env)
  assign("u2",data[,2],penden.env)
  
  tilde.Psi.d <-  array(NA, dim=c(dim(data)[1],ddb,p))
  index.b <- matrix(0:K)

  for (j in 1:p)
    {
      obj <- paste("u",j,sep="")
      if(base=="Bernstein") {
        if(obj==diff) tilde.Psi.d[,,j] <-  apply(index.b,1,bernstein,get(paste("u",j,sep=""),penden.env), n=K)
        if(obj!=diff) tilde.Psi.d[,,j] <-  int.bernstein(penden.env,Y=get(paste("u",j,sep=""),penden.env))
      }
      if(base=="B-spline") {
        if(obj==diff) tilde.Psi.d[,,j] <-  my.bspline(y=get(paste("u",j,sep=""),penden.env),K=K+1,q=q)$base.den
        if(obj!=diff) tilde.Psi.d[,,j] <-  int.bspline2(penden.env,Y=get(paste("u",j,sep=""),penden.env))
      }
    }
  assign("tilde.Psi.d",tilde.Psi.d,penden.env)
  assign("tilde.PSI.d.D",tilde.Psi.d[,Index.basis.D[,1],1],penden.env)
  
  for (j in 2:p)
    {
      assign("tilde.PSI.d.D",get("tilde.PSI.d.D",penden.env) * get("tilde.Psi.d",penden.env)[,Index.basis.D[,j],j],penden.env)
    }

  return(get("tilde.PSI.d.D",penden.env)%*%coef)
}
