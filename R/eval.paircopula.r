eval.paircopula <- function(val,K,int=FALSE,Index.basis.D,ck.val) {
  if(!is.matrix(val)) {
    if(is.data.frame(val)) val <- as.matrix(val) else stop("val has to be a data.frame or a matrix")
  }
  tilde.Psi.d <-  array(NA, dim=c((length(val)/2),K+1,2))
  val <- matrix(val,(length(val)/2),2)
  penden.env <- new.env()
  assign("K",K,penden.env)
  int.bernstein.help(penden.env)
  index.b <- matrix(0:K)
  for (j in 1:2)
    {
      if(int) tilde.Psi.d[,,j] <-  int.bernstein(x,Y=val[,j])
      else tilde.Psi.d[,,j] <- apply(index.b,1,bernstein,x=val[,j],n=K)
    }
  tilde.PSI.d.D <- tilde.Psi.d[,Index.basis.D[,1],1]
  
  for (j in 2:2)
    {  
      tilde.PSI.d.D <- tilde.PSI.d.D * tilde.Psi.d[,Index.basis.D[,j],j]
    }
  return(tilde.PSI.d.D%*%ck.val)
}
