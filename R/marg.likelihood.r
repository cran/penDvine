marg.likelihood <- function(penden.env,pen.likelihood,temp=FALSE) {
  help <- eigen(get("DDD.sum",penden.env),symmetric=TRUE)
  index <- which(help$values>1e-8)
  evalues <- help$values[index]
  Utilde <- help$vectors[,index]
  k1 <- 0.5*sum(log(get("lambda",penden.env)*evalues))
  k2 <- pen.likelihood
  if(!temp) eneu <- eigen(t(Utilde)%*%-get("Derv2.pen",penden.env)%*%Utilde)$values
  else eneu <- eigen(t(Utilde)%*%-get("Derv2.pen.temp",penden.env)%*%Utilde)$values
  if(any(eneu<=0)) print("eneu <=0")
  k3 <- -0.5*sum(log(eneu))
  if(!temp) assign("marg.log.like",k1+k2+k3,penden.env)
  else assign("marg.log.like.temp",k1+k2+k3,penden.env)
}
