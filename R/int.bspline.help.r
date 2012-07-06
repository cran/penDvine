int.bspline.help <- function(penden.env) {
  len.k <- get("K",penden.env)
  x.help <- seq(0,1,length=1001)
  int.base <- distr.func.help(base,knots=seq(0,1,length=(get("K",penden.env)+1)),penden.env,q=get("q",penden.env),y=seq(0,1,length=1001))
  assign("int.base",int.base,penden.env)
  assign("x.help",x.help,penden.env)
}
