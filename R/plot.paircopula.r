plot.paircopula <- function(x,val=NULL,marg=TRUE,plot=TRUE,int=FALSE,main.txt=NULL,
                    sub.txt=NULL,contour=FALSE,cuts=20,cex=1,cex.axes=1,
                    xlab=NULL,ylab=NULL,zlab=NULL,...) {
  if(!class(x)=="paircopula") stop("obj has to be of class paircopula")
  library(lattice)
  env <- list()
  p <- get("p",x)
  q <- 0
  d <- get("d",x)
  ddb <- get("ddb",x)
  Index.basis.D <- get("Index.basis.D",x)
  ck <- get("ck.val",x)
  D <- get("D",x)
  alpha <- 0
  base <- get("base",x)
  index.b <- matrix(0:get("dd",x))
  #if(int & base=="Bernstein") stop("Distribution function is not supported for Bernstein polynomials")

  if(is.null(val))
    {
      if(p!=2) stop("geht nicht!")
      else {
        x.grid <- seq(0,1,length=21)
        grid <- expand.grid(y1=x.grid, y2=x.grid)  
        tilde.Psi.d <-  array(NA, dim=c(dim(grid)[1],get("ddb",x),p))
        for (j in 1:p)
          {
            if(base=="Bernstein") {
              if(int) tilde.Psi.d[,,j] <-  int.bernstein(x,Y=grid[,j])
              else tilde.Psi.d[,,j] <-  apply(index.b,1,bernstein,x=grid[,j],n=get("dd",x))
            }
            if(base=="B-spline") {
              tilde.Psi.d[,,j] <-  my.bspline(y=grid[,j],K=get("K",x)+1,q=get("q",x))$base.den
            }  
          }
        tilde.PSI.d.D <- tilde.Psi.d[,Index.basis.D[,1],1]
        
        for (j in 2:p)
          {        
            tilde.PSI.d.D <- tilde.PSI.d.D * tilde.Psi.d[,Index.basis.D[,j],j]
          }
        grid[["plot"]] <- tilde.PSI.d.D%*%ck
        lam1 <- get("lambda",x)[1]               
        if(is.null(main.txt)) {
          main.txt <- substitute("K="*a*","*c*"="*d,list(a=D,c=parse(text="lambda")[[1]],d=lam1))
          main.txt <- as.expression(main.txt)
        }
          
        k <- dim(x$liste)[1]
        log.like <- round(get("log.like",x),3)
        pen.log.like <- round(get("pen.log.like",x),3)
        if(is.null(sub.txt)) sub.txt <- paste("log like=",log.like,", pen. log like= ",pen.log.like,", cAIC=",round(get("cAIC",x),3),sep="")
        if(int) z.txt <- "distribution" else z.txt <- "density"

        hh <- c("y1","y2")
        values <- as.formula(paste("plot~",paste(hh,collapse="*"),sep=""))
               
        if(!contour) obj1 <- wireframe(values,data=grid,outer=TRUE,sub=sub.txt,zlab=list(label=z.txt,cex=cex.axes),xlab=list(cex=cex.axes,label=xlab),
                                       ylab=list(cex=cex.axes,label=ylab),scales=list(arrows=FALSE,col="black",font=3,x=list(cex=cex),y=list(cex=cex),z=list(cex=cex)),
                                       main=main.txt,shade=TRUE,zlim=c(0,max(grid)), par.settings = list(axis.line = list(col = "transparent")),
                                       par.box = c(col = "black"))
        else obj1 <- contourplot(values, data=grid,outer=TRUE,sub=sub.txt,zlab=list(label=z.txt,cex=cex.axes),xlab=list(label=xlab,cex=cex.axes),
                                 ylab=list(label=ylab,cex=cex.axes),scales=list(arrows=FALSE,col="black",font=3,cex=cex),
                                 zlim=c(0,max(grid)),main=main.txt,shade=TRUE,cuts=cuts)
        
        if(marg) {
          T.marg <- get("T.marg",x)
          base <- get("tilde.Psi.knots.d",x)
          
          xx <- rep(seq(0,1,length=get("ddb",x)),p)
          density <- c()
          fac <- c()
          
          for(j in 1:p)
            {
              density <- c(density,round(base%*%(T.marg[,,j]%*%ck),5))
              fac <- c(fac,rep(j,get("ddb",x)))
            }
          datafr <- data.frame(xx,density,fac)
          graph.sets <-list(superpose.line=list(col=c(1:p),superpose.symbol = list(col = c(1:p))))
          
          obj2 <- xyplot(density~xx|fac,type="l",auto.key=list(space="right",title="marginal densities",sort=FALSE),par.settings=graph.sets)
          if(plot) {
            print(obj1,position=c(0,0.35,1,1),more=TRUE)
            print(obj2,position=c(0,0,1,0.35))
          }
          else return(list(density=obj1,marg.density=obj2))
        }
        else if(plot) {
          print(obj1)
          check <- any(grid$plot<0)
          print(check)
          if(check) print(grid[grid$plot<0,])
        }
        else return(list(density=obj1,grid=grid))
      }
    }
  else
    {
      if(!is.matrix(val)) {
        if(is.data.frame(val)) val <- as.matrix(val) else stop("val has to be a data.frame or a matrix")
      }
        tilde.Psi.d <-  array(NA, dim=c((length(val)/p),get("dd",x)+1,p))
        val <- matrix(val,(length(val)/p),p)

        for (j in 1:p)
          {
            if(int) tilde.Psi.d[,,j] <-  int.bernstein(x,Y=val[,j])
            else tilde.Psi.d[,,j] <- apply(index.b,1,bernstein,x=val[,j],n=get("dd",x))
          }
        tilde.PSI.d.D <- tilde.Psi.d[,Index.basis.D[,1],1]
   
        for (j in 2:p)
          {  
            tilde.PSI.d.D <- tilde.PSI.d.D * tilde.Psi.d[,Index.basis.D[,j],j]
          }
        datafr <- data.frame(val,tilde.PSI.d.D%*%ck)
        colnames(datafr)[p+1] <- "fit"
        return(datafr)
      }
}
