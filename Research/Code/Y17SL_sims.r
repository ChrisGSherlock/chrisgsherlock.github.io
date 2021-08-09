### Simulations for DAVB paper
## The code that actually produces the plots is inside an
## if (F) { }
## So that individual plots can be produced one at a time
## ( a couple of them take a few minutes to run )


norm<-function(x) {
    return(sqrt(sum(x*x)))
}

log.density.MVN<-function(x,mu,Sigma) {
    a=(x-mu)
    b=t(a)%*%solve(Sigma)%*%a
    return(-0.5*b[1,1])
}

lp.power<-function(x,a) {
    return(-norm(x)^a/a)
}
glp.power<-function(x,a) {
    return(-x*norm(x)^(a-2))
}

##||X||^a /a ~ Gam(d/a,1)
r.power<-function(n,d,a) {
    trans=rgamma(n,d/a)
    return((a*trans)^(1/a))
}

quantile.power<-function(q,d,a,...) {
    trans=qgamma(q,d/a,...)
    return((a*trans)^(1/a))    
}

## RWM ##

## lambda = scaling
## V = variance matrix (up to scaling)
## ... = parameters for lpfn

darwm.fixed.V.t2c<-function(nits,x0,lambda,V,thin=1,print=0,condfn,pi.pow,pihat.pow,pihat.sc) {
  nx<-length(x0)
  X.out<-matrix(nrow=floor(nits/thin)+1,ncol=nx+1)
  n.accept<-0
  x.curr<-x0
  lp.curr<-lp.power(x.curr,pi.pow)
  lphat.curr<-lp.power(x.curr/pihat.sc,pihat.pow)
  X.out[1,]<-c(x.curr,lp.curr)
  k<-2

  Vsqrt<-t(chol(V))
  d<-dim(V)
  
  for (i in 1:nits) {
      x.prop<-x.curr+lambda*Vsqrt%*% rnorm(nx)
      lp.prop<-lp.power(x.prop,pi.pow)
      lphat.prop<-lp.power(x.prop/pihat.sc,pihat.pow)
      R<-min(0,lphat.prop-lphat.curr)
      R<-R+min(0,lp.prop+lphat.curr-lp.curr-lphat.prop)
      if (log(runif(1))<R) {
          x.curr<-x.prop
          lp.curr<-lp.prop; lphat.curr=lphat.prop
          n.accept<-n.accept+1
      }
      if (i/thin == floor(i/thin)) {
          X.out[k,]<-c(x.curr,lp.curr)
          k<-k+1
      }
      if (print>0) {
          if (i/print == floor(i/print)) {
              print(paste("Iteration=",i,". alpha=",n.accept/i,sep=""))
              print(c(x.curr,lp.curr))
          }
      }
      if (condfn(x.curr)) {
          break;
      }
  }

  if (condfn(x0)) {
      condit=0
  }
  else if (condfn(x.curr)) {
      condit=i
  }
  else {
      condit=-1
  }
  print(paste("Condit=",condit," Acceptance rate=",n.accept/i,sep=""))
  return(list(acc=n.accept/i,condit=condit,X=X.out))
}

## IS ##

## ...= parameters for lp
## rqfn - simulate from q
## lpqfn - log density of q

daind.samp.t2c<-function(nits,x0,thin=1,print=0,condfn,q.pow,q.sc,pi.pow,pihat.pow,pihat.sc) {
  nx<-length(x0)
  X.out<-matrix(nrow=floor(nits/thin)+1,ncol=nx+1)
  n.accept<-0
  x.curr<-x0

  lp.curr<-lp.power(x.curr,pi.pow)
  lq.curr<-lp.power(x.curr/q.sc,q.pow)
  lphat.curr<-lp.power(x.curr/pihat.sc,pihat.pow)
  
  X.out[1,]<-c(x.curr,lp.curr)
  k<-2

  Vsqrt<-t(chol(V))
  
  for (i in 1:nits) {
      x.prop<-r.power(1,nx,q.pow)*q.sc
      lp.prop<-lp.power(x.prop,pi.pow)
      lq.prop<-lp.power(x.prop/q.sc,q.pow)
      lphat.prop<-lp.power(x.prop/pihat.sc,pihat.pow)

      R<-min(0,lphat.prop-lphat.curr + lq.curr-lq.prop)
      R<- R+min(0,lp.prop+lphat.curr-lp.curr-lphat.prop)
      if (log(runif(1))<R) {
          x.curr<-x.prop
          lp.curr<-lp.prop
          lq.curr<-lq.prop
          lphat.curr<-lphat.prop
          n.accept<-n.accept+1
      }
      if (i/thin == floor(i/thin)) {
          X.out[k,]<-c(x.curr,lp.curr)
          k<-k+1
      }
      if (print>0) {
          if (i/print == floor(i/print)) {
              print(paste("Iteration=",i,". alpha=",n.accept/i,sep=""))
              print(c(x.curr,lp.curr))
          }
      }
      if (condfn(x.curr)) {
          break;
      }
  }
  
  if (condfn(x0)) {
      condit=0
  }
  else if (condfn(x.curr)) {
      condit=i
  }
  else {
      condit=-1
  }
  print(paste("Condit=",condit," Acceptance rate=",n.accept/i,sep=""))
  return(list(acc=n.accept/i,condit=condit,X=X.out))
}


## MALA ##

# lambda = scaling
# obs = data structure containing obs
# V = variance matrix (up to scaling)
## mag>0 is for truncated MALA
damala.fixed.V.t2c<-function(nits,x0,lambda,V,mag=-1,thin=1,print=0,condfn,pi.pow,pihat.pow,pihat.sc) {
  nx<-length(x0)
  X.out<-matrix(nrow=floor(nits/thin)+1,ncol=nx+1)
  n.accept<-0; n.chop<-0
  x.curr<-x0
  lp.curr<-lp.power(x.curr,pi.pow)
  lphat.curr<-lp.power(x.curr/pihat.sc,pihat.pow)
  glp.curr<-glp.power(x.curr,pi.pow)
  if ((mag>0) && (norm(glp.curr)>mag)) {
      glp.curr=glp.curr*mag/norm(glp.curr)
  }
  X.out[1,]<-c(x.curr,lp.curr)
  k<-2

  Vsqrt<-t(chol(V))
  d<-dim(V)
  
  for (i in 1:nits) {
      a=lambda^2/2*V %*% glp.curr
      x.prop<-x.curr+a+lambda*Vsqrt%*% rnorm(nx)
      lp.prop<-lp.power(x.prop,pi.pow)
      lphat.prop<-lp.power(x.prop/pihat.sc,pihat.pow)
      glp.prop<-glp.power(x.prop,pi.pow)
      if ((mag>0) && (norm(glp.prop)>mag)) {
          glp.prop=glp.prop*mag/norm(glp.prop); n.chop=n.chop+1
      }

      q.curr2prop<-log.density.MVN(x.prop,x.curr+a,lambda^2*V)
      q.prop2curr<-log.density.MVN(x.curr,x.prop+lambda^2/2*glp.prop,lambda^2*V)

      R1<-min(0,lphat.prop-lphat.curr + q.prop2curr-q.curr2prop)
      R<- R1+min(0,lp.prop+lphat.curr-lp.curr-lphat.prop)
      if (log(runif(1))<R) {
          x.curr<-x.prop
          lp.curr<-lp.prop
          glp.curr<-glp.prop
          lphat.curr<-lphat.prop
          n.accept<-n.accept+1
      }
      if (i/thin == floor(i/thin)) {
          X.out[k,]<-c(x.curr,lp.curr)
          k<-k+1
      }
      if (print>0) {
          if (i/print == floor(i/print)) {
              print(paste("Iteration=",i,". alpha=",n.accept/i,sep=""))
              print(c(x.curr,lp.curr))
          }
      }
      if (condfn(x.curr)) {
          break;
      }
  }
  
  if (condfn(x0)) {
      condit=0
  }
  else if (condfn(x.curr)) {
      condit=i
  }
  else {
      condit=-1
  }
  print(paste("Condit=",condit," Acceptance rate=",n.accept/i,", Chop rate=",n.chop/i,sep=""))
  return(list(acc=n.accept/i,condit=condit,X=X.out))
}

## MALA using gradient of the approximation ##

# lambda = scaling
# obs = data structure containing obs
# V = variance matrix (up to scaling)
## mag>0 is for truncated MALA
damalahat.fixed.V.t2c<-function(nits,x0,lambda,V,mag=-1,thin=1,print=0,condfn,pi.pow,pihat.pow,pihat.sc) {
  nx<-length(x0)
  X.out<-matrix(nrow=floor(nits/thin)+1,ncol=nx+1)
  n.accept<-0
  x.curr<-x0
  lp.curr<-lp.power(x.curr,pi.pow)
  lphat.curr<-lp.power(x.curr/pihat.sc,pihat.pow)
  glp.curr<-glp.power(x.curr/pihat.sc,pihat.pow)/pihat.sc
##  print("Initial grad pre-shrink=")
##  print(glp.curr)
  if ((mag>0) && (norm(glp.curr)>mag)) {
      glp.curr=glp.curr*mag/norm(glp.curr)
  }
##  print("Initial=")
##  print(x.curr)
##  print("norm=")
##  print(norm(glp.curr))
##  print("grad=")
##  print(glp.curr)
  X.out[1,]<-c(x.curr,lp.curr)
  k<-2

  Vsqrt<-t(chol(V))
  d<-dim(V)
  
  for (i in 1:nits) {
      a=lambda^2/2*V %*% glp.curr
      x.prop<-x.curr+a+lambda*Vsqrt%*% rnorm(nx)
      lp.prop<-lp.power(x.prop,pi.pow)
      lphat.prop<-lp.power(x.prop/pihat.sc,pihat.pow)
      glp.prop<-glp.power(x.prop/pihat.sc,pihat.pow)/pihat.sc
##            print("grad Prop-pre-squish=")
##print(glp.prop)
      if ((mag>0) && (norm(glp.prop)>mag)) {
          glp.prop=glp.prop*mag/norm(glp.prop)
      }
##      print("Prop=")
##      print(x.prop)
##      print("norm=")
##      print(norm(glp.prop))
##      print("grad=")
##      print(glp.prop)
##      print("mag=")
##      print(mag)

      q.curr2prop<-log.density.MVN(x.prop,x.curr+a,lambda^2*V)
      q.prop2curr<-log.density.MVN(x.curr,x.prop+lambda^2/2*glp.prop,lambda^2*V)

      R1<-min(0,lphat.prop-lphat.curr + q.prop2curr-q.curr2prop)
      R<- R1+min(0,lp.prop+lphat.curr-lp.curr-lphat.prop)
##      print(c(R1,R-R1))
      if (log(runif(1))<R) {
          x.curr<-x.prop
          lp.curr<-lp.prop
          glp.curr<-glp.prop
          lphat.curr<-lphat.prop
          n.accept<-n.accept+1
      }
      if (i/thin == floor(i/thin)) {
          X.out[k,]<-c(x.curr,lp.curr)
          k<-k+1
      }
      if (print>0) {
          if (i/print == floor(i/print)) {
              print(paste("Iteration=",i,". alpha=",n.accept/i,sep=""))
              print(c(x.curr,lp.curr))
          }
      }
      if (condfn(x.curr)) {
          break;
      }
  }
  
  if (condfn(x0)) {
      condit=0
  }
  else if (condfn(x.curr)) {
      condit=i
  }
  else {
      condit=-1
  }
  print(paste("Condit=",condit," Acceptance rate=",n.accept/i,sep=""))
  return(list(acc=n.accept/i,condit=condit,X=X.out))
}


## Assumes d=5 and a=2, median
Conda2d5<-function(x) {
    return(norm(x)<2.086)
}
Conda1p5d5<-function(x) {
    return(norm(x)<2.7297)
}
Conda1d5<-function(x) {
    return(norm(x)<4.671)
}

if (F) {
    d=5
    V=diag(rep(1,d))
    xlab=expression("-log"[10]~"p"[0])    
    nrep=20
    qstarts=10^(-(1:6))

    pchrwmGd=3; pchrwmBd=4; colrwmGd=4; colrwmBd=2 ## blue +, red x
    pchisGd=1; pchisBd=2; colisGd=5; colisBd=6 ## cyan o, magenta triangle
    pchtmGd=1; pchtmBd=2; coltmGd=5; coltmBd=6 ## cyan o, magenta triangle
    pchmGd=20; pchmBd=17; colmGd=3; colmBd=1 ## green solid o, black striangle
    
    ## Figure 1  - ~ Examples 1 and 2
    rstartsa2=quantile.power(qstarts,5,2,lower.tail=FALSE)
    nstart=length(qstarts)
    Fig1scBd=matrix(nrow=nstart,ncol=nrep)
    Fig1scGd=matrix(nrow=nstart,ncol=nrep)
    pi.pow=2
    pihat.pow=2
    maxits=10000
    set.seed(1234567)
    
    for (i in 1:nstart) {
        for (j in 1:nrep) {
            s0=rnorm(d)
            s0=s0/norm(s0)
            x0=s0*rstartsa2[i]
            
            a=darwm.fixed.V.t2c(maxits,x0,2.4/sqrt(d),V,thin=1,print=0,Conda2d5,pi.pow,pihat.pow,2.0)
            Fig1scGd[i,j]=a$condit
            
            a=darwm.fixed.V.t2c(maxits,x0,2.4/sqrt(d),V,thin=1,print=0,Conda2d5,pi.pow,pihat.pow,0.5) 
            Fig1scBd[i,j]=a$condit

        }        
    }
    Fig1scBd[Fig1scBd==-1]=maxits+1
    Fig1scBd=log10(Fig1scBd); Fig1scGd=log10(Fig1scGd)
    ymx=max(max(Fig1scGd),max(Fig1scBd))
    ymn=min(min(Fig1scGd),min(Fig1scBd))

    xax=-log10(qstarts)
    xlim=c(min(xax)-.5,max(xax)+.5)
    plot(jitter(xax),Fig1scBd[,1],xlim=xlim,ylim=c(ymn,ymx),col=colrwmBd,xlab=xlab,ylab=expression("log"[10]~"Iterations until ||x|| < median"),pch=pchrwmBd,main=expression("log "*pi*"(x) = -||x||"^2*"/2,   log "*hat(pi)*"(x) = -||x||"^2*"/(2"*lambda^2*" )"))
    points(jitter(xax),Fig1scGd[,1],col=colrwmGd,pch=pchrwmGd)
    
    for (i in 2:nrep) {
        points(jitter(xax),Fig1scBd[,i],col=colrwmBd,pch=pchrwmBd)
        points(jitter(xax),Fig1scGd[,i],col=colrwmGd,pch=pchrwmGd)
    }
    abline(h=log10(maxits),lwd=2,lty=2)
   dev.print(pdf,file="Fig1.pdf")
 
    ## Figure 2 - examples 3 and 4

    set.seed(1234567)
    maxits=100000
    qstarts=10^(-(1:6))
    rstartsa1=quantile.power(qstarts,5,1,lower.tail=FALSE)
    nstart=length(rstartsa1)
    Fig2RWMscBd=matrix(nrow=nstart,ncol=nrep)
    Fig2RWMscGd=matrix(nrow=nstart,ncol=nrep)
    Fig2ISscBd=matrix(nrow=nstart,ncol=nrep)
    Fig2ISscGd=matrix(nrow=nstart,ncol=nrep)
    pi.pow=1
    pihat.pow=1
    q.pow=1
    
    for (i in 1:nstart) {
        for (j in 1:nrep) {
            s0=rnorm(d)
            s0=s0/norm(s0)
            x0=s0*rstartsa1[i]
            
            a=darwm.fixed.V.t2c(maxits,x0,2.4/sqrt(d),V,thin=1,print=0,Conda1d5,pi.pow,pihat.pow,2.0)
            Fig2RWMscGd[i,j]=a$condit
            
            a=darwm.fixed.V.t2c(maxits,x0,2.4/sqrt(d),V,thin=1,print=0,Conda1d5,pi.pow,pihat.pow,0.5) 
            Fig2RWMscBd[i,j]=a$condit

            a=daind.samp.t2c(maxits,x0,thin=1,print=0,Conda1d5,q.pow,1.0,pi.pow,pihat.pow,2.0)
            Fig2ISscGd[i,j]=a$condit
            a=daind.samp.t2c(maxits,x0,thin=1,print=0,Conda1d5,q.pow,1.0,pi.pow,pihat.pow,0.5)
            Fig2ISscBd[i,j]=a$condit

        }        
    }
    Fig2RWMscBd[Fig2RWMscBd==-1]=maxits+1
    Fig2ISscBd[Fig2ISscBd==-1]=maxits+1

    
    Fig2RWMscGd=log10(Fig2RWMscGd); Fig2RWMscBd=log10(Fig2RWMscBd)
    Fig2ISscGd=log10(Fig2ISscGd); Fig2ISscBd=log10(Fig2ISscBd)
    
    ymx=max(max(Fig2RWMscGd),max(Fig2RWMscBd),max(Fig2ISscGd),max(Fig2ISscBd))
    ymn=min(min(Fig2RWMscGd),min(Fig2RWMscBd),min(Fig2ISscGd),min(Fig2ISscBd))

    xax=-log10(qstarts)
    xlim=c(min(xax)-.5,max(xax)+.5)
    plot(jitter(xax),Fig2RWMscBd[,1],ylim=c(ymn,ymx),xlim=xlim,col=colrwmBd,xlab=xlab,ylab=expression("log"[10]~"Iterations until ||x|| < median"),pch=pchrwmBd,main=expression("log "*pi*"(x) = -||x|| ,   log "*hat(pi)*"(x) = -||x||/"*lambda))
    points(jitter(xax),Fig2RWMscGd[,1],col=colrwmGd,pch=pchrwmGd)
    points(jitter(xax),Fig2ISscBd[,1],col=colisBd,pch=pchisBd)
    points(jitter(xax),Fig2ISscGd[,1],col=colisGd,pch=pchisGd)
    
    for (i in 2:nrep) {
        points(jitter(xax),Fig2RWMscBd[,i],col=colrwmBd,pch=pchrwmBd)
        points(jitter(xax),Fig2RWMscGd[,i],col=colrwmGd,pch=pchrwmGd)
        points(jitter(xax),Fig2ISscBd[,i],col=colisBd,pch=pchisBd)
        points(jitter(xax),Fig2ISscGd[,i],col=colisGd,pch=pchisGd)
    }
    abline(h=log10(maxits),lwd=2,lty=2)
    dev.print(pdf,file="Fig2.pdf")

    



    ## Figure 3 - examples 5 and 6

    pi.pow=1.5
    pihat.pow.Gd=1.2
    pihat.pow.Bd=1.8
    pihat.sc=1
    maxits=10000
    nrep=10
    
    qstarts=10^(-(1:6)*4)
    rstartsa1p5=quantile.power(qstarts,5,pi.pow,lower.tail=FALSE)
##    rstartsa1p5[7]=25
    nstart=length(rstartsa1p5)
    Fig3RWMpBd=matrix(nrow=nstart,ncol=nrep)
    Fig3RWMpGd=matrix(nrow=nstart,ncol=nrep)
    Fig3TMpBd=matrix(nrow=nstart,ncol=nrep)
    Fig3TMpGd=matrix(nrow=nstart,ncol=nrep)
    Fig3MpBd=matrix(nrow=nstart,ncol=nrep)
    Fig3MpGd=matrix(nrow=nstart,ncol=nrep)
    set.seed(1234567)
    
    for (i in 1:nstart) {
        print(c("quantile=",i))
        for (j in 1:nrep) {
            s0=rnorm(d)
            s0=s0/norm(s0)
            x0=s0*rstartsa1p5[i]
            
            a=darwm.fixed.V.t2c(maxits,x0,1.4,V,thin=1,print=0,Conda1p5d5,pi.pow,pihat.pow.Gd,pihat.sc)
            Fig3RWMpGd[i,j]=a$condit
            
            a=darwm.fixed.V.t2c(maxits,x0,1.4,V,thin=1,print=0,Conda1p5d5,pi.pow,pihat.pow.Bd,pihat.sc) 
            Fig3RWMpBd[i,j]=a$condit

            a=damala.fixed.V.t2c(maxits,x0,1.6,V,2.4,thin=1,print=0,Conda1p5d5,pi.pow,pihat.pow.Gd,pihat.sc)
            Fig3TMpGd[i,j]=a$condit
            a=damala.fixed.V.t2c(maxits,x0,1.6,V,2.4,thin=1,print=0,Conda1p5d5,pi.pow,pihat.pow.Bd,pihat.sc)
            Fig3TMpBd[i,j]=a$condit

            a=damala.fixed.V.t2c(maxits,x0,1.6,V,-1,thin=1,print=0,Conda1p5d5,pi.pow,pihat.pow.Gd,pihat.sc)
            Fig3MpGd[i,j]=a$condit
            a=damala.fixed.V.t2c(maxits,x0,1.6,V,-1,thin=1,print=0,Conda1p5d5,pi.pow,pihat.pow.Bd,pihat.sc)
            Fig3MpBd[i,j]=a$condit
        }        
    }
    Fig3RWMpBd[Fig3RWMpBd==-1]=maxits+1; Fig3RWMpGd[Fig3RWMpGd==-1]=maxits+1
    Fig3TMpBd[Fig3TMpBd==-1]=maxits+1; Fig3TMpGd[Fig3TMpGd==-1]=maxits+1;
    Fig3MpBd[Fig3MpBd==-1]=maxits+1; Fig3MpGd[Fig3MpGd==-1]=maxits+1;
    
    Fig3RWMpGd=log10(Fig3RWMpGd); Fig3RWMpBd=log10(Fig3RWMpBd)
    Fig3TMpGd=log10(Fig3TMpGd); Fig3TMpBd=log10(Fig3TMpBd)
    Fig3MpGd=log10(Fig3MpGd); Fig3MpBd=log10(Fig3MpBd)

    ymx=max(max(Fig3RWMpGd),max(Fig3RWMpBd),max(Fig3TMpGd),max(Fig3TMpBd))
    ymx=max(ymx,max(Fig3MpGd),max(Fig3MpBd))
    ymn=min(min(Fig3RWMpGd),min(Fig3RWMpBd),min(Fig3TMpGd),min(Fig3TMpBd))
    ymn=min(ymn,min(Fig3MpGd),min(Fig3MpBd))
    
    xax=-log10(qstarts)
    xlim=c(min(xax)-.5,max(xax)+.5)
    plot(jitter(xax),Fig3RWMpBd[,1],ylim=c(ymn,ymx),col=colrwmBd,xlab=xlab,ylab=expression("log"[10]~"Iterations until ||x|| < median"),pch=pchrwmBd,xlim=xlim,main=expression("log "*pi*"(x) = -||x||"^beta*" ,   log "*hat(pi)*"(x) = -||x||"^gamma))
    points(jitter(xax),Fig3RWMpGd[,1],col=colrwmGd,pch=pchrwmGd)
    points(jitter(xax),Fig3TMpBd[,1],col=coltmBd,pch=pchtmBd)
    points(jitter(xax),Fig3TMpGd[,1],col=coltmGd,pch=pchtmGd)
    points(jitter(xax),Fig3MpBd[,1],col=colmBd,pch=pchmBd)
    points(jitter(xax),Fig3MpGd[,1],col=colmGd,pch=pchmGd)
    
    for (i in 2:nrep) {
        points(jitter(xax),Fig3RWMpBd[,i],col=colrwmBd,pch=pchrwmBd)
        points(jitter(xax),Fig3RWMpGd[,i],col=colrwmGd,pch=pchrwmGd)
        points(jitter(xax),Fig3TMpBd[,i],col=coltmBd,pch=pchtmBd)
        points(jitter(xax),Fig3TMpGd[,i],col=coltmGd,pch=pchtmGd)
        points(jitter(xax),Fig3MpBd[,i],col=colmBd,pch=pchmBd)
        points(jitter(xax),Fig3MpGd[,i],col=colmGd,pch=pchmGd)
    }
    abline(h=log10(maxits),lwd=2,lty=2)
    dev.print(pdf,file="Fig3.pdf")




    ## Figure 4 - examples 7 and 8

    pi.pow=1.5
    pihat.pow.Gd=1.2
    pihat.pow.Bd=1.8
    pihat.sc=1
    nrep=20

    qstarts=10^(-(1:6)*4)
    rstartsa1p5=quantile.power(qstarts,5,pi.pow,lower.tail=FALSE)
##    rstartsa1p5[7]=25
    nstart=length(rstartsa1p5)
    Fig4TMhpBd=matrix(nrow=nstart,ncol=nrep)
    Fig4TMhpGd=matrix(nrow=nstart,ncol=nrep)
    Fig4MhpGd=matrix(nrow=nstart,ncol=nrep)
    Fig4MhpBd=matrix(nrow=nstart,ncol=nrep)
    maxits=10000
    set.seed(1234567)
    
    for (i in 1:nstart) {
        print(i)
        for (j in 1:nrep) {
            s0=rnorm(d)
            s0=s0/norm(s0)
            x0=s0*rstartsa1p5[i]
            
            a=damalahat.fixed.V.t2c(maxits,x0,1.6,V,2.5,thin=1,print=0,Conda1p5d5,pi.pow,pihat.pow.Gd,pihat.sc)
            Fig4TMhpGd[i,j]=a$condit
            a=damalahat.fixed.V.t2c(maxits,x0,1.6,V,2.5,thin=1,print=0,Conda1p5d5,pi.pow,pihat.pow.Bd,pihat.sc)
            Fig4TMhpBd[i,j]=a$condit

            a=damalahat.fixed.V.t2c(maxits,x0,1.6,V,-1,thin=1,print=0,Conda1p5d5,pi.pow,pihat.pow.Gd,pihat.sc)
            Fig4MhpGd[i,j]=a$condit
            a=damalahat.fixed.V.t2c(maxits,x0,1.6,V,-1,thin=1,print=0,Conda1p5d5,pi.pow,pihat.pow.Bd,pihat.sc)
            Fig4MhpBd[i,j]=a$condit
        }        
    }
    Fig4TMhpBd[Fig4TMhpBd==-1]=maxits+1; Fig4TMhpGd[Fig4TMhpGd==-1]=maxits+1;
    Fig4MhpGd[Fig4MhpGd==-1]=maxits+1; Fig4MhpBd[Fig4MhpBd==-1]=maxits+1
    
    Fig4TMhpGd=log10(Fig4TMhpGd); Fig4TMhpBd=log10(Fig4TMhpBd)
    Fig4MhpGd=log10(Fig4MhpGd); Fig4MhpBd=log10(Fig4MhpBd)

    ymx=max(max(Fig4TMhpGd),max(Fig4TMhpBd))
    ymx=max(ymx,max(Fig4MhpGd),max(Fig4MhpBd))
    ymn=min(min(Fig4TMhpGd),min(Fig4TMhpBd))
    ymn=min(ymn,min(Fig4MhpGd),min(Fig4MhpBd))

    xax=-log10(qstarts)
    xlim=c(min(xax)-.5,max(xax)+.5)
    plot(jitter(xax),Fig4TMhpBd[,1],xlim=xlim,ylim=c(ymn,ymx),col=coltmBd,xlab=xlab,ylab=expression("log"[10]~"Iterations until ||x|| < median"),pch=pchtmBd,,main=expression("log "*pi*"(x) = -||x||"^beta*" ,   log "*hat(pi)*"(x) = -||x||"^gamma*", q uses "*nabla*" log "*hat(pi)))
    points(jitter(xax),Fig4TMhpGd[,1],col=coltmGd,pch=pchtmGd)
    points(jitter(xax),Fig4MhpBd[,1],col=colmBd,pch=pchmBd)
    points(jitter(xax),Fig4MhpGd[,1],col=colmGd,pch=pchmGd)
    
    for (i in 2:nrep) {
        points(jitter(xax),Fig4TMhpBd[,i],col=coltmBd,pch=pchtmBd)
        points(jitter(xax),Fig4TMhpGd[,i],col=coltmGd,pch=pchtmGd)
        points(jitter(xax),Fig4MhpBd[,i],col=colmBd,pch=pchmBd)
        points(jitter(xax),Fig4MhpGd[,i],col=colmGd,pch=pchmGd)
    }
    abline(h=log10(maxits),lwd=2,lty=2)
    dev.print(pdf,file="Fig4.pdf")

}

