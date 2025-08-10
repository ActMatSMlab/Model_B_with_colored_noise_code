c  this program simulates phase ordering in the 2-d nematic model.
      character*8 name2(5000)
	parameter(k=512,k2=k**2,khalf=k/2,n=2)

        real x(k,k,2),af(n,k,k),xn(n,k,k,n),b(k,k),cur(2,k,k)

        
	integer m(k),p(k)
        real data(8388608),datc(8388608)
	
	integer nn(2)
	real s(-khalf:khalf-1,-khalf:khalf-1),scatt(100000,0:2000)
	real c(-khalf:khalf-1,-khalf:khalf-1),corr(100000,0:2000)

!        open(unit=95,file='eta1.inp')
!        read(95,*) eta1    
     
	dt=0.001
	dx=1.0
	dx2=2*dx
	dxx=dx**2
	par=dt/dx**2
        phia=0.11
        af0 = 0.0

cccccccccccccccc==============================
      
c       anisotropic diffusion coeff

c       non-eqm coupling coeff
         at= 1.0

        alam=7.0
c        alam= phia*65
        diff=1.0
        bf = sqrt(10.0)

        frho=1.0
        rave=-0.7

        ad=0.0

	nens=1
	frac=0.5
	drun=1.0

        

        ntime=60000
        ist = 1
	nspa= 100
	nscatt=5000


        npro=1

	iseed=-123495 

	kcut=khalf*frac
	runmax=sqrt(2.0)*kcut
	kmax=runmax
	nout=runmax/drun
	pi=2*asin(1.0)
	nn(1)=k
	nn(2)=k
        open(111, file='fnameapol')
        read(111,*)(name2(i),i=1,5000)
c        open(112, file='fnamerhol')
c        read(112,*)(name2(i),i=1,900)


c   here, we initialize the neighbour table

        do i=1,k
          p(i)=i+1
          if(i.eq.k)p(i)=1
          m(i)=i-1
          if(i.eq.1)m(i)=k
        enddo
	
c  here, we start the iteration process


	
    
c  here, we set the initial conditions for Q


        sub=0.0
        do j=1,k
        do i=1,k
        den=0.1*(0.5-ran2(iseed))
        x(i,j,1)=rave+den
        sub=sub+den
        enddo
        enddo

        sub=sub/k2
       sumden=0.0
        do j=1,k
        do i=1,k
        x(i,j,1)=x(i,j,1)-sub
        sumden=sumden+x(i,j,1)
        enddo
        enddo
        sumden=sumden/k2

        ! sp_temp noise initialization

        do ic =1,n
          do j=1,k
        do i=1,k

        xn(ic,i,j,1) = 0.2*(0.5-ran2(iseed))

        enddo
        enddo
        enddo

c  here, we start the iteration process

	iold=1
	inew=2


        do itime=1,ntime

         do ic=1,n
        do j=1,k
        do i=1,k

           af(ic,i,j)= bf*gasdev(iseed)

           enddo
           enddo
           enddo
 
        do j=1,k
        do i=1,k
         
        px=(x(p(i),j,iold)-2*x(i,j,iold)+x(m(i),j,iold))/dxx
        py=(x(i,p(j),iold)-2*x(i,j,iold)+x(i,m(j),iold))/dxx

        rho=(px+py)  

                    
        b(i,j)=-rho+(x(i,j,iold)**3-x(i,j,iold))
        
        axx=(xn(1,p(i),j,iold)+xn(1,m(i),j,iold)-2*xn(1,i,j,iold))/dxx
       ayy=(xn(1,i,p(j),iold)+xn(1,i,m(j),iold)-2*xn(1,i,j,iold))/dxx

       bxx=(xn(2,p(i),j,iold)+xn(2,m(i),j,iold)-2*xn(2,i,j,iold))/dxx
       byy=(xn(2,i,p(j),iold)+xn(2,i,m(j),iold)-2*xn(2,i,j,iold))/dxx

        xn(1,i,j,inew)= xn(1,i,j,iold)
     1   - (at**(-1))*(xn(1,i,j,iold)
     1   -(alam**2)*(axx+ayy))*dt
     1   +(at**(-1))*af(1,i,j)*sqrt(dt)

         xn(2,i,j,inew)= xn(2,i,j,iold)
     1   - (at**(-1))*(xn(2,i,j,iold)
     1   -(alam**2)*(bxx+byy))*dt
     1   +(at**(-1))*af(2,i,j)*sqrt(dt)
       
        enddo
        enddo
        
cc......calculating current
	do j=1,k
	do i=1,k
	cx=(b(p(i),j)-b(m(i),j))/dx2
	cy=(b(i,p(j))-b(i,m(j)))/dx2
	
	cur(1,i,j)=cx
	cur(2,i,j)=cy
	end do
	end do
	
        sums=0.0
        do j=1,k
        do i=1,k

	curx_x=(cur(1,p(i),j)-cur(1,m(i),j))/dx2
        cury_y=(cur(2,i,p(j))-cur(2,i,m(j)))/dx2
        
        afxx=(xn(1,p(i),j,inew)-xn(1,m(i),j,inew))/dx2
        afyy=(xn(2,i,p(j),inew)-xn(2,i,m(j),inew))/dx2
     

        rhoxx=(b(p(i),j)-2*b(i,j)+b(m(i),j))/dxx
        rhoyy=(b(i,p(j))-2*b(i,j)+b(i,m(j)))/dxx     
    
        add1=rhoxx+rhoyy

        x(i,j,inew)=x(i,j,iold)+dt*add1*diff
     1              +frho*(afxx+afyy)*sqrt(dt*2*diff)

  
        sums=sums+x(i,j,inew)     

	enddo
	enddo
                
        sums=sums/k2

        sumxmod=sums
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        



	if(itime > ist)then
	if (mod(itime,20).eq.0)then
       write(*,*)itime,sums,iens,alam,at,diff
c	isee=itime/nspa
c        write(*,*)itime,iens,sumden,sumxmod
       write(53,*)itime,sumxmod
c        if(iens.eq.1) then

       npro=npro+1
        open(300,file =name2(npro))

c        open(200, file =name1(npro))
        do j=1,k
        do i=1,k
        
        write(300,*)i,j,x(i,j,inew)
c        endif
        enddo
        write(300,*)""
        enddo
        endif
	endif

c	enddo
	
c	endif

c   here, we rewrite the new configuration into the old one

	itemp=iold
	iold=inew
	inew=itemp
        

c  here, the loop on the time ends

	enddo
	stop
	end

c***********************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      real data(8388608)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     *k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END

c***********************************************************
        FUNCTION ran2(idum)
       INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
       REAL ran2,AM,EPS,RNMX
       PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     * IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     * NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
       INTEGER idum2,j,k,iv(NTAB),iy
       SAVE iv,iy,idum2
       DATA idum2/123456789/, iv/NTAB*0/, iy/0/
       if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
           k=idum/IQ1
           idum=IA1*(idum-k*IQ1)-k*IR1
           if (idum.lt.0) idum=idum+IM1
           if (j.le.NTAB) iv(j)=idum
11       continue
         iy=iv(1)
       endif
       k=idum/IQ1
       idum=IA1*(idum-k*IQ1)-k*IR1
       if (idum.lt.0) idum=idum+IM1
       k=idum2/IQ2
       idum2=IA2*(idum2-k*IQ2)-k*IR2
       if (idum2.lt.0) idum2=idum2+IM2
       j=1+iy/NDIV
       iy=iv(j)-idum2
       iv(j)=idum
       if(iy.lt.1)iy=iy+IMM1
       ran2=min(AM*iy,RNMX)

       return
       END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          FUNCTION gasdev(idum)
        INTEGER idum
        REAL gasdev
C USES ran1
CReturns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
cas the source of uniform deviates.
        INTEGER iset
        REAL fac,gset,rsq,v1,v2,ran2
        SAVE iset,gset
        DATA iset/0/
        if (idum.lt.0) iset=0           !Reinitialize.
        if (iset.eq.0) then             !We don’t have an extra deviate handy, so
1        v1=2.*ran2(idum)-1.0 !pick two uniform numbers in the square extend
         v2=2.*ran2(idum)-1.0   !g from -1 to +1 in each direction,
        rsq=v1**2+v2**2                 !see if they are in the unit circle,
        if(rsq.ge.1..or.rsq.eq.0.)goto 1        !and if they are not, try again.
        fac=sqrt(-2.*log(rsq)/rsq) !Now make the Box-Muller transformation to get
      !  two normal deviates. Return one and save
!the other for next time.
        gset=v1*fac
        gasdev=v2*fac
        iset=1          !Set flag.
        else !We have an extra deviate handy,
        gasdev=gset !so return it,
        iset=0 !and unset the flag.
        endif
        return
        END
