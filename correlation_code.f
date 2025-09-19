c  this program simulates modelb with color noise

	parameter(k=2048,k2=k**2,khalf=k/2,n=2)

        real x(k,k,2),af(n,k,k),xn(n,k,k,n),b(k,k)

        
	integer m(k),p(k)
        real data(8388608),datc(8388608)
	
	integer nn(2)
	real s(-khalf:khalf-1,-khalf:khalf-1),scatt(100000,0:2000)
	real c(-khalf:khalf-1,-khalf:khalf-1),corr(100000,0:2000)  
     
	dt=0.005	!time step
	dx=1.0		!grid size
	dx2=2*dx
	dxx=dx**2
	par=dt/dx**2	!stability check
 
cccccccccccccccc==============================
        at= 1.0 	         !tau in color noise
        alam=7.0	         ! lambda in color noise
        diff=1.0	         ! diffusivity
        bf = sqrt(10.0)	 ! strength of color noise

        frho=1.0		 ! noise strength
        rave=-0.8	         ! mean density of passive		


	nens=25		 ! number of realization
	frac=0.5
	drun=1.0

    ntime=8000		! simulation time
	nspa=10
	nscatt=400

	iseed=-123497 

	kcut=khalf*frac
	runmax=sqrt(2.0)*kcut
	kmax=runmax
	nout=runmax/drun
	pi=2*asin(1.0)
	nn(1)=k
	nn(2)=k

	open(unit=1,file='lens.dat',status='unknown')		! file for length from structure factor
	open(unit=2,file='lenc.dat',status='unknown')		! file for length from correlatio function

        open(unit=11,file='sf.1',status='unknown')		! files for storing structure factor
        open(unit=12,file='sf.2',status='unknown')
        open(unit=13,file='sf.3',status='unknown')
        open(unit=14,file='sf.4',status='unknown')
        open(unit=15,file='sf.5',status='unknown')
        open(unit=16,file='sf.6',status='unknown')
        open(unit=17,file='sf.7',status='unknown')
        open(unit=18,file='sf.8',status='unknown')
        open(unit=19,file='sf.9',status='unknown')
        open(unit=20,file='sf.10',status='unknown')
        open(unit=21,file='sf.11',status='unknown')
        open(unit=22,file='sf.12',status='unknown')
        open(unit=23,file='sf.13',status='unknown')
        open(unit=24,file='sf.14',status='unknown')
        open(unit=25,file='sf.15',status='unknown')
        open(unit=26,file='sf.16',status='unknown')
        open(unit=27,file='sf.17',status='unknown')
        open(unit=28,file='sf.18',status='unknown')
        open(unit=29,file='sf.19',status='unknown')
        open(unit=30,file='sf.20',status='unknown')
        ns=10

        open(unit=31,file='cf.1',status='unknown')		! files for storing correlation function
        open(unit=32,file='cf.2',status='unknown')
        open(unit=33,file='cf.3',status='unknown')
        open(unit=34,file='cf.4',status='unknown')
        open(unit=35,file='cf.5',status='unknown')
        open(unit=36,file='cf.6',status='unknown')
        open(unit=37,file='cf.7',status='unknown')
        open(unit=38,file='cf.8',status='unknown')
        open(unit=39,file='cf.9',status='unknown')
        open(unit=40,file='cf.10',status='unknown')
        open(unit=41,file='cf.11',status='unknown')
        open(unit=42,file='cf.12',status='unknown')
        open(unit=43,file='cf.13',status='unknown')
        open(unit=44,file='cf.14',status='unknown')
        open(unit=45,file='cf.15',status='unknown')
        open(unit=46,file='cf.16',status='unknown')
        open(unit=47,file='cf.17',status='unknown')
        open(unit=48,file='cf.18',status='unknown')
        open(unit=49,file='cf.19',status='unknown')
        open(unit=50,file='cf.20',status='unknown')
        nc=30

c   here, we initialize the structure factor and
c   the correlation function to zero

	do i=1,100000
	do j=0,2000
	
	scatt(i,j)=0.0
	corr(i,j)=0.0
	
	enddo
	enddo

c   here, we initialize the neighbour table

        do i=1,k
        p(i)=i+1
        if(i.eq.k)p(i)=1
        m(i)=i-1
        if(i.eq.1)m(i)=k
        enddo
	
c  here, we start the iteration process


	do iens=1,nens
    
c  here, we set the initial conditions for phi


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

        af(ic,i,j)= bf*gasdev(iseed)		!Gaussian white noise

        enddo
        enddo
        enddo
 
        do j=1,k
        do i=1,k
        
        ! terms in ModelB 
        
        px=(x(p(i),j,iold)-2*x(i,j,iold)+x(m(i),j,iold))/dxx
        py=(x(i,p(j),iold)-2*x(i,j,iold)+x(i,m(j),iold))/dxx

        rho=(px+py)  

        ! chemical potential
                    
        b(i,j)=-rho+(x(i,j,iold)**3-x(i,j,iold))
        
        axx=(xn(1,p(i),j,iold)+xn(1,m(i),j,iold)-2*xn(1,i,j,iold))/dxx
        ayy=(xn(1,i,p(j),iold)+xn(1,i,m(j),iold)-2*xn(1,i,j,iold))/dxx

        bxx=(xn(2,p(i),j,iold)+xn(2,m(i),j,iold)-2*xn(2,i,j,iold))/dxx
        byy=(xn(2,i,p(j),iold)+xn(2,i,m(j),iold)-2*xn(2,i,j,iold))/dxx
        
	! Update equation for color noise
	
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

        sums=0.0
        do j=1,k
        do i=1,k

        afxx=(xn(1,p(i),j,inew)-xn(1,m(i),j,inew))/dx2
        afyy=(xn(2,i,p(j),inew)-xn(2,i,m(j),inew))/dx2
     

        rhoxx=(b(p(i),j)-2*b(i,j)+b(m(i),j))/dxx
        rhoyy=(b(i,p(j))-2*b(i,j)+b(i,m(j)))/dxx     
    
        add1=rhoxx+rhoyy
	
	! Update equation for phi (Model B + color noise)
	
        x(i,j,inew)=x(i,j,iold)+dt*add1*diff
     1              +frho*(afxx+afyy)*sqrt(dt*2*diff)

  
        sums=sums+x(i,j,inew)     

	enddo
	enddo
                
        sums=sums/k2

        sumxmod=sums
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        


c   here, we calculate the contribution to the scattering function
c   and the correlation function

        if((itime/nspa)*nspa.eq.itime)then

	isee=itime/nspa

c        write(53,*)itime,iens,sumden,sumxmod

c  here, we calculate the statistical properties of the fields
	
	do j=-khalf,khalf-1
	do i=-khalf,khalf-1
	s(i,j)=0.0
	c(i,j)=0.0
	enddo
	enddo

	do i=1,2*k2
	datc(i)=0.0
	enddo
	
c  structure factor for phi

	sum=0.0
	l=1
	do j=1,k
	do i=1,k

	adel=x(i,j,inew)-sumxmod
	data(l)=sign(1.0,adel)		! hardening data
	data(l+1)=0.0
	l=l+2
	enddo
	enddo

	call fourn(data,nn,2,1)
	l=1
	do j=0,khalf-1
	do i=0,khalf-1
	sc=data(l)**2+data(l+1)**2
	s(i,j)=s(i,j)+sc
	datc(l)=datc(l)+sc
	sum=sum+sc
	l=l+2
	enddo
	do i=-khalf,-1
	sc=data(l)**2+data(l+1)**2
	s(i,j)=s(i,j)+sc
	datc(l)=datc(l)+sc
	sum=sum+sc
	l=l+2
	enddo
	enddo
	do j=-khalf,-1
	do i=0,khalf-1
	sc=data(l)**2+data(l+1)**2
	s(i,j)=s(i,j)+sc
	datc(l)=datc(l)+sc
	sum=sum+sc
	l=l+2
	enddo
	do i=-khalf,-1
	sc=data(l)**2+data(l+1)**2
	s(i,j)=s(i,j)+sc
	datc(l)=datc(l)+sc
	sum=sum+sc
	l=l+2
	enddo
	enddo
	sum=sum/k2

c  correlation function for phi

	call fourn(datc,nn,2,-1)

	l=1
	do j=0,khalf-1
	do i=0,khalf-1
	c(i,j)=datc(l)
	l=l+2
	enddo
	do i=-khalf,-1
	c(i,j)=datc(l)
	l=l+2
	enddo
	enddo
	do j=-khalf,-1
	do i=0,khalf-1
	c(i,j)=datc(l)
	l=l+2
	enddo
	do i=-khalf,-1
	c(i,j)=datc(l)
	l=l+2
	enddo
	enddo


c  here, we spherically average quantities for output

	run=0.0
	runup=run+drun
	do iout=0,nout
	if(iout.eq.0)then
	
	scatt(isee,iout)=scatt(isee,iout)+s(0,0)/sum/nens
	corr(isee,iout)=corr(isee,iout)+c(0,0)/nens
    
	
	go to 1
	else
	icont=0
	conts=0.0
	contc=0.0
	do i=0,kmax
	do j=0,kmax
	check=sqrt(i**2+j**2+0.0)
	if((check.gt.run).and.(check.le.runup))then
	icont=icont+4
	
	conts=conts+s(i,j)+s(-i,j)+s(i,-j)+s(-i,-j)
	contc=contc+c(i,j)+c(-i,j)+c(i,-j)+c(-i,-j)


	endif
	enddo
	enddo
	
	conts=conts/icont
	contc=contc/icont


	scatt(isee,iout)=scatt(isee,iout)+conts/sum/nens
	corr(isee,iout)=corr(isee,iout)+contc/nens
       

	endif	
1	run=run+drun
	runup=run+drun
	enddo
	
	endif

c   here, we rewrite the new configuration into the old one

	itemp=iold
	iold=inew
	inew=itemp
        

c  here, the loop on the time ends

	enddo

c  here, the loop on the ensemble ends

	enddo
	
c  here, we output the length scales and scaled quantities
	

	nsee=ntime/nspa
	nskip=nscatt/nspa

	do isee=1,nsee
	
c  length from structure factor for Q
	
	r1=0.0
	r2=0.0
	do i=0,nout
	if(i.eq.0)then
	run=0.0
	else
	run=(i+0.5)*drun*(2*pi/k/dx)
	endif
	r1=r1+run*scatt(isee,i)
	r2=r2+scatt(isee,i)
	enddo
	r=r1/r2
	
	write(1,*)isee*nspa*dt,1/r
	
	
c  length from correlation function for Q

	do i=0,nout
	ord1=corr(isee,i)/corr(isee,0)
	ord2=corr(isee,i+1)/corr(isee,0)

	if((ord1.ge.0.5).and.(ord2.le.0.5))then
	
	xlow=i*dx
	xup=(i+1)*dx
	top=0.5*dx+(xlow*ord2-xup*ord1)
	x0=top/(ord2-ord1)
	write(2,*)isee*nspa*dt,x0
	go to 2
	
	endif
	
	enddo
	
2       dum=0.0

	
c  scaled structure factors and correlation functions for Q

	if((isee/nskip)*nskip.eq.isee)then
	ns=ns+1
	nc=nc+1
	
	do i=0,nout
	if(i.eq.0)then
	run=0.0
	dist=0.0
	else
	run=(i+0.5)*drun*(2*pi/k/dx)
	dist=(i+0.5)*drun*dx
	endif
        write(ns,*)run,scatt(isee,i)
        write(nc,*)dist,corr(isee,i)/corr(isee,0)

	enddo
	endif

c  scaled structure factors and correlation functions for rho


        enddo

	stop
	end

c***********************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine for FFT

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
! Subroutine for random number genereation

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
! Subroutine for random number generation

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
