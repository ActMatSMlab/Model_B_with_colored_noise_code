c  this program simulates modelB with color noise

        character*8 name2(5000)
	parameter(k=512,k2=k**2,khalf=k/2,n=2)

        real x(k,k,2),af(n,k,k),xn(n,k,k,n),b(k,k),cur(2,k,k)

        
	integer m(k),p(k)
	integer nn(2)
	   
	dt=0.005		! time step
	dx=1.0			! grid size
	dx2=2*dx
	dxx=dx**2		
	par=dt/dx**2		! stability check
       
        af0 = 0.0

cccccccccccccccc==============================
        at= 1.0		! tau in color noise

        alam=7.0		! lambda in color noise
        diff=1.0		! diffusivity
        bf = sqrt(10.0)	! strength of color noise

        frho=1.0		! noise strength
        rave=-0.7		! mean density of passive

        ad=0.0

	nens=1			! number of realization
	frac=0.5
	drun=1.0

        

        ntime=60000		! total simulation time
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
        open(111, file='fnameapol')		! storing outputs
        read(111,*)(name2(i),i=1,5000)



c   here, we initialize the neighbour table

        do i=1,k
        p(i)=i+1
        if(i.eq.k)p(i)=1
        m(i)=i-1
        if(i.eq.1)m(i)=k
        enddo
	    
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

        af(ic,i,j)= bf*gasdev(iseed)

        enddo
        enddo
        enddo
 
        do j=1,k
        do i=1,k
         
        ! Terms in modelB
         
        px=(x(p(i),j,iold)-2*x(i,j,iold)+x(m(i),j,iold))/dxx
        py=(x(i,p(j),iold)-2*x(i,j,iold)+x(i,m(j),iold))/dxx

        rho=(px+py)  

	! Chemical potential
	                    
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

! Update equation for phi

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
        write(53,*)itime,sumxmod
	
	! storing output
	
        npro=npro+1
        open(300,file =name2(npro))

        do j=1,k
        do i=1,k
        
        write(300,*)i,j,x(i,j,inew)

        enddo
        write(300,*)""
        enddo
        endif
	endif


c   here, we rewrite the new configuration into the old one

	itemp=iold
	iold=inew
	inew=itemp
        

c  here, the loop on the time ends

	enddo
	stop
	end

c***********************************************************
! subroutine for random number generation

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
