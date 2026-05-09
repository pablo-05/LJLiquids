      PROGRAM MAIN
      IMPLICIT NONE
      REAL step, rmax, alpha, pi, rho, error,limit,beta,PHI
      INTEGER i,N,contador,contmax
      PARAMETER (rmax=4.096, alpha=0.01, beta=1.0/1.15, rho=0.75)
      PARAMETER (pi=4*atan(1.0), step=0.001, N=int(rmax/step+0.99))
      PARAMETER (limit=0.01,contmax=10000)
      real fold(N), fnew(N), c(N), h(N), g(N), k(N)
      real dk

      ! Initialize momentum space frequency vector
      dk = pi / rmax
      do i=1, N
        k(i) = (i-1) * dk
      end do

      ! Initialize fold array to zero value
      fold = 0.0

      ! ============================================
      ! =        1. g(r) using PY Closure          =
      ! ============================================
      ! Main loop of the iterative method
      error = 1.0 ! This value is just to enter the loop
	  contador=0
      do while((error.gt.limit).and.(contador.lt.contmax)) ! Condition for convergence
        error = 0.0 ! We use error as a cumulative variable, so we reset to zero
		contador=contador+1
        do i=1, N
          ! Calculate r*h(r) — sinft needs r*f(r) for the 3D radial transform
          h(i)=(i-1)*step*((fold(i)+1)*exp(-beta*PHI((i-1)*step))-1.0)

          ! Calculate r*c(r) = r*h(r) - r*f(r)
          c(i) = h(i) - (i-1)*step*fold(i)
        end do

        ! We perform both Fourier transforms
        call sinft(h,N)
        call sinft(c,N)

        ! k-space OZ: k*f_hat(k) = (4pi)^2 * step^2 * sinft(rh)*sinft(rc) / k
        do i=2, N
          fnew(i) = rho * h(i) * c(i) * (4.0 * pi * step)**2 / k(i)
        end do  
		fnew(1) =fnew(2)
        ! Inverse sinft gives r*f(r); divide by normalization to recover gamma(r)
        call sinft(fnew,N)

        ! We have to normalize the inverse values
        do i=2, N
          fnew(i) = fnew(i) / (2.0 * pi * (i-1) * step * rmax)
        end do
		fnew(1) =fnew(2)
        ! Calculate the total divergence (before mixing to avoid artificial limit bypassing)
        do i=1, N
          error = error + abs(fnew(i)-fold(i))
        end do

        do i=1, N
          fnew(i) = alpha * fnew(i) + (1.0-alpha) * fold(i)
        end do

        ! Set the arrays for a new iteration
        do i=1, N
          fold(i)=fnew(i)
        end do
      end do

      ! If convergence was achieved, calculate g(r)
      do i=1, N
          g(i) = exp( -beta * PHI((i-1)*step) ) * ( fold(i) + 1 )
      end do

      ! Save the results to memory
      open(11, File="PY.dat")
      do i=1, N
        write(11, *) (i-1)*step, g(i)
      end do
      close(11)
	  
      ! ============================================
      ! =        2. g(r) using HNC Closure          =
      ! ============================================
      ! We are iterating over h(r) instead of f(r), but to save space
      ! we are going to use the already declared fold to store positions there

      ! Get the last h(r) PY value
      do i=1, N
        ! Calculate r*h(r) — sinft needs r*f(r) for the 3D radial transform
        h(i) = 0.0 ! i*step * ((fold(i)+1)*exp(-beta*PHI(i*step)) - 1.0)
      end do

      ! MAIN LOOP of the iterative method
      error = 1.0 ! This value is just to enter the loop
	  contador=0
      do while((error.gt.limit).and.(contador.lt.contmax)) ! Condition for convergence
        error = 0.0 ! We use error as a cumulative variable, so we reset to zero
		contador=contador+1
        ! We keep the rh(r) values at fold as they are going to be overwritten
        ! Set the arrays for a new iteration
        do i=1, N
          fold(i) = h(i)
        end do
        ! We perform rh(r) Fourier transform to obtain qh(q)
        call sinft(h,N)

        ! We get h(q) and then use it to calculate qc(q)
        do i=2, N
          h(i) = h(i) * 4 * pi * step / k(i)
          c(i) = k(i) * h(i) / ( 1 + rho * h(i) ) ! We calculate qc(q)
        end do
		h(1) = h(2)
		c(1) = c(2)
        call sinft(c,N) ! Inverse Fourier transform to obtain rc(r)

        ! We have to normalize the inverse values to obtain c(r)
        ! We also perform on the spot the closure relation calculation
        do i=2, N
          c(i) = c(i) / (2.0 * pi * (i-1) * step * rmax)
          h(i)=exp(-beta*PHI((i-1)*step)+(fold(i)/((i-1)*step))-c(i) )
          h(i) = (i-1)*step * ( h(i) - 1.0 )
          
          error = error + abs(h(i)-fold(i)) ! Divergence before mixing
          
          h(i) = alpha * h(i) + (1.0-alpha) * fold(i)
        end do
		h(1) = h(2)
		c(1) = c(2)
		
      end do

      ! Save the results to memory
      open(11, File="HNC.dat")
	  write(11, *) 0., 0.
      do i=2, N
        write(11, *) (i-1)*step, (h(i)/((i-1)*step)) + 1.0
      end do
      close(11)

      end program


      REAL FUNCTION PHI(r)
      IMPLICIT NONE
      REAL r, rc, vc
      
      rc = 2.5
      
      if (r .lt. 0.3) then
        PHI = 1.e30 ! Hard-core repulsion, effectively infinite
      else if (r .lt. rc) then
        vc = 4.0 * (1.0/(rc**12) - 1.0/(rc**6))
        PHI = 4.0 * (1.0/(r**12) - 1.0/(r**6)) - vc
      else
        PHI = 0.0
      end if
      RETURN

      END FUNCTION
	  
c ---------------------------------------------------------------------
      SUBROUTINE sinft(y,n)
      INTEGER n
      REAL y(n)
CU    USES realft
      INTEGER j
      REAL sum,y1,y2
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n)
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      y(1)=0.0
      do 11 j=1,n/2
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
        y1=wi*(y(j+1)+y(n-j+1))
        y2=0.5*(y(j+1)-y(n-j+1))
        y(j+1)=y1+y2
        y(n-j+1)=y1-y2
11    continue
      call realft(y,n,+1)
      sum=0.0
      y(1)=0.5*y(1)
      y(2)=0.0
      do 12 j=1,n-1,2
        sum=sum+y(j)
        y(j)=y(j+1)
        y(j+1)=sum
12    continue
      return
      END

      SUBROUTINE realft(data,n,isign)
      INTEGER isign,n
      REAL data(n)
CU    USES four1
      INTEGER i,i1,i2,i3,i4,n2p3
      REAL c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if (isign.eq.1) then
        c2=-0.5
        call four1(data,n/2,+1)
      else
        c2=0.5
        theta=-theta
      endif
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n/2,-1)
      endif
      return
      END

      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
	  