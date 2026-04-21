      program MAIN
      implicit NONE
      REAL step, rmax, alpha, pi, rho, error,limit,betae
      INTEGER i,N
      PARAMETER (rmax=4, alpha=0.1, betae=1, rho=0.490)
      PARAMETER (pi=4*atan(real(1)), N=2**9, step=rmax/N)
      PARAMETER (limit=0.01)
      real fo(N), fn(N), c(N), h(N),g(N)
      print *, "the maximum valid k_r (in spherical coordinates) is "
     &, pi/step
      do i=1, N
        fo(i)=0
      end do
      error=1
      do while(error.gt.limit)
        error=0
        do i=1, N
          !calculamos rh(r) en el espacio de momentos
          h(i)=step*i*((fo(i)+1)*
     &exp(-4*betae*((1/(i*step)**12)-(1/(i*step)**6)))-1)
          !calculamos rc(r) en el espacio de momentos
          c(i)=h(i)-step*i*fo(i)
        end do
        call sinft(h,N)
        call sinft(c,N)
        do i= 1, N
        !calculamos qfnew(q)=(4pi)**2 (dr)**2 sinft(rh)sinft(rc)/q con q=pi i/rmax
          fn(i) = rho * h(i) * c(i) * 16 * pi *step**2 / i * rmax
        end do
        call sinft(fn,N)
        do i=1, N
        !obtenemos f(r)=2/N * dq/(r *2pi**2) sinft(qfnew(q)) con dq=pi/rmax y r = i*step
        fn(i) = fn(i) / (N * pi * i * step * rmax)
        end do
        fn= alpha * fn + (1-alpha) * fo
        do i =1, N
          error = error + abs(fn(i)-fo(i))
        end do
        do i=1, N
          fo(i)=fn(i)
        end do
      end do
      do i=1, N
          g(i)=(fn(i)+1)*
     &exp(-4*betae*((1/(i*step)**12)-(1/(i*step)**6)))
      end do
      open(unit=11, file='theogofr.dat', status='replace')
      do i=1, N
        write(11, *) i*step, g(i)
      end do
      close(11)
      end program