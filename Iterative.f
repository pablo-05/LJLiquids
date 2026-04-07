      program MAIN
      implicit NONE
      REAL*8 step, rmax, alpha, pi
      INTEGER i,j,k,N
      PARAMETER (rmax=8.d0, alpha=0.1d0)
      PARAMETER (pi=4.d0*datan(1.d0), N=2**11, step=rmax/N)
      real fo(N), fn(N), c(N), h(N)
      print *, "the maximum valid k_r (in spherical coordinates) is "
     &, pi/step
      do i=1, N
        fo(i)=exp(-100 * i*step)
      end do
      call sinft(fo,N)
      open(unit=11, file='fourier.dat', status='replace')
      do i=1, N
        write(11, *) i*pi/rmax, fo(i)*step
      end do
      close(11)
      end program