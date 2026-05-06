      PROGRAM MAIN
      IMPLICIT NONE
      REAL step, rmax, alpha, pi, rho, error,limit,beta,PHI
      INTEGER i,N
      PARAMETER (rmax=4.096, alpha=0.01, beta=1.0/0.95, rho=0.75)
      PARAMETER (pi=4*atan(1.0), step=0.001, N=int(rmax/step+0.99))
      PARAMETER (limit=0.01)
      real fold(N), fnew(N), c(N), h(N), g(N), k(N)
      real dk

      ! Initialize momentum space frequency vector
      dk = pi / rmax
      do i=1, N
        k(i) = i * dk
      end do

      ! Initialize fold array to zero value
      fold = 0.0

      ! ============================================
      ! =        1. g(r) using PY Closure          =
      ! ============================================
      ! Main loop of the iterative method
      error = 1.0 ! This value is just to enter the loop
      do while(error.gt.limit) ! Condition for convergence
        error = 0.0 ! We use error as a cumulative variable, so we reset to zero
        do i=1, N
          ! Calculate r*h(r) — sinft needs r*f(r) for the 3D radial transform
          h(i) = i*step * ((fold(i)+1.0)*exp(-beta*PHI(i*step)) - 1.0)

          ! Calculate r*c(r) = r*h(r) - r*f(r)
          c(i) = h(i) - i*step*fold(i)
        end do

        ! We perform both Fourier transforms
        call sinft(h,N)
        call sinft(c,N)

        ! k-space OZ: k*f_hat(k) = (4pi)^2 * step^2 * sinft(rh)*sinft(rc) / k
        do i=1, N
          fnew(i) = rho * h(i) * c(i) * (4.0 * pi * step)**2 / k(i)
        end do  

        ! Inverse sinft gives r*f(r); divide by normalization to recover gamma(r)
        call sinft(fnew,N)

        ! We have to normalize the inverse values
        do i=1, N
          fnew(i) = fnew(i) / (2.0 * pi * i * step * rmax)
        end do

        ! Calculate the total divergence (before mixing to avoid artificial limit bypassing)
        do i=1, N
          error = error + abs(fnew(i)-fold(i))
        end do

        do i=1, N
          fnew(i) = alpha * fnew(i) + ((1.0-alpha) * fold(i))
        end do

        ! Set the arrays for a new iteration
        do i=1, N
          fold(i)=fnew(i)
        end do
      end do

      ! If convergence was achieved, calculate g(r)
      do i=1, N
        g(i) = exp( -beta * PHI(i*step) ) * ( fold(i) + 1.0 )
      end do

      ! Save the results to memory
      open(11, File="PY.dat")
      do i=1, N
        write(11, *) i*step, g(i)
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
        h(i) = 0.0 ! i*step * ((fold(i)+1.0)*exp(-beta*PHI(i*step)) - 1.0)
      end do

      ! MAIN LOOP of the iterative method
      error = 1.0 ! This value is just to enter the loop
      do while(error.gt.limit) ! Condition for convergence
        error = 0.0 ! We use error as a cumulative variable, so we reset to zero

        ! We keep the rh(r) values at fold as they are going to be overwritten
        ! Set the arrays for a new iteration
        do i=1, N
          fold(i) = h(i)
        end do
        ! We perform rh(r) Fourier transform to obtain qh(q)
        call sinft(h,N)

        ! We get h(q) and then use it to calculate qc(q)
        do i=1, N
          h(i) = h(i) * 4 * pi * step / k(i)
          c(i) = k(i) * h(i) / ( 1 + rho * h(i) ) ! We calculate qc(q)
        end do

        call sinft(c,N) ! Inverse Fourier transform to obtain rc(r)

        ! We have to normalize the inverse values to obtain c(r)
        ! We also perform on the spot the closure relation calculation
        do i=1, N
          c(i) = c(i) / (2.0 * pi * i * step * rmax)
          h(i) = exp( -beta*PHI(i*step) + (fold(i)/(i*step)) - c(i) )
          h(i) = i*step * ( h(i) - 1.0 )
          
          error = error + abs(h(i)-fold(i)) ! Divergence before mixing
          
          h(i) = alpha * h(i) + (1.0-alpha) * fold(i)
        end do
      end do

      ! Save the results to memory
      open(11, File="HNC.dat")
      do i=1, N
        write(11, *) i*step, (h(i)/(i*step)) + 1.0
      end do
      close(11)

      end program


      REAL FUNCTION PHI(r)
      IMPLICIT NONE
      REAL r

      PHI = 4.0 * ( 1.0/(r**12) - 1.0/(r**6) )
      RETURN

      END FUNCTION