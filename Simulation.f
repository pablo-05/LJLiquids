      program MAIN
      implicit NONE
      
      ! Parameters and Variables
      real*8 L, ran, deltamax, pi, dist, rbin, rho, r_in, r_out
      real*8 r_centro, V_capa, g_r, vol_frac
      integer Np, Nmc, Nmov, i, j, k, nbin, Nbines, Nmeas, Ntarg, n, m
      logical sol
      
      ! Parameters setup
      parameter (L=14.d0, Np=256, Nmc=100000, deltamax=0.2d0)
      parameter (pi=4*datan(1.d0), rbin=0.1d0)
      
      ! Arrays
      real*8 r(Np, 3), rnew(3), bines(0:500)
      
      ! Calculate Volume Fraction (phi) for info
      vol_frac = (Np * (4.d0/3.d0) * pi * (1.d0**3)) / (L**3)
      print *, "Volume Fraction (phi): ", vol_frac
      
      ! 1. Initialization: Simple Cubic Lattice
      ! This ensures particles start without overlaps
      nbin = nint(Np**(1.d0/3.d0))
      k = 1
      do i=0, nbin-1
        do j=0, nbin-1
          do m=0, nbin-1
            if (k .le. Np) then
              r(k,1) = i * (L/nbin) + 0.5d0
              r(k,2) = j * (L/nbin) + 0.5d0
              r(k,3) = m * (L/nbin) + 0.5d0
              k = k + 1
            end if
          end do
        end do
      end do

      ! 2. Thermalization / Equilibration
      Nmov = 0
      do k=1, Nmc
        call random_number(ran)
        i = int(ran * Np) + 1
        
        do j=1,3
          call random_number(ran)
          rnew(j) = r(i,j) + deltamax * (2.d0 * ran - 1.d0)
          ! Periodic Boundary Conditions
          if (rnew(j) .ge. L) rnew(j) = rnew(j) - L
          if (rnew(j) .lt. 0.d0) rnew(j) = rnew(j) + L
        end do
        
        call solap(Np, i, r, rnew, L, sol)
        if (sol) then
          Nmov = Nmov + 1
          r(i,1:3) = rnew(1:3)
        end if
      end do
      print *, "Equilibration Acceptance: ", dble(Nmov)/Nmc

      ! 3. Production and g(r) Measurement
      bines = 0.d0
      Nmeas = 0
      Ntarg = 2000 ! Number of samples to take
      
      do k=1, Ntarg * 20 ! Sample every 20 MC steps
        call random_number(ran)
        i = int(ran * Np) + 1
        do j=1,3
          call random_number(ran)
          rnew(j) = r(i,j) + deltamax * (2.d0 * ran - 1.d0)
          if (rnew(j) .ge. L) rnew(j) = rnew(j) - L
          if (rnew(j) .lt. 0.d0) rnew(j) = rnew(j) + L
        end do
        
        call solap(Np, i, r, rnew, L, sol)
        if (sol) r(i,1:3) = rnew(1:3)

        ! Take Measurement
        if (mod(k, 20) .eq. 0) then
          Nmeas = Nmeas + 1
          do n = 1, Np - 1
            do m = n + 1, Np
              call distance(r(n,1:3), r(m,1:3), L, dist)
              if (dist .lt. L/2.d0) then
                i = int(dist / rbin)
                bines(i) = bines(i) + 2.d0
              end if
            end do
          end do
        end if
      end do

      ! 4. Normalization and Output
      rho = dble(Np) / (L**3)
      Nbines = int((L/2.d0) / rbin)
      open(unit=10, file='gofr.dat', status='replace')

      do i = 1, Nbines - 1
        r_in = dble(i) * rbin
        r_out = r_in + rbin
        r_centro = r_in + (rbin / 2.d0)
        
        ! Shell Volume: V = 4/3 * pi * (r_out^3 - r_in^3)
        V_capa = (4.d0/3.d0) * pi * (r_out**3 - r_in**3)
        
        ! g(r) = <n(r)> / (rho * V_shell)
        g_r = bines(i) / (dble(Nmeas) * dble(Np) * rho * V_capa)
        
        write(10, '(F10.4, F12.6)') r_centro, g_r
      end do
      close(10)
      print *, "Simulation finished. Results in gofr.dat"

      end program

      ! Subroutine: Check for Hard Sphere overlaps (radius = 1.0, sigma = 2.0)
      subroutine solap(Np, i, r, rnew, L, sol)
      implicit none
      integer i, Np, k
      real*8 r(Np, 3), rnew(3), L, dist
      logical sol

      sol = .true.
      do k=1, Np
        if (k .ne. i) then
          call distance(r(k,1:3), rnew, L, dist)
          if (dist .lt. 2.d0) then ! Overlap if distance < sigma
            sol = .false.
            return
          end if
        end if
      end do
      end subroutine

      ! Subroutine: Minimum Image Convention Distance
      subroutine distance(r1, r2, L, dist)
      implicit none
      real*8 r1(3), r2(3), L, d(3), dist
      integer k

      do k=1, 3
        d(k) = r1(k) - r2(k)
        ! Apply nearest image
        if (d(k) .gt. L/2.d0) d(k) = d(k) - L
        if (d(k) .lt. -L/2.d0) d(k) = d(k) + L
      end do
      dist = dsqrt(d(1)**2 + d(2)**2 + d(3)**2)
      end subroutine

