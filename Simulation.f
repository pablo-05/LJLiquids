      program MAIN
      implicit NONE
      
      ! Parameters and Variables
      REAL*8 L, ran, deltamax, pi, dist, rbin, rho, r_in, r_out
      REAL*8 r_centro, V_capa, g_r, vol_frac, radius
      INTEGER Np, STEPS, Nmov, Ntarg, n, initCubeSide, Nbins, Nmeas
      INTEGER i, j, k, m ! Indices for the program. m for particles,
      ! i, j, k for coordinates
      LOGICAL nOverlap
      
      ! Parameters setup
      PARAMETER (L=14.d0, Np=256, STEPS=1000, deltamax=0.2d0)
      PARAMETER (pi=4*datan(1.d0), rbin=0.1d0, radius=1.d0)
      
      ! Arrays
      REAL*8 r(Np, 3), rnew(3), bins(0:500) ! Positions of every part.,
      ! new generated move for one and bins to store g(r)
      
      ! Calculate Fraction of occupied volume with cubeside L and Np s.
      vol_frac = (Np * (4.d0/3.d0) * pi * (radius**3)) / (L**3)
      print *, "Occupied Volume Fraction: ", vol_frac
      

      ! ============================================
      ! =           1. Initialization              =
      ! ============================================
      ! This ensures particles start without overlaps. We are
      ! initializing to a simple cubic lattice with equal separation
      initCubeSide = int(Np**(1.d0/3.d0)) ! Max cube side needed
      m = 1
      do i=0, initCubeSide
        do j=0, initCubeSide
          do k=0, initCubeSide
            if (m .le. Np) then
              r(m,1) = i * (L/initCubeSide) + 0.5d0
              r(m,2) = j * (L/initCubeSide) + 0.5d0
              r(m,3) = k * (L/initCubeSide) + 0.5d0
              m = m + 1
              ! If we have initialized all particles, exit loops
              if (m.eq.Np+1) goto 999
            end if
          end do
        end do
      end do
999   continue


      ! ============================================
      ! =    2. Thermalization / Equilibration     =
      ! ============================================
      Nmov = 0 ! Number of accepted, total performed moves
      do k=1, STEPS    ! We are going to perform STEPS MC steps
        do j=1, Np     ! For each MC step, iterate over every particle
998       call random_number(ran) ! We do a random particle selection
          m = int(ran * Np) + 1 ! int() produces from 0 to Np-1, thus +1
          
          ! We generate a random movement in each direction from 0 to
          ! deltamax for the selected particle.
          do i=1,3 ! Three coordinates
            call random_number(ran) ! Generate random variable
            ! Store the new pos in rnew to later evaluate it
            ! We are making the random be in (-1, 1) so it is in each dir.
            rnew(i) = r(m,i) + deltamax * (2.d0 * ran - 1.d0)

            ! Apply Periodic Boundary Conditions
            if (rnew(i) .ge. L) rnew(i) = rnew(i) - L
            if (rnew(i) .lt. 0.d0) rnew(i) = rnew(i) + L
          end do
          
          ! Subroutine to check for overlapping
          call overlapping(Np, m, r, rnew, L, nOverlap, radius)
          if (nOverlap) then
            r(i,1:3) = rnew(1:3)
          else
            Nmov = Nmov + 1
            goto 998
          end if
        end do
      end do
      print *, "Acceptance rate: ",STEPS*Np/dble(Nmov+STEPS*Np)



      ! ============================================
      ! =   3. g(r) Production and Measurement     =
      ! ============================================
      bins = 0.d0
      Nmeas = 0
      Ntarg = 2000 ! Number of samples to take
      
      do k=1, Ntarg * 20 ! Sample every 20 MC steps
        CALL random_number(ran)
        i = int(ran * Np) + 1
        do j=1,3
          CALL random_number(ran)
          rnew(j) = r(i,j) + deltamax * (2.d0 * ran - 1.d0)
          if (rnew(j) .ge. L) rnew(j) = rnew(j) - L
          if (rnew(j) .lt. 0.d0) rnew(j) = rnew(j) + L
        end do
        
        CALL OVERLAPPING(Np, i, r, rnew, L, nOverlap, radius)
        if (nOverlap) r(i,1:3) = rnew(1:3)

        ! Take Measurement
        if (mod(k, 20) .eq. 0) then
          Nmeas = Nmeas + 1
          do n = 1, Np - 1
            do m = n + 1, Np
              CALL DISTANCE(r(n,1:3), r(m,1:3), L, dist)
              if (dist .lt. L/2.d0) then
                i = int(dist / rbin)
                bins(i) = bins(i) + 2.d0
              end if
            end do
          end do
        end if
      end do

      ! 4. Normalization and Output
      rho = dble(Np) / (L**3)
      Nbins = int((L/2.d0) / rbin)
      open(unit=10, file='gofr.dat', status='replace')

      do i = 1, Nbins - 1
        r_in = dble(i) * rbin
        r_out = r_in + rbin
        r_centro = r_in + (rbin / 2.d0)
        
        ! Shell Volume: V = 4/3 * pi * (r_out^3 - r_in^3)
        V_capa = (4.d0/3.d0) * pi * (r_out**3 - r_in**3)
        
        ! g(r) = <n(r)> / (rho * V_shell)
        g_r = bins(i) / (dble(Nmeas) * dble(Np) * rho * V_capa)
        
        write(10, '(F10.4, F12.6)') r_centro, g_r
      end do
      close(10)
      print *, "Simulation finished. Results in gofr.dat"

      end program



      ! Subroutine: Check for Hard Sphere overlaps
      SUBROUTINE OVERLAPPING(Np, m, r, rnew, L, nOverlap, radius)
      IMPLICIT NONE
      INTEGER m, Np, i
      REAL*8 r(Np, 3), rnew(3), L, dist, radius
      LOGICAL nOverlap

      ! We start by saying there is no overlap
      nOverlap = .true.
      do i=1, Np ! We check against every other particle
        if (i .ne. m) then ! except ourselves
          ! We calculate the distance between the new pos. and the part.
          CALL DISTANCE(r(i,1:3), rnew, L, dist)
          if (dist .lt. 2.d0*radius) then ! Overlap if distance < sigma
            nOverlap = .false.
            return
          end if
        end if
      end do
      end subroutine

      ! Subroutine: Minimum Image Convention Distance
      SUBROUTINE DISTANCE(r1, r2, L, dist)
      IMPLICIT NONE
      REAL*8 r1(3), r2(3), L, d(3), dist
      INTEGER k

      do k=1, 3
        d(k) = r1(k) - r2(k)
        ! Apply nearest image
        if (d(k) .gt. L/2.d0) d(k) = d(k) - L
        if (d(k) .lt. -L/2.d0) d(k) = d(k) + L
      end do
      dist = dsqrt(d(1)**2 + d(2)**2 + d(3)**2)
      end subroutine

