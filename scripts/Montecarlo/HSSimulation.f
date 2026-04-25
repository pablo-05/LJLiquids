      program MAIN
      implicit NONE
      
      ! Parameters and Variables
      REAL*8 L, ran, deltamax, pi, dist, rbin, rho, r_in, r_out
      REAL*8 r_middle, V_shell, g_r, vol_frac, radius
      INTEGER Np, STEPS, Nmov, n, initCubeSide, Nbins, Nmeas
      INTEGER i, j, k, t, m ! Indices for the program. m for particles,
      ! i, j, k for coordinates
      LOGICAL nOverlap
      
      ! Parameters setup
      PARAMETER (L=16.d0, Np=343, STEPS=10000, deltamax=0.2d0)
      PARAMETER (pi=4*datan(1.d0), rbin=0.01d0, radius=1.d0)
      PARAMETER (Nbins=int((L/2.d0)/rbin))
      
      ! Arrays
      REAL*8 r(Np, 3), rnew(3), bins(0:Nbins-1) ! Positions of every part.,
      ! new generated move for one and bins to store g(r)
      
      ! Calculate Fraction of occupied volume with cubeside L and Np s.
      vol_frac = (Np * (4.d0/3.d0) * pi * (radius**3)) / (L**3)
      print *, "Occupied Volume Fraction: ", vol_frac
      

      ! ============================================
      ! =           1. Initialization              =
      ! ============================================
      ! This ensures particles start without overlaps. We are
      ! initializing to a simple cubic lattice with equal separation
      initCubeSide = int(dble(Np)**(1.d0/3.d0) + 0.1d0) ! Max cube side needed
      m = 1
      do i=0, initCubeSide - 1
        do j=0, initCubeSide - 1
          do k=0, initCubeSide - 1
            if (m .le. Np) then
              r(m,1) = i * (L/initCubeSide) + radius
              r(m,2) = j * (L/initCubeSide) + radius
              r(m,3) = k * (L/initCubeSide) + radius
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
      bins = 0.d0
      Nmeas = 0
      t = 1
      do k=1, STEPS    ! We are going to perform STEPS MC steps
        !print*,k
        do j=1, Np     ! For each MC step, iterate over every particle
          call random_number(ran) ! We do a random particle selection
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
            r(m,1:3) = rnew(1:3)
            Nmov = Nmov + 1
          end if        
        end do
        ! ============================================
        ! =   3. g(r) Production and Measurement     =
        ! ============================================
        ! Take Measurement every 20 STEPS after the 0.4*STEPS
        if (t .ge. 20 .and. k .ge. 4000) then
          t = 0
          Nmeas = Nmeas + 1
          do n = 1, Np - 1
            do m = n + 1, Np
              CALL DISTANCE(r(n,1), r(n,2), r(n,3), r(m,1), r(m,2),
     +                      r(m,3), L, dist)
              if (dist .lt. L/2.d0) then
                i = int(dist / rbin)
                bins(i) = bins(i) + 2.d0
              end if
            end do
          end do
        end if
        t = t+1
      end do
      print *, "Acceptance rate: ", Nmov/dble(STEPS*Np)


      ! ============================================
      ! =       4. Normalization and Output        =
      ! ============================================
      rho = dble(Np) / (L**3)
      open(unit=10, file='gofr.dat', status='replace')

      do i = 0, Nbins - 1
        r_in = dble(i) * rbin
        r_out = r_in + rbin
        r_middle = r_in + (rbin / 2.d0)
        
        ! Shell Volume: V = 4/3 * pi * (r_out^3 - r_in^3)
        V_shell = (4.d0/3.d0) * pi * (r_out**3 - r_in**3)
        
        ! g(r) = <n(r)> / (rho * V_shell)
        g_r = bins(i) / (dble(Nmeas) * dble(Np) * rho * V_shell)
        
        write(10, '(F10.4, F12.6)') r_middle, g_r
      end do
      close(10)
      print *, "Simulation finished. Results in gofr.dat"

      END PROGRAM

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
          CALL DISTANCE(r(i,1), r(i,2), r(i,3), rnew(1), rnew(2),
     +                  rnew(3), L, dist)
          if (dist .lt. 2.d0*radius) then ! Overlap if distance < sigma
            nOverlap = .false.
            return
          end if
        end if
      end do
      END SUBROUTINE

      ! Subroutine: Minimum Image Convention Distance
      SUBROUTINE DISTANCE(x1, y1, z1, x2, y2, z2, L, dist)
      IMPLICIT NONE
      REAL*8 x1, y1, z1, x2, y2, z2, L, dist, dx, dy, dz

      dx = x1 - x2
      ! Apply nearest image
      if (dx .gt. L/2.d0) dx = dx - L
      if (dx .lt. -L/2.d0) dx = dx + L
      
      dy = y1 - y2
      if (dy .gt. L/2.d0) dy = dy - L
      if (dy .lt. -L/2.d0) dy = dy + L
      
      dz = z1 - z2
      if (dz .gt. L/2.d0) dz = dz - L
      if (dz .lt. -L/2.d0) dz = dz + L
      
      dist = dsqrt(dx**2 + dy**2 + dz**2)
      END SUBROUTINE

