      program MAIN
      implicit NONE
      
      ! Parameters and Variables
      REAL*8 L, ran, deltamax, pi, dist, rbin, rho, r_in, r_out
      REAL*8 r_middle, V_shell, g_r
      INTEGER Np, STEPS, Nmov, n, initCubeSide, Nbins, Nmeas
      INTEGER i, j, k, t, m ! Indices for the program. m for particles,
      ! i, j, k for coordinates
      LOGICAL nAccepted
      
      ! Parameters setup
      ! For density rho = 0.75, L = (Np/0.75)^(1/3). With Np=343, L ~ 7.7046d0
      PARAMETER (L=7.7046d0, Np=343, STEPS=50000, deltamax=0.2d0)
      PARAMETER (pi=4*datan(1.d0), rbin=0.1d0)
      PARAMETER (Nbins=int(4.1d0/rbin))
      
      ! Arrays
      REAL*8 r(Np, 3), rnew(3), bins(0:Nbins-1) ! Positions of every part.,
      ! new generated move for one and bins to store g(r)
      
      ! In Lennard-Jones reduced units (sigma=1), density is just N/V
      rho = dble(Np) / (L**3)
      print *, "Reduced Density (rho*): ", rho
      

      ! ============================================
      ! =           1. Initialization              =
      ! ============================================
      ! This ensures particles start without overlaps. We are
      ! initializing to a simple cubic lattice with equal separation
      ! Using +0.999d0 ensures 343^(1/3) doesn't truncate to 6 due to FP math
      initCubeSide = int(dble(Np)**(1.d0/3.d0) + 0.999d0) ! Max cube side needed
      m = 1
      do i=0, initCubeSide - 1
        do j=0, initCubeSide - 1
          do k=0, initCubeSide - 1
            if (m .le. Np) then
              r(m,1) = i * (L/initCubeSide)
              r(m,2) = j * (L/initCubeSide)
              r(m,3) = k * (L/initCubeSide)
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
        call PROGRESS_BAR(k, STEPS)
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
          call ENERGY(Np, m, r, rnew, L, nAccepted)
          if (nAccepted) then
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
              if (dist .lt. 4.1d0) then
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
      SUBROUTINE ENERGY(Np, m, r, rnew, L, nAccepted)
      IMPLICIT NONE
      INTEGER m, Np, i
      REAL*8 r(Np, 3), rnew(3), L, dist, newE, oldE, POTENTIAL,
     +      RAN, beta
      LOGICAL nAccepted

      ! We start by saying we accept it
      nAccepted = .true.
      newE = 0.d0
      oldE = 0.d0
      
      ! We are at T=0.95, so beta=1/T
      beta = 1.d0 / 1.15d0
      
      do i=1, Np ! We check against every other particle
        if (i .ne. m) then ! except ourselves
          ! We calculate the distance between the new pos. and the part.
          CALL DISTANCE(r(i,1), r(i,2), r(i,3), rnew(1), rnew(2),
     +                  rnew(3), L, dist)
          newE = newE + POTENTIAL(dist)
          CALL DISTANCE(r(i,1), r(i,2), r(i,3), r(m,1), r(m,2),
     +                  r(m,3), L, dist)
          oldE = oldE + POTENTIAL(dist)
        end if
      end do

      ! We accept the move if the energy is lower, or with a Boltzmann
      ! probability if it is higher. We are at T=1, so kT=1.
      if (newE .gt. oldE) then
        CALL RANDOM_NUMBER(RAN)
        if (EXP(-beta*(newE-oldE)) .lt. RAN) then
          nAccepted = .false.
        end if
      end if
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

      REAL*8 FUNCTION POTENTIAL(dist)
      IMPLICIT NONE
      REAL*8 dist

      POTENTIAL = 4.d0 * (1.d0/(dist**12) - 1.d0/(dist**6))

      END FUNCTION

      ! Subroutine: Prints a progress bar
      SUBROUTINE PROGRESS_BAR(step, total_steps)
      IMPLICIT NONE
      INTEGER step, total_steps, percent, i
      CHARACTER(len=50) bar

      percent = int((dble(step)/dble(total_steps))*100.d0)

      do i = 1, 50
        if (i .le. (percent/2)) then
          bar(i:i) = '#'
        else
          bar(i:i) = '-'
        end if
      end do

      write(*, '(A, A, A, I3, A)', advance='no') 
     +    char(13)//'[', bar, '] ', percent, '%'
      
      if (step .eq. total_steps) then
        write(*,*)
      end if

      END SUBROUTINE

