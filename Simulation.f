      program MAIN
      Implicit NONE
      real*8 L, ran, deltamax, pi, dist, rbin, r1(3), r2(3)
      integer Np, Nmc, Nes, Nmeas, i, j, k, nbin
      logical sol
      
      ! Usamos .d0 para asegurar que son reales de doble precisión
      parameter (L=14.d0, Np=256, Nmc=1000000, deltamax=L/100.d0
     &, pi=4.d0*datan(1.d0), rbin=1.d0, 
     & nbin=int(l/rbin*dsqrt(3.d0)/2.d0))
      real*8 r(Np, 3), rnew(3), bines(0:nbin)
      print *, Np*4.d0/3*pi/(L**3)

      Nmeas = int(Np**(1.d0/3))
      Nes = Np

      ! 1. Inicialización de posiciones
      ! Usamos dble() para evitar el error de la división entera en Fortran (cambiar inicialización)
      do i=0,Nmeas
        do j=0,Nmeas
          do k= 0,Nmeas
            r(Nes,1) = 1+(2+(L-Nmeas*2)/(Nmeas))*i
            r(Nes,2) = 1+(2+(L-Nmeas*2)/(Nmeas))*j
            r(Nes,3) = 1+(2+(L-Nmeas*2)/(Nmeas))*k

            Nes = Nes - 1
            if (Nes.eq.0) goto 999
          end do
        end do
      end do
999   continue

      Nmeas = 0
      Nes = 0
      
      ! 2. Bucle principal de Monte Carlo
      do k=1, Nmc
        ! Elegir partícula al azar
        call random_number(ran)
        i = int(ran * Np) + 1
        
        ! Proponer movimiento
        do j=1,3
          call random_number(ran)
          ! Movimiento centrado en 0 usando reales (-1.0 a 1.0)
          rnew(j) = r(i,j) + deltamax * (2.d0 * ran - 1.d0)
          
          ! Aplicar Condiciones de Contorno Periódicas (Pac-Man)
          if (rnew(j) .ge. L) then
            rnew(j) = rnew(j) - L
          else if (rnew(j) .lt. 0.d0) then
            rnew(j) = rnew(j) + L
          end if
        end do
        
        Nes = Nes + 1
        
        ! Comprobar solapamientos con las demás partículas
        call solap(Np, i, r, rnew, L, sol)
        
        ! Aceptar o rechazar
        if (sol) then
          Nmeas = Nmeas + 1
          do j=1,3
            r(i,j) = rnew(j)
          end do
        end if
      end do

      ! Imprimir tasa de aceptación final
      print *, "Tasa de aceptacion: ", Nmeas * 1.d0 / Nes
      
      do i=1,Np
        do j=1,i
          do k=1,3
            r1(k)=r(i,k)
            r2(k)=r(j,k)
          end do
          call distance(r1, r2, L, dist)
          bines(int(dist/rbin)) = bines(int(dist/rbin)) + 2 ! Pendiente inicializar a 0 bines
        end do  
      end do

      end program

      ! 3. Subrutina para calcular la distancia e imagen mínima
      SUBROUTINE solap(Np, i, r, rnew, L, sol)
      implicit NONE
      integer i, Np, j, k
      real*8 r(Np, 3), rnew(3), L, dist, r1(3)
      logical sol

      sol = .true.
      do k=1, Np
        if (k .ne. i) then
          do j=1,3
            r1(j)=r(k,j)
          end do

          call distance(r1,rnew, L, dist)
          
          ! Comprobar solapamiento
          if (dist .lt. 2.d0) then
            sol = .false.
            EXIT

          end if
        end if
      end do
      END SUBROUTINE

      SUBROUTINE distance( r1, r2, L, dist)
      implicit NONE
      integer k
      real*8 r1(3), r2(3), L, d(3), dist

      do k=1, 3
        d(k) = dabs(r1(k) - r2(k))
        if (d(k) .gt. L / 2.d0) then
          d(k) = L - d(k)
        end if
      end do
      dist = dsqrt(d(1)**2 + d(2)**2 + d(3)**2)

      END SUBROUTINE

