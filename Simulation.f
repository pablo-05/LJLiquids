      program MAIN
      Implicit NONE
      real*8 L, ran, deltamax, pi
      integer Np, Nmc, Nes, Nmeas, i, j, k
      logical sol
      
      ! Usamos .d0 para asegurar que son reales de doble precisión
      parameter (L=10.d0, Np=100, Nmc=100000, deltamax=L/100.d0
     &, pi=4.d0*datan(1.d0))
      real*8 r(Np, 3), rnew(3)
      print *, Np*4.d0/3*pi/(L**3)

      Nmeas = 0
      Nes = 0

      ! 1. Inicialización de posiciones
      ! Usamos dble() para evitar el error de la división entera en Fortran (cambiar inicialización)
      do i=1,Np
        do j=1,3
          r(i,j) = L * dble(i) / dble(Np)
        end do
      end do
      
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
        call distance(Np, i, r, rnew, L, sol)
        
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
      
      end program

      ! 3. Subrutina para calcular la distancia e imagen mínima
      SUBROUTINE distance(Np, i, r, rnew, L, sol)
      implicit NONE
      integer i, Np, j, k
      real*8 r(Np, 3), rnew(3), L, d(3)
      logical sol

      sol = .true.
      do k=1, Np
        if (k .ne. i) then
          do j=1, 3
            d(j) = dabs(rnew(j) - r(k,j))
            if (d(j) .gt. L / 2.d0) then
              d(j) = L - d(j)
            end if
          end do
          
          ! Comprobar solapamiento
          if (dsqrt(d(1)**2 + d(2)**2 + d(3)**2) .lt. 2.d0) then
            sol = .false.
            EXIT
          end if
        end if
      end do
      END SUBROUTINE

