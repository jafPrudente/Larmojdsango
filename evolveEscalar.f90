subroutine evolveEscalar

   !===============================================!
   !      Esta subrutina evoluciona el Sistema     !
   !               escalar con ICN                 !
   !===============================================!

   !------------------------------------------------
   ! Usamos el módulo 'arrays' para declarar arreglos
   ! y el módulo 'vars' para las variables de entrada.
   use vars
   use arrays

   implicit none

   !------------------------------------------------
   ! Definimos variables importantes.
   integer j, k, n

   real(8) dtw, dif
   real(8) Srr_p, Srr

   !------------------------------------------------
   ! Alojamos en memoria los arreglos.
   allocate( phi(Nr), vphi(Nr) )
   allocate( phi1(Nr), phi2(Nr) )
   allocate( phi1_f(Nr), phi2_f(Nr) )
   allocate( phi1_p(Nr), phi2_p(Nr) )

   allocate( pi1(Nr), pi2(Nr), psi1(Nr), psi2(Nr) )
   allocate( pi1_f(Nr), pi2_f(Nr), psi1_f(Nr), psi2_f(Nr) )
   allocate( pi1_p(Nr), pi2_p(Nr), psi1_p(Nr), psi2_p(Nr) )

   !------------------------------------------------
   ! Llamamos a las condiciones iniciales.
   call initialEscalar

   !------------------------------------------------
   ! Calculamos adot en dadas las condiciones iniciales.
   Srr = r(j)*alpha(j)*( psi1(j)*pi1(j) + psi2(j)*pi2(j) )
   adot(j) = -Srr

   !------------------------------------------------
   ! Si solo queremos la condición inicial ya la tendríamos
   ! y finalizamos esta subrutina con el 'return'.
   if (evolveType == 1) then
      return
   end if

   !------------------------------------------------
   ! Calculamos la perturbación y calculamos el valor de las funciones métricas ya perturbadas.
   write(*,*) '--------'
   write(*,*) 'm_i --> ', (rmax/dos)*( uno - (uno/a(Nr)**2) )

   if ( delta /= 0 ) then
      call perturbacionEscalar
      phi1 = phi1 + deltaEscalar

      do j=2, Nr-1
         psi1(j) = (phi1(j+1) - phi1(j-1))*medio/dr
      end do

      call metricEscalar
   end if
   write(*,*) 'm_p --> ', (rmax/dos)*( uno - (uno/a(Nr)**2) )

   !------------------------------------------------
   ! Imprimimos información inicial a pantalla.
   write(*,*) 'Comenzamos a evolucionar...'
   print *,'----------------------------'
   print *,'|  Iteracion  | Tiempo Fis |'
   print *,'----------------------------'

   write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',0,'    | ',t,'  |'

   !------------------------------------------------
   ! Guardamos información al tiempo cero.

   ! Abrimos los archivos.
   open(10, file = './' // trim(dirname) // '/phi1.rt')
   open(11, file = './' // trim(dirname) // '/phi2.rt')
   open(12, file = './' // trim(dirname) // '/phi.rt')
   open(13, file = './' // trim(dirname) // '/phi1.t')
   open(14, file = './' // trim(dirname) // '/phi2.t')

   open(20, file = './' // trim(dirname) // '/a.rt')
   open(21, file = './' // trim(dirname) // '/alpha.rt')
   open(22, file = './' // trim(dirname) // '/metric.rt')
   open(23, file = './' // trim(dirname) // '/adot.rt')
   open(24, file = './' // trim(dirname) // '/error.t')

   open(30, file = './' // trim(dirname) // '/masa.rt')
   open(31, file = './' // trim(dirname) // '/masa.t')

   ! Escribimos la información en sí.
   call wardaditoStatic
   call wardaditoDinamic

   !================================================
   ! Comenzamos con el ciclo principal de evolución.
   do n=1, Nt

      !------------------------------------------------
      ! Imprimimos información en pantalla.
      if (mod(n,savedataT) == 0) then
         write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',n,'    | ',t,'  |'
      end if

      !------------------------------------------------
      ! Avanzamos en el tiempo.
      t = t + dt

      !------------------------------------------------
      ! Guardamos información previa.
      phi1_p = phi1
      phi2_p = phi2

      psi1_p = psi1
      psi2_p = psi2

      pi1_p = pi1
      pi2_p = pi2

      a_p     = a
      alpha_p = alpha

      !================================================
      ! Iniciamos el ciclo interno del ICN.
      do k=1, 3

         !------------------------------------------------
         ! Calculamos los lados derechos.
         call fuenteEscalar

         !------------------------------------------------
         ! Calculamos el tamaño de paso temporal en función de k.
         if (k < 3) then
            dtw = medio*dt
         else
            dtw = dt
         end if

         !------------------------------------------------
         ! Actualizamos las funciones.
         do j=1, Nr
            phi1(j) = phi1_p(j) + dtw*phi1_f(j)
            phi2(j) = phi2_p(j) + dtw*phi2_f(j)
         end do

         do j=2, Nr-1
            psi1(j) = psi1_p(j) + dtw*psi1_f(j)
            psi2(j) = psi2_p(j) + dtw*psi2_f(j)
         end do

         do j=2, Nr-1
            pi1(j)  = pi1_p(j)  + dtw*pi1_f(j)
            pi2(j)  = pi2_p(j)  + dtw*pi2_f(j)
         end do

         !------------------------------------------------
         ! Definimos condiciones de frontera.

         ! Condiciones en el origen (simetrías).
         psi1(1) = - psi1(2)
         psi2(1) = - psi2(2)

         pi1(1)  = + pi1(2)
         pi2(1)  = + pi2(2)

         ! Condiciones en la frontera exterior.
         dif = uno + dtw/dr + dtw/r(Nr)

         pi1(Nr) = (pi1(Nr-1)*(dtw/dr - uno) + pi1_p(Nr)*(uno - dtw/dr - dtw/r(Nr)) &
            + pi1_p(Nr-1)*(dtw/dr + uno))/dif
         pi2(Nr) = (pi2(Nr-1)*(dtw/dr - uno) + pi2_p(Nr)*(uno - dtw/dr - dtw/r(Nr)) &
            + pi2_p(Nr-1)*(dtw/dr + uno))/dif

         psi1(Nr) = -pi1(Nr) - phi1(Nr)/r(Nr)
         psi2(Nr) = -pi2(Nr) - phi2(Nr)/r(Nr)

         !------------------------------------------------
         ! Evolucionamos la métrica.
         call metricEscalar

         !================================================
         ! Termina el ciclo interno del ICN.
      end do

      !------------------------------------------------
      ! Calculamos adot (12 de Argelia/Siddhartha).
      do j=1, Nr
         Srr_p = r(j)*alpha_p(j)*( psi1_p(j)*pi1_p(j) + psi2_p(j)*pi2_p(j) )
         Srr = r(j)*alpha(j)*( psi1(j)*pi1(j) + psi2(j)*pi2(j) )

         adot(j) = (a(j) - a_p(j))/dt &
            - medio*(Srr_p + Srr)
      end do

      !------------------------------------------------
      ! Guardamos información.

      ! Primero los arreglos puramente temporales.
      if ( mod(n,savedataT) == 0 ) then
         call wardaditoStatic
      end if

      ! Ahora los arreglos que dependen de 'r' y 't'.
      if ( mod(n,savedataT) == 0 ) then
         call wardaditoDinamic
      end if

      !================================================
      ! Finaliza el ciclo principal de evolución.
   end do
   write(*,*) 'Hemos terminado, Camarada UwUr'

   !------------------------------------------------
   ! Cerramos las puertas que abrimos.
   close(10)
   close(11)
   close(12)
   close(12)
   close(14)

   close(20)
   close(21)
   close(22)
   close(23)
   close(24)

end subroutine
