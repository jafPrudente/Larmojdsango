subroutine metricEscalar

   !===============================================!
   !      Esta subrutina calcula la métrica        !
   !              del sistema escalar              !
   !===============================================!
   ! Solucionamos (10) y (11) de Argelia/Siddhartha.

   !------------------------------------------------
   ! Usamos el módulo 'arrays' para declarar arreglos
   ! y el módulo 'vars' para las variables de entrada.
   use arrays
   use vars

   implicit none

   !------------------------------------------------
   ! Declaramos variables útiles.
   integer j

   real(8) rho, geo
   real(8) aux, aux_p
   real(8) k1, k2
   real(8) x, y

   !------------------------------------------------
   ! Solucionaremos usando Runge-Kutta a orden 2 (regla de dos tercios).

   !------------------------------------------------
   ! Primero para a.
   a(2) = uno

   do j=2, Nr-1

      rho = (r(j)*medio)*( pi1(j)**2 + pi2(j)**2 + psi1(j)**2 + psi2(j)**2 &
         + (a(j)**2)*(phi1(j)**2 + phi2(j)**2) )

      geo = ( a(j)**2 - uno )/( dos*r(j) )

      aux_p = rho - geo

      k1  = a(j)*aux_p

      x   = r(j) + dr*medio
      y   = a(j) + k1*dr*medio

      !------------------------------------------------
      rho = (x*medio)*( pi1(j+1)**2 + pi2(j+1)**2 + psi1(j+1)**2 + psi2(j+1)**2 &
         + y**2*(phi1(j+1)**2 + phi2(j+1)**2) )

      geo = ( y**2 - uno )/( dos*x )

      aux = rho - geo

      k2  = y*medio*( aux + aux_p )

      !------------------------------------------------
      a(j+1) = a(j) + k2*dr

   end do

   a(1) = a(2)

   !------------------------------------------------
   ! Ahora solucionamos para alpha.
   ! Necesitamos solucionar del final hacia el inicio porque
   ! solo conocemos el valor en la frontera, no en el inicio.
   alpha(Nr) = uno/a(Nr) !-a(Nr)*(phi1(Nr)/r(Nr) + psi1(Nr))/pi1(Nr)

   do j=Nr-1, 2, -1

      rho = (r(j+1)*medio)*( pi1(j+1)**2 + pi2(j+1)**2 + psi1(j+1)**2 + psi2(j+1)**2 &
         - a(j+1)**2*(phi1(j+1)**2 + phi2(j+1)**2) )

      geo = ( a(j+1)**2 - uno )/( dos*r(j+1) )

      aux_p = geo + rho

      k1 = alpha(j+1)*aux_p

      x  = r(j+1) - dr*medio
      y  = alpha(j+1) - k1*dr*medio

      !------------------------------------------------
      rho = (x*medio)*( pi1(j)**2 + pi2(j)**2 + psi1(j)**2 + psi2(j)**2 &
         - a(j)**2*(phi1(j)**2 + phi2(j)**2) )

      geo = ( a(j)**2 - uno )/( dos*x )

      aux = geo + rho

      k2  = y*medio*( aux + aux_p )

      !------------------------------------------------
      alpha(j) = alpha(j+1) - k2*dr

   end do

   alpha(1) = alpha(2)

end subroutine
