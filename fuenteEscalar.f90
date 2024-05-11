subroutine fuenteEscalar

   !===============================================!
   !      Esta subrutina calcula las fuentes       !
   !                    escalares                  !
   !===============================================!

   !------------------------------------------------
   ! Usamos el módulo 'arrays' para declarar arreglos
   ! y el módulo 'vars' para las variables de entrada.
   use arrays
   use vars

   implicit none

   !------------------------------------------------
   ! Declaramos variables útiles.
   integer j

   ! Definimos el rhs para \phi, son las ecs (9) de Argelia/Siddhartha.
   do j=1, Nr
      phi1_f(j) = ( alpha(j)/a(j) )*pi1(j)
      phi2_f(j) = ( alpha(j)/a(j) )*pi2(j)
   end do

   ! Definimos el rhs para las \psi, son las ecs (9) de Argelia/Siddhartha.
   do j=2, Nr-1
      psi1_f(j) = ( (alpha(j+1)/a(j+1))*pi1(j+1) - (alpha(j-1)/a(j-1))*pi1(j-1) )/(dos*dr)
      psi2_f(j) = ( (alpha(j+1)/a(j+1))*pi2(j+1) - (alpha(j-1)/a(j-1))*pi2(j-1) )/(dos*dr)

      pi1_f(j) = (uno/r(j)**2)*( r(j+1)**2*(alpha(j+1)/a(j+1))*psi1(j+1) &
         - r(j-1)**2*(alpha(j-1)/a(j-1))*psi1(j-1) )/(dos*dr) &
         -a(j)*alpha(j)*phi1(j)

      pi2_f(j) = (uno/r(j)**2)*( r(j+1)**2*(alpha(j+1)/a(j+1))*psi2(j+1) &
         - r(j-1)**2*(alpha(j-1)/a(j-1))*psi2(j-1) )/(dos*dr) &
         -a(j)*alpha(j)*phi2(j)
   end do

end subroutine
