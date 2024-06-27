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
   integer i
   real(8) idrh
   idrh = medio/dr
   
   !------------------------------------------------
   ! Llenamos las fuentes de phi.
     do i=1, Nr
        phi1_f(i) = alpha(i)*pi1(i)/a(i)
        phi2_f(i) = alpha(i)*pi2(i)/a(i)
     end do
   
   !------------------------------------------------
   ! Llenamos las fuentes de psi.
     do i=2, Nr-1
        psi1_f(i) = idrh*(alpha(i+1)*pi1(i+1)/a(i+1) &
                - alpha(i-1)*pi1(i-1)/a(i-1))
        psi2_f(i) = idrh*(alpha(i+1)*pi2(i+1)/a(i+1) &
                - alpha(i-1)*pi2(i-1)/a(i-1))
     end do
   
   !------------------------------------------------
   ! Llenamos las fuentes de pi.
     do i=2, 10
        pi1_f(i) = 3.0D0*((alpha(i+1)*r(i+1)**2*psi1(i+1)/a(i+1) &
               - alpha(i-1)*r(i-1)**2*psi1(i-1)/a(i-1)) &
               /(r(i+1)**3 - r(i-1)**3)) &
               - a(i)*alpha(i)*phi1(i)
   
        pi2_f(i) = 3.0D0*((alpha(i+1)*r(i+1)**2*psi2(i+1)/a(i+1) &
               - alpha(i-1)*r(i-1)**2*psi2(i-1)/a(i-1)) &
               /(r(i+1)**3 - r(i-1)**3)) &
               - a(i)*alpha(i)*phi2(i)
     end do
   
     do i=11, Nr-1
        pi1_f(i) = idrh*((alpha(i+1)*r(i+1)**2*psi1(i+1)/a(i+1) &
               - alpha(i-1)*r(i-1)**2*psi1(i-1)/a(i-1)) &
               )/r(i)**2 &
               - a(i)*alpha(i)*phi1(i)
   
        pi2_f(i) = idrh*((alpha(i+1)*r(i+1)**2*psi2(i+1)/a(i+1) &
               - alpha(i-1)*r(i-1)**2*psi2(i-1)/a(i-1)) &
               )/r(i)**2 &
               - a(i)*alpha(i)*phi2(i)
     end do
   
end subroutine
