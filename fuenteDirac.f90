subroutine fuenteDirac

   !===============================================!
   !      Esta subrutina calcula las fuentes       !
   !                   de Dirac                    !
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

   !------------------------------------------------
   ! Definimos el rhs para F1, F2, G1 y G2 son (35) de Daka.

   do j=3, Nr-1
      F1_f(j) = -( alpha(j)/r(j) )*( G1(j) - r(j)*F2(j) ) &
         - dsqrt( alpha(j)/a(j) )*( &
         + dsqrt( alpha(j+1)/a(j+1) )*G1(j+1) &
         - dsqrt( alpha(j-1)/a(j-1) )*G1(j-1) )*medio/dr

      F2_f(j) = -( alpha(j)/r(j) )*( G2(j) + r(j)*F1(j) ) &
         - dsqrt( alpha(j)/a(j) )*( &
         + dsqrt( alpha(j+1)/a(j+1) )*G2(j+1) &
         - dsqrt( alpha(j-1)/a(j-1) )*G2(j-1) )*medio/dr

      G1_f(j) = ( alpha(j)/r(j) )*( F1(j) - r(j)*G2(j) ) &
         - dsqrt( alpha(j)/a(j) )*( &
         + dsqrt( alpha(j+1)/a(j+1) )*F1(j+1) &
         - dsqrt( alpha(j-1)/a(j-1) )*F1(j-1) )*medio/dr

      G2_f(j) = ( alpha(j)/r(j) )*( F2(j) + r(j)*G1(j) ) &
         - dsqrt( alpha(j)/a(j) )*( &
         + dsqrt( alpha(j+1)/a(j+1) )*F2(j+1) &
         - dsqrt( alpha(j-1)/a(j-1) )*F2(j-1) )*medio/dr
   end do

   F1_f(2) = F1_f(3)
   F2_f(2) = F2_f(3)
   G1_f(2) = G1_f(3)
   G2_f(2) = G2_f(3)

end subroutine
