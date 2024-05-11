subroutine metricDirac

   !===============================================!
   !      Esta subrutina calcula la métrica        !
   !             del sistema de Dirac              !
   !===============================================!
   ! Solucionamos (36) de Daka (con las convenciones de (37)).

   !------------------------------------------------
   ! Usamos el módulo 'arrays' para declarar arreglos
   ! y el módulo 'vars' para las variables de entrada.
   use arrays
   use vars

   implicit none

   !------------------------------------------------
   ! Declaramos variables útiles.
   integer j

   real(8) rho, Srr, geo
   real(8) aux, aux_p
   real(8) k1, k2
   real(8) x, y

   real(8), allocatable, dimension (:) :: xokas

   !------------------------------------------------
   ! Alojamos en memoria el arreglo.
   allocate( xokas(Nr) )

   !------------------------------------------------
   ! Solucionaremos usando Runge-Kutta a orden 2 (regla de dos tercios).

   !------------------------------------------------
   ! Definimos la parte de las derivadas.
   do j=3, Nr-1
      xokas(j) = F1(j)*(G2(j+1) - G2(j-1))*medio/dr &
         -F2(j)*(G1(j+1) - G1(j-1))*medio/dr &
         +G1(j)*(F2(j+1) - F2(j-1))*medio/dr &
         -G2(j)*(F1(j+1) - F1(j-1))*medio/dr
   end do

   xokas(2)  = F1(2)*(-G2(4) + 4*G2(3) - 3*G2(2))/(2*dr) &
      - F2(2)*(-G1(4) + 4*G1(3) - 3*G1(2))/(2*dr) &
      + G1(2)*(-F2(4) + 4*F2(3) - 3*F2(2))/(2*dr) &
      - G2(2)*(-F1(4) + 4*F1(3) - 3*F1(2))/(2*dr)

   xokas(Nr) = F1(Nr)*(3*G2(Nr) - 4*G2(Nr-1) + G2(Nr-2))/(2*dr) &
      - F2(Nr)*(3*G1(Nr) - 4*G1(Nr-1) + G1(Nr-2))/(2*dr) &
      + G1(Nr)*(3*F2(Nr) - 4*F2(Nr-1) + F2(Nr-2))/(2*dr) &
      - G2(Nr)*(3*F1(Nr) - 4*F1(Nr-1) + F1(Nr-2))/(2*dr)

   !------------------------------------------------
   ! Primero para a.
   a(2) = uno

   do j=2, Nr-1

      rho = ( dos/( (r(j)*a(j))**2 ) )*( a(j)*( F1(j)**2 + F2(j)**2 - G1(j)**2 - G2(j)**2 ) &
         + ( dos*a(j)/r(j) )*( F1(j)*G2(j) - F2(j)*G1(j) ) &
         + xokas(j) )

      geo = ( a(j)**2 - uno )/( dos*r(j) )

      aux_p = r(j)*( a(j)**2 )*rho - geo

      k1  = a(j)*aux_p

      x   = r(j) + dr*medio
      y   = a(j) + k1*dr*medio

      !------------------------------------------------
      rho = ( dos/( (x*y)**2 ) )*( y*( F1(j+1)**2 + F2(j+1)**2 - G1(j+1)**2 - G2(j+1)**2 ) &
         + ( dos*y/x )*( F1(j+1)*G2(j+1) - F2(j+1)*G1(j+1) ) &
         + xokas(j+1) )

      geo = ( y**2 - uno )/( dos*x )

      aux = x*( y**2 )*rho - geo

      k2  = y*( aux + aux_p )*medio

      !------------------------------------------------
      a(j+1) = a(j) + k2*dr

   end do

   a(1) = a(2)

   !------------------------------------------------
   ! Ahora solucionamos para alpha.
   ! Necesitamos solucionar del final hacia el inicio porque
   ! solo conocemos el valor en la frontera, no en el inicio.
   alpha(Nr) = uno/a(Nr)

   do j=Nr-1, 2, -1

      Srr = ( dos/( (r(j+1)*a(j+1))**2 ) )*xokas(j+1)

      geo = ( a(j+1)**2 - uno )/( dos*r(j+1) )

      aux_p = r(j+1)*(a(j+1)**2)*Srr + geo

      k1  = alpha(j+1)*aux_p

      x  = r(j+1) - medio*dr
      y  = alpha(j+1) - medio*k1*dr

      !------------------------------------------------
      Srr = ( dos/( (x*a(j))**2 ) )*xokas(j)

      geo = ( a(j)**2 - uno )/( dos*x )

      aux = x*(a(j)**2)*Srr + geo

      k2 = y*medio*( aux + aux_p )

      !------------------------------------------------
      alpha(j) = alpha(j+1) - k2*dr

   end do

   alpha(1) = alpha(2)

end subroutine
