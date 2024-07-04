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
   integer i

   real(8) idr
   real(8) aux, ssa, rk1, rhoa, rhoa1, rhoa2, rrhoa
   real(8) ssalpha, aux1, ssa2, aux2

   idr = uno/dr

   !------------------------------------------------
   ! Comenzamos solucionando para a.
   a(1) = uno
   a(2) = uno

   do i=3, Nr
      rhoa = medio*( (pi1(i-1)**2  + pi2(i-1)**2 &
                   + psi1(i-1)**2 + psi2(i-1)**2)/a(i-1)**2 &
                   + phi1(i-1)**2 + phi2(i-1)**2 )

      ssa = a(i-1)*( medio*(uno - a(i-1)**2)/r(i-1) &
            + r(i-1)*a(i-1)**2*rhoa )

      rk1 = a(i-1) + medio*dr*ssa

      !------------------------------------------------
      rhoa1 = medio*( (pi1(i-1)**2  + pi2(i-1)**2 &
                   + psi1(i-1)**2 + psi2(i-1)**2)/rk1**2 &
                   + phi1(i-1)**2 + phi2(i-1)**2 )

      rhoa2 = medio*( (pi1(i)**2  + pi2(i)**2 &
                   + psi1(i)**2 + psi2(i)**2)/rk1**2 &
                   + phi1(i)**2 + phi2(i)**2 )

      rrhoa = medio* (r(i-1)*rhoa1 + r(i)*rhoa2)

      ssa = rk1* ( (uno - rk1**2)/(r(i-1) + r(i)) + rk1**2*( rrhoa ))

      a(i) = a(i-1) + dr*(ssa)
   end do


   !------------------------------------------------
   ! Ahora solucionamos para a.
   alpha(Nr) = uno/a(Nr)

   do i=Nr-1, 2, -1
      rhoa = medio*( (pi1(i+1)**2  + pi2(i+1)**2 &
                   + psi1(i+1)**2 + psi2(i+1)**2)/a(i+1)**2 &
                   + phi1(i+1)**2 + phi2(i+1)**2 )

      ssa =  medio*(uno - a(i+1)**2)/r(i+1) &
            + r(i+1)*a(i+1)**2*rhoa 

      aux1 = ( a(i+1)**2 - uno )/r(i+1) +  ssa  & 
            - r(i+1)*a(i+1)**2*(phi1(i+1)**2 + phi2(i+1)**2)

      ssalpha = alpha(i+1)*aux1 

      aux = alpha(i+1) - medio*dr*ssalpha

      !------------------------------------------------
      rhoa = medio*( (pi1(i)**2  + pi2(i)**2 &
                   + psi1(i)**2 + psi2(i)**2)/a(i)**2 &
                   + phi1(i)**2 + phi2(i)**2 )

      ssa2 =  medio*(uno - a(i)**2)/r(i) &
            + r(i)*a(i)**2*rhoa 

      aux2 = ( a(i)**2 - uno )/r(i) +  ssa2  & 
            - r(i)*a(i)**2*(phi1(i)**2 + phi2(i)**2)

      ssalpha = medio*(aux1+aux2)
  
      alpha(i) = alpha(i+1) - dr*aux*(ssalpha)
  end do

  alpha(1) = alpha(2)

end subroutine
