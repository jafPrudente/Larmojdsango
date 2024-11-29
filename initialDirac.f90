subroutine initialDirac

   !===============================================!
   !      Esta subrutina soluciona el sistema      !
   !           estacionario de Dirac               !
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
   integer ward
   real(8) alpha0

   !------------------------------------------------
   ! Apicamos el shooting.
   write(*,*) 'Iniciamos el shooting...'
   call shootingDirac

   !------------------------------------------------
   ! Regresamos el reescalamiento.
   alpha0 = uno/( a(Nr)*alpha(Nr) )
   alpha  = alpha0*alpha

   !------------------------------------------------
   ! Imprimimos variables de interés.
   write(*,*) 'w   --> ', w
   write(*,*) 'wt  --> ', w*alpha0
   write(*,*) 's   --> ', a(Nr)*alpha(Nr)

   ! Guardamos en un archivo.
   open(ward, file = './' // trim(dirname) // '/sol.dat')

   do j=2, Nr, savedataR
      write(ward, "(5ES16.8)") r(j), a(j), alpha(j), F(j), G(j)
   end do

   close(ward)

   ! Llenamos la parte real e imaginaria de las funciones.
   do j=1, Nr
      F1(j) = F(j)
      F2(j) = cero
      G1(j) = cero
      G2(j) = G(j)
   end do

end subroutine initialDirac

!------------------------------------------------
! Definimos una subrutina para las derivadas.

subroutine derivadasDirac(w, r, a, alpha, F, G, da, dalpha, dF, dG)

   implicit none

   real(8) uno, dos
   real(8) w, r
   real(8) rho, srr, zet
   real(8) a, alpha, F, G
   real(8) da, dalpha, dF, dG

   uno    = 1.0D0
   dos    = 2.0D0

   rho = ( (dos*w)/(r**2 * a * alpha) )*(F**2 + G**2)
   srr = ( dos/(r**2 * a) )*( (w/alpha)*(F**2 + G**2) - (dos*F*G)/r - F**2 + G**2 )
   zet = ( r*a**2 )*(rho - srr) - (a**2 - uno)/r

   da = r*(a**3)*rho - a*(a**2 - uno)/(dos*r)
   dalpha = r*(a**2)*alpha*srr + alpha*(a**2 - uno)/(dos*r)
   dF = F*( zet/dos + a/r ) - G*( a + w*a/alpha )
   dG = G*( zet/dos - a/r ) - F*( a - w*a/alpha )

end subroutine derivadasDirac

!================================================
! Definimos una subrutina para solucionar.
subroutine solucionadorDirac

   !------------------------------------------------
   ! Usamos el módulo 'arrays' para declarar arreglos
   ! y el módulo 'vars' para las variables de entrada.
   use arrays
   use vars

   implicit none

   !------------------------------------------------
   ! Declaramos variables útiles.
   integer j
   real(8) F0

   ! Declaramos variables del método.
   real(8) k1a, k1alpha, k1F, k1G
   real(8) k2a, k2alpha, k2F, k2G
   real(8) k3a, k3alpha, k3F, k3G
   real(8) k4a, k4alpha, k4F, k4G

   ! Declaramos variables que guarden a las funciones derivadas.
   real(8) da, dalpha, dF, dG

   !------------------------------------------------
   ! Definimos las condiciones iniciales (con simetrías).
   F0   = campo0

   r(2) = dr
   a(2) = uno
   alpha(2) = uno
   F(2) = F0*dr
   G(2) = tercio*F0*(w - uno)*(dr**2)

   r(1) = -r(2)
   a(1) = a(2)
   alpha(1) = alpha(2)
   F(1) = -F(2)
   G(1) =  G(2)

   !------------------------------------------------
   ! Hacemos el método de Runge-Kutta en sí.
   do j=2, Nr-1

      !------------------------------------------------
      call derivadasDirac(w, r(j), a(j), alpha(j), F(j), G(j), &
         da, dalpha, dF, dG)

      k1a = da
      k1alpha = dalpha
      k1F = dF
      k1G = dG

      !------------------------------------------------
      call derivadasDirac(w, r(j) + medio*dr, a(j) + medio*k1a*dr, alpha(j) + medio*k1alpha*dr, F(j) + medio*k1F*dr, G(j) + medio*k1G*dr, &
         da, dalpha, dF, dG)

      k2a = da
      k2alpha = dalpha
      k2F = dF
      k2G = dG

      !------------------------------------------------
      call derivadasDirac(w, r(j) + medio*dr, a(j) + medio*k2a*dr, alpha(j) + medio*k2alpha*dr, F(j) + medio*k2F*dr, G(j) + medio*k2G*dr, &
         da, dalpha, dF, dG)

      k3a = da
      k3alpha = dalpha
      k3F = dF
      k3G = dG

      !------------------------------------------------
      call derivadasDirac(w, r(j) + dr, a(j) + k3a*dr, alpha(j) + k3alpha*dr, F(j) + k3F*dr, G(j) + k3G*dr, &
         da, dalpha, dF, dG)

      k4a = da
      k4alpha = dalpha
      k4F = dF
      k4G = dG

      !------------------------------------------------
      r(j+1) = r(j) + dr
      a(j+1) = a(j) + (k1a + dos*k2a + dos*k3a + k4a)*(dr*sexto)
      alpha(j+1) = alpha(j) + (k1alpha + dos*k2alpha + dos*k3alpha + k4alpha)*(dr*sexto)
      F(j+1) = F(j) + (k1F + dos*k2F + dos*k3F + k4F)*(dr*sexto)
      G(j+1) = G(j) + (k1G + dos*k2G + dos*k3G + k4G)*(dr*sexto)

   end do

end subroutine solucionadorDirac

!================================================
! Definimos una subrutina de shooting.
subroutine shootingDirac

   !------------------------------------------------
   ! Usamos el módulo 'arrays' para declarar arreglos
   ! y el módulo 'vars' para las variables de entrada.
   use arrays
   use vars

   implicit none

   integer counter

   real(8) dw, tolerancia

   ! Definimos la primer diferencial de 'w'.
   dw = 0.000001

   ! Definimos la tolerancia.
   tolerancia = dr**3

   ! Iniciamos el contador en ceros.
   counter = 0

   ! Solucionamos una vez.
   call solucionadorDirac

   !------------------------------------------------
   ! Iniciamos el sí el shooting.
   do while ( abs( F(Nr) ) > tolerancia )

      ! Tenemos dos opciones, que al final 'F' sea mayor o menor a cero.
      ! '>' para estados base;
      ! '<' para estados excitados.
      if ( F(Nr) > tolerancia ) then
         ! Actualizamos la w.
         w = w + dw

         ! Volvemos a solucionar.
         call solucionadorDirac

      else

         ! Aquí la 'w' ha sido muy grande, entonces debemos regresar.
         w = w - dw

         ! Ahora hacemos más pequeña la varaición y volvemos a avanzar.
         dw = medio*dw
         w = w + dw

         ! Volvemos a solucionar.
         call solucionadorDirac

      end if

      ! Actualizamos el contador.
      counter = counter + 1

   end do

   write(*,*) 'Iteraciones del shooting: ', counter

end subroutine shootingDirac
