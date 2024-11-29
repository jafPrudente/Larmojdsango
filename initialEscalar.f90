subroutine initialEscalar

   !===============================================!
   !      Esta subrutina soluciona el sistema      !
   !         estacionario de Klein-Gordon          !
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
   call shootingKG

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
      write(ward, "(4ES16.8)") r(j), a(j), alpha(j), phi(j)
   end do

   close(ward)

   !------------------------------------------------
   ! Calculamos las funciones auxiliares.

   ! Llenamos la parte real e imaginaria de las funciones.
   phi1 = phi
   phi2 = cero

   ! Para Psi.
   psi1 = vphi
   psi2 = cero

   ! Ahora para Pi.
   pi1 = cero
   pi2 = -(w*alpha0)*phi*( a/alpha )

end subroutine

!------------------------------------------------
! Definimos una subrutina para las derivadas.
subroutine derivadasKG(w, r, a, alpha, phi, vphi, da, dalpha, dphi, dvphi)

   implicit none

   real(8) w, r, a, alpha, phi, vphi
   real(8) da, dalpha, dphi, dvphi
   real(8) uno, dos, medio

   uno = 1.0D0
   dos = 2.0D0
   medio = uno/dos

   da     = a*( (uno - a**2)/(dos*r) + (r*medio)*(vphi**2 + ((a*phi)**2)*(uno + (w/alpha)**2)) )
   dalpha = alpha*( (a**2 - uno)/r + da/a - r*(a*phi)**2 )
   dphi   = vphi
   dvphi  = -(dos/r + dalpha/alpha - da/a)*vphi + (a**2)*(uno - (w/alpha)**2)*phi

end subroutine derivadasKG

subroutine solucionadorKG

   !------------------------------------------------
   ! Usamos el módulo 'arrays' para declarar arreglos
   ! y el módulo 'vars' para las variables de entrada.
   use arrays
   use vars

   implicit none

   !------------------------------------------------
   ! Declaramos variables útiles.
   integer j
   real(8) phi0

   ! Declaramos variables del método.
   real(8) k1a, k1alpha, k1phi, k1vphi
   real(8) k2a, k2alpha, k2phi, k2vphi
   real(8) k3a, k3alpha, k3phi, k3vphi
   real(8) k4a, k4alpha, k4phi, k4vphi

   ! Declaramos variables que guarden a las funciones derivadas.
   real(8) da, dalpha, dphi, dvphi

   !------------------------------------------------
   ! Definimos las condiciones iniciales (con simetrías).
   phi0     = campo0

   r(2)     = dr
   a(2)     = 1
   alpha(2) = 1
   phi(2)   = phi0 !phi0 + (phi0*medio)*(uno - w**2)*r(1)**2
   vphi(2)  = 0 !phi0*(uno - w**2)*r(1)

   r(1)     = -r(2)
   a(1)     = +a(2)
   alpha(1) = +alpha(2)
   phi(1)   = +phi(2)
   vphi(1)  = -vphi(2)

   !------------------------------------------------
   ! Hacemos el método de Runge-Kutta en sí.
   do j=2, Nr-1

      !------------------------------------------------
      call derivadasKG(w, r(j), a(j), alpha(j), phi(j), vphi(j), &
         da, dalpha, dphi, dvphi)

      k1a     = da
      k1alpha = dalpha
      k1phi   = dphi
      k1vphi  = dvphi

      !------------------------------------------------
      call derivadasKG(w, r(j) + medio*dr, a(j) + medio*k1a*dr, alpha(j) + medio*k1alpha*dr, phi(j) + medio*k1phi*dr, &
         vphi(j) + medio*k1vphi*dr, da, dalpha, dphi, dvphi)

      k2a     = da
      k2alpha = dalpha
      k2phi   = dphi
      k2vphi  = dvphi

      !------------------------------------------------
      call derivadasKG(w, r(j) + medio*dr, a(j) + medio*k2a*dr, alpha(j) + medio*k2alpha*dr, phi(j) + medio*k2phi*dr, &
         vphi(j) + medio*k2vphi*dr, da, dalpha, dphi, dvphi)

      k3a     = da
      k3alpha = dalpha
      k3phi   = dphi
      k3vphi  = dvphi

      !------------------------------------------------
      call derivadasKG(w, r(j) + dr, a(j) + k3a*dr, alpha(j) + k3alpha*dr, phi(j) + k3phi*dr, vphi(j) + k3vphi*dr, &
         da, dalpha, dphi, dvphi)

      k4a     = da
      k4alpha = dalpha
      k4phi   = dphi
      k4vphi  = dvphi

      !------------------------------------------------
      r(j+1)     = r(j) + dr
      a(j+1)     = a(j) + (k1a + dos*k2a + dos*k3a + k4a)*(dr*sexto)
      alpha(j+1) = alpha(j) + (k1alpha + dos*k2alpha + dos*k3alpha + k4alpha)*(dr*sexto)
      phi(j+1)   = phi(j) + (k1phi + dos*k2phi + dos*k3phi + k4phi)*(dr*sexto)
      vphi(j+1)  = vphi(j) + (k1vphi + dos*k2vphi + dos*k3vphi + k4vphi)*(dr*sexto)

   end do

end subroutine solucionadorKG

!================================================
! Definimos una subrutina de shooting.
subroutine shootingKG

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
   call solucionadorKG

   !------------------------------------------------
   ! Iniciamos el sí el shooting.
   do while ( abs( phi(Nr) ) > tolerancia )

      ! Tenemos dos opciones, que al final 'F' sea mayor o menor a cero.
      ! '>' para estados base;
      ! '<' para estads excitados.
      if ( phi(Nr) > tolerancia ) then
         ! Actualizamos la w.
         w = w + dw

         ! Volvemos a solucionar.
         call solucionadorKG

      else

         ! Aquí la 'w' ha sido muy grande, entonces debemos regresar.
         w = w - dw

         ! Ahora hacemos más pequeña la varaición y volvemos a avanzar.
         dw = medio*dw
         w = w + dw

         ! Volvemos a solucionar.
         call solucionadorKG

      end if

      ! Actualizamos el contador.
      counter = counter + 1

   end do

   write(*,*) 'Iteraciones del shooting: ', counter

end subroutine shootingKG
