subroutine initialProca

   !===============================================!
   !      Esta subrutina soluciona el sistema      !
   !           estacionario de Proca               !
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
   real(8) s0

   !------------------------------------------------
   ! Apicamos el shooting.
   write(*,*) 'Iniciamos el shooting...'
   call shootingProca

   !------------------------------------------------
   ! Regresamos el reescalamiento.
   s0 = 1/s(Nr)

   !------------------------------------------------
   ! Imprimimos variables de interés.
   write(*,*) 'w   --> ', w
   write(*,*) 'wt  --> ', w*s0
   write(*,*) 's   --> ', s(Nr)*s0

   ! Guardamos en un archivo.
   open(ward, file = './' // trim(dirname) // '/sol.dat')

   do j=2, Nr, savedataR
      write(ward, "(5ES16.8)") r(j), m(j), s(j), F(j), G(j)
   end do

   close(ward)

   ! Llenamos la parte real e imaginaria de las funciones.
   do j=1, Nr
      F1(j) = F(j)
      F2(j) = cero
      G1(j) = cero
      G2(j) = G(j)
   end do

   !------------------------------------------------
   ! Encontremos los coeficientes métricos.
   do j=2, Nr
      a(j)     = uno/dsqrt(uno - (dos*m(j)/r(j)) )
      alpha(j) = s0*s(j)*dsqrt(uno - (dos*m(j)/r(j)) )
   end do

   a(1)     = a(2)
   alpha(1) = alpha(2)

end subroutine initialProca

 !------------------------------------------------
 ! Definimos una subrutina para las derivadas.

subroutine derivadasProca(w, r, m, s, F, G, dm, ds, dF, dG)

   implicit none

   real(8) uno, dos, cuatro, medio
   real(8) w, r, m, s, F, G, N
   real(8) dN, dm, ds, dF, dG

   medio  = 1.0D0/2.0D0
   uno    = 1.0D0
   dos    = 2.0D0
   cuatro = 4.0D0

   N  = uno - dos*m/r
   dN = dos*(m/r**2) - (dos/r)*( (medio*r**2)*( ((s*N*G)/w)**2 + G**2*N + (F**2)/(N*s**2) ) )

   dm = (medio*r**2)*( ((s*N*G)/w)**2 + G**2*N + (F**2)/(N*s**2) )
   ds = r*s*( G**2 + (F/(N*s))**2 )
   dF = (w - (s**2*N)/w)*G
   dG = - ( ds/s + dos/r + dN/N )*G - (w/(s**2*N**2))*F

end subroutine derivadasProca

 !================================================
 ! Definimos una subrutina para solucionar.
subroutine solucionadorProca

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
   real(8) k1m, k1s, k1F, k1G
   real(8) k2m, k2s, k2F, k2G
   real(8) k3m, k3s, k3F, k3G
   real(8) k4m, k4s, k4F, k4G

   ! Declaramos variables que guarden a las funciones derivadas.
   real(8) dm, ds, dF, dG

   !------------------------------------------------
   ! Definimos las condiciones iniciales (con simetrías).
   F0   = campo0

   r(2) = dr
   m(2) = (sexto*F0**2)*r(2)**3
   s(2) = uno + (medio*F0**2)*r(2)**2
   F(2) = F0 + ( sexto*(F0*(uno - w**2)) )*r(2)**2
   G(2) = -tercio*F0*w*r(2)

   r(1) = -r(2)
   m(1) = -m(2)
   s(1) =  s(2)
   F(1) =  F(2)
   G(1) = -G(2)

   !------------------------------------------------
   ! Hacemos el método de Runge-Kutta en sí.
   do j=2, Nr-1

      !------------------------------------------------
      call derivadasProca(w, r(j), m(j), s(j), F(j), G(j), &
         dm, ds, dF, dG)

      k1m = dm
      k1s = ds
      k1F = dF
      k1G = dG

      !------------------------------------------------
      call derivadasProca(w, r(j) + medio*dr, m(j) + medio*k1m*dr, s(j) + medio*k1s*dr, F(j) + medio*k1F*dr, G(j) + medio*k1G*dr, &
         dm, ds, dF, dG)

      k2m = dm
      k2s = ds
      k2F = dF
      k2G = dG

      !------------------------------------------------
      call derivadasProca(w, r(j) + medio*dr, m(j) + medio*k2m*dr, s(j) + medio*k2s*dr, F(j) + medio*k2F*dr, G(j) + medio*k2G*dr, &
         dm, ds, dF, dG)

      k3m = dm
      k3s = ds
      k3F = dF
      k3G = dG

      !------------------------------------------------
      call derivadasProca(w, r(j) + dr, m(j) + k3m*dr, s(j) + k3s*dr, F(j) + k3F*dr, G(j) + k3G*dr, &
         dm, ds, dF, dG)

      k4m = dm
      k4s = ds
      k4F = dF
      k4G = dG

      !------------------------------------------------
      r(j+1) = r(j) + dr
      m(j+1) = m(j) + (k1m + dos*k2m + dos*k3m + k4m)*(dr*sexto)
      s(j+1) = s(j) + (k1s + dos*k2s + dos*k3s + k4s)*(dr*sexto)
      F(j+1) = F(j) + (k1F + dos*k2F + dos*k3F + k4F)*(dr*sexto)
      G(j+1) = G(j) + (k1G + dos*k2G + dos*k3G + k4G)*(dr*sexto)

   end do

end subroutine solucionadorProca

 !================================================
 ! Definimos una subrutina de shooting.
subroutine shootingProca

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
   call solucionadorProca

   !------------------------------------------------
   ! Iniciamos el sí el shooting.
   do while ( abs( F(Nr) ) > tolerancia )

      ! Tenemos dos opciones, que al final 'F' sea mayor o menor a cero.
      ! '>' para estados base;
      ! '<' para estads excitados.
      if ( F(Nr) > tolerancia ) then
         ! Actualizamos la w.
         w = w + dw

         ! Volvemos a solucionar.
         call solucionadorProca

      else

         ! Aquí la 'w' ha sido muy grande, entonces debemos regresar.
         w = w - dw

         ! Ahora hacemos más pequeña la varaición y volvemos a avanzar.
         dw = medio*dw
         w = w + dw

         ! Volvemos a solucionar.
         call solucionadorProca

      end if

      ! Actualizamos el contador.
      counter = counter + 1

   end do

   write(*,*) 'Iteraciones del shooting: ', counter

end subroutine shootingProca

