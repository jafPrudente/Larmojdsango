subroutine perturbacionEscalar

   !===============================================!
   !    Esta subrutina calcula la perturbación     !
   !               del sistema de KG               !
   !===============================================!

   !------------------------------------------------
   ! Usamos el módulo 'arrays' para declarar arreglos
   ! y el módulo 'vars' para las variables de entrada.
   use arrays
   use vars

   implicit none

   !------------------------------------------------
   ! Declaramos variables útiles.
   integer indicePerturbador

   real(8) amp, sigma, r0

   real(8) medio_maximo, FWHM
   integer j, indice_der
   logical encontrado_der

   !================================================
   ! Definimos la amplitud de la perturbación.
   amp = delta*phi1(2)

   !================================================
   ! Definimos la anchura a media altura (FWHM).

   ! Encontramos el valor de la mitad de la altura.
   medio_maximo = phi1(2)/2

   ! Iniciamos variables lógicas que dirán cuándo encontramos los valores de la anchura.
   encontrado_der = .false.

   ! Con un ciclo buscamos los índices de F1 donde se encuentra la media altura.
   do j=1, Nr
      if ( phi1(j) <= medio_maximo ) then
         indice_der = j
         encontrado_der = .true.
         exit
      end if
   end do

   ! Calculamos el valor de la FWHM en sí.
   if (encontrado_der) then
      FWHM = r(indice_der) - r(2)
   end if

   ! Calculamos la sigma en sí.
   sigma = deltaFWHM*FWHM

   !================================================
   ! Encontramos el valor de r0 donde centraremos la perturbación (será en el máximo).
   r0 = r(2)

   !================================================
   ! Definimos la perturbación en sí.
   deltaEscalar = amp*dexp( -(r - r0)**2/(dos*sigma**2) )

   ! Guardamos la perturbación en un archivo.
   open(indicePerturbador, file = './' // trim(dirname) // '/perturbacion.dat')

   do j=2, Nr, savedataR
      write(indicePerturbador, "(3ES16.8)") r(j), deltaEscalar(j), phi1(j) + deltaEscalar(j)
   end do

   close(indicePerturbador)

end subroutine
