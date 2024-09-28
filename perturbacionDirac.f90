subroutine perturbacionDirac

   !===============================================!
   !    Esta subrutina calcula la perturbación     !
   !             del sistema de Dirac              !
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
   integer j, indice_der, indice_izq
   logical encontrado_der, encontrado_izq

   integer indice_max

   !================================================
   ! Definimos la amplitud de la perturbación.
   amp = delta*maxval(F1)

   !================================================
   ! Definimos la anchura a media altura (FWHM).

   ! Encontramos el valor de la mitad de la altura.
   medio_maximo = maxval(F1)/2

   ! Iniciamos variables lógicas que dirán cuándo encontramos los valores de la anchura.
   encontrado_izq = .false.
   encontrado_der = .false.

   ! Con un ciclo buscamos los índices de F1 donde se encuentra la media altura.
   do j=1, Nr
      if (.not. encontrado_izq .and. F1(j) >= medio_maximo) then
         indice_izq = j
         encontrado_izq = .true.
      end if

      if (encontrado_izq .and. F1(j) <= medio_maximo) then
         indice_der = j
         encontrado_der = .true.
         exit
      end if
   end do

   ! Calculamos el valor de la FWHM en sí.
   if (encontrado_izq .and. encontrado_der) then
      FWHM = r(indice_der) - r(indice_izq)
   end if

   ! Calculamos la sigma en sí.
   sigma = deltaFWHM*FWHM

   !================================================
   ! Encontramos el valor de r0 donde centraremos la perturbación (será en el máximo).
   indice_max = maxloc(F1, dim=1)
   r0 = r(indice_max)

   !================================================
   ! Definimos la perturbación en sí.
   deltaDirac = amp*(r/r0)*dexp( -(r - r0)**2/(dos*sigma**2) )

   ! Guardamos la perturbación en un archivo.
   open(indicePerturbador, file = './' // trim(dirname) // '/perturbacion.dat')

   do j=2, Nr, savedataR
      write(indicePerturbador, "(3ES16.8)") r(j), deltaDirac(j), F1(j) + deltaDirac(j)
   end do

   close(indicePerturbador)

end subroutine
