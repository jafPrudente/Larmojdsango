program main

   ! Verde que te quiero verde.

   !===============================================!
   !      Este programa evoluciona estrellas       !
   !===============================================!

   !------------------------------------------------
   ! Usamos el módulo 'arrays' para declarar arreglos
   ! y el módulo 'vars' para las variables de entrada.
   use arrays
   use vars

   implicit none

   ! Llamamos la subrutina que lee la información de entrada.
   call lector

   !------------------------------------------------
   ! Allojamos en memoria los arreglos generales para comenzar la evolución.
   allocate( r(Nr) )
   allocate( m(Nr), s(Nr) )
   allocate( a(Nr), alpha(Nr) )
   allocate( a_p(Nr), alpha_p(Nr) )
   allocate( adot(Nr) )

   !------------------------------------------------
   ! Evolucionamos en sí: tenemos una opción
   ! para cada tipo de materia.
   if(      matterType == 1 ) then
      call evolveEscalar

   else if( matterType == 2 ) then
      call evolveDirac

   else if( matterType == 3 ) then
      call evolveProca
   end if


end program
