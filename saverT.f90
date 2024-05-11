subroutine wardaditoStatic

   !===============================================!
   !      Esta subrutina guarda los arreglos       !
   !                  generados                    !
   !===============================================!

   !------------------------------------------------
   ! Usamos el módulo 'arrays' para declarar arreglos
   ! y el módulo 'vars' para las variables de entrada.
   use arrays
   use vars

   implicit none

   !===================================================
   ! Dependiendo del tipo de materia guardaremos distintos datos.
   if(      matterType == 1 ) then
      !------------------------------------------------
      ! Para Klein-Gordon.
      write(13,"(3ES16.8)") t, phi1(2)
      write(14,"(3ES16.8)") t, phi2(2)
      write(24,"(3ES16.8)") t, maxval( abs(adot) )


      !===================================================
   else if( matterType == 2 ) then
      !------------------------------------------------
      ! Para Dirac.
      write(13,"(2ES16.8)") t, F1(586)
      write(78,"(2ES16.8)") t, 0.11704489*dcos(0.924*t + 0.15)
      write(14,"(2ES16.8)") t, maxval( F2 )
      write(23,"(2ES16.8)") t, maxval( G1 )
      write(24,"(2ES16.8)") t, maxval( G2 )

      write(35,"(2ES16.8)") t, maxval( abs(adot) )
      write(41,"(2ES16.8)") t, (rmax/dos)*( uno - (uno/a(Nr)**2) )


      !===================================================
   else if( matterType == 3 ) then
      call evolveProca
   end if

end subroutine
