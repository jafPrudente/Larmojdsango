subroutine wardaditoDinamic

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

   integer j

   !===================================================
   ! Dependiendo del tipo de materia guardaremos distintos datos.
   if(      matterType == 1 ) then
      !------------------------------------------------
      ! Para Klein-Gordon.
      do j=2, Nr, savedataR
         write(10,"(3ES16.8)") t, r(j), phi1(j)
         write(11,"(3ES16.8)") t, r(j), phi2(j)
         write(12,"(3ES16.8)") t, r(j), phi1(j)**2 + phi2(j)**2

         write(20,"(3ES16.8)") t, r(j), a(j)
         write(21,"(3ES16.8)") t, r(j), alpha(j)
         write(22,"(3ES16.8)") t, r(j), a(j)*alpha(j)
         write(23,"(3ES16.8)") t, r(j), adot(j)

         write(30,"(3ES16.8)") t, r(j), (r(j)/dos)*( uno - (uno/a(j)**2) )
      end do

      ! Dejamos dos espacios para las gráficas de gnuplot.
      write(10,*)
      write(10,*)
      write(11,*)
      write(11,*)
      write(12,*)
      write(12,*)

      write(20,*)
      write(20,*)
      write(21,*)
      write(21,*)
      write(22,*)
      write(22,*)
      write(23,*)
      write(23,*)

      write(30,*)
      write(30,*)


      !===================================================
   else if( matterType == 2 ) then
      !------------------------------------------------
      ! Para Dirac.
      do j=2, Nr, savedataR
         write(10,"(3ES16.8)") t, r(j), F1(j)
         write(11,"(3ES16.8)") t, r(j), F2(j)
         write(12,"(3ES16.8)") t, r(j), F(j)

         write(20,"(3ES16.8)") t, r(j), G1(j)
         write(21,"(3ES16.8)") t, r(j), G2(j)
         write(22,"(3ES16.8)") t, r(j), G(j)

         write(30,"(3ES16.8)") t, r(j), a(j)
         write(31,"(3ES16.8)") t, r(j), alpha(j)
         write(32,"(3ES16.8)") t, r(j), a(j)*alpha(j)
         write(33,"(3ES16.8)") t, r(j), adot(j)

         write(40,"(3ES16.8)") t, r(j), (r(j)/dos)*( uno - (uno/a(j)**2) )
      end do

      ! Dejamos dos espacios para las gráficas de gnuplot.
      write(10,*)
      write(10,*)
      write(11,*)
      write(11,*)
      write(12,*)
      write(12,*)

      write(20,*)
      write(20,*)
      write(21,*)
      write(21,*)
      write(22,*)
      write(22,*)

      write(30,*)
      write(30,*)
      write(31,*)
      write(31,*)
      write(32,*)
      write(32,*)
      write(33,*)
      write(33,*)

      write(40,*)
      write(40,*)


      !===================================================
   else if( matterType == 3 ) then
      call evolveProca
   end if

end subroutine
