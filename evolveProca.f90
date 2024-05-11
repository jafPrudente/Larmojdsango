subroutine evolveProca

   !===============================================!
   !      Esta subrutina evoluciona el Sistema     !
   !               de Proca con ICN                !
   !===============================================!

   !------------------------------------------------
   ! Usamos el módulo 'arrays' para declarar arreglos
   ! y el módulo 'vars' para las variables de entrada.
   use vars
   use arrays

   implicit none

   !------------------------------------------------
   ! Definimos variables importantes.
   integer j, k, n

   real(8) dtw

   !------------------------------------------------
   ! Alojamos en memoria los arreglos.
   allocate( F(Nr), G(Nr) )
   allocate( F1(Nr), F2(Nr), G1(Nr), G2(Nr) )

   !------------------------------------------------
   ! Abriremos archivos.

   !------------------------------------------------
   ! Llamamos a las condiciones iniciales.
   call initialProca
   adot = cero

   !------------------------------------------------
   ! Si solo queremos la condición inicial ya la tendríamos
   ! y finalizamos esta subrutina con el 'return'.
   if (evolveType == 1) then
      return
   end if

   !------------------------------------------------
   ! Imprimimos información inicial a pantalla.


   !------------------------------------------------
   ! Guardamos información al tiempo cero.


   !================================================
   ! Comenzamos con el ciclo principal de evolución.


   !================================================
   ! Finaliza el ciclo principal de evolución.

   !------------------------------------------------
   ! Cerramos las puertas que abrimos.

end subroutine
