subroutine evolveDirac

   !===============================================!
   !      Esta subrutina evoluciona el Sistema     !
   !               de Dirac con ICN                !
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
   real(8) Srr, Srr_p

   !------------------------------------------------
   ! Alojamos en memoria los arreglos.
   allocate( F(Nr), G(Nr) )
   allocate( F1(Nr), F2(Nr), G1(Nr), G2(Nr) )
   allocate( F1_f(Nr), F2_f(Nr), G1_f(Nr), G2_f(Nr) )
   allocate( F1_p(Nr), F2_p(Nr), G1_p(Nr), G2_p(Nr) )
   allocate( deltaDirac(Nr) )

   !------------------------------------------------
   ! Llamamos a las condiciones iniciales.
   call initialDirac

   !------------------------------------------------
   ! Calculamos adot en dadas las condiciones iniciales.
   Srr = (dos / (r(j)**2 * a(j)) )*( &
      + F1(j)*(F2(j+1) - F2(j-1)) &
      - F2(j)*(F1(j+1) - F1(j-1)) &
      + G1(j)*(G2(j+1) - G2(j-1)) &
      - G2(j)*(G1(j+1) - G1(j-1)) )*medio/dr

   adot(j) = r(j)*a(j)*alpha(j)*Srr

   !------------------------------------------------
   ! Si solo queremos la condición inicial ya la tendríamos
   ! y finalizamos esta subrutina con el 'return'.
   if (evolveType == 1) then
      return
   end if

   !------------------------------------------------
   ! Calculamos la perturbación y calculamos el valor de las funciones métricas ya perturbadas.
   write(*,*) '--------'
   write(*,*) 'm_i --> ', (rmax/dos)*( uno - (uno/a(Nr)**2) )
   call perturbacionDirac
   F1 = F1 + deltaDirac
   call metricDirac
   write(*,*) 'm_p --> ', (rmax/dos)*( uno - (uno/a(Nr)**2) )

   !------------------------------------------------
   ! Imprimimos información inicial a pantalla.
   write(*,*) 'Comenzamos a evolucionar...'
   print *,'----------------------------'
   print *,'|  Iteracion  | Tiempo Fis |'
   print *,'----------------------------'

   write(*,"(A4,I7,A6,ES9.2,A3)") ' |  ',0,'    | ',t,'  |'

   !------------------------------------------------
   ! Guardamos información al tiempo cero.

   ! Abrimos los archivos.
   open(10, file = './' // trim(dirname) // '/F1.rt')
   open(11, file = './' // trim(dirname) // '/F2.rt')
   open(12, file = './' // trim(dirname) // '/F.rt')
   open(13, file = './' // trim(dirname) // '/F1.t')
   open(14, file = './' // trim(dirname) // '/F2.t')
   open(78, file = './' // trim(dirname) // '/frec.t')

   open(20, file = './' // trim(dirname) // '/G1.rt')
   open(21, file = './' // trim(dirname) // '/G2.rt')
   open(22, file = './' // trim(dirname) // '/G.rt')
   open(23, file = './' // trim(dirname) // '/G1.t')
   open(24, file = './' // trim(dirname) // '/G2.t')

   open(30, file = './' // trim(dirname) // '/a.rt')
   open(31, file = './' // trim(dirname) // '/alpha.rt')
   open(32, file = './' // trim(dirname) // '/metric.rt')
   open(33, file = './' // trim(dirname) // '/adot.rt')
   open(35, file = './' // trim(dirname) // '/error.t')

   open(40, file = './' // trim(dirname) // '/masa.rt')
   open(41, file = './' // trim(dirname) // '/masa.t')

   ! Escribimos la información en sí.
   call wardaditoStatic
   call wardaditoDinamic


   !================================================
   ! Comenzamos con el ciclo principal de evolución.
   do n=1, Nt

      !------------------------------------------------
      ! Imprimimos información en pantalla.
      if (mod(n,savedataT) == 0) then
         write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',n,'    | ',t,'  |'
      end if

      !------------------------------------------------
      ! Avanzamos en el tiempo.
      t = t + dt

      !------------------------------------------------
      ! Guardamos información previa.
      do j = 1, Nr
         F1_p(j) = F1(j)
         F2_p(j) = F2(j)

         G1_p(j) = G1(j)
         G2_P(j) = G2(j)

         a_p(j)  = a(j)
         alpha_p(j) = alpha(j)
      end do

      !================================================
      ! Iniciamos el ciclo interno del ICN.
      do k=1, 3

         !------------------------------------------------
         ! Calculamos el tamaño de paso temporal en función de k.
         if (k < 3) then
            dtw = medio*dt
         else
            dtw = dt
         end if

         !------------------------------------------------
         ! Calculamos los lados derechos.
         call fuenteDirac

         !------------------------------------------------
         ! Actualizamos las funciones.
         do j=2, Nr-1
            F1(j) = F1_p(j) + dtw*F1_f(j)
            F2(j) = F2_p(j) + dtw*F2_f(j)

            G1(j) = G1_p(j) + dtw*G1_f(j)
            G2(j) = G2_p(j) + dtw*G2_f(j)
         end do

         !------------------------------------------------
         ! Definimos condiciones de frontera.

         ! Condiciones en el origen (simetrías).
         F1(1) = -F1(2)
         F2(1) = -F2(2)

         G1(1) = G1(2)
         G2(1) = G2(2)

         ! Condiciones en la frontera exterior.
         F1(Nr) = F1_p(Nr-1) + ((dr - dtw)/(dr + dtw))*( F1_p(Nr) - F1(Nr-1) )
         F2(Nr) = F2_p(Nr-1) + ((dr - dtw)/(dr + dtw))*( F2_p(Nr) - F2(Nr-1) )

         G1(Nr) = G1_p(Nr-1) + ((dr - dtw)/(dr + dtw))*( G1_p(Nr) - G1(Nr-1) )
         G2(Nr) = G2_p(Nr-1) + ((dr - dtw)/(dr + dtw))*( G2_p(Nr) - G2(Nr-1) )

         !------------------------------------------------
         ! Evolucionamos la métrica.
         call metricDirac

         !================================================
         ! Termina el ciclo interno del ICN.
      end do

      !------------------------------------------------
      ! Calculamos adot.
      adot(2) = cero
      adot(1) = cero

      do j=3, Nr-1
         Srr = (dos / (r(j)**2 * a(j)) )*( &
            + F1(j)*(F2(j+1) - F2(j-1)) &
            - F2(j)*(F1(j+1) - F1(j-1)) &
            + G1(j)*(G2(j+1) - G2(j-1)) &
            - G2(j)*(G1(j+1) - G1(j-1)) )*medio/dr

         Srr_p = (dos / (r(j)**2 * a_p(j)) )*( &
            + F1_p(j)*(F2_p(j+1) - F2_p(j-1)) &
            - F2_p(j)*(F1_p(j+1) - F1_p(j-1)) &
            + G1_p(j)*(G2_p(j+1) - G2_p(j-1)) &
            - G2_p(j)*(G1_p(j+1) - G1_p(j-1)) )*medio/dr

         adot(j) = ( a(j) - a_p(j) )/dt &
            + r(j)*medio*( a(j)*alpha(j)*Srr + a_p(j)*alpha_p(j)*Srr_p )
      end do

      adot(Nr) = cero

      !------------------------------------------------
      ! Guardamos información.

      ! Primero los arreglos puramente temporales.
      if ( mod(n,savedataT) == 0 ) then
         call wardaditoStatic
      end if

      ! Ahora los arreglos que dependen de 'r' y 't'.
      if ( mod(n,savedataT) == 0 ) then
         call wardaditoDinamic
      end if

      !================================================
      ! Finaliza el ciclo principal de evolución.
   end do
   write(*,*) 'Hemos terminado, Camarada UwUr'

   !------------------------------------------------
   ! Cerramos las puertas que abrimos.
   close(10)
   close(11)
   close(12)
   close(13)
   close(14)
   close(78)

   close(20)
   close(21)
   close(22)
   close(23)
   close(24)

   close(30)
   close(31)
   close(32)
   close(33)
   close(35)

   close(40)
   close(41)

end subroutine
