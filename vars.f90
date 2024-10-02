module vars

   !===============================================!
   !      En este módulo se solicitan datos        !
   !===============================================!

   implicit none

   ! Declaramos las variables que pediremos.
   integer evolveType, matterType
   integer Nt, Nr
   integer lineasT, puntosR
   integer savedataT, savedataR
   integer indMax

   real(8) campo0, w
   real(8) rmax, t, tf, courant
   real(8) dr, dt

   real(8) campoMax

   real(8) delta, deltaFWHM

   real(8) sexto, cuarto, tercio, medio
   real(8) uno, dos, tres, cuatro, cero

   character(len=10) dirname
   character(len=10) char_campoEntero, char_campoDecimal
   integer campoEntero

contains

   subroutine lector

      ! Llamamos a los arreglos.
      use arrays

      ! Definimos un entero para referirnos al 'input.par'.
      integer fid

      ! Abrimos el 'input.par'.
      open(newunit = fid, file = './input.par', action = 'read', position = 'rewind')

      !------------------------------------------------
      ! Tipo de evolución: 1 -> estática (solo c.i.); 2 -> Dinámica (evolución)
      read(fid,*) evolveType

      !------------------------------------------------
      ! Tipo de materia: 1 -> KG; 2 -> Dirac; 3 -> Proca
      read(fid,*) matterType

      !------------------------------------------------
      ! Parámetro inicial de la estrella, psi0, F0 o P0.
      read(fid,*) campo0

      !------------------------------------------------
      ! Parámetro de la estrella \omega.
      read(fid,*) w

      !------------------------------------------------
      ! Límite superior espacial.
      read(fid,*) rmax

      !------------------------------------------------
      ! Tiempo de integración final (Nt = tf*(Nr/(courant*rf))).
      read(fid,*) tf

      !------------------------------------------------
      ! Tamaño de paso en la malla radial.
      read(fid,*) dr

      !------------------------------------------------
      ! Cuántas líneas temporales se quieren guardar.
      read(fid,*) lineasT

      !------------------------------------------------
      ! Cuántos puntos en la malla radial se quieren.
      read(fid,*) puntosR

      !------------------------------------------------
      ! Parámetro de Courant.
      read(fid,*) courant

      !------------------------------------------------
      ! Amplitud de la perturbación en relación al máximo.
      read(fid,*) delta

      !------------------------------------------------
      ! Anchura de la perturbación en función a la FWHM.
      read(fid,*) deltaFWHM

      !------------------------------------------------
      ! Cerramos las puertas que abrimos.
      close(fid)

      !------------------------------------------------
      ! Definimos las variables que nos interesan.
      dt = courant*dr
      Nr = int( rmax/dr )
      Nt = int( (tf/dt) ) ! tf*( Nr/(courant*rmax) )

      savedataT = int( Nt/lineasT )
      savedataR = int( Nr/puntosR )

      !------------------------------------------------
      ! Iniciamos variables
      sexto  = 1.0D0/6.0D0
      cuarto = 1.0D0/4.0D0
      tercio = 1.0D0/3.0D0
      medio  = 1.0D0/2.0D0

      uno    = 1.0D0
      dos    = 2.0D0
      tres   = 3.0D0
      cuatro = 4.0D0
      cero   = 0.0D0

      t = cero

      dr = (rmax)/dble(Nr)
      dt = courant*dr

      indMax = 0
      campoMax = cero

      !------------------------------------------------
      ! Creamos la carpeta que guardará la información.
      write(char_campoDecimal, '(F0.3)') campo0

      campoEntero = int(campo0)
      write(char_campoEntero, '(I1)') campoEntero

      if (campoEntero == 0) then
         dirname = trim( trim(char_campoEntero) // trim(char_campoDecimal) )
      elseif (campoEntero > 0) then
         dirname = trim( char_campoDecimal )
      end if

      call execute_command_line ('mkdir -p ' // trim(dirname) )
      call execute_command_line ('cp ./input.par ' // './' // trim(dirname) // '/input.par')

   end subroutine

end module
