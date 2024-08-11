module arrays

   !===============================================!
   !       En este módulo se definen arreglos      !
   !===============================================!

   implicit none

   real(8), allocatable, dimension (:) :: r

   real(8), allocatable, dimension (:) :: m, s
   real(8), allocatable, dimension (:) :: a, alpha
   real(8), allocatable, dimension (:) :: a_p, alpha_p
   real(8), allocatable, dimension (:) :: adot

   ! Arreglos para el caso escalar.
   real(8), allocatable, dimension (:) :: phi, vphi
   real(8), allocatable, dimension (:) :: phi1, phi2
   real(8), allocatable, dimension (:) :: phi1_f, phi2_f
   real(8), allocatable, dimension (:) :: phi1_p, phi2_p

   real(8), allocatable, dimension (:) :: pi1, pi2, psi1, psi2
   real(8), allocatable, dimension (:) :: pi1_f, pi2_f, psi1_f, psi2_f
   real(8), allocatable, dimension (:) :: pi1_p, pi2_p, psi1_p, psi2_p

   real(8), allocatable, dimension (:) :: deltaEscalar

   ! Arreglos para el caso de Dirac.
   real(8), allocatable, dimension (:) :: F, G
   real(8), allocatable, dimension (:) :: F1, F2, G1, G2
   real(8), allocatable, dimension (:) :: F1_f, F2_f, G1_f, G2_f
   real(8), allocatable, dimension (:) :: F1_p, F2_p, G1_p, G2_P
   real(8), allocatable, dimension (:) :: deltaDirac

end module
