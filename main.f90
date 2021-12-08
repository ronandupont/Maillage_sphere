!===========
program main 
!===========
!
  use deftype
  use maillage
  use maillage_sphere

!
        implicit none
!
  integer                     :: nb_pts,i,j,n_s,n_f
  real(kind=kind(0.d0)), dimension(:,:),allocatable :: S,Sb,Se
  real ( kind = 8)            :: R
  integer, dimension (:,:), allocatable :: E
 

!=================================================================
!======= Nombre de points et rayon souhaités sur la sphère =======
!=================================================================

write(*,"('Entrez le nombre de points minimum désirés:')"); read(*,*) nb_pts
R=1.d0

CALL sphereIco(nb_pts,R,E,S,n_s,n_f)
print*,"Nombre de points de la sphère:",n_s

!============================================================
!======= Ecriture du tableau des sommest et des faces =======
!============================================================

OPEN ( UNIT =48 , FILE = 'S.dat')
DO i =1,n_s
	WRITE (48 ,*) ( S (i , j ) , j =1 , 3 )
ENDDO
CLOSE (48)


OPEN ( UNIT =48 , FILE = 'E.dat')
DO i =1,n_f
	WRITE (48 ,*) ( E (i , j ) , j =1 , 3 )
ENDDO
CLOSE (48)


!===============
end program main
!===============
