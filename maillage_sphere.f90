module maillage_sphere


contains

subroutine sphereOcta(nb_pts,R,Ef,Sf,n_s,n_f)

!=========================== DOCUMENTATION ============================= 
! 
!    AUTEUR 
! CROGUENNEC/DUPONT
!    DATE DE CREATION 
! Octobre 2020
!
!    ENTREES 
!  nb_pts  : nombre de points souhaité sur la sphère
!  R       : rayon de la sphère
!
!    SORTIES 
!   Ef        : Talbeau des faces (connectivité des 3 numéros de sommets associés)
!   Sf        : Tableau des sommets
!   n_s       : Nombre de sommets (taille de S)
!   n_f       : Nombre de faces (taille de E)
!
!
!     UTILISATION
!    Pour obtenir les coordonnées des sommets de la face i:
!     S(Ef(i,1),1) coord en x sur le sommet 1
!     S(Ef(i,1),2) coord en y sur le sommet 1
!     S(Ef(i,1),3) coord en z sur le sommet 1

!=========================== DEBUT DES DECLARATIONS ==================== 
  use deftype
  use maillage

!
        implicit none
!
  integer                     :: N,i, j, k,l,nimp,nk,taille_S,taille_S_bis,a1,a2,a3,a12,a13,a23,a4,a5,a6,n_f,n_s,s1,s2,s3,nb_pts
  real ( kind = 8)            :: x12,x13,x23,R
  real ( kind = 8)            :: y12,y13,y23
  real ( kind = 8)            :: z12,z13,z23
  real ( kind = 8)            :: nrm12,nrm13,nrm23,nrm
  real(kind=kind(0.d0)), dimension(:,:),allocatable :: S,Sb,Se,Sf
  integer, dimension (:,:), allocatable :: E2,Ef

  Type(chaineR3), pointer :: maille
  Type(chaineR2), pointer :: R2
  Type(chaineR1), pointer :: R1

  n=int(dlog(nb_pts-2.d0)/dlog(4.d0)-1d-12) ! nombre d'itérations
  n_f=2*4**(n+1) !nombre de faces à la n-ième itération
  n_s=int(4**(n+1)+2.d0) !nombre de sommets à la n-ième itération

!     allocations
  allocate(R1)      !rangeR1
  allocate(R2)      !rangeR2
  allocate(maille)  !rangeR3
  allocate(Ef(n_f,3)) ! Tableaux des entiers renvoyant les numéro de sommets (i=face ,j=sommets)
  allocate(E2(n_f,3)) ! Tableaux 2 des entiers renvoyant les numéro de sommets (i=face ,j=sommets)
  allocate(S(0:n_f,3))
  allocate(Sb(0:n_f,3)) ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)
  allocate(Se(0:n_f,3)) ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)
 ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)
  allocate(Sf(n_s,3)) ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)
 ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)


!=========================== DEBUT DU CODE EXECUTABLE ==================


!Initialisation de S0

!S_1
Sb(0,1)=2.d0**0.5/2.d0
Sb(0,2)=2.d0**0.5/2.d0
Sb(0,3)=0.d0

!S_2
Sb(1,1)=-2.d0**0.5/2.d0
Sb(1,2)=2.d0**0.5/2.d0
Sb(1,3)=0.d0

!S_3
Sb(2,1)=2.d0**0.5/2.d0
Sb(2,2)=-2.d0**0.5/2.d0
Sb(2,3)=0.d0

!S_4
Sb(3,1)=-2.d0**0.5/2.d0
Sb(3,2)=-2.d0**0.5/2.d0
Sb(3,3)=0.d0

!S_5
Sb(4,1)=0.d0
Sb(4,2)=0.d0
Sb(4,3)=1.d0

!S_6
Sb(5,1)=0.d0
Sb(5,2)=0.d0
Sb(5,3)=-1.d0


! Norme du rayon de la sphère
nrm=(Sb(1,1)**2+Sb(1,2)**2+Sb(1,3)**2)**0.5


! Calcul des points avec un rayon R
Sb(0,1)=Sb(0,1)/nrm*R
Sb(0,2)=Sb(0,2)/nrm*R
Sb(0,3)=Sb(0,3)/nrm*R

Sb(1,1)=Sb(1,1)/nrm*R
Sb(1,2)=Sb(1,2)/nrm*R
Sb(1,3)=Sb(1,3)/nrm*R

Sb(2,1)=Sb(2,1)/nrm*R
Sb(2,2)=Sb(2,2)/nrm*R
Sb(2,3)=Sb(2,3)/nrm*R

Sb(3,1)=Sb(3,1)/nrm*R
Sb(3,2)=Sb(3,2)/nrm*R
Sb(3,3)=Sb(3,3)/nrm*R

Sb(4,1)=Sb(4,1)/nrm*R
Sb(4,2)=Sb(4,2)/nrm*R
Sb(4,3)=Sb(4,3)/nrm*R

Sb(5,1)=Sb(5,1)/nrm*R
Sb(5,2)=Sb(5,2)/nrm*R
Sb(5,3)=Sb(5,3)/nrm*R

nrm=R

!Initialisation de Ef

!Face 1
Ef(1,1)=0
Ef(1,2)=1
Ef(1,3)=2

!Face 2
Ef(2,1)=1
Ef(2,2)=2
Ef(2,3)=5

!Face 3
Ef(3,1)=5
Ef(3,2)=2
Ef(3,3)=4

!Face 4
Ef(4,1)=0
Ef(4,2)=2
Ef(4,3)=4

!Face 5
Ef(5,1)=0
Ef(5,2)=1
Ef(5,3)=3

!Face 6
Ef(6,1)=1
Ef(6,2)=5
Ef(6,3)=3

!Face 7
Ef(7,1)=3
Ef(7,2)=4
Ef(7,3)=5

!Face 8
Ef(8,1)=0
Ef(8,2)=3
Ef(8,3)=4


!Initialisation des sommets initiaux dans rangeR

!Sommet 1
!R1
R1%z=Sb(0,3)
!R1%num=1

!R2
R2%y=Sb(0,2)
R2%cr1=>R1

!R3
maille%x=Sb(0,1)
maille%cr2=>R2

taille_S=1 !1er sommet rangé

! Ranger les 3 autres sommets
CALL rangeR3(Sb(1,1),Sb(1,2),Sb(1,3),maille,taille_S) !Sommet 2
CALL rangeR3(Sb(2,1),Sb(2,2),Sb(2,3),maille,taille_S) !Sommet 3
CALL rangeR3(Sb(3,1),Sb(3,2),Sb(3,3),maille,taille_S) !Sommet 4
CALL rangeR3(Sb(4,1),Sb(4,2),Sb(4,3),maille,taille_S) !Sommet 5
CALL rangeR3(Sb(5,1),Sb(5,2),Sb(5,3),maille,taille_S) !Sommet 6

CALL parcoursR3(maille) !numerotation

! Trouver les nums
CALL donne_num(Sb(0,1),Sb(0,2),Sb(0,3),maille,a1,nimp)
CALL donne_num(Sb(1,1),Sb(1,2),Sb(1,3),maille,a2,nimp)
CALL donne_num(Sb(2,1),Sb(2,2),Sb(2,3),maille,a3,nimp)
CALL donne_num(Sb(3,1),Sb(3,2),Sb(3,3),maille,a4,nimp)
CALL donne_num(Sb(4,1),Sb(4,2),Sb(4,3),maille,a5,nimp)
CALL donne_num(Sb(5,1),Sb(5,2),Sb(5,3),maille,a6,nimp)

! Affectation des premiers sommets à S
S(a1,1)=Sb(0,1)
S(a1,2)=Sb(0,2)
S(a1,3)=Sb(0,3)

S(a2,1)=Sb(1,1)
S(a2,2)=Sb(1,2)
S(a2,3)=Sb(1,3)

S(a3,1)=Sb(2,1)
S(a3,2)=Sb(2,2)
S(a3,3)=Sb(2,3)

S(a4,1)=Sb(3,1)
S(a4,2)=Sb(3,2)
S(a4,3)=Sb(3,3)

S(a5,1)=Sb(4,1)
S(a5,2)=Sb(4,2)
S(a5,3)=Sb(4,3)

S(a6,1)=Sb(5,1)
S(a6,2)=Sb(5,2)
S(a6,3)=Sb(5,3)

!Corps du programme
DO i=1,n
	Se=S
	taille_S_bis=taille_S
  DO j=1,2*4**i
    s1=Ef(j,1) !Sommet 1 de face j
    s2=Ef(j,2) !Sommet 2 de face j
    s3=Ef(j,3) !Sommet 3 de face j

    !Nouveaux coord en x non normalisés
    x12=(S(s1,1)+S(s2,1))/2.d0
    x13=(S(s1,1)+S(s3,1))/2.d0
    x23=(S(s2,1)+S(s3,1))/2.d0

    !Nouveaux coord en y non normalisés
    y12=(S(s1,2)+S(s2,2))/2.d0
    y13=(S(s1,2)+S(s3,2))/2.d0
    y23=(S(s2,2)+S(s3,2))/2.d0

    !Nouveaux coord en z non normalisés
    z12=(S(s1,3)+S(s2,3))/2.d0
    z13=(S(s1,3)+S(s3,3))/2.d0
    z23=(S(s2,3)+S(s3,3))/2.d0


    !Normes

    nrm12=(x12**2+y12**2+z12**2)**0.5
    nrm13=(x13**2+y13**2+z13**2)**0.5
    nrm23=(x23**2+y23**2+z23**2)**0.5

    !Normalisation
    x12=x12*nrm/nrm12
    y12=y12*nrm/nrm12
    z12=z12*nrm/nrm12

    x13=x13*nrm/nrm13
    y13=y13*nrm/nrm13
    z13=z13*nrm/nrm13

    x23=x23*nrm/nrm23
    y23=y23*nrm/nrm23
    z23=z23*nrm/nrm23

  !Ajout à maille(R3)
  CALL rangeR3(x12,y12,z12,maille,taille_S)
  CALL rangeR3(x13,y13,z13,maille,taille_S)
  CALL rangeR3(x23,y23,z23,maille,taille_S)

  ! Affectation à S des nouveaux sommets dans la liste des anciens:
  S(taille_S_bis+3*(j-1),1)=x12
  S(taille_S_bis+3*(j-1),2)=y12
  S(taille_S_bis+3*(j-1),3)=z12

  S(taille_S_bis+3*(j-1)+1,1)=x13
  S(taille_S_bis+3*(j-1)+1,2)=y13
  S(taille_S_bis+3*(j-1)+1,3)=z13

  S(taille_S_bis+3*(j-1)+2,1)=x23
  S(taille_S_bis+3*(j-1)+2,2)=y23
  S(taille_S_bis+3*(j-1)+2,3)=z23
  ENDDO

	!Numérotation  
	CALL parcoursR3(maille) 
	
  !Ajout des sommets
	DO k=0,taille_S_bis+2*4**i*3-1
		CALL donne_num(S(k,1),S(k,2),S(k,3),maille,a1,nimp)
		Sb(a1,1)=S(k,1)
		Sb(a1,2)=S(k,2)
		Sb(a1,3)=S(k,3)
	ENDDO
	S=Sb  

	DO j=1,2*4**i
    s1=Ef(j,1) !Sommet 1 de face j
    s2=Ef(j,2) !Sommet 2 de face j
    s3=Ef(j,3) !Sommet 3 de face j

    !Nouveaux coord en x non normalisés
    x12=(Se(s1,1)+Se(s2,1))/2.d0
    x13=(Se(s1,1)+Se(s3,1))/2.d0
    x23=(Se(s2,1)+Se(s3,1))/2.d0

    !Nouveaux coord en y non normalisés
    y12=(Se(s1,2)+Se(s2,2))/2.d0
    y13=(Se(s1,2)+Se(s3,2))/2.d0
    y23=(Se(s2,2)+Se(s3,2))/2.d0

    !Nouveaux coord en z non normalisés
    z12=(Se(s1,3)+Se(s2,3))/2.d0
    z13=(Se(s1,3)+Se(s3,3))/2.d0
    z23=(Se(s2,3)+Se(s3,3))/2.d0


    !Normes
    nrm12=(x12**2+y12**2+z12**2)**0.5
    nrm13=(x13**2+y13**2+z13**2)**0.5
    nrm23=(x23**2+y23**2+z23**2)**0.5

    !Normalisation
    x12=x12*nrm/nrm12
    y12=y12*nrm/nrm12
    z12=z12*nrm/nrm12

    x13=x13*nrm/nrm13
    y13=y13*nrm/nrm13
    z13=z13*nrm/nrm13

    x23=x23*nrm/nrm23
    y23=y23*nrm/nrm23
    z23=z23*nrm/nrm23

    ! Numérotations centre d'arretes		  
    CALL donne_num(x12,y12,z12,maille,a12,nimp)
    CALL donne_num(x13,y13,z13,maille,a13,nimp)
    CALL donne_num(x23,y23,z23,maille,a23,nimp)

    ! Numérotations anciens sommets		  
    CALL donne_num(Se(s1,1),Se(s1,2),Se(s1,3),maille,a1,nimp)
    CALL donne_num(Se(s2,1),Se(s2,2),Se(s2,3),maille,a2,nimp)
    CALL donne_num(Se(s3,1),Se(s3,2),Se(s3,3),maille,a3,nimp)

    !Face 1
    E2(1+4*(j-1),1)=a1
    E2(1+4*(j-1),2)=a12
    E2(1+4*(j-1),3)=a13

    !Face 2
    E2(2+4*(j-1),1)=a12
    E2(2+4*(j-1),2)=a23
    E2(2+4*(j-1),3)=a13

    !Face 3
    E2(3+4*(j-1),1)=a13
    E2(3+4*(j-1),2)=a3
    E2(3+4*(j-1),3)=a23

    !Face 4
    E2(4+4*(j-1),1)=a12
    E2(4+4*(j-1),2)=a23
    E2(4+4*(j-1),3)=a2
	ENDDO
	Ef=E2
ENDDO
Sf(1:taille_S,:)=S(0:taille_S-1,:)

end subroutine sphereOcta

subroutine sphereTetra(nb_pts,R,Ef,Sf,n_s,n_f)

!=========================== DOCUMENTATION ============================= 
! 
!    AUTEUR 
! CROGUENNEC/DUPONT
!    DATE DE CREATION 
! Octobre 2020
!
!    ENTREES 
!  nb_pts  : nombre de points souhaité sur la sphère
!  R       : rayon de la sphère
!
!    SORTIES 
!   Ef        : Talbeau des faces (connectivité des 3 numéros de sommets associés)
!   Sf        : Tableau des sommets
!   n_s       : Nombre de sommets (taille de S)
!   n_f       : Nombre de faces (taille de E)
!
!
!     UTILISATION
!    Pour obtenir les coordonnées des sommets de la face i:
!     S(Ef(i,1),1) coord en x sur le sommet 1
!     S(Ef(i,1),2) coord en y sur le sommet 1
!     S(Ef(i,1),3) coord en z sur le sommet 1

!=========================== DEBUT DES DECLARATIONS ==================== 
  use deftype
  use maillage

!
        implicit none
!
  integer                     :: N,i, j, k,l,nimp,nk,taille_S,taille_S_bis,a1,a2,a3,a12,a13,a23,a4,n_f,n_s,s1,s2,s3,nb_pts,taille_E
  real ( kind = 8)            :: x12,x13,x23,R
  real ( kind = 8)            :: y12,y13,y23
  real ( kind = 8)            :: z12,z13,z23
  real ( kind = 8)            :: nrm12,nrm13,nrm23,nrm
  real(kind=kind(0.d0)), dimension(:,:),allocatable :: S,Sb,Se,Sf
  integer, dimension (:,:), allocatable :: E2,Ef

  Type(chaineR3), pointer :: maille
  Type(chaineR2), pointer :: R2
  Type(chaineR1), pointer :: R1

  n=int((dlog((nb_pts-2.5)*2.d0+1.d0)/dlog(4.d0))-1d-12) ! nombre d'itérations
  n_f=4**(n+1) !nombre de faces à la n-ième itération
  n_s=int(2.5+(4**(n+1)-1.d0)/2.d0) !nombre de sommets à la n-ième itération
  taille_E=n_f

!     allocations
  allocate(R1)      !rangeR1
  allocate(R2)      !rangeR2
  allocate(maille)  !rangeR3
  allocate(Ef(n_f,3)) ! Tableaux des entiers renvoyant les numéro de sommets (i=face ,j=sommets)
  allocate(E2(n_f,3)) ! Tableaux 2 des entiers renvoyant les numéro de sommets (i=face ,j=sommets)
  allocate(S(0:n_f,3))
  allocate(Sb(0:n_f,3)) ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)
  allocate(Se(0:n_f,3)) ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)
 ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)
  allocate(Sf(n_s,3)) ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)
 ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)


!=========================== DEBUT DU CODE EXECUTABLE ==================


!S_1
Sb(0,1)=-0.5
Sb(0,2)=-((0.75)**0.5)/3.d0
Sb(0,3)=-((2.d0/3.d0)**0.5)/4.d0

!S_2
Sb(1,1)=0.5
Sb(1,2)=-((0.75)**0.5)/3.d0
Sb(1,3)=-((2.d0/3.d0)**0.5)/4.d0

!S_3
Sb(2,1)=0
Sb(2,2)=2*((0.75)**0.5)/3.d0
Sb(2,3)=-((2.d0/3.d0)**0.5)/4.d0

!S_4
Sb(3,1)=0.d0
Sb(3,2)=0.d0
Sb(3,3)=3.d0*((2.d0/3.d0)**0.5)/4.d0

! Norme du rayon de la sphère
nrm=(Sb(1,1)**2+Sb(1,2)**2+Sb(1,3)**2)**0.5


! Calcul des points avec un rayon R
Sb(0,1)=Sb(0,1)/nrm*R
Sb(0,2)=Sb(0,2)/nrm*R
Sb(0,3)=Sb(0,3)/nrm*R

Sb(1,1)=Sb(1,1)/nrm*R
Sb(1,2)=Sb(1,2)/nrm*R
Sb(1,3)=Sb(1,3)/nrm*R

Sb(2,1)=Sb(2,1)/nrm*R
Sb(2,2)=Sb(2,2)/nrm*R
Sb(2,3)=Sb(2,3)/nrm*R

Sb(3,1)=Sb(3,1)/nrm*R
Sb(3,2)=Sb(3,2)/nrm*R
Sb(3,3)=Sb(3,3)/nrm*R

nrm=R

!Initialisation de Ef

!Face 1
Ef(1,1)=1
Ef(1,2)=2
Ef(1,3)=3

!Face 2
Ef(2,1)=1
Ef(2,2)=2
Ef(2,3)=0

!Face 3
Ef(3,1)=1
Ef(3,2)=3
Ef(3,3)=0

!Face 4
Ef(4,1)=2
Ef(4,2)=3
Ef(4,3)=0


!Initialisation des sommets initiaux dans rangeR

!Sommet 1
!R1
R1%z=Sb(0,3)
!R1%num=1

!R2
R2%y=Sb(0,2)
R2%cr1=>R1

!R3
maille%x=Sb(0,1)
maille%cr2=>R2

taille_S=1 !1er sommet rangé

! Ranger les 3 autres sommets
CALL rangeR3(Sb(1,1),Sb(1,2),Sb(1,3),maille,taille_S) !Sommet 2
CALL rangeR3(Sb(2,1),Sb(2,2),Sb(2,3),maille,taille_S) !Sommet 3
CALL rangeR3(Sb(3,1),Sb(3,2),Sb(3,3),maille,taille_S) !Sommet 4
CALL rangeR3(Sb(3,1),Sb(3,2),Sb(3,3),maille,taille_S) !Essai de remettre le sommet 4 en doublon

CALL parcoursR3(maille) !numerotation

! Trouver les nums
CALL donne_num(Sb(0,1),Sb(0,2),Sb(0,3),maille,a1,nimp)
CALL donne_num(Sb(1,1),Sb(1,2),Sb(1,3),maille,a2,nimp)
CALL donne_num(Sb(2,1),Sb(2,2),Sb(2,3),maille,a3,nimp)
CALL donne_num(Sb(3,1),Sb(3,2),Sb(3,3),maille,a4,nimp)

! Affectation des premiers sommets à S
S(a1,1)=Sb(0,1)
S(a1,2)=Sb(0,2)
S(a1,3)=Sb(0,3)

S(a2,1)=Sb(1,1)
S(a2,2)=Sb(1,2)
S(a2,3)=Sb(1,3)

S(a3,1)=Sb(2,1)
S(a3,2)=Sb(2,2)
S(a3,3)=Sb(2,3)

S(a4,1)=Sb(3,1)
S(a4,2)=Sb(3,2)
S(a4,3)=Sb(3,3)

!Corps du programme
DO i=1,n
	Se=S
	taille_S_bis=taille_S
  DO j=1,4**i
    s1=Ef(j,1) !Sommet 1 de face j
    s2=Ef(j,2) !Sommet 2 de face j
    s3=Ef(j,3) !Sommet 3 de face j

    !Nouveaux coord en x non normalisés
    x12=(S(s1,1)+S(s2,1))/2.d0
    x13=(S(s1,1)+S(s3,1))/2.d0
    x23=(S(s2,1)+S(s3,1))/2.d0

    !Nouveaux coord en y non normalisés
    y12=(S(s1,2)+S(s2,2))/2.d0
    y13=(S(s1,2)+S(s3,2))/2.d0
    y23=(S(s2,2)+S(s3,2))/2.d0

    !Nouveaux coord en z non normalisés
    z12=(S(s1,3)+S(s2,3))/2.d0
    z13=(S(s1,3)+S(s3,3))/2.d0
    z23=(S(s2,3)+S(s3,3))/2.d0


    !Normes

    nrm12=(x12**2+y12**2+z12**2)**0.5
    nrm13=(x13**2+y13**2+z13**2)**0.5
    nrm23=(x23**2+y23**2+z23**2)**0.5

    !Normalisation
    x12=x12*nrm/nrm12
    y12=y12*nrm/nrm12
    z12=z12*nrm/nrm12

    x13=x13*nrm/nrm13
    y13=y13*nrm/nrm13
    z13=z13*nrm/nrm13

    x23=x23*nrm/nrm23
    y23=y23*nrm/nrm23
    z23=z23*nrm/nrm23

  !Ajout à maille(R3)
  CALL rangeR3(x12,y12,z12,maille,taille_S)
  CALL rangeR3(x13,y13,z13,maille,taille_S)
  CALL rangeR3(x23,y23,z23,maille,taille_S)

  ! Affectation à S des nouveaux sommets dans la liste des anciens:
  S(taille_S_bis+3*(j-1),1)=x12
  S(taille_S_bis+3*(j-1),2)=y12
  S(taille_S_bis+3*(j-1),3)=z12

  S(taille_S_bis+3*(j-1)+1,1)=x13
  S(taille_S_bis+3*(j-1)+1,2)=y13
  S(taille_S_bis+3*(j-1)+1,3)=z13

  S(taille_S_bis+3*(j-1)+2,1)=x23
  S(taille_S_bis+3*(j-1)+2,2)=y23
  S(taille_S_bis+3*(j-1)+2,3)=z23
  ENDDO

	!Numérotation  
	CALL parcoursR3(maille) 
	
  !Ajout des sommets
	DO k=0,taille_S_bis+4**i*3-1
		CALL donne_num(S(k,1),S(k,2),S(k,3),maille,a1,nimp)
		Sb(a1,1)=S(k,1)
		Sb(a1,2)=S(k,2)
		Sb(a1,3)=S(k,3)
	ENDDO
	S=Sb  

	DO j=1,4**i
    s1=Ef(j,1) !Sommet 1 de face j
    s2=Ef(j,2) !Sommet 2 de face j
    s3=Ef(j,3) !Sommet 3 de face j

    !Nouveaux coord en x non normalisés
    x12=(Se(s1,1)+Se(s2,1))/2.d0
    x13=(Se(s1,1)+Se(s3,1))/2.d0
    x23=(Se(s2,1)+Se(s3,1))/2.d0

    !Nouveaux coord en y non normalisés
    y12=(Se(s1,2)+Se(s2,2))/2.d0
    y13=(Se(s1,2)+Se(s3,2))/2.d0
    y23=(Se(s2,2)+Se(s3,2))/2.d0

    !Nouveaux coord en z non normalisés
    z12=(Se(s1,3)+Se(s2,3))/2.d0
    z13=(Se(s1,3)+Se(s3,3))/2.d0
    z23=(Se(s2,3)+Se(s3,3))/2.d0


    !Normes
    nrm12=(x12**2+y12**2+z12**2)**0.5
    nrm13=(x13**2+y13**2+z13**2)**0.5
    nrm23=(x23**2+y23**2+z23**2)**0.5

    !Normalisation
    x12=x12*nrm/nrm12
    y12=y12*nrm/nrm12
    z12=z12*nrm/nrm12

    x13=x13*nrm/nrm13
    y13=y13*nrm/nrm13
    z13=z13*nrm/nrm13

    x23=x23*nrm/nrm23
    y23=y23*nrm/nrm23
    z23=z23*nrm/nrm23

    ! Numérotations centre d'arretes		  
    CALL donne_num(x12,y12,z12,maille,a12,nimp)
    CALL donne_num(x13,y13,z13,maille,a13,nimp)
    CALL donne_num(x23,y23,z23,maille,a23,nimp)

    ! Numérotations anciens sommets		  
    CALL donne_num(Se(s1,1),Se(s1,2),Se(s1,3),maille,a1,nimp)
    CALL donne_num(Se(s2,1),Se(s2,2),Se(s2,3),maille,a2,nimp)
    CALL donne_num(Se(s3,1),Se(s3,2),Se(s3,3),maille,a3,nimp)

    !Face 1
    E2(1+4*(j-1),1)=a1
    E2(1+4*(j-1),2)=a12
    E2(1+4*(j-1),3)=a13

    !Face 2
    E2(2+4*(j-1),1)=a12
    E2(2+4*(j-1),2)=a23
    E2(2+4*(j-1),3)=a13

    !Face 3
    E2(3+4*(j-1),1)=a13
    E2(3+4*(j-1),2)=a3
    E2(3+4*(j-1),3)=a23

    !Face 4
    E2(4+4*(j-1),1)=a12
    E2(4+4*(j-1),2)=a23
    E2(4+4*(j-1),3)=a2
	ENDDO
	Ef=E2
ENDDO
Sf(1:taille_S,:)=S(0:taille_S-1,:)

end subroutine sphereTetra


subroutine sphereIco(nb_pts,R,Ef,Sf,n_s,n_f)

!=========================== DOCUMENTATION ============================= 
! 
!    AUTEUR 
! CROGUENNEC/DUPONT
!    DATE DE CREATION 
! Octobre 2020
!
!    ENTREES 
!  nb_pts  : nombre de points souhaité sur la sphère
!  R       : rayon de la sphère
!
!    SORTIES 
!   Ef        : Talbeau des faces (connectivité des 3 numéros de sommets associés)
!   Sf        : Tableau des sommets
!   n_s       : Nombre de sommets (taille de S)
!   n_f       : Nombre de faces (taille de E)
!
!
!     UTILISATION
!    Pour obtenir les coordonnées des sommets de la face i:
!     S(Ef(i,1),1) coord en x sur le sommet 1
!     S(Ef(i,1),2) coord en y sur le sommet 1
!     S(Ef(i,1),3) coord en z sur le sommet 1

!=========================== DEBUT DES DECLARATIONS ==================== 
  use deftype
  use maillage

!
        implicit none
!
  integer                     :: N,i, j, k,l,nimp,nk,taille_S,taille_S_bis,a1,a2,a3,a12,a13,a23,a4,a5,a6,a7,a8,a9,a10,a11,a12_,n_f
  integer                     :: n_s,s1,s2,s3,nb_pts
  real ( kind = 8)            :: x12,x13,x23,R,phi
  real ( kind = 8)            :: y12,y13,y23
  real ( kind = 8)            :: z12,z13,z23
  real ( kind = 8)            :: nrm12,nrm13,nrm23,nrm
  real(kind=kind(0.d0)), dimension(:,:),allocatable :: S,Sb,Se,Sf
  integer, dimension (:,:), allocatable :: E2,Ef

  Type(chaineR3), pointer :: maille
  Type(chaineR2), pointer :: R2
  Type(chaineR1), pointer :: R1

  n=int((dlog((nb_pts-2.d0)/10.d0)/dlog(4.d0))+1.d0-1d-12) ! nombre d'itérations
  n_f=5*4**(n+1) !nombre de faces à la n-ième itération
  n_s=int(10*4**n+2.d0) !nombre de sommets à la n-ième itération
  phi=(1.d0+5**0.5)/2.d0
!     allocations
  allocate(R1)      !rangeR1
  allocate(R2)      !rangeR2
  allocate(maille)  !rangeR3
  allocate(Ef(n_f,3)) ! Tableaux des entiers renvoyant les numéro de sommets (i=face ,j=sommets)
  allocate(E2(n_f,3)) ! Tableaux 2 des entiers renvoyant les numéro de sommets (i=face ,j=sommets)
  allocate(S(0:n_f,3))
  allocate(Sb(0:n_f,3)) ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)
  allocate(Se(0:n_f,3)) ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)
 ! Tableau des réels renvoyant position (i=num sommet, j=1,2,3 - x y z)
 


!=========================== DEBUT DU CODE EXECUTABLE ==================


!S_1
Sb(0,1)=phi
Sb(0,2)=0.d0
Sb(0,3)=1.d0

!S_2
Sb(1,1)=phi
Sb(1,2)=0.d0
Sb(1,3)=-1.d0

!S_3
Sb(2,1)=-phi
Sb(2,2)=0.d0
Sb(2,3)=1.d0

!S_4
Sb(3,1)=-phi
Sb(3,2)=0.d0
Sb(3,3)=-1.d0

!S_5
Sb(4,1)=0.d0
Sb(4,2)=1.d0
Sb(4,3)=phi

!S_6
Sb(5,1)=0.d0
Sb(5,2)=-1.d0
Sb(5,3)=phi

!S_7
Sb(6,1)=0.d0
Sb(6,2)=1.d0
Sb(6,3)=-phi

!S_8
Sb(7,1)=0.d0
Sb(7,2)=-1.d0
Sb(7,3)=-phi

!S_9
Sb(8,1)=1.d0
Sb(8,2)=phi
Sb(8,3)=0.d0

!S_10
Sb(9,1)=-1.d0
Sb(9,2)=phi
Sb(9,3)=0.d0

!S_11
Sb(10,1)=1.d0
Sb(10,2)=-phi
Sb(10,3)=0.d0

!S_12
Sb(11,1)=-1.d0
Sb(11,2)=-phi
Sb(11,3)=0.d0

! Norme du rayon de la sphère
nrm=(Sb(1,1)**2+Sb(1,2)**2+Sb(1,3)**2)**0.5

! Calcul des points avec un rayon R
Sb(0,1)=Sb(0,1)/nrm*R
Sb(0,2)=Sb(0,2)/nrm*R
Sb(0,3)=Sb(0,3)/nrm*R

Sb(1,1)=Sb(1,1)/nrm*R
Sb(1,2)=Sb(1,2)/nrm*R
Sb(1,3)=Sb(1,3)/nrm*R

Sb(2,1)=Sb(2,1)/nrm*R
Sb(2,2)=Sb(2,2)/nrm*R
Sb(2,3)=Sb(2,3)/nrm*R

Sb(3,1)=Sb(3,1)/nrm*R
Sb(3,2)=Sb(3,2)/nrm*R
Sb(3,3)=Sb(3,3)/nrm*R

Sb(4,1)=Sb(4,1)/nrm*R
Sb(4,2)=Sb(4,2)/nrm*R
Sb(4,3)=Sb(4,3)/nrm*R

Sb(5,1)=Sb(5,1)/nrm*R
Sb(5,2)=Sb(5,2)/nrm*R
Sb(5,3)=Sb(5,3)/nrm*R

Sb(6,1)=Sb(6,1)/nrm*R
Sb(6,2)=Sb(6,2)/nrm*R
Sb(6,3)=Sb(6,3)/nrm*R

Sb(7,1)=Sb(7,1)/nrm*R
Sb(7,2)=Sb(7,2)/nrm*R
Sb(7,3)=Sb(7,3)/nrm*R

Sb(8,1)=Sb(8,1)/nrm*R
Sb(8,2)=Sb(8,2)/nrm*R
Sb(8,3)=Sb(8,3)/nrm*R

Sb(9,1)=Sb(9,1)/nrm*R
Sb(9,2)=Sb(9,2)/nrm*R
Sb(9,3)=Sb(9,3)/nrm*R

Sb(10,1)=Sb(10,1)/nrm*R
Sb(10,2)=Sb(10,2)/nrm*R
Sb(10,3)=Sb(10,3)/nrm*R

Sb(11,1)=Sb(11,1)/nrm*R
Sb(11,2)=Sb(11,2)/nrm*R
Sb(11,3)=Sb(11,3)/nrm*R


nrm=R

!Initialisation de Ef

!Face 1
Ef(1,1)=0
Ef(1,2)=1
Ef(1,3)=2

!Face 2
Ef(2,1)=0
Ef(2,2)=1
Ef(2,3)=3

!Face 3
Ef(3,1)=0
Ef(3,2)=2
Ef(3,3)=4

!Face 4
Ef(4,1)=0
Ef(4,2)=3
Ef(4,3)=6

!Face 5
Ef(5,1)=0
Ef(5,2)=4
Ef(5,3)=6

!Face 6
Ef(6,1)=1
Ef(6,2)=2
Ef(6,3)=5

!Face 7
Ef(7,1)=1
Ef(7,2)=3
Ef(7,3)=7

!Face 8
Ef(8,1)=1
Ef(8,2)=5
Ef(8,3)=7

!Face 9
Ef(9,1)=2
Ef(9,2)=4
Ef(9,3)=8

!Face 10
Ef(10,1)=2
Ef(10,2)=5
Ef(10,3)=8

!Face 11
Ef(11,1)=3
Ef(11,2)=6
Ef(11,3)=9

!Face 12
Ef(12,1)=3
Ef(12,2)=7
Ef(12,3)=9

!Face 13
Ef(13,1)=4
Ef(13,2)=6
Ef(13,3)=10

!Face 14
Ef(14,1)=4
Ef(14,2)=8
Ef(14,3)=10

!Face 15
Ef(15,1)=5
Ef(15,2)=7
Ef(15,3)=11

!Face 16
Ef(16,1)=5
Ef(16,2)=8
Ef(16,3)=11

!Face 17
Ef(17,1)=6
Ef(17,2)=9
Ef(17,3)=10

!Face 18
Ef(18,1)=7
Ef(18,2)=9
Ef(18,3)=11

!Face 19
Ef(19,1)=8
Ef(19,2)=10
Ef(19,3)=11

!Face 20
Ef(20,1)=9
Ef(20,2)=10
Ef(20,3)=11



! Ecriture de xk dans un fichier sous forme de matrice
OPEN ( UNIT =48 , FILE = 'Efi.dat')
DO i =1,20
	WRITE (48 ,*) ( Ef (i , j ) , j =1 , 3 )
ENDDO
CLOSE (48)

!Initialisation des sommets initiaux dans rangeR

!Sommet 1
!R1
R1%z=Sb(0,3)
!R1%num=1

!R2
R2%y=Sb(0,2)
R2%cr1=>R1

!R3
maille%x=Sb(0,1)
maille%cr2=>R2

taille_S=1 !1er sommet rangé

! Ranger les 3 autres sommets
CALL rangeR3(Sb(1,1),Sb(1,2),Sb(1,3),maille,taille_S) !Sommet 2
CALL rangeR3(Sb(2,1),Sb(2,2),Sb(2,3),maille,taille_S) !Sommet 3
CALL rangeR3(Sb(3,1),Sb(3,2),Sb(3,3),maille,taille_S) !Sommet 4
CALL rangeR3(Sb(4,1),Sb(4,2),Sb(4,3),maille,taille_S) !Sommet 5
CALL rangeR3(Sb(5,1),Sb(5,2),Sb(5,3),maille,taille_S) !Sommet 6
CALL rangeR3(Sb(6,1),Sb(6,2),Sb(6,3),maille,taille_S) !Sommet 7
CALL rangeR3(Sb(7,1),Sb(7,2),Sb(7,3),maille,taille_S) !Sommet 8
CALL rangeR3(Sb(8,1),Sb(8,2),Sb(8,3),maille,taille_S) !Sommet 9
CALL rangeR3(Sb(9,1),Sb(9,2),Sb(9,3),maille,taille_S) !Sommet 10
CALL rangeR3(Sb(10,1),Sb(10,2),Sb(10,3),maille,taille_S) !Sommet 11
CALL rangeR3(Sb(11,1),Sb(11,2),Sb(11,3),maille,taille_S) !Sommet 12


CALL parcoursR3(maille) !numerotation

! Trouver les nums
CALL donne_num(Sb(0,1),Sb(0,2),Sb(0,3),maille,a1,nimp)
CALL donne_num(Sb(1,1),Sb(1,2),Sb(1,3),maille,a2,nimp)
CALL donne_num(Sb(2,1),Sb(2,2),Sb(2,3),maille,a3,nimp)
CALL donne_num(Sb(3,1),Sb(3,2),Sb(3,3),maille,a4,nimp)
CALL donne_num(Sb(4,1),Sb(4,2),Sb(4,3),maille,a5,nimp)
CALL donne_num(Sb(5,1),Sb(5,2),Sb(5,3),maille,a6,nimp)
CALL donne_num(Sb(6,1),Sb(6,2),Sb(6,3),maille,a7,nimp)
CALL donne_num(Sb(7,1),Sb(7,2),Sb(7,3),maille,a8,nimp)
CALL donne_num(Sb(8,1),Sb(8,2),Sb(8,3),maille,a9,nimp)
CALL donne_num(Sb(9,1),Sb(9,2),Sb(9,3),maille,a10,nimp)
CALL donne_num(Sb(10,1),Sb(10,2),Sb(10,3),maille,a11,nimp)
CALL donne_num(Sb(11,1),Sb(11,2),Sb(11,3),maille,a12_,nimp)
! Affectation des premiers sommets à S
S(a1,1)=Sb(0,1)
S(a1,2)=Sb(0,2)
S(a1,3)=Sb(0,3)

S(a2,1)=Sb(1,1)
S(a2,2)=Sb(1,2)
S(a2,3)=Sb(1,3)

S(a3,1)=Sb(2,1)
S(a3,2)=Sb(2,2)
S(a3,3)=Sb(2,3)

S(a4,1)=Sb(3,1)
S(a4,2)=Sb(3,2)
S(a4,3)=Sb(3,3)

S(a5,1)=Sb(4,1)
S(a5,2)=Sb(4,2)
S(a5,3)=Sb(4,3)

S(a6,1)=Sb(5,1)
S(a6,2)=Sb(5,2)
S(a6,3)=Sb(5,3)

S(a7,1)=Sb(6,1)
S(a7,2)=Sb(6,2)
S(a7,3)=Sb(6,3)

S(a8,1)=Sb(7,1)
S(a8,2)=Sb(7,2)
S(a8,3)=Sb(7,3)

S(a9,1)=Sb(8,1)
S(a9,2)=Sb(8,2)
S(a9,3)=Sb(8,3)

S(a10,1)=Sb(9,1)
S(a10,2)=Sb(9,2)
S(a10,3)=Sb(9,3)

S(a11,1)=Sb(10,1)
S(a11,2)=Sb(10,2)
S(a11,3)=Sb(10,3)

S(a12_,1)=Sb(11,1)
S(a12_,2)=Sb(11,2)
S(a12_,3)=Sb(11,3)

! Ecriture de xk dans un fichier sous forme de matrice
OPEN ( UNIT =48 , FILE = 'Sbi.dat')
DO i =0 ,11
	WRITE (48 ,*) ( S (i , j ) , j =1 , 3 )
ENDDO
CLOSE (48)




!Corps du programme
DO i=1,n
	Se=S
	taille_S_bis=taille_S
  DO j=1,5*4**i
    s1=Ef(j,1) !Sommet 1 de face j
    s2=Ef(j,2) !Sommet 2 de face j
    s3=Ef(j,3) !Sommet 3 de face j
 !   print*,s1,s2,s3
!print*,S(s1,1),S(s1,2),S(s1,3)
!print*,S(s2,1),S(s2,2),S(s2,3)
!print*,S(s3,1),S(s3,2),S(s3,3)

    !Nouveaux coord en x non normalisés
    x12=(S(s1,1)+S(s2,1))/2.d0
    x13=(S(s1,1)+S(s3,1))/2.d0
    x23=(S(s2,1)+S(s3,1))/2.d0

    !Nouveaux coord en y non normalisés
    y12=(S(s1,2)+S(s2,2))/2.d0
    y13=(S(s1,2)+S(s3,2))/2.d0
    y23=(S(s2,2)+S(s3,2))/2.d0

    !Nouveaux coord en z non normalisés
    z12=(S(s1,3)+S(s2,3))/2.d0
    z13=(S(s1,3)+S(s3,3))/2.d0
    z23=(S(s2,3)+S(s3,3))/2.d0


    !Normes

    nrm12=(x12**2+y12**2+z12**2)**0.5
    nrm13=(x13**2+y13**2+z13**2)**0.5
    nrm23=(x23**2+y23**2+z23**2)**0.5

!  print*,nrm23
    !Normalisation
    x12=x12*nrm/nrm12
    y12=y12*nrm/nrm12
    z12=z12*nrm/nrm12

    x13=x13*nrm/nrm13
    y13=y13*nrm/nrm13
    z13=z13*nrm/nrm13

    x23=x23*nrm/nrm23
    y23=y23*nrm/nrm23
    z23=z23*nrm/nrm23

!    print*,x23,y23,z23
  !Ajout à maille(R3)
  CALL rangeR3(x12,y12,z12,maille,taille_S)
  CALL rangeR3(x13,y13,z13,maille,taille_S)
  CALL rangeR3(x23,y23,z23,maille,taille_S)

  ! Affectation à S des nouveaux sommets dans la liste des anciens:

  S(taille_S_bis+3*(j-1),1)=x12
  S(taille_S_bis+3*(j-1),2)=y12
  S(taille_S_bis+3*(j-1),3)=z12

  S(taille_S_bis+3*(j-1)+1,1)=x13
  S(taille_S_bis+3*(j-1)+1,2)=y13
  S(taille_S_bis+3*(j-1)+1,3)=z13

  S(taille_S_bis+3*(j-1)+2,1)=x23
  S(taille_S_bis+3*(j-1)+2,2)=y23
  S(taille_S_bis+3*(j-1)+2,3)=z23
  ENDDO

	!Numérotation  
	CALL parcoursR3(maille) 

  !Ajout des sommets
	DO k=0,taille_S_bis+5*4**i*3-1
		CALL donne_num(S(k,1),S(k,2),S(k,3),maille,a1,nimp)
		Sb(a1,1)=S(k,1)
		Sb(a1,2)=S(k,2)
		Sb(a1,3)=S(k,3)
	ENDDO
	S=Sb  

	DO j=1,5*4**i
    s1=Ef(j,1) !Sommet 1 de face j
    s2=Ef(j,2) !Sommet 2 de face j
    s3=Ef(j,3) !Sommet 3 de face j

    !Nouveaux coord en x non normalisés
    x12=(Se(s1,1)+Se(s2,1))/2.d0
    x13=(Se(s1,1)+Se(s3,1))/2.d0
    x23=(Se(s2,1)+Se(s3,1))/2.d0

    !Nouveaux coord en y non normalisés
    y12=(Se(s1,2)+Se(s2,2))/2.d0
    y13=(Se(s1,2)+Se(s3,2))/2.d0
    y23=(Se(s2,2)+Se(s3,2))/2.d0

    !Nouveaux coord en z non normalisés
    z12=(Se(s1,3)+Se(s2,3))/2.d0
    z13=(Se(s1,3)+Se(s3,3))/2.d0
    z23=(Se(s2,3)+Se(s3,3))/2.d0


    !Normes
    nrm12=(x12**2+y12**2+z12**2)**0.5
    nrm13=(x13**2+y13**2+z13**2)**0.5
    nrm23=(x23**2+y23**2+z23**2)**0.5

    !Normalisation
    x12=x12*nrm/nrm12
    y12=y12*nrm/nrm12
    z12=z12*nrm/nrm12

    x13=x13*nrm/nrm13
    y13=y13*nrm/nrm13
    z13=z13*nrm/nrm13

    x23=x23*nrm/nrm23
    y23=y23*nrm/nrm23
    z23=z23*nrm/nrm23

    ! Numérotations centre d'arretes		  
    CALL donne_num(x12,y12,z12,maille,a12,nimp)
    CALL donne_num(x13,y13,z13,maille,a13,nimp)
    CALL donne_num(x23,y23,z23,maille,a23,nimp)

    ! Numérotations anciens sommets		  
    CALL donne_num(Se(s1,1),Se(s1,2),Se(s1,3),maille,a1,nimp)
    CALL donne_num(Se(s2,1),Se(s2,2),Se(s2,3),maille,a2,nimp)
    CALL donne_num(Se(s3,1),Se(s3,2),Se(s3,3),maille,a3,nimp)

    !Face 1
    E2(1+4*(j-1),1)=a1
    E2(1+4*(j-1),2)=a12
    E2(1+4*(j-1),3)=a13

    !Face 2
    E2(2+4*(j-1),1)=a12
    E2(2+4*(j-1),2)=a23
    E2(2+4*(j-1),3)=a13

    !Face 3
    E2(3+4*(j-1),1)=a13
    E2(3+4*(j-1),2)=a3
    E2(3+4*(j-1),3)=a23

    !Face 4
    E2(4+4*(j-1),1)=a12
    E2(4+4*(j-1),2)=a23
    E2(4+4*(j-1),3)=a2
	ENDDO
	Ef=E2
ENDDO
n_S=taille_S
print*,n_S
allocate(Sf(n_S,3))
Sf(1:taille_S,:)=S(0:taille_S-1,:)


end subroutine sphereIco



end module maillage_sphere
