!==================
module maillage
!==================
!
contains

subroutine construct_maillage_local(locl,maille,phil,xl,yl,zl,x,y,z,phi,loc,delta,dx,dy,dz,m)

use deftype
implicit none

  
  integer, dimension (:,:,:), allocatable, intent(out) :: locl    !tableau qui renseigne les indices des voisins
  Type(chaineR3), pointer :: maille ! maillage sous forme liste chainee ordonnee tridimensionnelle coorespondant a locl
  real (kind = 8), dimension (:), allocatable, intent(out) :: phil,xl,yl,zl !valeur de phi et position sur le petit maillage construit
  real (kind = 8), dimension (:), intent(in) :: x,y,z,phi!donnee de valeur de phi et position du maillage depuis lequel on extrait le maillage (delta vosinage de phi=0)
integer, dimension (:,:,:), intent(in) :: loc!donnee du maillage par tableau d'indice du maillage total (grossier)
  real (kind = 8), intent(in) :: delta,dx,dy,dz  ! taille du voisinage et pas de la discretisation de maille
  integer, dimension (:), intent(in) :: m !niveau de raffinement  dans les trois directions a appliquer 
        
  integer :: i,il, np,i1,i2,i3
  real (kind = 8) :: m1,m2,m3,dxphi,dyphi,dzphi
  il=0
  do i=1,size(x)
!parcours des points a un delta voisinage de l interface
     if (abs(phi(i))<delta) then
        m1=0.5d0/float(m(1))
        m2=0.5d0/float(m(2))
        m3=0.5d0/float(m(3))
        dxphi=phi(i)-phi(loc(i,4,1))
        if (dxphi*phi(i)<0) dxphi=phi(loc(i,1,1))-phi(i)
        dxphi=dxphi/dx
!
        dyphi=phi(i)-phi(loc(i,5,1))
        if (dyphi*phi(i)<0) dyphi=phi(loc(i,2,1))-phi(i)
        dyphi=dyphi/dy
!
        dzphi=phi(i)-phi(loc(i,6,1))
        if (dxphi*phi(i)<0) dxphi=phi(loc(i,3,1))-phi(i)
        dzphi=dzphi/dz
        do i1=-m(1),m(1)
           do i2=-m(2),m(2)
              do i3=-m(3),m(3)
                 call rangeR3 (x(i)+i1*m1,y(i)+i2*m2,z(i)+i3*m3, maille,il)
                 maille%cr2%cr1%phi=phi(i)+dxphi*i1*m1+dyphi*i2*m2+dzphi*i3*m3 ! valeur de phi interpole
              end do
           end do
        end do
     endif
  end do
  allocate(locl(il,7,2),phil(il),xl(il),yl(il),zl(il))
  call numerotebis(maille,locl,phil,xl,yl,zl) !numerotation et construction de locl , phil, xl,yl,zl

end subroutine construct_maillage_local


subroutine correction_maillage_local(loc,taille_loc,phil,xl,yl,zl,phimax,dx,dy,dz)
! objectif: retirer de loc les points de maillages eloignes et rajouter ceux manquants du phimax voisinage
use deftype
implicit none

  
  integer, dimension (:,:,:), intent(inout) :: loc
  integer, intent(inout) :: taille_loc
  real (kind = 8), dimension (:), intent(inout) :: phil,xl,yl,zl
  real (kind = 8), intent (in) :: phimax,dx,dy,dz

        
  integer :: i,il, np, taille_m,inew, num, j,jj,interieur
  real (kind = 8) :: xi,yi,zi,valj 
  logical :: i_bord
 Type(chaineR3), pointer :: maille
    type (chainej), Pointer   :: liste_lib !indices a liberer
    type (chainej), Pointer   :: liste_add,liste_temp !liste des points du bord qui ont un voisin -au moins- a inserer et classe par ordre de phi


!intialisation maille

  Type(chaineR2), pointer :: R2
  Type(chaineR1), pointer :: R1
!

allocate(R1,R2,maille)
R1%z=0.d0
Nullify(R1%suiv) ; Nullify(R1%prec)
R2%y=0.d0
R2%cr1=>R1
nullify(R2%suiv) ; Nullify(R2%prec)
maille%x=-1.d100
maille%cr2=>R2
nullify(maille%suiv) ; Nullify(maille%prec)


! parcours du tableau loc: marquage des indices à liberer "liste_lib" dans une liste, corrections des indices voisins
    np=0
    allocate(liste_lib)
    do il=1,taille_loc
       if (abs(phil(il))>=phimax) exit
    end do
    liste_lib%j=il
           liste_lib%val=0.d0
           nullify(liste_lib%suiv)
           nullify(liste_lib%prec)
    do i=il+1,taille_loc
       if (abs(phil(i))>=phimax) then
          call ajoutval(i,0.d0,liste_lib,np)
          do j=1,3
             loc(loc(i,j,1),j+3,1)=loc(i,j,1)
          end do
          do j=4,6
             loc(loc(i,j,1),j-3,1)=loc(i,j,1)
          end do
       end if
    end do
! ajout d'indice complementaire
 do i=taille_loc,size(loc,1)
    call ajoutval(i,0.d0,liste_lib,np)
 end do
! parcours: marquage des points au bord avec voisin à creer dans une liste classee par taille de abs(phi) liste_add,
! ainsi qu une chaineR3 de ces points 

    do il=1,taille_loc
       if (abs(phil(il))<phimax) then
          do j=1,6
             i_bord=.true.
             if (loc(il,j,1)==i) exit
             i_bord=.false.
          end do
          if (i_bord) then
             call rangeR3 (xl(il),yl(il),zl(il), maille,taille_m)
             maille%cr2%cr1%num=il
             np=0
             allocate(liste_add)
             liste_add%j=il
             liste_add%val=0.d0
             nullify(liste_add%suiv)
             nullify(liste_add%prec)
             exit
          endif
       end if
    end do
    do i=il,taille_loc
       if (abs(phil(i))<phimax) then
          do j=1,6
             i_bord=.true.
             if (loc(i,j,1)==i) exit
             i_bord=.false.
          end do
          if (i_bord) then
             call rangeR3 (xl(i),yl(i),zl(i), maille,taille_m)
             maille%cr2%cr1%num=i
             maille%cr2%cr1%interieur=2
             call class_val(i,abs(phil(i)),liste_add,np)
          endif
       end if
    end do


!creation des voisins par parcours de la liste dans l'ordre et ajout à la liste et a chaineR3, avec utilisation des indices liberes de liste_lib choisis par proximite.
!ne pas oublier de corriger les voisins dans loc(:,1:6,1) en s'appuyant sur la petite chaineR3


    do while (associated(liste_add%prec))
       liste_add=>liste_add%prec
    end do
    liste_temp=>liste_add
    do while (associated(liste_temp))
       i=liste_temp%j
       do j=1,6
          if (loc(i,j,1)==i) then
             xi=xl(i)
             yi=yl(i)
             zi=zl(i)
             if (j==1) xi=xi+1.d0
             if (j==4) xi=xi-1.d0
             if (j==2) yi=yi+1.d0
             if (j==5) yi=yi-1.d0
             if (j==3) zi=zi+1.d0
             if (j==6) zi=zi-1.d0

             call rangeR3(xi,yi,zi, maille,taille_m)
             call donne_val_phi(xi,yi,zi,j,dx,dy,dz,maille,phil,valj)
             if (abs(valj)<phimax) then
                call donne_indice(i,liste_lib,inew)    

!ajouter a liste_add les eventuels nouveaux points ?????

                taille_loc=max(taille_loc,inew)
                maille%cr2%cr1%num=inew
                loc(inew,7,1)=inew
                phil(inew)=valj
                loc(i,j,1)=inew
                loc(inew,mod(j+2,6)+1,1)=i
                do jj=j+1,6
                   xi=xl(i)
                   yi=yl(i)
                   zi=zl(i)
                   if (jj==4) xi=xi-1.d0
                   if (jj==2) yi=yi+1.d0
                   if (jj==5) yi=yi-1.d0
                   if (jj==3) zi=zi+1.d0
                   if (jj==6) zi=zi-1.d0
                   call donne_num (xi,yi,zi, maille,num,interieur)
                   loc(i,jj,1)=num
                   loc(num,mod(jj+2,6)+1,1)=i
                end do
                do jj=1,j-1
                   xi=xl(i)
                   yi=yl(i)
                   zi=zl(i)
                   if (jj==1) xi=xi+1.d0
                   if (jj==4) xi=xi-1.d0
                   if (jj==2) yi=yi+1.d0
                   if (jj==5) yi=yi-1.d0
                   if (jj==3) zi=zi+1.d0
                   call donne_num (xi,yi,zi, maille,num,interieur)
                   loc(i,jj,1)=num
                   loc(num,mod(jj+2,6)+1,1)=i
                end do
             end if
       
          end if

       end do
    end do


end subroutine correction_maillage_local


subroutine donne_val_phi(xi,yi,zi,j,dx,dy,dz,maille,phi,valj)

use deftype
implicit none

real (kind = 8),  intent(in) :: xi,yi,zi
integer, intent(in) :: j
real (kind = 8), intent(in) :: dx,dy,dz
Type(chaineR3), pointer :: maille
real (kind = 8), dimension (:), intent(in) :: phi
real (kind = 8), intent(out) :: valj
        
  integer :: i1,i2,interieur
  real (kind = 8) :: xj,yj,zj,phixm,phiym,phizm
  

xj=xi
if (j==1) xj=xj+dx
if (j==4) xj=xj-dx
yj=yi
if (j==2) yj=yj+dy
if (j==5) yj=yj-dy
zj=zi
if (j==3) zj=zj+dz
if (j==6) zj=zj-dz
i1=0
i2=0
call donne_num (xj+dx,yj,zj, maille,i1,interieur)
call donne_num (xj-dx,yj,zj, maille,i2,interieur)
phixm=1.d20
if (i1>0) phixm=min(abs(phi(i1)),phixm)
if (i2>0) phixm=min(abs(phi(i2)),phixm)
 !
i1=0
i2=0
call donne_num (xj,yj+dy,zj, maille,i1,interieur)
call donne_num (xj,yj-dy,zj, maille,i2,interieur)
phiym=1.d20
if (i1>0) phiym=min(abs(phi(i1)),phiym)
if (i2>0) phiym=min(abs(phi(i2)),phiym)
!
i1=0
i2=0
call donne_num (xj,yj,zj, maille,i1,interieur)
call donne_num (xj,yj,zj, maille,i2,interieur)
phizm=1.d20
if (i1>0) phizm=min(abs(phi(i1)),phizm)
if (i2>0) phizm=min(abs(phi(i2)),phizm)
!calcul base sur resolution eikonale

valj=0.d0 ! a faire!!!!
end subroutine donne_val_phi


subroutine class_val(i,r,liste_add,np)

use deftype
implicit none

integer, intent(in) :: i
real (kind = 8), intent(in) :: r
type (chainej), Pointer   :: liste_add
integer, intent(inout) :: np



end subroutine class_val

subroutine donne_indice(i,liste_lib,inew)
use deftype
implicit none

integer, intent(in) :: i
type (chainej), Pointer   :: liste_lib
integer, intent(inout) :: inew

inew=i ! a faire



end subroutine donne_indice



subroutine construct_maillage(phi,dx,dist_secure,localise_u,localise_v,localise_w,&
localise_p,phi_u,phi_v,phi_w,phi_p,x,y,z,vit_vois_p)

use deftype
implicit none
        real (kind = 8), dimension (:,:,:), intent(in) :: phi
        real (kind = 8) :: dx,dist_secure
        integer, dimension (:,:,:), allocatable, intent(out) :: localise_u,localise_v,localise_w,localise_p
        real (kind = 8), dimension (:), allocatable, intent(out) :: phi_u,phi_v,phi_w,phi_p,x,y,z
        integer, dimension (:,:), allocatable, intent(out) :: vit_vois_p

        Type(chaineR3), pointer :: maille_px,maille_py,maille_pz,maille_u,maille_v,maille_w
        integer :: taille

taille = 0
      call defmaille_pression_x(maille_px,phi,taille,3.d0*dx,dist_secure)
print*,'taille pressionx',taille
      call defmaille_vitesse(maille_px,maille_u,taille,1)
      allocate (localise_u(taille,7,2),phi_u(taille))
print*,'taille u',taille
      call defmaille_pression_y(maille_py,phi,taille,3.d0*dx,dist_secure)
print*,'taille pressiony',taille
      call defmaille_vitesse(maille_py,maille_v,taille,2)
      allocate (localise_v(taille,7,2),phi_v(taille))
print*,'taille v',taille
      call defmaille_pression_z(maille_pz,phi,taille,3.d0*dx,dist_secure)
      allocate (localise_p(taille,7,2),phi_p(taille),x(taille),y(taille),z(taille), vit_vois_p(taille,6))
print*,'taille pressionz',taille
      call defmaille_vitesse(maille_pz,maille_w,taille,3)
      allocate (localise_w(taille,7,2),phi_w(taille))
print*,'taille w',taille
      call parcoursR3(maille_pz)
      call parcoursR3(maille_u)
      call parcoursR3(maille_v)
      call parcoursR3(maille_w)
!
print*,'apres parcours numerotant'
      call numerote(maille_u,localise_u,phi_u)
      call numerote(maille_v,localise_v,phi_v)
      call numerote(maille_w,localise_w,phi_w)
      call numerotebis(maille_pz,localise_p,phi_p,x,y,z)
      call vision_mac(localise_p,x,y,z,maille_u,maille_v,maille_w,vit_vois_p)

!penser a desallouer les maillages

    end subroutine construct_maillage

subroutine numerotemorsebis(lch,k)
    use deftype
    implicit none
    Type(listechaine), dimension(:), intent(in) :: lch
    integer, intent(out) ::k
    type (chainej), Pointer   :: chj
    integer :: i
    
    k=0
    do i=1, size(lch)
       chj=>lch(i)%liste
       do while (associated(chj%prec))
          chj=>chj%prec
       end do
       do while (associated(chj))
          k=k+1
          chj%k=k
          chj=>chj%suiv
       end do
    end do
!print*,'verification de la taille', k
end subroutine numerotemorsebis
!==================================================
!=======================================================================
   subroutine newbiskfonctionij(i,j,k,lch)
	use deftype
    implicit none
    integer, intent(in) :: i,j
    integer, intent(out) :: k
    type(listechaine), dimension(:), intent(in)  :: lch
    type(chainej), Pointer :: chj
    integer ::l
!
    chj=>lch(i)%liste
 do while (associated(chj))
       if (j>chj%j) chj=>chj%suiv
       if (j<chj%j) chj=>chj%prec
       if (j==chj%j) exit
    end do
    k=chj%k
end subroutine newbiskfonctionij
!=======================================================================
!=======================================================================
 subroutine test_symetrie(lch)
    use deftype
    implicit none
    Type(listechaine), dimension(:), intent(in)       :: lch

! 
    type (chainej), Pointer            :: chj,chjb
    integer :: i,k,tailleA
    real (kind=8) :: delta
    print*,'test symetrie'
  
    do k=1,size(lch)
       chj=>lch(k)%liste
       do while (associated(chj%prec))
          chj=>chj%prec
       end do
       do while (associated(chj))
          chjb=>lch(chj%j)%liste
          do while (associated(chjb%prec))
             chjb=>chjb%prec
          end do
          do while (associated(chjb))
             if (chjb%j==k) then
                delta=chjb%val-chj%val
!print*,'i,j,jb,val,valj',k,chj%j,chjb%j,chj%val,chjb%val
                if (abs(delta)>1.d-12) print*,'non symetrie'
             end if
             chjb=>chjb%suiv
          end do
         
          chj=>chj%suiv
       end do
    end do
  
  end subroutine test_symetrie
!=======================================================================
!=======================================================================
 subroutine produit(lch,u,v)
    use deftype
    implicit none
    Type(listechaine), dimension(:), intent(in)       :: lch
    real(kind=8),dimension(:), intent(in )  :: u
    real(kind=8),dimension(:), intent(out)  :: v
! 
    type (chainej), Pointer            :: chj
    integer :: i,k,tailleA
    
    v=0.d0
    do i=1,size(lch)
       chj=>lch(i)%liste
       do while (associated(chj%prec))
          chj=>chj%prec
       end do
       do while (associated(chj))
          v(i)=v(i)+chj%val*u(chj%j)
          chj=>chj%suiv
       end do
    end do
  end subroutine produit
!=======================================================================
!=======================================================================
 subroutine conversionchainemorsebis(lch,A)
    use deftype
    implicit none
    Type(listechaine), dimension(:), intent(in)       :: lch
    type(morse),dimension(:), intent(out)  :: A
! 
    type (chainej), Pointer            :: chj
    integer :: i,k,tailleA
    
    k=0
    do i=1,size(lch)
       chj=>lch(i)%liste
       do while (associated(chj%prec))
          chj=>chj%prec
       end do
       do while (associated(chj))
          k=k+1
   !       print*,'k',k,chj%j
          A(k)%i=i
          A(k)%j=chj%j
          A(k)%val=chj%val
          chj=>chj%suiv
       end do
    end do
tailleA=k
!print*,'verification de la taille', tailleA
end subroutine conversionchainemorsebis
!=======================================================================
 subroutine desallouelistebis(ch)
    use deftype
    implicit none
    Type(listechaine), dimension(:), intent(inout)         :: ch
! 
    type (chainej), Pointer            :: chj,chjb,chjbb,chjbef
    integer :: i
    
! deallocation 
! deallocation partie droite de sous chaine
    do i=1,size(ch)
       chj=>ch(i)%liste
       if (associated(chj)) then
          chjb=>chj
          chjbb=>chj%prec
          do while (associated(chjb))
             chjbef=>chjb
             chjb=>chjb%suiv
             deallocate(chjbef)
          end do
          do while (associated(chjbb))
             chjbef=>chjbb
             chjbb=>chjbb%prec
             deallocate(chjbef)
          end do
       endif

    end do
!
  end subroutine desallouelistebis
!=======================================================================
subroutine annuleval_liste(chj)
    use deftype
    implicit none
    Type(chainej), Pointer        :: chj
! 
    type (chainej), Pointer            :: chjb,chjbt
    
! 
!  partie droite de sous chaine

       chjb=>chj
       do while (associated(chjb))
          chjbt=>chjb
          chjb=>chjb%suiv
          chjbt%val=0.d0
       end do
!  partie gauche de sous chaine
       chjb=>chj
       do while (associated(chjb))
          chjbt=>chjb
          chjb=>chjb%prec
          chjbt%val=0.d0
       end do
!  point de depart de sous chaine
       if (associated (chj)) chj%val=0.d0
!
end subroutine annuleval_liste
!!$ subroutine ajout(i,j,val,ch,t)
!!$use deftype
!!$    implicit none
!!$ 
!!$    integer, intent(in)                  :: i,j
!!$    Type(dblechaine), Pointer            :: ch
!!$    real*8, intent(in)                   :: val
!!$    integer, intent(inout)               :: t
!!$
!!$    logical                         :: droite, droitedep
!!$!    integer                         :: k,valntemp
!!$!    real*8                          :: value
!!$
!!$   Type(dblechaine), Pointer             :: chloc,temp
!!$   Type(chainej), Pointer                :: loclj 
!!$
!!$    ! 
!!$    !
!!$    !
!!$    ! inrterdit de rentrer dans cette routine si la chaine pointe sur null
!!$    droite=.false.
!!$    if (i>=ch%i) droite=.true.
!!$    do while (associated(ch))
!!$          temp=>ch
!!$          if (i>ch%i) then
!!$             if (droite) then
!!$                ch=>ch%suiv
!!$             else
!!$                t=t+1
!!$                call insertdroite(i,j,val,ch)
!!$                return
!!$             end if
!!$          else
!!$             if (i==ch%i) then 
!!$                call ajoutval(j,val,ch%lj,t)
!!$                return
!!$             else
!!$                if (droite) then 
!!$                   t=t+1
!!$                   call insertgauche(i,j,val,ch)
!!$                   return
!!$                else
!!$                   ch=>ch%prec
!!$                end if
!!$             end if
!!$          end if
!!$       end do
!!$!cas du bout de liste droite ou gauche
!!$       allocate (loclj)
!!$       t=t+1
!!$       loclj%j=j
!!$       loclj%val=val
!!$       nullify(loclj%suiv)
!!$       nullify(loclj%prec)
!!$       allocate(chloc)
!!$       chloc%i=i
!!$       chloc%lj=>loclj
!!$       if (droite) then
!!$!ajout à droite de temp
!!$          chloc%prec=>temp
!!$          nullify(chloc%suiv)
!!$          temp%suiv=>chloc
!!$       else
!!$!ajout à gauche  de temp
!!$          chloc%suiv=>temp
!!$          nullify(chloc%prec)
!!$          temp%prec=>chloc
!!$       end if
!!$       ch=>temp
!!$       return
!!$  end subroutine ajout

!=======================================================================
  subroutine ajoutval(j,val,ch,t)
    use deftype
    implicit none
    integer, intent(in)               :: j
    Type(chainej), Pointer            :: ch
    real(kind=kind(0.d0)), intent(in) :: val
    integer, intent(inout)            :: t
    logical                           :: droite
    Type(chainej), Pointer            :: temp, loclj 
    !
    !
!print*,'entree ajout val',t,j
       droite=.false.
       if (j>=ch%j) droite=.true.
       do while (associated(ch))
          temp=>ch
          if (j>ch%j) then
             if (droite) then
                ch=>ch%suiv
             else
                t=t+1
                call insertdroitelj(j,val,ch)
                return
             end if
          else
             if (j==ch%j) then 
                ch%val=ch%val+val
                !             print*,'ajout au point declare',t,j
                return
             else
                if (droite) then 
                   t=t+1
                   call insertgauchelj(j,val,ch)
                   return
                else
                   ch=>ch%prec
                end if
             end if
          end if
       end do
       !cas du bout de liste droite ou gauche
       allocate (loclj)
       t=t+1
       loclj%j=j
       loclj%val=val
       if (droite) then
!ajout a droite de temp
          loclj%prec=>temp
          nullify(loclj%suiv)
          temp%suiv=>loclj
       else
!ajout a gauche de temp
          loclj%suiv=>temp
          nullify(loclj%prec)
          temp%prec=>loclj
       end if
       ch=>temp
       return
  end subroutine ajoutval
!=======================================================================


!=======================================================================
  subroutine ecraseval(j,val,ch,t)
    use deftype
    implicit none
    integer, intent(in)               :: j
    Type(chainej), Pointer            :: ch
    real(kind=kind(0.d0)), intent(in) :: val
    integer, intent(inout)            :: t
    logical                           :: droite
    Type(chainej), Pointer            :: temp, loclj 
    !
    !
!print*,'entree ajout val',t,j
       droite=.false.
       if (j>=ch%j) droite=.true.
       do while (associated(ch))
          temp=>ch
          if (j>ch%j) then
             if (droite) then
                ch=>ch%suiv
             else
                t=t+1
print*,'probleme dans ecrase val'
                call insertdroitelj(j,val,ch)
                return
             end if
          else
             if (j==ch%j) then 
                ch%val=val
                !             print*,'ajout au point declare',t,j
                return
             else
                if (droite) then 
                   t=t+1
print*,'probleme dans ecrase val'
                   call insertgauchelj(j,val,ch)
                   return
                else
                   ch=>ch%prec
                end if
             end if
          end if
       end do
       !cas du bout de liste droite ou gauche
print*,'probleme dans ecrase val'
       allocate (loclj)
       t=t+1
       loclj%j=j
       loclj%val=val
       if (droite) then
!ajout a droite de temp
          loclj%prec=>temp
          nullify(loclj%suiv)
          temp%suiv=>loclj
       else
!ajout a gauche de temp
          loclj%suiv=>temp
          nullify(loclj%prec)
          temp%prec=>loclj
       end if
       ch=>temp
       return
     end subroutine ecraseval
!=======================================================================

!=======================================================================
 subroutine insertgauchelj(j,val,ch)
  use deftype
    implicit none
    integer, intent(in)               :: j
    Type(chainej), Pointer            :: ch
    real(kind=kind(0.d0)), intent(in) :: val
    Type(chainej), Pointer            :: temp, loclj
    !
!cas du bout de liste droite ou gauche
    temp=>ch%prec
       allocate (loclj)
       loclj%j=j
       loclj%val=val
       loclj%suiv=>ch
       ch%prec=>loclj
       loclj%prec=>temp
       temp%suiv=>loclj
       return
  end subroutine insertgauchelj
!=======================================================================

!=======================================================================
 subroutine insertdroitelj(j,val,ch)
  use deftype
    implicit none
    integer, intent(in)               :: j
    Type(chainej), Pointer            :: ch
    real(kind=kind(0.d0)), intent(in) :: val
    Type(chainej), Pointer            :: temp, loclj
    !
!cas du bout de liste droite ou gauche
    temp=>ch%suiv
       allocate (loclj)
       loclj%j=j
       loclj%val=val
       loclj%suiv=>temp
       temp%prec=>loclj
       loclj%prec=>ch
       ch%suiv=>loclj
       return
  end subroutine insertdroitelj
!=======================================================================



subroutine defmaille_vitesse(P3,U3,taille,ind)
use deftype
  implicit none

  !-------------------------------------------------
  !          declaration des variables
  !-------------------------------------------------
!
   Type(chaineR3), pointer :: P3,U3,dblech

!  

    type (chaineR2), Pointer            :: ch2, R2
    type (chaineR1), Pointer            :: ch1, R1  
    integer :: taille,ind
!

allocate(R1,R2,U3)
R1%z=0.d0
Nullify(R1%suiv) ; Nullify(R1%prec)
R2%y=0.d0
R2%cr1=>R1
nullify(R2%suiv) ; Nullify(R2%prec)
U3%x=-1.d100
U3%cr2=>R2
nullify(U3%suiv) ; Nullify(U3%prec)
!
!
! il suffit de parcours les faces
! pour chaque points (x,y, z) issu des maillages
! il reste à les ranger
!
!
!

 


  taille =0  

    dblech=>P3
    do while (associated(dblech%prec))
       dblech=>dblech%prec
    end do
    do while (associated(dblech))
       ch2=>dblech%cr2
       do while (associated(ch2%prec))
          ch2=>ch2%prec
       end do
       do while(associated(ch2))
          ch1=>ch2%cr1
          do while (associated(ch1%prec))
             ch1=>ch1%prec
          end do
          do while(associated(ch1).and.associated(ch1%suiv))
             if (ind==3) call  rangeR3 (dblech%x,ch2%y,(ch1%z+ch1%suiv%z)*0.5d0, U3,taille)
             if (ind==2) call  rangeR3 (ch2%y,(ch1%z+ch1%suiv%z)*0.5d0,dblech%x, U3,taille)
             if (ind==1) call  rangeR3 ((ch1%z+ch1%suiv%z)*0.5d0, dblech%x,ch2%y, U3,taille)
             U3%cr2%cr1%phi=(ch1%phi+ch1%suiv%phi)*0.5d0
             U3%cr2%cr1%interieur=0    ! exterieur sauf si cas suivant
           !  if (ch1%phi+ch1%suiv%phi<0.d0) U3%cr2%cr1%interieur=2 !bord 
             if (ch1%interieur==1.and.ch1%suiv%interieur==1) U3%cr2%cr1%interieur=1 !interieur
           !  if (ch1%interieur+ch1%suiv%interieur==1) U3%cr2%cr1%interieur=2 !bord
            ! if (ch1%interieur==0.and.ch1%suiv%interieur==0) U3%cr2%cr1%interieur=0 !exterieur
!!$! ajout  des bords inferieurs: en fait inutile
!!$            if (associated(ch1).and.associated(ch1%prec)) then
!!$                if (ind==3) call  rangeR3 (dblech%x,ch2%y,(ch1%z+ch1%prec%z)*0.5d0, U3,taille)
!!$                if (ind==2) call  rangeR3 (ch2%y,(ch1%z+ch1%prec%z)*0.5d0,dblech%x, U3,taille)
!!$                if (ind==1) call  rangeR3 ((ch1%z+ch1%prec%z)*0.5d0, dblech%x,ch2%y, U3,taille)
!!$                U3%cr2%cr1%phi=(ch1%phi+ch1%prec%phi)*0.5d0
!!$                U3%cr2%cr1%interieur=0    ! exterieur sauf si cas suivant
!!$                if (ch1%interieur==1.and.ch1%prec%interieur==1) U3%cr2%cr1%interieur=1 !interieur
!!$             end if
!!$!fin d'ajout
             ch1=>ch1%suiv
          end do
          ch2=>ch2%suiv
       end do
       dblech=>dblech%suiv
    end do

end subroutine defmaille_vitesse
  !-----------------------------------------------



subroutine defmaille_pression_x(R3,phi,taille,dist, delta)
use deftype
  implicit none

  !-------------------------------------------------
  !          declaration des variables
  !-------------------------------------------------
!
   Type(chaineR3), pointer :: R3
   real(kind=8), dimension (:,:,:), intent(in)  :: phi
  integer, intent(out)                              :: taille
  real(kind=8), intent(in)                              :: dist, delta
!
  integer                                           :: i,j,k,t
!  
  Type(chaineR2), pointer :: R2
  Type(chaineR1), pointer :: R1
!
 taille=0
!print*,'???'
allocate(R1,R2,R3)
R1%z=0.d0
Nullify(R1%suiv) ; Nullify(R1%prec)
R2%y=0.d0
R2%cr1=>R1
nullify(R2%suiv) ; Nullify(R2%prec)
R3%x=-1.d100
R3%cr2=>R2
nullify(R3%suiv) ; Nullify(R3%prec)
!
!
! il suffit de parcours les faces
! pour chaque points (x,y, z) issu des maillages
! il reste à les ranger
!
!
!
t =0
do k=2,size(phi,3)-1
  do j=2,size(phi,2)-1
     do i=2,size(phi,1)-1
!print*,'???',size(phi,3),size(phi,2),size(phi,1)
        if (phi(i,j,k)<dist) then
          ! print*,'avant ranger3',j
           call rangeR3 (j*1.d0,k*1.d0,i*1.d0, R3,taille)
          ! print*,'apres ranger3',R3%x
           R3%cr2%cr1%interieur=0
           R3%cr2%cr1%phi=phi(i,j,k)
           if (phi(i,j,k)+phi(i-1,j,k)<delta.and.phi(i,j,k)+phi(i+1,j,k)<delta.and.&
                phi(i,j,k)+phi(i,j-1,k)<delta.and.phi(i,j,k)+phi(i,j+1,k)<delta.and.&
                phi(i,j,k)+phi(i,j,k-1)<delta.and.phi(i,j,k)+phi(i,j,k+1)<delta&
            !    .and.k>2.and.k<size(phi,3)-1.and.j>2.and.j<size(phi,2)-1.and.&
            !    i>2.and.i<size(phi,1)-1)
                ) R3%cr2%cr1%interieur=1
           t=t+1
        end if
!print*,t
     end do
  end do
end do

end subroutine defmaille_pression_x
subroutine defmaille_pression_y(R3,phi,taille,dist, delta)
use deftype
  implicit none

  !-------------------------------------------------
  !          declaration des variables
  !-------------------------------------------------
!
   Type(chaineR3), pointer :: R3
   real(kind=8), dimension (:,:,:), intent(in)  :: phi
  integer, intent(out)                              :: taille
  real(kind=8), intent(in)                              :: dist, delta
!
  integer                                           :: i,j,k,t
!  
  Type(chaineR2), pointer :: R2
  Type(chaineR1), pointer :: R1
!
 taille=0
allocate(R1,R2,R3)
R1%z=0.d0
Nullify(R1%suiv) ; Nullify(R1%prec)
R2%y=0.d0
R2%cr1=>R1
nullify(R2%suiv) ; Nullify(R2%prec)
R3%x=-1.d100
R3%cr2=>R2
nullify(R3%suiv) ; Nullify(R3%prec)
!
!
! il suffit de parcours les faces
! pour chaque points (x,y, z) issu des maillages
! il reste à les ranger
!
!
!
t =0
do k=2,size(phi,3)-1
  do j=2,size(phi,2)-1
     do i=2,size(phi,1)-1
        if (phi(i,j,k)<dist) then
           call rangeR3 (k*1.d0,i*1.d0,j*1.d0, R3,taille)
           R3%cr2%cr1%phi=phi(i,j,k)
           R3%cr2%cr1%interieur=0
          ! print*,'i,jk',i,k,j,R3%x,R3%cr2%y,R3%cr2%cr1%z
           if (phi(i,j,k)+phi(i-1,j,k)< delta.and.phi(i,j,k)+phi(i+1,j,k)< delta.and.&
                phi(i,j,k)+phi(i,j-1,k)< delta.and.phi(i,j,k)+phi(i,j+1,k)< delta.and.&
                phi(i,j,k)+phi(i,j,k-1)< delta.and.phi(i,j,k)+phi(i,j,k+1)< delta) R3%cr2%cr1%interieur=1
           t=t+1
        end if
     end do
  end do
end do

end subroutine defmaille_pression_y
subroutine defmaille_pression_z(R3,phi,taille,dist,delta)
use deftype
  implicit none

  !-------------------------------------------------
  !          declaration des variables
  !-------------------------------------------------
!
   Type(chaineR3), pointer :: R3
   real(kind=8), dimension (:,:,:), intent(in)  :: phi
  integer, intent(out)                              :: taille
  real(kind=8), intent(in)                              :: dist, delta
!
  integer                                           :: i,j,k,t
!  
  Type(chaineR2), pointer :: R2
  Type(chaineR1), pointer :: R1
!
 taille=0
allocate(R1,R2,R3)
R1%z=0.d0
Nullify(R1%suiv) ; Nullify(R1%prec)
R2%y=0.d0
R2%cr1=>R1
nullify(R2%suiv) ; Nullify(R2%prec)
R3%x=-1.d100
R3%cr2=>R2
nullify(R3%suiv) ; Nullify(R3%prec)
!
!
! il suffit de parcours les faces
! pour chaque points (x,y, z) issu des maillages
! il reste à les ranger
!
!
!
t =0
do k=2,size(phi,3)-1
  do j=2,size(phi,2)-1
     do i=2,size(phi,1)-1
        if (phi(i,j,k)<dist) then
           call rangeR3 (i*1.d0,j*1.d0,k*1.d0, R3,taille)
       !    print*,'i,jk',i,j,k,R3%x,R3%cr2%y,R3%cr2%cr1%z
           R3%cr2%cr1%phi=phi(i,j,k)
           R3%cr2%cr1%interieur=0
           if (phi(i,j,k)+phi(i-1,j,k)<delta.and.phi(i,j,k)+phi(i+1,j,k)<delta.and.&
                phi(i,j,k)+phi(i,j-1,k)<delta.and.phi(i,j,k)+phi(i,j+1,k)<delta.and.&
                phi(i,j,k)+phi(i,j,k-1)<delta.and.phi(i,j,k)+phi(i,j,k+1)<delta) R3%cr2%cr1%interieur=1
           t=t+1
        end if
     end do
  end do
end do

end subroutine defmaille_pression_z
  !-----------------------------------------------



 subroutine parcoursR3(ch)
use deftype
    implicit none
    Type(chaineR3), Pointer               :: ch
! 
   Type(chaineR3), Pointer               :: dblech
    type (chaineR2), Pointer            :: ch2
    type (chaineR1), Pointer            :: ch1
    integer :: num


    
   num=0
    dblech=>ch
    do while (associated(dblech%prec))
       dblech=>dblech%prec
    end do
    do while (associated(dblech))
       ch2=>dblech%cr2
       do while (associated(ch2%prec))
          ch2=>ch2%prec
       end do
       do while(associated(ch2))
          ch1=>ch2%cr1
          do while (associated(ch1%prec))
             ch1=>ch1%prec
          end do
          do while(associated(ch1))
             ch1%num=num
             num=num+1
             !print*,'x y z num',dblech%x,ch2%y,ch1%z,ch1%num
             
             ch1=>ch1%suiv
          end do
          ch2=>ch2%suiv
       end do
       dblech=>dblech%suiv
    end do
!print*,'dernier num dans parcours',num-1
end subroutine parcoursR3
!

 subroutine numerote(ch,tab,phi)
use deftype
    implicit none
 
    Type(chaineR3), Pointer            :: ch
    integer, dimension(:,:,:), intent(out)    :: tab 
    real (kind = 8), dimension (:), intent(out)  :: phi
!
   Type(chaineR3), Pointer               :: dblech,tempR3
    type (chaineR2), Pointer            :: ch2
    type (chaineR1), Pointer            :: ch1
    integer :: num,interieur


    

    dblech=>ch
    do while (associated(dblech%prec))
       dblech=>dblech%prec
    end do
    dblech=>dblech%suiv ! on evite le point d initialisation
    do while (associated(dblech))
       ch2=>dblech%cr2
       do while (associated(ch2%prec))
          ch2=>ch2%prec
       end do
       do while(associated(ch2))
          ch1=>ch2%cr1
          do while (associated(ch1%prec))
             ch1=>ch1%prec
          end do
          do while(associated(ch1))
!             print*,'toto'
!             tab(ch1%num,:)=ch1%num
             num=ch1%num
             phi(ch1%num)=ch1%phi
             tab(ch1%num,7,1)=num
             tab(ch1%num,7,2)=ch1%interieur
             if (num==0) print*,'pb de premier noeud'
!             print*,'num',num
             tempR3=>dblech
!print*,'num pos pointeur',num!,dblech%x+1.d0,ch2%y,ch1%z
             call donne_num (dblech%x+1.d0,ch2%y,ch1%z, tempR3,num,interieur)
!print*,'num à droite en x',num
             tab(ch1%num,1,1)=num
             tab(ch1%num,1,2)=interieur
             tempR3=>dblech
!print*,'num pos pointeur',num,tempR3%x,tempR3%cr2%y,tempR3%cr2%cr1%z
             num=ch1%num
             call donne_num (dblech%x ,ch2%y+1.d0,ch1%z , tempR3,num,interieur)
    !                print*,'num à droite en y',num
             tab(ch1%num,2,1)=num
             tab(ch1%num,2,2)=interieur
             tempR3=>dblech
             num=ch1%num
             call donne_num (dblech%x,ch2%y , ch1%z+1.d0, tempR3,num,interieur)
   !          print*,'num à droite en z',num
             tab(ch1%num,3,1)=num
             tab(ch1%num,3,2)=interieur
             tempR3=>dblech
             num=ch1%num
             call donne_num (dblech%x-1.d0,ch2%y,ch1%z, tempR3,num,interieur)
  !           print*,'num à gauche en x',num
             if (num==0) then
                tab(ch1%num,4,1)=ch1%num
                else
                tab(ch1%num,4,1)=num
             end if
             tab(ch1%num,4,2)=interieur
             tempR3=>dblech
             num=ch1%num
             call donne_num (dblech%x ,ch2%y-1.d0,ch1%z , tempR3,num,interieur)
      !       print*,'num à gauche en y',num,ch1%num,size(tab,1)
             tab(ch1%num,5,1)=num
             tab(ch1%num,5,2)=interieur
!print*,'toto'
             tempR3=>dblech
!print*,'avant le 6',num,tempR3%x,tempR3%cr2%y,tempR3%cr2%cr1%z
             num=ch1%num
             call donne_num (dblech%x,ch2%y , ch1%z-1.d0, tempR3,num,interieur)
 !            print*,'num à gauche en z',num
             tab(ch1%num,6,1)=num
             tab(ch1%num,6,2)=interieur
             ch1=>ch1%suiv
!print*,'totofin1'
          end do
          ch2=>ch2%suiv
!print*,'totofi2'
       end do
       dblech=>dblech%suiv
!print*,'totofin3'
    end do    
  end subroutine numerote
!================================================================================================================
 subroutine numerotebis(ch,tab,phi,x,y,z)
use deftype
    implicit none
 
    Type(chaineR3), Pointer            :: ch
    integer, dimension(:,:,:), intent(out)    :: tab 
    real (kind = 8), dimension (:), intent(out)  :: phi,x,y,z
!
   Type(chaineR3), Pointer               :: dblech,tempR3
    type (chaineR2), Pointer            :: ch2
    type (chaineR1), Pointer            :: ch1
    integer :: num,interieur


    

    dblech=>ch
    do while (associated(dblech%prec))
       dblech=>dblech%prec
    end do
    dblech=>dblech%suiv ! on evite le point d initialisation
    do while (associated(dblech))
       ch2=>dblech%cr2
       do while (associated(ch2%prec))
          ch2=>ch2%prec
       end do
       do while(associated(ch2))
          ch1=>ch2%cr1
          do while (associated(ch1%prec))
             ch1=>ch1%prec
          end do
          do while(associated(ch1))
!             print*,'toto'
!             tab(ch1%num,:)=ch1%num
             num=ch1%num
             phi(ch1%num)=ch1%phi
             x(ch1%num)=dblech%x
             y(ch1%num)=ch2%y
             z(ch1%num)=ch1%z
             tab(ch1%num,7,1)=num
             tab(ch1%num,7,2)=ch1%interieur
             if (num==0) print*,'pb de premier noeud'
!             print*,'num',num
             tempR3=>dblech
!print*,'num pos pointeur',num!,dblech%x+1.d0,ch2%y,ch1%z
             call donne_num (dblech%x+1.d0,ch2%y,ch1%z, tempR3,num,interieur)
!print*,'num à droite en x',num
            ! if (interieur==0) num=ch1%num
             tab(ch1%num,1,1)=num
             tab(ch1%num,1,2)=interieur
             tempR3=>dblech
             num=ch1%num
!print*,'num pos pointeur',num,tempR3%x,tempR3%cr2%y,tempR3%cr2%cr1%z
             call donne_num (dblech%x ,ch2%y+1.d0,ch1%z , tempR3,num,interieur)
    !                print*,'num à droite en y',num
            ! if (interieur==0) num=ch1%num
             tab(ch1%num,2,1)=num
             tab(ch1%num,2,2)=interieur
             tempR3=>dblech
             num=ch1%num
             call donne_num (dblech%x,ch2%y , ch1%z+1.d0, tempR3,num,interieur)
   !          print*,'num à droite en z',num
            ! if (interieur==0) num=ch1%num
             tab(ch1%num,3,1)=num
             tab(ch1%num,3,2)=interieur
             tempR3=>dblech
             num=ch1%num
             call donne_num (dblech%x-1.d0,ch2%y,ch1%z, tempR3,num,interieur)
  !           print*,'num à gauche en x',num
            ! if (interieur==0) num=ch1%num
              if (num==0) then
                tab(ch1%num,4,1)=ch1%num
                else
                tab(ch1%num,4,1)=num
             end if  
             tab(ch1%num,4,2)=interieur
             tempR3=>dblech
             num=ch1%num
             call donne_num (dblech%x ,ch2%y-1.d0,ch1%z , tempR3,num,interieur)
      !       print*,'num à gauche en y',num,ch1%num,size(tab,1)
            ! if (interieur==0) num=ch1%num
             tab(ch1%num,5,1)=num
             tab(ch1%num,5,2)=interieur
!print*,'toto'
             tempR3=>dblech
!print*,'avant le 6',num,tempR3%x,tempR3%cr2%y,tempR3%cr2%cr1%z
             num=ch1%num
             call donne_num (dblech%x,ch2%y , ch1%z-1.d0, tempR3,num,interieur)
 !            print*,'num à gauche en z',num
            ! if (interieur==0) num=ch1%num
             tab(ch1%num,6,1)=num
             tab(ch1%num,6,2)=interieur
             ch1=>ch1%suiv
!print*,'totofin1'
          end do
          ch2=>ch2%suiv
!print*,'totofi2'
       end do
       dblech=>dblech%suiv
!print*,'totofin3'
    end do    
  end subroutine numerotebis

!================================================================================================================
 subroutine centrage(vois,u,v,w,up,vp,wp)
use deftype
    implicit none
    integer, dimension(:,:), intent(in)    :: vois
    real (kind = 8), dimension (:), intent(in)  :: u,v,w
    real (kind = 8), dimension (:), intent(out)  :: up,vp,wp
    integer ::i
up=0.d0
vp=0.d0
wp=0.d0
    do i=1,size(vois,1)
!print*,'i',i,vois(i,:)
       if (vois(i,1)>0.and.vois(i,4)>0) up(i)=(u(vois(i,1))+u(vois(i,4)))*0.5d0
       if (vois(i,2)>0.and.vois(i,5)>0) vp(i)=(v(vois(i,2))+v(vois(i,5)))*0.5d0
       if (vois(i,3)>0.and.vois(i,6)>0) wp(i)=(w(vois(i,3))+w(vois(i,6)))*0.5d0
    end do

  end subroutine centrage

!================================================================================================================
 subroutine vision_mac(loc,x,y,z,ch_u,ch_v,ch_w,vois)
use deftype
    implicit none
 
    integer, dimension(:,:,:), intent(in)    :: loc
    real (kind = 8), dimension (:), intent(in)  :: x,y,z
    Type(chaineR3), Pointer            :: ch_u,ch_v,ch_w
    integer, dimension(:,:), intent(out)    :: vois

!
    Type(chaineR3), Pointer               :: tempR3
    integer :: i,num,interieur
    vois=1
    do i=1,size(loc,1)
      ! if (loc(i,7,2)==1) then 
! critere de maille pression entouree de maille pression: les vitesses ont donc du sens
    !   if (loc(i,1,1)/=loc(i,7,1).and.loc(i,2,1)/=loc(i,7,1).and.loc(i,3,1)/=loc(i,7,1).and.loc(i,4,1)/=loc(i,7,1)&
    !        .and.loc(i,5,1)/=loc(i,7,1).and.loc(i,6,1)/=loc(i,7,1)) then
          tempR3=>ch_u
          num=0
          call donne_num (x(i)+0.5d0,y(i),z(i), tempR3,num,interieur)
          if (num==0) then 
             call donne_num (x(i)-0.5d0,y(i),z(i), tempR3,num,interieur)
            ! print*,'probleme1',ch_u%cr2%cr1%num
             vois(i,1)=num!tempR3%cr2%cr1%num
          else
             vois(i,1)=num
          endif
          num=0
          call donne_num (x(i)-0.5d0,y(i),z(i), tempR3,num,interieur)
          if (num==0) then 
             call donne_num (x(i)+0.5d0,y(i),z(i), tempR3,num,interieur)
!             print*,'probleme4'
             vois(i,4)=num!tempR3%cr2%cr1%num
          else
             vois(i,4)=num
          endif
          tempR3=>ch_v
          num=0
          call donne_num (x(i),y(i)+0.5d0,z(i), tempR3,num,interieur)
          if (num==0) then 
             call donne_num (x(i),y(i)-0.5d0,z(i), tempR3,num,interieur)
           !  print*,'probleme2'
             vois(i,2)=num!tempR3%cr2%cr1%num
          else
             vois(i,2)=num
          endif
          num=0
          call donne_num (x(i),y(i)-0.5d0,z(i), tempR3,num,interieur)
          if (num==0) then 
             call donne_num (x(i),y(i)+0.5d0,z(i), tempR3,num,interieur)
            ! print*,'probleme5'
             vois(i,5)=num!tempR3%cr2%cr1%num
          else
             vois(i,5)=num
          endif
          tempR3=>ch_w
          num=0
          call donne_num (x(i),y(i),z(i)+0.5d0, tempR3,num,interieur)
          if (num==0) then 
             call donne_num (x(i),y(i),z(i)-0.5d0, tempR3,num,interieur)
             vois(i,3)=num!tempR3%cr2%cr1%num
             !print*,'probleme3',tempR3%cr2%cr1%num,tempR3%x-x(i),tempR3%cr2%y-y(i),tempR3%cr2%cr1%z-z(i)
          else
             vois(i,3)=num
          endif
          call donne_num (x(i),y(i),z(i)-0.5d0, tempR3,num,interieur)
          if (num==0) then 
             call donne_num (x(i),y(i),z(i)+0.5d0, tempR3,num,interieur)
             !print*,'probleme6',tempR3%cr2%cr1%num,tempR3%x-x(i),tempR3%cr2%y-y(i),tempR3%cr2%cr1%z-z(i)
             vois(i,6)=num!tempR3%cr2%cr1%num
          else
             vois(i,6)=num
          endif
 !      end if
    end do
  end subroutine vision_mac
!================================================================================================================
 subroutine donne_num(x,y,z,R3,num,interieur)
use deftype
    implicit none
 
    Type(chaineR3), Pointer            :: R3
    real*8, intent(in)                   :: x,y,z
    integer, intent(inout)               :: num,interieur

    logical                         :: droite, droitedep

   Type(chaineR3), Pointer             :: chloc,temp
   Type(chaineR2), Pointer                :: locR2 
   Type(chaineR1), Pointer                :: locR1

    ! 
    !
    !
    ! inrterdit de rentrer dans cette routine si la chaine pointe sur null
    droite=.false.
    if (x>=R3%x) droite=.true.
    do while (associated(R3))
          temp=>R3
          if (x>R3%x+1.d-12) then
             if (droite) then
                R3=>R3%suiv
             else
      !     print*,'probleme'
                R3=>temp
                return
             end if
          else
             if (abs(x-R3%x)<1.d-12) then 
 !               print*,'avant R2'
                call donne_numR2(y,z,R3%cr2,num,interieur)
                return
             else
                if (droite) then 
  !      print*,'probleme'
                   R3=>temp
                    return
                else
                   R3=>R3%prec
                end if
             end if
          end if
       end do
       R3=>temp
!print*,'probleme'
       return
     end subroutine donne_num
!-----------------------------------------------
!-----------------------------------------------
 subroutine donne_numR2(y,z,R2,t,in)
use deftype
    implicit none
 
    Type(chaineR2), Pointer            :: R2
    real*8, intent(in)                   :: y,z
    integer, intent(inout)               :: t,in

    logical                         :: droite, droitedep

   Type(chaineR2), Pointer             :: chloc,temp
   Type(chaineR1), Pointer                :: locR1

    ! 
    !
    !
    ! inrterdit de rentrer dans cette routine si la chaine pointe sur null
    droite=.false.
!print*,'coucou',R2%y
    if (y>=R2%y) droite=.true.
!print*,'droite',droite
    do while (associated(R2))
          temp=>R2
          if (y>R2%y+1.d-12) then
             if (droite) then
                R2=>R2%suiv
             else
                R2=>temp
       !      print*,'probleme'
                return
             end if
          else
             if (abs(y-R2%y)<1.d-12) then 
            !    print*,'avant R1'
                call donne_numR1(z,R2%cr1,t,in)
                return
             else
                if (droite) then 
       !     print*,'probleme'
                   R2=>temp
                   return
                else
                   R2=>R2%prec
                end if
             end if
          end if
!print*,'que passa',R2%y
       end do
       R2=>temp
!print*,'probleme'
       return
     end subroutine donne_numR2
!--------------------------
!-----------------------------------------------
 subroutine donne_numR1(z,R1,t,in)
use deftype
    implicit none
 
    Type(chaineR1), Pointer            :: R1
    real*8, intent(in)                   :: z
    integer, intent(inout)               :: t,in

    logical                         :: droite, droitedep

   Type(chaineR1), Pointer             :: chloc,temp
    ! 
    !
    !
    ! inrterdit de rentrer dans cette routine si la chaine pointe sur null
    droite=.false.
    if (z>=R1%z) droite=.true.
    do while (associated(R1))
          temp=>R1
          if (z>R1%z+1.d-12) then
             if (droite) then
                R1=>R1%suiv
             else
       !       print*,'probleme'
                R1=>temp
                return
             end if
          else
             if (abs(z-R1%z)<1.d-12) then 
                t=R1%num
                in=R1%interieur
                return
             else
                if (droite) then 
        !      print*,'probleme'
                   R1=>temp
                   return
                else
                   R1=>R1%prec
                end if
             end if
          end if
       end do
       R1=>temp
!print*,'probleme'
       return
  end subroutine donne_numR1
!---------------------------
 subroutine rangeR3(x,y,z,R3,t)
use deftype
    implicit none
 
    Type(chaineR3), Pointer            :: R3
    real*8, intent(in)                   :: x,y,z
    integer, intent(inout)               :: t

    logical                         :: droite, droitedep

   Type(chaineR3), Pointer             :: chloc,temp
   Type(chaineR2), Pointer                :: locR2 
   Type(chaineR1), Pointer                :: locR1

    ! 
    !
    !
    ! inrterdit de rentrer dans cette routine si la chaine pointe sur null
    droite=.false.
    if (x>=R3%x) droite=.true.
    do while (associated(R3))
          temp=>R3
          if (x>R3%x+1.d-12) then
             if (droite) then
                R3=>R3%suiv
             else
                t=t+1
                call insertdroitex(x,y,z,R3)
                return
             end if
          else
             if (abs(x-R3%x)<1e-12) then 
                call rangeR2(y,z,R3%cr2,t)
                return
             else
                if (droite) then 
                   t=t+1
                   call   insertgauchex(x,y,z,R3)
                   return
                else
                   R3=>R3%prec
                end if
             end if
          end if
       end do
!cas du bout de liste droite ou gauche
       allocate (locR1)
       locR1%z=z
       Nullify(locR1%suiv) ; Nullify(locR1%prec)
       allocate (locR2)
       t=t+1
       locR2%y=y
       locR2%cr1=>locR1
       nullify(locR2%suiv)
       nullify(locR2%prec)
       allocate(chloc)
       chloc%x=x
       chloc%cr2=>locR2
       if (droite) then
!ajout à droite de temp
          chloc%prec=>temp
          nullify(chloc%suiv)
          temp%suiv=>chloc
       else
!ajout à gauche  de temp
          chloc%suiv=>temp
          nullify(chloc%prec)
          temp%prec=>chloc
       end if
!!!modif
          R3=>chloc
!       R3=>temp
!       print*,'dans range',R3%x
       return
  end subroutine rangeR3
!-----------------------------------------------
 subroutine rangeR2(y,z,R2,t)
use deftype
    implicit none
 
    Type(chaineR2), Pointer            :: R2
    real*8, intent(in)                   :: y,z
    integer, intent(inout)               :: t

    logical                         :: droite, droitedep

   Type(chaineR2), Pointer             :: chloc,temp
   Type(chaineR1), Pointer                :: locR1

    ! 
    !
    !
    ! inrterdit de rentrer dans cette routine si la chaine pointe sur null
    droite=.false.
    if (y>=R2%y) droite=.true.
    do while (associated(R2))
          temp=>R2
          if (y>R2%y+1.d-12) then
             if (droite) then
                R2=>R2%suiv
             else
                t=t+1
                call insertdroitey(y,z,R2)
                return
             end if
          else
             if (abs(y-R2%y)<1e-12) then 
                call rangeR1(z,R2%cr1,t)
                return
             else
                if (droite) then 
                   t=t+1
                   call insertgauchey(y,z,R2)
                   return
                else
                   R2=>R2%prec
                end if
             end if
          end if
       end do
!cas du bout de liste droite ou gauche
       allocate (locR1)
       locR1%z=z
Nullify(locR1%suiv) ; Nullify(locR1%prec)
       t=t+1
       allocate(chloc)
       chloc%y=y
       chloc%cr1=>locR1
       if (droite) then
!ajout à droite de temp
          chloc%prec=>temp
          nullify(chloc%suiv)
          temp%suiv=>chloc
       else
!ajout à gauche  de temp
          chloc%suiv=>temp
          nullify(chloc%prec)
          temp%prec=>chloc
       end if
!modif
       R2=>chloc!temp
       return
  end subroutine rangeR2

!-----------------------------------------------
 subroutine rangeR1(z,R1,t)
use deftype
    implicit none
 
    Type(chaineR1), Pointer            :: R1
    real*8, intent(in)                   :: z
    integer, intent(inout)               :: t

    logical                         :: droite, droitedep

   Type(chaineR1), Pointer             :: chloc,temp
    ! 
    !
    !
    ! inrterdit de rentrer dans cette routine si la chaine pointe sur null
    droite=.false.
    if (z>=R1%z) droite=.true.
    do while (associated(R1))
          temp=>R1
          if (z>R1%z+1.d-12) then
             if (droite) then
                R1=>R1%suiv
             else
                t=t+1
                call insertdroitez(z,R1)
                return
             end if
          else
             if (abs(z-R1%z)<1e-12) then 
                return
             else
                if (droite) then 
                   t=t+1
                   call insertgauchez(z,R1)
                   return
                else
                   R1=>R1%prec
                end if
             end if
          end if
       end do
!cas du bout de liste droite ou gauche
       t=t+1
       allocate(chloc)
       chloc%z=z
       if (droite) then
!ajout à droite de temp
          chloc%prec=>temp
          nullify(chloc%suiv)
          temp%suiv=>chloc
       else
!ajout à gauche  de temp
          chloc%suiv=>temp
          nullify(chloc%prec)
          temp%prec=>chloc
       end if
       R1=>chloc!temp
       return
  end subroutine rangeR1

subroutine insertgauchex(x,y,z,ch)
use deftype
    implicit none
    real*8, intent(in)                   :: x,y,z
    Type(chaineR3), Pointer               :: ch
!
    Type(chaineR3), Pointer             :: temp, locch
    type(chaineR2), Pointer            :: locR2
    type(chaineR1), Pointer            :: locR1
    ! 
    !
    !
!cas du bout de liste droite ou gauche
    temp=>ch%prec
       allocate (locch)
       locch%x=x
       allocate (locR2)
       locR2%y=y
       nullify(locR2%suiv)
       nullify(locR2%prec)
       allocate (locR1)
       locR1%z=z
       nullify(locR1%suiv)
       nullify(locR1%prec)
       locR2%cr1=>locR1
       locch%cr2=>locR2
       locch%suiv=>ch
       ch%prec=>locch
       locch%prec=>temp
       temp%suiv=>locch
!modif
       ch=>locch
       return
  end subroutine insertgauchex

 subroutine insertdroitex(x,y,z,ch)
use deftype
    implicit none
    real*8, intent(in)                   :: x,y,z
    Type(chaineR3), Pointer               :: ch
!
    Type(chaineR3), Pointer             :: temp, locch
    type(chaineR2), Pointer            :: locR2
    type(chaineR1), Pointer            :: locR1
    ! 
    !
    !
!cas du bout de liste droite ou gauche
    temp=>ch%suiv
       allocate (locch)
       locch%x=x
       allocate (locR2)
       locR2%y=y
       nullify(locR2%suiv)
       nullify(locR2%prec)
       allocate (locR1)
       locR1%z=z
       nullify(locR1%suiv)
       nullify(locR1%prec)
       locR2%cr1=>locR1
       locch%cr2=>locR2
       locch%suiv=>temp
       temp%prec=>locch
       locch%prec=>ch
       ch%suiv=>locch
!modif
       ch=>locch
       return
  end subroutine insertdroitex

!-----------------------------------------------
subroutine insertgauchey(y,z,ch)
use deftype
    implicit none
    real*8, intent(in)                   :: y,z
    Type(chaineR2), Pointer               :: ch
!
    Type(chaineR2), Pointer             :: temp, locch
    type(chaineR1), Pointer            :: locR1
    ! 
    !
    !
!cas du bout de liste droite ou gauche
    temp=>ch%prec
       allocate (locch)
       locch%y=y
       allocate (locR1)
       locR1%z=z
       nullify(locR1%suiv)
       nullify(locR1%prec)
       locch%cr1=>locR1
       locch%suiv=>ch
       ch%prec=>locch
       locch%prec=>temp
       temp%suiv=>locch
!modif
       ch=>locch
       return
  end subroutine insertgauchey

 subroutine insertdroitey(y,z,ch)
use deftype
    implicit none
    real*8, intent(in)                   :: y,z
    Type(chaineR2), Pointer               :: ch
!
    Type(chaineR2), Pointer             :: temp, locch
    type(chaineR1), Pointer            :: locR1
    ! 
    !
    !
!cas du bout de liste droite ou gauche
    temp=>ch%suiv
       allocate (locch)
       locch%y=y
       allocate (locR1)
       locR1%z=z
       nullify(locR1%suiv)
       nullify(locR1%prec)
       locch%cr1=>locR1
       locch%suiv=>temp
       temp%prec=>locch
       locch%prec=>ch
       ch%suiv=>locch
!modif
       ch=>locch
       return
  end subroutine insertdroitey

!-----------------------------------------------
!-----------------------------------------------
subroutine insertgauchez(z,ch)
use deftype
    implicit none
    real*8, intent(in)                   :: z
    Type(chaineR1), Pointer               :: ch
!
    Type(chaineR1), Pointer             :: temp, locch
    ! 
    !
    !
!cas du bout de liste droite ou gauche
    temp=>ch%prec
       allocate (locch)
       locch%z=z
       locch%suiv=>ch
       ch%prec=>locch
       locch%prec=>temp
       temp%suiv=>locch
!modif
       ch=>locch
       return
  end subroutine insertgauchez

 subroutine insertdroitez(z,ch)
use deftype
    implicit none
    real*8, intent(in)                   :: z
    Type(chaineR1), Pointer               :: ch
!
    Type(chaineR1), Pointer             :: temp, locch
    ! 
    !
    !
!cas du bout de liste droite ou gauche
      temp=>ch%suiv
       allocate (locch)
       locch%z=z
       locch%suiv=>temp
       temp%prec=>locch
       locch%prec=>ch
       ch%suiv=>locch
!modif
       ch=>locch
       return
  end subroutine insertdroitez




!======================
end module maillage
!======================
