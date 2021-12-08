module operateur


contains




!================================================================================================================
  subroutine laplacien ( lch,dx,dy,dz,Nx,Ny,Nz,nk )
!================================================================================================================
use deftype

implicit none
! pour inserer dans la liste

type (listechaine), dimension(:), intent(inout) :: lch
        real (kind=8), intent(in) :: dx,dy,dz
        integer, intent(in) :: Nx,Ny,Nz
        integer, intent(inout) :: nk

        integer :: i,j,k,l

print*,'debut de initialisation'
       do i=1, Nx*Ny*Nz
          allocate (lch(i)%liste)
           lch(i)%liste%j=i
           lch(i)%liste%val=2/dx**2+2/dy**2+2/dz**2
           nullify(lch(i)%liste%suiv)
           nullify(lch(i)%liste%prec)
        end do
        nk=Nx*Ny*Nz
        l=0
        do k=1,Nz
           do j=1,Ny
              do i=1,Nx
                 l=l+1!l=i+(j-1)*Nx+(k-1)*Nx*Ny
                 if (i<Nx) call ajoutval(l+1,-1/dx**2,lch(l)%liste,nk)
                 if (i>1) call ajoutval(l-1,-1/dx**2,lch(l)%liste,nk)
                 if (j<Ny) call ajoutval(l+Nx,-1/dy**2,lch(l)%liste,nk)
                 if (j>1) call ajoutval(l-Nx,-1/dy**2,lch(l)%liste,nk)
                 if (k<Nz) call ajoutval(l+Nx*Ny,-1/dz**2,lch(l)%liste,nk)
                 if (k>1) call ajoutval(l-Nx*Ny,-1/dz**2,lch(l)%liste,nk)
              end do
           end do
        end do
 
  end subroutine laplacien
!================================================================================================================
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


    end do
!
  end subroutine desallouelistebis
!=======================================================================







end module operateur
