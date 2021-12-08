module deftype

  implicit none



  Type chaines
     
     integer             :: i
     Type(chaines), Pointer :: prec,suiv

  End Type chaines
 Type chainej
     
     integer             :: j,k
     real*8              :: val
     Type(chainej), Pointer :: prec,suiv
  End Type chainej
!
  Type dblechaine
     
     integer             :: i
     Type(chainej), Pointer       :: lj
     Type(dblechaine), Pointer :: prec,suiv

  End Type dblechaine
!
  Type chaineR3
     real*8             :: x
     Type(chaineR2), Pointer       :: cr2
     Type(chaineR3), Pointer :: prec,suiv

  End Type chaineR3
 !
  Type chaineR2
     
     real*8             :: y
     Type(chaineR1), Pointer       :: cr1
     Type(chaineR2), Pointer :: prec,suiv

  End Type chaineR2
 !
  Type chaineR1
     integer            :: num, interieur
     real*8             :: z,phi
     Type(chaineR1), Pointer :: prec,suiv

  End Type chaineR1 
 Type morse
  
    integer                          :: i,j
    real(kind = 8)                   :: val

  End Type morse 
!
!----- Debut Type ------------------------------------------------------
  Type listechaine     
     Type(chainej), Pointer    :: liste
  End Type listechaine
!----- Fin Type --------------------------------------------------------
end module deftype
