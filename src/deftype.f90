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
