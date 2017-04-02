 MODULE intList
 use stepPotentials
 
 type intData
    integer patnum, simp ! patnum: num atom de la parella, simp: index en llista de la parella
    type(stepPotInt) :: stepPt
    real :: xsum
    real*8 :: timp
    real :: deltak    
 end type intData
 
 type intpList
    integer nats
    type(intData), pointer :: iData(:)
 end type intpList
 
 type closeones
     integer patnum, simp
     real dmin
 endtype closeones
 
 type overlapscheck
     integer nats
     type(closeOnes), pointer :: vecinos(:)
 end type overlapscheck
     
 
 CONTAINS
 
 function allocateintPList(natom, ioerr) result (pl)
 integer, intent(IN) :: natom
 integer, intent(OUT) :: ioerr
 type (intpList) pl
  allocate (pl%iData(natom), stat=ioerr)
  pl%nats=0
 end function allocateintPList
 
  function allocateCList(natom, ioerr) result (pl)
 integer, intent(IN) :: natom
 integer, intent(OUT) :: ioerr
 type (overlapscheck) pl
  allocate (pl%vecinos(natom), stat=ioerr)
  pl%nats=0
 end function allocateCList
 
 END MODULE intList
