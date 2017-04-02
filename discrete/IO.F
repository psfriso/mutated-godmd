!===============================================================================   
 subroutine readCommLine (files, unit_i, unit_o, NFILES)
 USE commLine
  integer, intent(IN) :: NFILES
  integer, intent(OUT) :: unit_i, unit_o
  type(commLineOption), intent(INOUT) :: files(NFILES) 
!
  call inputArgs(files)
  unit_i = openFn(files, '-i')
  if (fileName(files,'-o').ne.'log') then
        unit_o = openFn(files, '-o')
  else
        unit_o = 6
  endif
 end subroutine readCommLine
!===============================================================================   
 subroutine programHeader(unit_o)
   integer unit_o
#ifdef OMP
   integer omp_get_max_threads
#endif
   write (unit_o, *) "================================================="
   write (unit_o, *) "=                                               ="
   write (unit_o, *) "=                    GOdMD                      ="
   write (unit_o, *) "=                                               ="
   write (unit_o, *) "= P Sfriso, A. Emperador, JL. Gelpi, M.Orozco   =" 
   write (unit_o, *) "=                                               ="
   write (unit_o, *) "=                  (c) 2013                     =" 
   write (unit_o, *) "================================================="
   write (unit_o, *)
#ifdef OMP
   write (unit_o, '(" Running OpenMP version on ", i3, " processors")') omp_get_max_threads()
   write (unit_o, *)
#endif
 end subroutine programHeader
!===============================================================================   
 subroutine errorAllocmem (ioerr, text)
  integer ioerr
  character(len=*) text
  if (ioerr.ne.0) then
    write (0, '("Error allocating memory",a30)') text
    stop 1
  endif
 end subroutine errorAllocmem
!===============================================================================   
 subroutine writeSnapshot (unit_traj, ibloc, r, atom, res, chain, rnum, natom)
 use geometryDP
  integer, intent(IN) :: natom, unit_traj, ibloc
  type(pointDP), intent(IN) :: r(natom)
  character(len=4), intent(IN) :: atom(natom), res(natom)
  character(len=1), intent(IN) :: chain(natom)
  character(len=5), intent(IN) :: rnum(natom)
  integer i, j
     write (unit_traj, '("MODEL",8X,I4)') ibloc
     do i = 1,natom
      j=i
      write (unit_traj, '("ATOM",2X,I5,1X,A4,1X,A3,1X,A1,A5,3X,3F8.3)') i, atom(j), res(j), chain(j), rnum(j), r(i)
     enddo
     write (unit_traj, '("ENDMDL")')
 end subroutine writeSnapshot
!===============================================================================   
! subroutine writeRestartFiles(files, r, v, RSTV, natom, NFILES)
! use geometry
! use geometryDP
! use commline
!   integer, intent(IN) :: natom, RSTV, NFILES
!   type(pointDP), intent (IN) :: r(natom), v(natom)
!   type(commLineOption), intent(IN) :: files(NFILES)
!
!   type(point), allocatable :: rsp(:) ! Single precision version of r for binary I/O
!   integer unit_rst, unit_rstv, i
!  
!   allocate (rsp(natom))
!   do i = 1, natom
!        rsp(i) = DPtoSP(r(i))
!   enddo
!   unit_rst = openFn(files, '-rst')
!   write (unit_rst) natom
!   write (unit_rst) (rsp(i), i = 1, natom)
!   close(unit_rst)
!!
 !  if (RSTV .eq. 1) then
 !       unit_rstv = openFn(files, '-rstv')
 !       write (unit_rstv) natom
 !       write (unit_rstv) (v(i), i = 1, natom)
 !       close(unit_rstv)
 !   endif
 !   deallocate(rsp)
 !   !
!end subroutine writeRestartFiles
!===============================================================================   
subroutine writeEnergies(unit_ener, unit_o, temps, epotgo, ekin0, natom,fractionDone,error,&
    aux_acc,xbeta,saco,discardedEner,wrong,encero,WENER, WLOG,ierr)
 ! EN MAIN:
 !           call writeEnergies(unit_ener, unit_o, iev, ierr, temps,&
 !                 epotgo, epot, ekin0, natom,fractionDone,error,aux_acc,xbeta ,.true.,.true.)
 !   
    use constants
    integer, intent(IN) :: unit_ener, unit_o, natom, saco,encero
    INTEGER, INTENT(IN) :: wrong,ierr
    real, intent(IN) ::  epotgo, ekin0, xbeta,aux_acc,fractionDone
    real*8,Intent(IN):: error
    real*8, intent(IN) :: temps
    real, INTENT(IN) :: discardedEner
    logical, intent(IN) :: WENER, WLOG
    real :: Kb= 8.314 ! J/(mol K)

    character(len = 92), parameter :: &
    fmtEner = '("Time ",f8.2, " Epot",f10.3, " Temp(K)",f8.2 ," Dener ",f12.2," Salto ",I4," Error ",I3)'
    character(len = 105), parameter :: fmtLog = &
    '(" Time(ps): ",f8.2, " Epot(kcal/mol): ",f10.3, "  Done(%): ",f4.2," RMSd(A):",f8.3,f8.2,f5.2,X,I3,X,2I4)'
      !    
      ! FALTA Kb
    if (WENER) &
    write (unit_ener, fmtEner)temps * 1.e-3, epotgo, ekin0/1.5/natom/Kb ,discardedEner,saco,wrong
    if (WLOG) & 
    write (unit_o, fmtLog) temps * 1.e-3, epotgo,fractionDone,error,xbeta,aux_acc,saco,encero,ierr
    
end subroutine writeEnergies
!===============================================================================   
 SUBROUTINE writePairList(unit_o,natom,contarBondedInt,stepPts)
      use stepPotentials
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: unit_o
      INTEGER, INTENT(IN) :: natom
      INTEGER, INTENT(IN) :: contarBondedInt
      TYPE( stepPotInt), INTENT (IN) :: stepPts(natom,natom)
      !
    write (unit_o, '(" Initial Pairs list")')
    write (unit_o, '(" ==================")')
    write (unit_o, '(" Total:              ",i9)') natom * (natom - 1) / 2
    write (unit_o, '(" Bonded:             ",i9)') contarBondedInt
    write (unit_o, '(" Non Bonded:         ",i9)') count(stepPts % active)
    write (unit_o, '("   GO 1 SQ WELL      ",i9)') count(stepPts % nstep .eq. 2 .and.stepPts%active.and.stepPts%tipInt==SSEC)
    write (unit_o, '("   GO 2 SQ WELL      ",i9)') count(stepPts % nstep .eq. 4 .and.stepPts%active)
    write (unit_o, '(" Fake:     ",i9)') count(stepPts % tipInt .eq. 3 .and.stepPts%active)
    !
    write (unit_o, *)
    write (unit_o, '(" System setup completed")')
    write (unit_o, *)
  END SUBROUTINE  writePairList
 !=================================================================================