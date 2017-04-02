    !
! Authors : Josep Ll. Gelpi, Agusti Emperador
! Subject : Structure module
! Revision: $id$
!
MODULE Structure
    USE Constants
    USE geometry
    USE stepPotentials
    !
    type pairlist
        integer npairs
        integer, pointer :: list(:,:)
    end type pairlist

    type Fragment
        integer ini, fin
    end type Fragment

    type atomData
        character(len = 4) atomId
        character(len = 4) resIdbyAtom
        character(len = 5) realResID !CON INSERCIONES
        integer resNum ! ref to residue num per atom, unique values
        integer resNumPDB ! ref to PDB residue num per atom
        integer molNum ! ref to mol num per atom
        character(len = 4) atType
        character(len = 1) chainId
        type(point) r
    end type atomData

    type residueData
        character(len = 4) resId
        integer resNumPDB
        character(len = 1) chainId
        integer ini, fin, ica, molres ! atom pointers per residue
    end type residueData
    
    
    type struc
        integer natoms, nres, nmol
        type(atomData), pointer :: ats(:)
        type(residueData), pointer :: res(:)
        integer, pointer :: bonds(:,:)
        real, pointer :: distat2(:,:)
        type (pairList) cov
        type (Fragment), pointer :: mols(:) ! ini and fin RESIDUE per fragment
        integer moltype
        type(intDataTop), pointer :: iData(:,:)
    end type struc
    
!===============================================================================
!   TIPOS DE DATOS HEREDADOS 
!===============================================================================
   type atom
     character(LEN=4) :: atomId
     character(LEN=4) :: atType
     character(LEN=1) :: element
     integer resNum
     integer ind
   end type atom
!
 type residue
   character(LEN=4) :: resId
   character(LEN=30) :: resName
   type (atom) atoms(MAXATPERRES)
   integer :: bonds (MAXATPERRES, MAXATPERRES)
   integer nAtoms
   integer ind
 end type residue
!
  type residueLibrary
    type (residue) :: residues(MAXRES)
    integer nRes
  end type residueLibrary
!
 type atType
   character(len=4) potId
   integer ind
   real qq, gfree, vol, evdw, rvdw, rhc, mas, xlamb
 end type atType
 
 type ffprm
   integer ntypes 
   type(atType) :: types(MAXATYPES)
 end type ffprm

 character*30, parameter :: potFmt = '(a4,7f8.3)'
 type(atType), parameter :: atTypeNull = atType('',0,0.,0.,0.,0.,0.,0.,0.,3.5)
!  

!===============================================================================
CONTAINS
!===============================================================================
!  FUNCIONES HEREDADAS
!===============================================================================
!================================================================================
 function emptyResidue () result (rr)
   type(residue) rr
   rr%resId=''
   rr%resName=''
   rr%atoms = atom('','','',0,0)
   rr%bonds = NULL
 end function emptyResidue
!================================================================================
 function getResidue (rl, resId) result (rr)
 use txtUtils
  type(residue) rr
  type(residueLibrary) rl
  character(LEN=*) resId
  integer i
  
  rr=emptyResidue()
  i=1
  do while (i.lt.rl%nRes.and..not.eqId(rl%residues(i)%resId,resId))
   i=i+1
  end do
  if (eqId(rl%residues(i)%resId,resId)) then
   rr = rl%residues(i)
  endif
  
 end function getResidue
!===============================================================================
 function getAtomFromResidue (rres,atomId) result(at)
 use txtUtils
  type(atom) at
  type(residue), intent(IN) :: rres
  character(LEN=*) atomId
  integer i
  
  i=1
  do while (i.le.rres%nAtoms.and..not.eqId(rres%atoms(i)%atomId,atomId)) 
   i=i+1
  enddo
  if (eqId(rres%atoms(i)%atomId,atomId)) then
   at = rres%atoms(i)
  else
   at = atom('','','',0,0)
  endif
  
 end function getAtomFromResidue
!================================================================================
 function getAtomFromLib (rl, resId, atomId) result(at)
 use txtUtils
  type(atom) at
  type(residue) rr
  type(residueLibrary) rl
  character(len=*) resId,atomId
  rr = getResidue(rl, resId)
  if (eqId(rr%resId,resId)) then
   at = getAtomFromResidue(rr,atomId)
  else
   at = atom('','','',0,0)
  endif
 end function getAtomFromLib
!================================================================================
 function getAtTypeFromLib (rl, resId, atomId) result(ty)
  type(residueLibrary) rl
  type(atom) at
  character(*) resId, atomId
  character(len=4) ty
  at = getAtomFromLib(rl, resId, atomId)
  ty = at%atType
 end function getAtTypeFromLib
!===============================================================================
!

!===============================================================================
    function findAtominRes(str, ires, atomId) result (p)
        use txtUtils
        integer p
        type(struc), intent(IN) :: str
        integer, intent(IN) :: ires
        character(*), intent(IN) :: atomId
        !
        p = str % res(ires) % ini
        do while (p .lt. str % res(ires) % fin.and..not.eqId(str % ats(p) % atomId, atomId))
            p = p + 1
        enddo
        if (.not.eqId(str % ats(p) % atomId, atomId)) p = 0

    end function findAtominRes
    !===============================================================================
    function allocatePairList(npairs) result (pl)
        type(pairList) pl
        integer, intent(IN) :: npairs
        allocate (pl % list(npairs, 2))
        pl % npairs = npairs
        pl % list = 0
    end function allocatePairList
    !===============================================================================
    function allocateStructureAtoms(natoms) result (st)
        type(struc) st
        integer, intent(IN) :: natoms
        allocate (st % ats(natoms))
        allocate (st % bonds(natoms, natoms), st % distat2(natoms, natoms), st % iData(natoms, natoms))
        st % ats % atomId = ''
        st % ats % resIdbyAtom = ''
        st % ats % chainId = ''
        st % ats % realResID = ''
        st % ats % resnum = 0
        st % ats % resnumPDB = 0
        st % ats % molnum = 0
        st % ats % r = point(0, 0, 0)
        st % ats % atType = ''
        st % bonds = NULL
        st % distat2 = 0.
        st % moltype = PROT
        st % natoms = natoms
        st % iData = intDataTop(0,NULL,0.,0.)
    end function allocateStructureAtoms
    !===============================================================================
    subroutine allocateStructureResidues(st, nres)
        type(struc) st
        integer, intent(IN) :: nres
        allocate(st % res(nres))
        st % res % resId = ''
        st % res % resNumPDB = 0
        st % res % chainId = ''
        st % res % molres = 0
        st % res % ini = 99999
        st % res % fin = 0
        st % res % ica = 0
        st % nres = nres
    end subroutine allocateStructureResidues
    !===============================================================================
    subroutine allocateStructureMols(st, nmol)
        type(struc) st
        integer, intent(IN) :: nmol
        allocate(st % mols(nmol))
        st % nmol = nmol
        st % mols % ini = 99999
        st % mols % fin = 0
    end subroutine allocateStructureMols
 
    !===============================================================================
    function loadStructurePDB(unt) result (st)
        use txtUtils
        type(struc) st
        integer, intent(IN) :: unt
        integer i, j, nmol, nres, natoms
        character(len = 1) refch
        character(len= 5):: refres
        character(len = 80) str
        real:: aux
        !
        natoms = 0
        nmol = 0
     !   insertion=' '
        !
        refch = ''
        10 read(unt, '(a80)', end = 20) str
        ! ELIMINEM EL H per posicio, pendent repasar canvis de format
        if (validAtomLine(str)) natoms = natoms + 1
        if (terStr(str).or.(str(22:22) .ne. refch).or.nmol .eq. 0) then
            nmol = nmol + 1
            refch = str(22:22)
        endif
        goto 10
        20 continue
        st = allocateStructureAtoms(natoms)
        call allocateStructureMols(st, nmol)
        !
        rewind(unt)
        !
        i = 0
        natoms = 0
        nmol = 0
        nres = 0
        refch = ''
        refres = ''
        30 read(unt, '(a80)', end = 40) str
        if (terStr(str).and.nmol .gt. 0) then
            if (st % mols(nmol) % fin .ne. 0) nmol = nmol + 1
        endif
        if (validAtomLine(str)) then
            !       write (6,*) str
            i = i + 1
            ! (13X,A4,1X,A3,1X,A1,1X,A4,4X,3F8.3,2f6.2)
           ! ( A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2 )!
            read (str,'(12X,A4,1X,A3,1X,A1,A5,3X,3F8.3,2f6.2)', end = 40) st % ats(i) % atomId, st % ats(i) % resIdbyAtom, &
            st % ats(i) % chainId,  st % ats(i) % realResID, &
            st % ats(i) % r % x, st % ats(i) % r % y, st % ats(i) % r % z,aux
           ! call removesp(st % ats(i) % realResID)
           ! st % ats(i) % realResID=st % ats(i) % realResID//insertion
           ! Write (*,'(13X,A3,1X,A4,A1,1X,A4,4X,3F8.3,2f6.2)') st % ats(i) % atomId, st % ats(i) % resIdbyAtom, &
           ! st % ats(i) % chainId, st % ats(i) % realResID, &
           ! st % ats(i) % r % x, st % ats(i) % r % y, st % ats(i) % r % z,aux
            if (refch .ne. st % ats(i) % chainId.or.nmol .eq. 0) then ! salvem PDBs sense cadena
                if (nmol .eq. 0) then
                    nmol = nmol + 1
                elseif (st % mols(nmol) % fin .ne. 0) then
                    nmol = nmol + 1
                endif
            endif
      
            !
            st % ats(i) % molNum = nmol
            if (refch .ne. st % ats(i) % chainId.or.refres .ne. st % ats(i) % realResID) nres = nres + 1
            refch = st % ats(i) % chainId
            refres = st % ats(i) % realResID
            st % ats(i) % resNum = nres
            st % mols(nmol) % ini = min(st % ats(i) % resNum, st % mols(nmol) % ini)
            st % mols(nmol) % fin = max(st % ats(i) % resNum, st % mols(nmol) % fin)
            natoms = natoms + 1
        endif
        goto 30
        40 continue
        ! cas de TER doblat al final
        do while (st % mols(nmol) % fin .eq. 0)
            nmol = nmol - 1
        enddo
        st % nmol = nmol
        do i = 1, st % natoms - 1
            do j = i + 1, st % natoms
                st % distat2(i, j) = calcDist2(st % ats(i) % r, st % ats(j) % r)
                st % distat2(j, i) = st % distat2(i, j)
            enddo
        enddo
        !
        call allocateStructureResidues(st, nres)
        !
        do i = 1, st % natoms
            st % res(st % ats(i) % resNum) % resId = st % ats(i) % resIdbyAtom
            st % res(st % ats(i) % resNum) % resNumPDB = st % ats(i) % resNumPDB
            st % res(st % ats(i) % resNum) % chainId = st % ats(i) % chainId
            st % res(st % ats(i) % resNum) % ini = min(i, st % res(st % ats(i) % resNum) % ini)
            st % res(st % ats(i) % resNum) % fin = max(i, st % res(st % ats(i) % resNum) % fin)
            st % res(st % ats(i) % resNum) % molres = st % ats(i) % molNum
            if (eqId(st % ats(i) % atomId, 'CA')) st % res(st % ats(i) % resNum) % ica = i
        enddo
        ! PENDENT determinar tipus de molecula
        st % molType = PROT
    end function loadStructurePDB
 !===============================================================================
 ! prepares newStr as two succesive loadStructurePDB
    function mergeStructures(str1, str2 ) result (newStr)
        use geometry
        type(struc), intent(IN) :: str1, str2
        type(struc) newStr
        integer i, j
        !
        newStr = allocateStructureAtoms(str1 % natoms + str2 % natoms)
        !
        newStr % ats(1:str1 % natoms) = str1 % ats(1:str1 % natoms)
        !
        newStr % distat2(1:str1 % natoms, 1:str1 % natoms) = str1 % distat2(1:str1 % natoms, 1:str1 % natoms)
        !
        newStr % ats(1 + str1 % natoms:newStr % natoms) = str2 % ats(1:str2 % natoms)

        newStr % ats(1 + str1 % natoms:newStr % natoms) % resnum = str2 % ats(1:str2 % natoms) % resnum + str1 % nres
        newStr % ats(1 + str1 % natoms:newStr % natoms) % molnum = str2 % ats(1:str2 % natoms) % molnum + str1 % nmol
        !
        newStr % distat2(1 + str1 % natoms:newStr % natoms, 1 + str1 % natoms:newStr % natoms) = &
        str2 % distat2(1:str2 % natoms, 1:str2 % natoms)
        !
        do i = 1, str1 % natoms
            do j = 1, str2 % natoms
                newStr % distat2(i, j + str1 % natoms) = 1.E10
                newStr % distat2(j + str1 % natoms, i) = newStr % distat2(i, j + str1 % natoms)
            enddo
        enddo
        !
        call allocateStructureResidues(newStr, str1 % nres + str2 % nres)
         newStr % res(1:str1 % nres) = str1 % res(1:str1 % nres)
  !      
         newStr % res(1 + str1 % nres:newStr % nres) = str2 % res(1:str2 % nres)
        !
        newStr % res(1 + str1 % nres:newStr % nres) % ini = str2 % res(1:str2 % nres) % ini + str1 % natoms
        newStr % res(1 + str1 % nres:newStr % nres) % fin = str2 % res(1:str2 % nres) % fin + str1 % natoms
        newStr % res(1 + str1 % nres:newStr % nres) % ica = str2 % res(1:str2 % nres) % ica + str1 % natoms
        newStr % res(1 + str1 % nres:newStr % nres) % molres = str2 % res(1:str2 % nres) % molres + str1 % nmol
        !
        call allocateStructureMols(newStr, str1 % nmol + str2 % nmol)
        newStr % mols(1:str1 % nmol) = str1 % mols(1:str1 % nmol)
        !
        newStr % mols(1 + str1 % nmol:newStr % nmol) % ini = str2 % mols(1:str2 % nmol) % ini + str1 % nres
        newStr % mols(1 + str1 % nmol:newStr % nmol) % fin = str2 % mols(1:str2 % nmol) % fin + str1 % nres
        ! PENDENT determinar tipus de molecula
        newStr % molType = COMPLEX
    end function mergeStructures
    !===============================================================================
    !===============================================================================
   logical function validAtomLine(str)
       character(*), intent(IN) :: str
       character(len=3), DIMENSION(33)::aa
     !!  aa=(/'THR','SER','PRO','GLU','ASP','CYS','HIS','CYG','HIE',&
      !      'LYS','ARG','ALA','VAL','ILE','ASN','GLN',&
      !      'PHE','GLY','TRP','MET','LEU','TYR','CYX','HIP',&
      !      'CYN','CPR','LYP','GLH','ASH','NAL','MSE','YCM','MHO'/)
      !      validAtomLine=.FALSE.
       ! unIfaqui(str(1:4).eq.'ATOM'.and.ANY(aa == str(18:21))) THEN  
      ! IF(ANY(aa == str(18:21))) THEN        
            validAtomLine = ((str(14:16).eq.'CA ').or.(str(13:15).eq.'CA '))
      ! ENDIF
    end function validAtomLine
    !==============================================================================
    logical function terStr(str)
        character(*), intent(IN) :: str
        terStr = (str(1:3) .eq. 'TER')
    end function terStr
    !===============================================================================
    subroutine assignAtType(str, resLib)
!        use resLibrary
        use txtUtils
        type(struc), intent(INOUT) :: str
        type(residueLibrary), intent(IN) :: resLib
        type(atom) att
        integer i

        do i = 1, str % natoms
            att = getAtomFromLib(resLib, str % ats(i) % resIdbyAtom, str % ats(i) % atomId)
        enddo
    end subroutine assignAtType
    !===============================================================================
    subroutine setBonds(str)
!        use resLibrary
        type(struc), intent(INOUT) :: str
        integer i, j

        if (str % moltype .eq. PROT.or.str % molType .eq. COMPLEX) then
            call setPeptideBonds(str) ! pendent fer crossbonds opcional
        endif

        call prepBondList(str)

        do i = 1, str % natoms - 1
            do j = i + 1, str % natoms
                str % bonds(j, i) = str % bonds(i, j)
            enddo
        enddo

    end subroutine setBonds
    !===============================================================================
    subroutine setPeptideBonds(str)
        type(struc), intent(INOUT) :: str
        integer nr
        !do nr = 1, str % nres - 1
        ! !   WRITE (*,*) "HOLA +1",str % res(nr) % molres,str % res(nr + 1) % molres
        !    if (str % res(nr) % molres == str % res(nr + 1) % molres) then
        !        ! Esto significa que pertenecen a la misma molecula (aun no esta activado
        !        ! el crossb )
        !        str % bonds(str % res(nr) % ica, str % res(nr + 1) % ica) = COVB ! CA-CA
        !  !      write(*,*) " ENLACE",covb,nr
        !    endif
        !enddo
        
        DO nr = 1, str %natoms-1
            str % bonds(nr, nr+1) = COVB
        END DO
        
    end subroutine setPeptideBonds
  !===============================================================================
    subroutine prepBondList(str)
        type(struc), intent(INOUT) :: str
        integer i, j
        ! PEDRO
        str % cov = allocatePairList(count(str % bonds .eq. COVB))
        str % cov % npairs = 0
        do i = 1, str % natoms - 1
            do j = i + 1, str % natoms
                if (str % bonds(i, j) .eq. COVB) then
                    str % cov % npairs = str % cov % npairs + 1
                    str % cov % list(str % cov % npairs, 1) = i
                    str % cov % list(str % cov % npairs, 2) = j
                endif
            enddo
        enddo
        
   !    Write(*,*)"NPAIRS ",  str % cov % npairs, count(str % bonds .eq. COVB)
        
    end subroutine prepBondList
    !===============================================================================
    subroutine setStepPotentials(str)
        !use potentials
        type (struc), intent(INOUT) :: str
        integer :: i,j
        call setStepPotCOV(str)
        call setStepPotSS(str)
        
         do i= 2, str % natoms
            do j = 1, i-1
                str % idata (j,i) = str % idata (i,j)
              !  write (*,*) "*****", i, j, str % iData(i, j)
            enddo
        enddo
        
    end subroutine setStepPotentials
    !===============================================================================
    subroutine setStepPotCOV(str)
        type (struc), intent(INOUT) :: str
        integer i, j, k
        !
        ! CAMBIO DISTANCIA
        !
        ! Para todos los enlaces covalentes
        do i = 1, str % cov % npairs
            j = str % cov % list(i, 1)
            k = str % cov % list(i, 2)
           ! str % iData(j, k) = intDataTop(getstepPotDefIndex(stL, 'COVB'), COVB, sqrt(str % distat2(j, k))*0.95, 0.)
           str % iData(j, k) = intDataTop(1, COVB, sqrt(str % distat2(j, k))*0.95, 0.)
 
            !      write (*,*) "+++++", i, j, str % iData(j, k)
            str % iData(k, j) = str % iData(j, k)
        enddo
    end subroutine setStepPotCOV
    !===============================================================================
    subroutine setStepPotSS(str)
        type (struc), intent(INOUT) :: str
        integer i, j
        
        do i= 2, str % natoms
            do j = 1, i-1
                IF ( str % iData(i, j) %id /= 1 ) &
                str % iData(i, j) = intDataTop(2, SSEC, sqrt(str % distat2(i, j)), 0.)
        !        write (*,*) "*****", i, j, str % iData(i, j)
            enddo
        enddo
    end subroutine setStepPotSS

    function writeAtomId(at) result (txt)
        type(atomData) at
        character(len = 15) txt
        write (txt, '(a4,1x,a1,i4,1x,a4)') at % resIdByAtom, at % chainId, at % resNumPDB, at % atomId
    end function writeAtomId
    !===============================================================================
    function writeResidueId(rr) result (txt)
        type(residueData) rr
        character(len = 10) txt
        write (txt, '(a4,1x,a1,i4)') rr % resId, rr % chainId, rr % resNumPDB
    end function writeResidueId
    !===============================================================================
END MODULE Structure
