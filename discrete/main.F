!===========================================================================
!
!  GOdMD, mutating residues
!  Pedro Sfriso
!
!==========================================================================
program discrete

    use commLine
    use geometry
    use geometryDP
    use txtUtils
    use stepPotentials
    use intList
    use paramSet
    use dims_utils
    use ls_rmsd
    use least_sq_fit
    use Structure
    use ran_mod
    use doubleDMD
    use ANM

    !
    IMPLICIT NONE
    !
#include "vars.h"


  call readCommLine(files, unit_i, unit_o, NFILES)
  call readInputParamSet(unit_i)
 !
 ! SetUp
 !
  str = loadStructurePDB(openFn(files, '-pdbin'))
  strTarg= loadStructurePDB(openFn(files, '-pdbtarg'))

 !
 call setBonds(str)
 !
 call  setStepPotentials(str)
 natom=str%natoms
 !
 ! Check if the system is a protein (minimo numero de atomos)
 if (natom < natomMin ) THEN
     write(unit_o,*)" Only systems above ", natomMin, " residues are accepted"
     stop
 end if
 !
 !
 write (unit_o, '(" Molecules: "i4," Residues:   ",i4," Atoms:  ",i6)') str%nmol, str%nres, natom
 write (unit_o, '(" Molecules: "i4," Residues:   ",i4," Atoms:  ",i6)') strTarg%nmol, strTarg%nres, strTarg%natoms
 !
 allocate (atoms(natom), PDBrnum(natom), res(natom), atp(natom), chain(natom),molnum(natom), &
           stepPtsDef(natom,natom),PDBTargetID(strTarg%natoms),stat = ioerr)
 call errorAllocmem(ioerr, 'NAtom')

 ! ==== ENTRA unit pairs, sale oid, natomCommon  ====
 !
 ! El allocate de oid hace que sea mejor dejarlo en main.
 !
 !  Lectura de la estructura final
 !
     unit_pairs= openFn(files, '-p1')
    !
    !  Vamos a leer la lista de pares
    !
    k=0
     readloop: DO
       READ(unit_pairs,*,IOSTAT=status)resOrig
       IF ( status /= 0 ) EXIT
       k=k+1
    END DO readloop
    natomCommon=k
    WRITE(unit_o,*)"NATOM COMMON", natomCommon
    ALLOCATE( oid(natomCommon),commonResidue(natomCommon))
    REWIND(unit_pairs)
    k=1
    readloop2: DO
       READ(unit_pairs,*,IOSTAT=status)resOrig
       IF ( status /= 0 ) EXIT
       call removesp(resOrig)
       commonResidue(k)=resOrig
       k=k+1
    END DO readloop2
    CLOSE(unit_pairs)
 ! Como guardar ya en oid el entero...
! ======================================================
    !
 allocate (rtarg(natomCommon), stat = ioerr)
 call errorAllocmem(ioerr, 'Coordinates')
 !
 ! PROBLEMA CON KB
 TEMP=TEMP
 !
 atoms=str%ats(:) %atomId
 res= str %ats(:) %resIdByAtom
 chain= str %ats(:) %chainId
 atp= str %ats(:) %atType
 molnum= str %ats(:) %molnum
 stepPtsDef=str %idata
 !
 DO i=1,natom
     call removesp(str %ats(i) %realResID)
     PDBrnum(i)= str %ats(i) %realResID
     call removesp(atoms(i) )
     call removesp(res(i) )
     call removesp(chain(i) )
     call removesp(atp(i) )
 END DO
 !
 oid=0
 k=1
   DO i=1,natom
        IF(ANY(commonResidue==PDBrnum(i)) ) THEN
            oid(k)=i
            k=k+1
        END IF
    END DO
    !
 !  Ahora decider que residuos del target son los comunes
 unit_pairs2= openFn(files, '-p2')
 k=1
   readloopTarget: DO
       READ(unit_pairs2,*,IOSTAT=status)resOrig
       IF ( status /= 0 ) EXIT
       call removesp(resOrig)
       commonResidue(k)=resOrig
       k=k+1
    END DO readloopTarget
 CLOSE(unit_pairs2)
 !
 DO i=1,strTarg%natoms
     call removesp(strTarg %ats(i) %realResID)
     PDBTargetID(i)= strTarg %ats(i) %realResID
 END DO
 !
 k=1
  DO i=1,strTarg%natoms
       IF(ANY(commonResidue==PDBTargetID(i)) ) THEN
      !     rspTarg(k)= strTarg % ats(i) %r
           rtarg(k) = SPtoDP(strTarg % ats(i) %r)
      !     write(unit_o,*)PDBTargetID(i)
           k=k+1
       END IF
   END DO


 ! Coordinates
 allocate ( r(natom), stat = ioerr)
 call errorAllocmem(ioerr, 'Coordinates')

!
    do i = 1, natom
        r(i) = SPtoDP(str % ats(i) % r)
    enddo
    !
 deallocate(commonResidue,PDBTargetID)
 !
 !END SETUP
 !
    call cpu_time(tinit)
    call readCommLine(files, unit_i, unit_o, NFILES)
    call programHeader(unit_o)
    call printFileSummary(files, unit_o)
    !
    call writeInputParamSet(unit_o)
    !
    ! Todos los atomos son CA es una asignacion a una lista entera
    evdw = 0.15
    rvdw = 2.0
    rhc  = 1.9 !1.6
    xm   = 0.012
    xsum = 1.0/xm !0.5 * (1./xm + 1./xm)
    distHC=(rhc + rhc)**2
    xmassa = xm*natom !sum(xm(1:natom))
    !
    !
    unit_ener = openFn(files, '-ener')
 !   STOP
    !
    ! Necesito tambien las distancias ( en A)
    allocate(dist_INI(natomCommon,natomCommon),dist_TARG(natomCommon,natomCommon))
    dist_INI=0.0
    dist_TARG=0.0
    ! CALL ERROR!!!
    DO j=1,natomCommon
        DO i=1,natomCommon
            dist_INI(i,j)=calcDistDP(r(oid(i)),r(oid(j)) )
            !
        END DO
    END DO
    !
      DO j=1,natomCommon
        DO i=1,natomCommon
            !
            dist_TARG(i,j)=calcDistDP(rtarg(i),rtarg(j))
        END DO
    END DO
    !
   call decideRCut4GO(natomCommon,r(oid(:)),rtarg,xm,rgyrPercent,mxRcut4go,rgyr,rcut2GO4)
   Write(unit_o,*) "PNS rgyr  ", rgyr, " CutOff ", sqrt(rcut2GO4)
   !
   IF(natom > 750 .or. real(natomCommon)/real(natom)<0.75 ) THEN
       errorAcceptable=1.5
       WRITE(unit_o,*)" Adjust minimum RMSd to ",errorAcceptable," A"
   END IF
    !
    !
    allocate( error_evo(error_evo_size,2 ), stat=ioerr)
    ! CALL ERROR!!
    error_evo=0.0
    !
    !
    allocate( ener_evo(ener_evo_size,2 ) , stat=ioerr)
    ! CALL ERROR!!
    ener_evo=0.0
    !
    !
    allocate (v(natom), stepPts(natom, natom),stat = ioerr)
    call errorAllocmem(ioerr, 'Setup/Natom')
    !
    !
   ! Estructura de la proteina
     do j = 2, natomCommon
        do i = 1, j-1
                ! Las particulas interaccionan
           stepPts(oid(i),oid(j)) = getStepPotFromDEF(GOENER,&
                sigmago,&
                dtol, &
                wellWidth,&
                dist_INI(i,j),&
                dist_TARG(i,j),&
                stepPtsDef(oid(i),oid(j))%tipint )
                stepPts(oid(i),oid(j)) % dmin = distHC
                stepPts(oid(i),oid(j)) % dego = 0.0
                stepPts(oid(i),oid(j)) % saltos =10
             !
              stepPts(oid(j),oid(i))=stepPts(oid(i),oid(j))
        enddo
    enddo
    !
    ! Vamos a fijar los potenciales de los residuos sin equivalente en el target
    !
    do j = 2, natom
        do i = 1, j-1
            IF(ANY(oid==i).and.ANY(oid==j)) THEN
                CYCLE
            ELSE
                  IF (stepPtsDef(i,j)%tipint == COVB) THEN
                     stepPts(i,j) = getStepPotForDummy(GOENER,sigmago,calcDistDP(r(i),r(j)),COVB)
                     stepPts(i,j) % dmin = distHC
                  ELSE
                     IF( calcDistDP(r(i),r(j)) < radioInteraccion ) THEN
                        ! Las particulas interaccionan
                        stepPts(i,j) = getStepPotForDummy(GOENER,sigmaFake,calcDistDP(r(i),r(j) ),FAKE)
                        stepPts(i,j) % dmin = distHC
                     END IF
                  END IF
                  stepPts(j,i)=stepPts(i,j)
            END IF
        enddo
    enddo

    deallocate (stepPtsDef)
    !
    !
    ! a por los cambios en ca dist i i+1
    DO j=1,natomCommon-1
          i=j+1
             IF ( stepPts(oid(j),oid(i)) % tipInt == COVB) THEN
                IF( (ABS(dist_INI(i,j)-dist_TARG(i,j)))/MIN(dist_INI(i,j),dist_TARG(i,j)) > 0.02_DBL) THEN
                  WRITE(unit_o,*) "UNO QUE SE MUEVE" , oid(i),oid(j),dist_INI(i,j),dist_TARG(i,j)
                  stepPts(oid(j), oid(i) ) % step(1) %r = REAL(MIN(dist_INI(i,j) ,dist_TARG(i,j)))*0.98
                  stepPts(oid(i), oid(j) ) % step(1) %r = REAL(MIN(dist_INI(i,j),dist_TARG(i,j)))*0.98
                  stepPts(oid(j), oid(i) ) % step(2) %r = REAL(MAX(dist_INI(i,j),dist_TARG(i,j)))*1.02
                  stepPts(oid(i), oid(j) ) % step(2) %r = REAL(MAX(dist_INI(i,j),dist_TARG(i,j)))*1.02
                END IF
             END IF
     END DO

    !
    call fixDistances(natom, natomCommon,oid,dist_INI,dist_TARG,ssectol,stepPts)
    write(unit_o,*) " He podido fijar las distancias"
    !
    allocate(discardedEner(natom,natom),stat = ioerr)
    ! CALL ERROR
    discardedEner=0.0
    !
    ! NM
    !
    ALLOCATE(evec(nevecs,natom),evals(nevecs),stat = ioerr)
    DO i=1,nevecs
       evec(i,:)=pointDP(0._DBL,0._DBL,0._DBL)
    ENDDO

    call ANM1(natom,r,nevecs,evec,evals)
    !
    DO i=1,nevecs
        WRITE (unit_o, '( "EIGENVALUE ",i3, " Valor ", f10.5  )')i,evals(i)
    END DO
    !
    call modulateDEgo(unit_o,natom,natomCommon,oid,evec,evals,r,dist_INI,stepPts)
    !
    DEALLOCATE(evec,evals)
    !
    Write(unit_o,*) " NMA-Metadynamica OK"

    ! precalcul parells
    allocate(blist(natom), nblist(natom),clist(natom), stat = ioerr)
    call errorAllocMem(ioerr, 'Pair lists')

    do i = 1, natom
        blist(i) = allocateintPList(5, ioerr)
        call errorAllocMem(ioerr, 'Bonded List')
    enddo

    do j = 2, natom
        do i = 1, j - 1
            if (stepPts(i, j) % tipInt .eq. COVB ) then
                stepPts(i, j) % active = .false.
                stepPts(j, i) % active = .false.
! j > i only
                IF(calcDistDP(r(i),r(j)) >20.0d0 ) THEN
                    write(*,*) " uuups", i,j
                !
                ! no poner ninguna restriccion entre particulas que estan muy lejos
                !
                    stepPts(i, j) % tipInt=0
                    stepPts(j, i) % tipInt=0
                ELSE
                blist(i) % nats = blist(i) % nats + 1
                blist(i) % iData(blist(i) % nats) = intData(j, blist(j) % nats, stepPts(i, j), xsum, 1.e15, 0.)
              !  write(unit_o,*) i,j
                contarBondedInt=contarBondedInt+1
                ENDIF
             endif
          enddo
    enddo
    write(unit_o,*)"bonded list OK"

    !write(unit_o,'(A4)') res(oid(:) )

   ! STOP
    do i = 1, natom
        nblist(i) = allocateintPList(INT(0.6*natom), ioerr)
        call errorAllocMem(ioerr, 'Non Bonded List')
    enddo

    nblist%nats = 0

    do i = 1, natom
        clist(i) = allocateCList(INT(natom), ioerr)
        call errorAllocMem(ioerr, 'contacts List')
    enddo
!
! ENTRA clist, natom, r, distHC, rcut :::: sale clist

    call activateClist(natom, r, distHC, rcontacts,clist)
    !
    write(unit_o,*) " contact list OK ", sum(clist(:)%nats)


    write(unit_o,*)"non-bonded list OK"

    call activateStepPot(stepPts, r,rtarg,oid ,rcut1,rcut2GO4, natom,natomCommon, nblist, xsum,radioInteraccion)
    !
    DEALLOCATE(dist_TARG,dist_INI)
    write(unit_o,*)"step-potentials activated"

    call writePairList(unit_o,natom,contarBondedInt,stepPts)

    deallocate (atp)
    !
    call cpu_time(tsetup)
    ! suma l'energia potencial de la conformacio inicial
    call calcEpot(natom,natom, r, stepPts, epotgo, epotfis, epot)

    rcm = calcCM(natom, r, xm)
    !
    unit_traj = openFn(files, '-trj')

    do i = 1, natom
        r(i) = r(i) - SPtoDP(rcm)
    enddo

    iev = 0
    call thermalize(seed, iev, natom, temp, xmassa, v, xm, ekin)
    !DIMS
    !
    allocate( w(natomCommon) ,rprev(natom), stat=ioerr)
    call errorAllocMem(ioerr, 'DIMS')

    call rmsd(natomCommon, r(oid(:)), rtarg(:), 1, U, center2, center1, error0, .FALSE.)
    call pdbsuperpos(natom, U , center1, center2,r)
    write(unit_o,*)"RMSD INICIAL ", error0

    call MCWeigth(natom,oid,r, rtarg,w,natomCommon)
    score = distanceCA(natom, oid,r,rtarg,w,natomCommon)
    write(unit_o,*)"SCORE ", score
    auxiliar=seed+iev
    scoreprev=1.E10
    rprev=r
    iacc=0

    call removesp(fileref)
    open(55,file=fileref)
    call writeSnapshot(55, 0, r, atoms, res, chain, PDBrnum, natom)
    close(55)
    WRITE(unit_traj,*)"# GODMD CRD file generated for ",fileref
    !


    ekin0 = ekin
    !
    temps = 0.
    allocate (tpart(natom), ipart(natom), toUpdate(natom), stat = ioerr)
    call errorAllocMem(ioerr, ' Colision times')

    Do i=1,Ener_evo_size
        Ener_evo(i,1)=epotgo
        Ener_evo(i,2)=i
    end do
    !
    saco=0
    slope=1.E10
    error_evo=1.E-10
    accepted=0
    totalTrials=0
    stopping=.FALSE.
    DRMSd=5.5E-5
    aux_acc=0.
    ibloc=0
    contarTodos=count(stepPts % nstep == 4 .and.stepPts % active)
    EnerCont=0
    enchufaMeta=.FALSE.
    slopeEner=1.
    contarSaltos=0
    aux_acc=1.0
    calcAuxAcc=.FALSE.
    nwritten=1
    contarSaltos=0
    unForceIt=.FALSE.
    cuantosForceIt=0
    finishMD=.FALSE.
    !
    !
    !----------------------------------------------------------------------------
    do while (.not.stopping)

    !
        axran=-ibloc-seed
        initiation=ran1(axran)
        ibloc=ibloc+1
        tacum = 0.
    !
        EnerCont=EnerCont+1
        IF (EnerCont > Ener_evo_size ) EnerCont=1
        Ener_evo(EnerCont,2) = epotgo
        Ener_evo(EnerCont,1)= real(ibloc)

        slopeEner=lreg(Ener_evo_size,Ener_evo)

        IF (slopeEner < 1.0E-1 .and. abs(slopeEner) < 2.5E-2) enchufaMeta=.TRUE.

        if( enchufaMeta ) THEN
            startMD=.TRUE.
            contarSaltos=0
            call metaDMD(natom, r,stepPts,discardedEner)
            call actualizeStepPot(stepPts,natom, nblist)
            enchufaMeta=.FALSE.
            saco = saco +1
            Ener_evo(:,2)=0.0
              do j = 2, natom - 1
                 do i = 1, j-1
                    IF (stepPts(i,j) %active) then
                        IF ( stepPts(i,j) %nstep==4 ) then
                            IF ( stepPts(i,j) % saltos==0 )then
                                contarSaltos=contarSaltos+1
                            end if
                        end if
                    end if
                end do
            end do

        END IF
        !
        !


        !----------------------------------------------------------------------------
        do while (tacum .lt. TSNAP)
            !
            call activateClist(natom, r, distHC, rcontacts,clist)
            !
            tacrect = 0
            do while (tacrect .lt. TRECT)
                taccorr = 0.
                iev = 0
                ierr = 0

                !----------------------------------------------------------------------------
                do while (taccorr .lt. TCORR)
                    ekin0 = calcEkin(v, xm, natom)
                    call colisioBond(blist, r, v, xm, natom)
                    call colisioNonBond(temps, nblist,clist, r, v, xm, xsum, ierr,TMIN, natom)
                    call DMDIntLoop(nblist, r, v, xm, tacact, tpart, ipart, TACT, TMIN, temps, iev, natom, toUpdate)
                    taccorr = taccorr + tacact
                enddo
                ! end do while(taccorr.lt.tcorr)-----------------------------------------------
                call thermalize(seed, iev, natom, temp, xmassa, v, xm, ekin)
                tacrect = tacrect + taccorr
            enddo
        !
         if ( StartMD.and..not.finishMD ) then ! Restraints do not apply up to 3 (arbitrary) time step to allow relaxation in case
         ! of a terrible structure
            score = distanceCA(natom,oid,r,rtarg,w,natomCommon)
            auxiliar=seed+iev
            totalTrials=totalTrials+1
            if (MCCheck(error,auxiliar, xbeta , scoreprev, score)) then
               r=rprev
            else
               rprev=r
               iacc=iacc+1
               scoreprev=score
               tacum=tacum+tacrect
               accepted=accepted+1
            endif
         else
               iacc=iacc+1
               tacum=tacum+tacrect
         endif

        write(unit_traj,'(10f8.3)')r
               call writeEnergies(unit_ener, unit_o, temps, epotgo, ekin0, natom,fractionDone,error,&
              aux_acc,xbeta,saco,sum(discardedEner(:,:)),count(stepPts(:,:) % saltos <0), &
              contarSaltos,.TRUE., .TRUE.,ierr)

        enddo
        ! end DO while(tacum.lt.tsnap) -------------------------------------------------

        IF(calcAuxAcc) aux_acc= accepted*1./totalTrials
      ! SELF_STOP

         errcont=errcont+1
         IF (errcont > error_evo_size ) errcont=1
         error_evo(errcont,2)=error
         error_evo(errcont,1)=real(ibloc)
         ! DECIDE TO STOP
         slope=abs(lreg(error_evo_size,error_evo))

         if ( (MOD(real(ibloc),5.0) == 0 .and.StartMD) .or.aux_acc<0.2) then
             xbeta= updateXbeta(xbeta,aux_acc,acceptance)
             accepted=1
             totalTrials=1
             calcAuxAcc=.TRUE.
         END IF

         !
         ! NO MD
         call AreInB(natom,stepPts,r,numEnB)

         fractionDone=Real(numEnB)/Real(contarTodos)

         call rmsd(natomCommon, r(oid(:) ),rtarg(:), 1, U, center2, center1, error,.FALSE.)

         call pdbsuperpos(natom, U , center1, center2,r)
         ! ways to stop
         !
         IF ( ibloc > NBLOC ) THEN
             stopping=.TRUE.
             exitStatus=1
         END IF

         !
         IF (slope < DRMSD .and. error < 0.4*error0 .and. .not.stopping) THEN
           IF(ibloc > error_evo_size .and. error > errorAcceptable .and. cuantosForceIt < 6) THEN
             IF(.not.unForceIt) error_evo(errcont,2)=1.E10
             call ForceIT(natom,natomCommon,oid,r, rtarg,w)
             xbeta=2*xbeta
             scoreprev=1.E10
             unForceIt=.TRUE.
             cuantosForceIt=cuantosForceIt+1
           END IF
         END IF


       IF (mod(ibloc,nskip)==0) THEN
      !    call writeSnapShot(unit_traj, nwritten, r, atoms, res, chain, rnum, natom)
      !
       !  write(unit_o,*) " errores ", ierr
         call calcEpot(natom,natom, r, stepPts, epotgo, epotfis, epot)
         !
         call writeEnergies(unit_ener, unit_o, temps, epotgo, ekin0, natom,fractionDone,error,&
              aux_acc,xbeta,saco,sum(discardedEner(:,:)),count(stepPts(:,:) % saltos <0), &
              contarSaltos,.TRUE., .TRUE.,ierr)
          write(unit_traj,'(10f8.3)')r
          nwritten=nwritten+1
      END IF

    enddo
    ! end DO iblock-----------------------------------------------------------------
    close(5)
    call cpu_time(tfin)
    write (unit_o, *)
    write (unit_o, '(" T I M I N G S ")')
    write (unit_o, '(" ============= ")')
    write (unit_o, '(" Setup: " f10.2," s")') tsetup - tinit
    write (unit_o, '(" Traj:  ",f10.2," s (",f10.2," ns/h, ",f7.2," ns/day)")') tfin - tsetup, &
    temps / (tfin - tsetup) * 3.6 / 1000., temps / (tfin - tsetup) * 3.6 / 1000. * 24
    write (unit_o, *)
    write (unit_o, *)
 !===============================================================================
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation Information                                   |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Timming                                                  |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Number of iterations              | IBLOC  |",I12  ," |")') IBLOC
   write (unit_o, '(" | Number of particles               | natom  |",I12  ," |")') NATOM
   write (unit_o, '(" | Total Time (minutes)             |tfin-tini|",f12.3," |")')  (tfin - tinit)/60.
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Structure                                                |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Number of particles in common  |natomCommon|",I12  ," |")') NATOMCOMMON
   write (unit_o, '(" | Changing distances cut-off        |rcut2GO4|",f12.3," |")') rcut2go4
   write (unit_o, '(" | Number of double min wells     |contarTodos|",I12," |")') contarTodos
   write (unit_o, '(" | Initial RMSD (A)                  | error0 |",f12.3," |")') error0
   write (unit_o, '(" | Final RMSD (A)                    |  error |",f12.3," |")') error
   write (unit_o, '(" | Fraction native contacts don  |fractionDone|",f12.3," |")') fractionDone
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Exit status                                |",I12  ," |")') exitStatus
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, *)
   write (unit_o,'("CHECK",2f8.3,X,I1,X,I4)') (error0-error)/error0, fractionDone, exitStatus, ibloc
!===============================================================================

end
!=============================================================================
!===============================================================================
 subroutine writeInputParamSet (unit_o)
     use paramSet
   integer, intent(IN) :: unit_o
   ! Copiem l'arxiu de parametres. Pendent format
   write (unit_o, *)
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | CALCULATION PARAMETERS                                   |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation settings                                      |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation Time (ps) (Nbloc x TSnap)       |",f12.3," |")') NBLOC * TSNAP / 1.e3
   write (unit_o, '(" | Output structure (fs)             | TSnap  |",f12.3," |")') TSNAP
   write (unit_o, '(" | Re-scoring target (fs)            | Trect  |",f12.3," |")') TRECT
   write (unit_o, '(" | Update velocities (fs)            | Tcorr  |",f12.3," |")') TCORR
   write (unit_o, '(" | Update Lists, collision times (fs)| Tact   |",f12.3," |")') TACT
   write (unit_o, '(" | Min. accepted colision time (fs)  | TMin   |",f12.8," |")') TMIN
   write (unit_o, '(" | Temperature (K)                   | Temp   |",6X,f6.2, " |")') TEMP
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | DIMS                                                     |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" |Acceptance                          |ACCEPTANCE|",f12.8," |")') acceptance
   write (unit_o, '(" |GO Energy   (kcal)                  |GOENER    |",f12.8," |")') GOENER
   write (unit_o, '(" |Changing distances well width(A)    |sigmago   |",f12.8," |")') sigmago
   write (unit_o, '(" |Tolerance between changes(A)        |dtol      |",f12.8," |")') dtol
   write (unit_o, '(" |Not-Changing distances well width(A)|wellwidth |",f12.8," |")') wellwidth
   write (unit_o, '(" |Dummy tolerance (A)                 |sigmaFake |",f12.8," |")') sigmaFAke
   write (unit_o, '(" |Factor scaling NMA coupling         |scalFactor|",f12.8," |")') scalFactor
   write (unit_o, '(" |Cut-off in Radius Gyration        | rgyrPercent|",f12.8," |")') rgyrpercent
   write (unit_o, '(" |Cut-off not changing GO interact.(A)|rcut1     |",f12.8," |")') rcut1
   write (unit_o, '(" |Minimum Energy of any well (kcal) | minEnerWell|",f12.8," |")') minEnerWell
   write (unit_o, '(" |Min time to eq. energy(ps)        |Ener_evo_size |"I5"    |")') Ener_evo_size
   write (unit_o, '(" |Number of ENM Eigenvectors          |NEVECS    |"I5  "    |")') nevecs
   write (unit_o, '(" |Fraction of energy to add at MD     |DEGO_BASE |",f12.8," |")') DEGO_BASE
   write (unit_o, '(" |Minimum acceptable RMSD (A)    |errorAcceptable|",f12.8," |")') errorAcceptable
   write (unit_o, '(" |Number of frames to skip            |NSKIP     |" I8    " |")') nskip
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Other                                                    |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Random generator seed                      |",7X,i5  " |")') seed
   write (unit_o, '(" | IDAB, IGOAB, IWR, ISOLV                    |",4X,4i2," |")') IDAB, IGOAB, IWR, ISOLV
   write (unit_o, '(" ------------------------------------------------------------")')
 end subroutine writeInputParamSet
!===============================================================================
