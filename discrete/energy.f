! ===============================================
 SUBROUTINE activateClist(natom, r, distHC, rcontacts,clist)
        USE geometryDP
        USE intList
        !
        IMPLICIT NONE
        !
        TYPE(overlapscheck), DIMENSION(natom), INTENT(INOUT) :: clist
        INTEGER, INTENT(IN) :: natom
        TYPE (pointDP) , DIMENSION(natom), INTENT(IN) :: r
        REAL*8, INTENT(IN) :: rcontacts
        REAL, INTENT(IN) :: distHC
        !
        integer :: i, j
        !
        clist%nats = 0
        do i=1,natom
            clist(i)%vecinos=closeOnes(0,0,0.0)
        end do
        !
        ! vamos a generar la lista de vecinos para sacar posibles overlaps
        !
        do j = 3, natom
            do i = 1, j - 2 ! el j-1 ya lo tiene en cuenta el bonded term
                IF(calcDist2DP(r(i),r(j)) < rcontacts ) THEN
                    clist(i) % nats = clist(i) % nats + 1
                    clist(j) % nats = clist(j) % nats + 1
                    !
                    clist(i) % vecinos(clist(i) % nats ) = closeOnes(j,clist(i) % nats,distHC)
                    clist(j) % vecinos(clist(j) % nats ) = closeOnes(i,clist(j) % nats,distHC)
                end if
            end do
        end do
    END SUBROUTINE activateClist

!===============================================================================
   subroutine actualizeStepPot(stepPts,natom, nblist)
   use stepPotentials
   use intList
   integer, intent(IN) :: natom
   type(stepPotInt), intent(IN) :: stepPts(natom,natom)
   type(intpList), intent(INOUT) :: nblist(natom)
   integer i, k,j
  
   do i = 1,natom
   
       do k=1,nblist(i)%nats 
       
           j=nblist(i)%idata(k)%patnum

         if ( stepPts(i,j)%active .and. stepPts(i,j)%nstep ==4 ) THEN

             nblist(i)%iData(k) %stepPt =  stepPts(i,j)
           
       END IF
  
      enddo
  enddo
  
  
   end subroutine actualizeStepPot
!===============================================================================
!===============================================================================
   subroutine activateStepPot(stepPts,r,rtarg,oid,rcut2GO2,rcut2GO4,natom,natomCommon,nblist, xsum,radioInteraccion)
   use stepPotentials
   use geometryDP
   use intList
   integer, intent(IN) :: natom,natomCommon
   real*8, intent(IN) :: rcut2GO2,rcut2GO4
   type(pointDP), intent(IN) :: r(natom),rtarg(natomCommon)
   type(stepPotInt), intent(INOUT) :: stepPts(natom,natom)
   INTEGER, DIMENSION(natomCommon),INTENT(IN)::oid
   type(intpList), intent(INOUT) :: nblist(natom)
   real, intent(IN) :: xsum
   real*8, intent(IN) :: radioInteraccion
   integer i, j,ii,jj
   real*8 :: rij_INI,rij_TARG
   logical :: IsActive
!
   where (stepPts%active)
      stepPts%active = .false.
   end where
      
   nblist%nats = 0
   IsActive=.FALSE.
   
 
   do jj = 2,natomCommon
    !  rj = r(j) ! intentem millorar cache hits a r(i)
      do ii = 1,jj-1
          i=oid(ii)
          j=oid(jj)
      secstruct:  if (stepPts(i,j)%tipInt == SSEC ) THEN
! inline 
         rij_INI = calcDist2DP(r(i),r(j))
         rij_TARG =calcDist2DP(rtarg(ii),rtarg(jj))
       nsteps:  IF( stepPts(i,j)%nstep == 2) THEN
             
           isactive1:  IF(rij_INI< rcut2GO2) THEN
                stepPts(i,j)%active = .true.
            !    write(*,*) "ACTIVE 1", i,j
                nblist(i)%nats = nblist(i)%nats + 1
                nblist(j)%nats = nblist(j)%nats + 1
                nblist(i)%iData(nblist(i)%nats) = intData(j, nblist(j)%nats, stepPts(i,j), xsum, 1.e15, 0.)
                nblist(j)%iData(nblist(j)%nats) = intData(i, nblist(i)%nats, stepPts(i,j), xsum, 1.e15, 0.)
             END IF isactive1
   
           ELSEIF  ( stepPts(i,j)%nstep == 4) THEN
     
    isactive21:   IF (MOD(real(i+j),2.0)==0 ) THEN
                 
    isactive22:     IF ((rij_INI<rcut2GO4 ).OR.(rij_TARG<rcut2GO4)) THEN
                                
                     stepPts(i,j)%active = .true.
               !      write(*,*) "ACTIVE 2", i,j
                     nblist(i)%nats = nblist(i)%nats + 1
                     nblist(j)%nats = nblist(j)%nats + 1
                     nblist(i)%iData(nblist(i)%nats) = intData(j, nblist(j)%nats, stepPts(i,j), xsum, 1.e15, 0.)
                     nblist(j)%iData(nblist(j)%nats) = intData(i, nblist(i)%nats, stepPts(i,j), xsum, 1.e15, 0.)
                END IF isactive22
        
                 END IF isactive21
                 
                
         ENDIF nsteps

       ENDIF secstruct
     enddo
   enddo
             
  

  
   do j = 2,natom
    !  rj = r(j) ! intentem millorar cache hits a r(i)
      do i = 1,j-1
  
          if (stepPts(i,j)%tipInt == FAKE ) THEN
              if(calcDistDP(r(i),r(j))<radioInteraccion) THEN
               stepPts(i,j)%active = .true.
        !       write(*,*) "ACTIVE ", i,j
               nblist(i)%nats = nblist(i)%nats + 1
               nblist(j)%nats = nblist(j)%nats + 1
               nblist(i)%iData(nblist(i)%nats) = intData(j, nblist(j)%nats, stepPts(i,j), xsum, 1.e15, 0.)
               nblist(j)%iData(nblist(j)%nats) = intData(i, nblist(i)%nats, stepPts(i,j), xsum, 1.e15, 0.)
              end if
          ENDIF
      enddo
   enddo
  
 end subroutine activateStepPot
!===============================================================================
 
 subroutine thermalize(seed, iev, natom, TEMP, xmassa, v, xm, ekin)
 use geometryDP
 use ran_mod
   integer, intent(IN) :: iev, seed, natom
   real, intent(IN) ::  xmassa, temp, xm
   type(pointDP), intent(OUT) :: v(natom)
   real, intent(OUT) :: ekin
   integer i, kk,NegSeed
   real calcEkin
   type(pointDP) vcm
   real*8 fi, sto
   real :: Kb= 8.314 ! J/(mol K)

!
   negSeed=seed
   fi=ran1(negSeed)
  ! fi=ran1(seed)
   !   write(*,*)"iev desde thermalize", iev
   vcm = pointDP(0., 0., 0.)
   do i = 1,natom
      kk = (seed + 3 * i + 1 + iev)
!   write(*,*)"kk desde thermalize", kk
      fi=ran1(kk)
      v(i)%x = fi - 0.5
!   write(*,*)"kk desde thermalize", kk
       kk = (seed + 3 * i + 1 + iev)
      fi=ran1(kk)
      v(i)%y = fi - 0.5
!   write(*,*)"kk desde thermalize", kk
      kk = (seed + 3 * i + 1 + iev)
      fi=ran1(kk)
      v(i)%z = fi - 0.5
!
      vcm = vcm + xm*1.d0 * v(i)
   enddo
   vcm = (1.d0/xmassa) * vcm
 !  write(*,*)" vel CM", vcm
   do i = 1,natom
      v(i) = v(i) - vcm
   enddo
   sto = sqrt(1.5 * natom * TEMP * kB / calcEkin(v, xm, natom))
 !  write(*,*) sto
   do i = 1,natom
      v(i) =  sto * v(i)
 !     write(*,*) v(i)

   enddo

   ekin = 1.5 * natom * TEMP *kb! No cal recalcularla 
   
 end subroutine thermalize
!========================================================================
!========================================================================
   subroutine calcEpot(natom,natomBare, r, stepPts, epotgo, epotfis, epot)
        use geometryDP
        use stepPotentials
        use constants

        integer, intent(IN) :: natom,natomBare
        type(pointDP), intent(IN) :: r(natom)
        type(stepPotInt), intent(IN) :: stepPts(natom, natom)
        real, intent(OUT) :: epotgo, epotfis
        real epot(MAXTIPINT)
        real dist
        integer i, j, k
        !PENDENT TREBALLAR AMB NBLIST
        epotgo = 0.
        epotfis = 0.
        epot = 0.
        do j = 2, natomBare
            do i = 1, j - 1
                if (stepPts(i, j) % active) then
                    dist = sqrt(real(calcDist2DP(r(i), r(j))))
                    k = stepPts(i, j) % nstep
                    do while ((k .gt. 1).and.dist .lt. stepPts(i, j) % step(k) % r)
                        epot(stepPts(i, j) % tipInt) = epot(stepPts(i, j) % tipInt) - stepPts(i, j) % step(k) % e / FACTE
                        k = k - 1
                    enddo
                    if (dist .lt. stepPts(i, j) % step(k) % r) &
                    epot(stepPts(i, j) % tipInt) = epot(stepPts(i, j) % tipInt) - stepPts(i, j) % step(k) % e / FACTE
                endif
            enddo
       enddo
       epotgo = epot(SSEC)
       epotfis = sum(epot) - epotgo
   end subroutine calcEpot
    !========================================================================
 pure function  calcEkin (v, xm, natom) result (ekin)
 use geometryDP
   integer, intent(IN) :: natom
   real ekin
   type(pointDP), intent(IN) :: v(natom)
   real, intent(IN) :: xm
   real*8, parameter :: a2 = 1.d-20*1.d30
   integer i
   ekin = 0.
   do i = 1,natom
      ekin = ekin + 0.5 * xm* real(A2 * dotDP(v(i), v(i)))
   enddo
 end function calcEkin
!===============================================================================

