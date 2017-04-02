!
! File:   DMD.f
! Author: gelpi
!
! Created on 14 de marzo de 2012, 16:54
!
subroutine dmdIntloop(nblist, r, v, nxm, tacact, tpart, ipart, TACT, TMIN, temps, iev, natom, toUpdate)
use intList
use geometryDP
   integer, intent(IN) :: natom
   real, intent(IN) :: nxm
   real*8, intent(IN) :: TACT, TMIN
   type(pointDP), intent(INOUT) :: r(natom), v(natom)
   type(intpList), intent(INOUT) :: nblist(natom)
   integer, intent(INOUT) :: iev
   real*8, intent(INOUT) :: tacact, temps
   logical toUpdate(natom)
!
   real*8 tpart(natom)
   real*8 tevent, tevent1
   integer ipart(natom)
   integer mem1, mem2, npair1, npair2, i,j
!
   tacact = 0.
   toUpdate = .true.
! evolucio temporal
!------------------------------------------------------------------------------
   do while (tacact .lt. TACT)
! busca quina es la propera colisio
       do i = 1, natom
            if (toUpdate(i)) then
                tpart(i) = 1.d15               
                do j = 1, nblist(i)%nats
                    if (nblist(i)%iData(j)%timp.lt.tpart(i)) then
                        tpart(i) = nblist(i)%iData(j)%timp
                        ipart(i) = j ! ipart recull index NO num d'atom
                    endif
                    toUpdate(i) = .false.
                enddo  
            endif
        enddo
        tevent = 1.e15
        call nextCol(mem1, mem2, npair1, npair2, tevent, ipart, tpart, natom, nblist)
        !
        tevent1 = tevent - temps
        temps = tevent
        tacact = tacact + tevent1
        iev = iev + 1
        ! translacio i variacio dels temps
        do i = 1, natom ! mantenim la versio inline degut a la barreja de tipus de real !!
            r(i)%x = r(i)%x + v(i)%x * tevent1
            r(i)%y = r(i)%y + v(i)%y * tevent1
            r(i)%z = r(i)%z + v(i)%z * tevent1
        enddo
        call updateV(r, v, nblist(mem1)%iData(npair1)%deltak, nxm, nblist(mem1)%iData(npair1)%xsum, &
        mem1, mem2, natom)
        ! Assegurem que trespassem la barrera en la direcció correcta
        r(mem1)%x = r(mem1)%x + v(mem1)%x * tevent1/DBLE(100000.)
        r(mem1)%y = r(mem1)%y + v(mem1)%y * tevent1/DBLE(100000.)
        r(mem1)%z = r(mem1)%z + v(mem1)%z * tevent1/DBLE(100000.)
        r(mem2)%x = r(mem2)%x + v(mem2)%x * tevent1/DBLE(100000.)
        r(mem2)%y = r(mem2)%y + v(mem2)%y * tevent1/DBLE(100000.)
        r(mem2)%z = r(mem2)%z + v(mem2)%z * tevent1/DBLE(100000.)
        !
        ! ara actualitza els temps de colisio per a les dues particules que han xocat
        ! anulem la colisio per a que no es repeteixi
        nblist(mem1)%iData(npair1)%timp = 1.d15
        nblist(mem2)%iData(npair2)%timp = 1.d15
        call updateXocPart(mem1, nblist, temps, r, v, TMIN, natom, toUpdate)
        call updateXocPart(mem2, nblist, temps, r, v, TMIN, natom, toUpdate)
    enddo
!    deallocate (tpart, ipart)
    ! end do while (tacact.lt.tact)------------------------------------------------
end subroutine dmdIntloop

