program calcintsDUO
    use intCalcs
    implicit none

    integer :: i,j,k,l,m,n
    integer :: iuninp,iunout,iunit,ios
    integer :: npoints, nvib, nj, nvals
    integer :: numvibN, numvibC, numJN, numJC, nOmN, nOmC, nLN, nLC
    integer :: ndictN, ndictC, ndys, ns1, ns2, nstates, cntmask

    real(kind=8) :: Tvib, Trot, ZPEN, ZPEC, EN, EC, rmin, rmax
    real(kind=8) :: nvib_r, nj_r, multN, multC, DJ, aux

    real(kind=8), allocatable :: vibmat(:,:), vibmatN(:,:), vibmatC(:,:)
    real(kind=8), allocatable :: JvalsN(:,:,:), JvalsC(:,:,:)
    real(kind=8), allocatable :: indexN(:,:,:,:), indexC(:,:,:,:)
    real(kind=8), allocatable :: keyN(:,:,:,:), keyC(:,:,:,:)
    real(kind=8), allocatable :: rvals(:), DJvals(:), rcomp(:)
    real(kind=8), allocatable :: evals(:,:,:,:,:,:), relInt(:,:,:,:,:,:)
    
    type(coef_dictt), allocatable :: coef_dictN(:), coef_dictC(:)
    type(dys_dictt), allocatable :: dysdict(:)
    type(state_map), allocatable :: s1(:), s2(:)
    
    character(len=256) :: fname, fnameN, fnameC, fileout, path, line
    character(len=256) :: fNvib,fCvib,fNJval,fCJval
    character(len=32) :: key, stt1, stt2
    character(len=2) :: energy_unit
    
    logical, allocatable :: mask(:)

    open(newunit=iuninp,file='fort.inp',action='read')

    allocate(s1(1), s2(1))

    do
        read(iuninp,'(A)',iostat=ios) line
        if (ios /= 0) exit

        read(line,*,iostat=ios) key
        if (ios /= 0) cycle

        if (len_trim(line) == 0) cycle
        if (line(1:1) == "!") cycle

        select case(trim(key))
            case("TempRV")
                read(line,*,iostat=ios) key, Trot, Tvib
            case("LimR")
                read(line,*,iostat=ios) key, rmin, rmax
            case("NVJ")
                read(line,*,iostat=ios) key, nvib_r, nj_r
                nvib = int(nvib_r) ; nj = int(nj_r)
            case("DJ")
                nvals = count_values(line) - 1
                allocate(DJvals(nvals))
                read(line,*,iostat=ios) key, DJvals
            case("npts")
                read(line,*,iostat=ios) key, npoints
            case("Name")
                read(line,*,iostat=ios) key, stt1, stt2
                stt1 = trim(stt1)
                stt2 = trim(stt2)
            case("pSys")
                line = trim(line)
                read(line,*,iostat=ios) key, fnameN, fnameC
                print*, fnameN, fnameC
            case("ZPE")
                read(line,*,iostat=ios) key, ZPEN, ZPEC
            case("Etot")
                read(line,*,iostat=ios) key, EN, EC
            case("Numv")
                read(line,*,iostat=ios) key, numvibN, numvibC
            case("NumJ")
                read(line,*,iostat=ios) key, numJN, numJC
            case("Mult")
                read(line,*,iostat=ios) key, multN, multC
            case("Lval")
                read(line,*,iostat=ios) key, nLN, nLC
            case("Dyson")
                read(line,*,iostat=ios) key, path
            case("IndSts")
                nOmN = int(multN)*(nLN+1) ; nOmC = int(multC)*(nLC+1)
                allocate(s1(1)%idx(nOmN),s2(1)%idx(nOmC))
                read(line,*,iostat=ios) key, s1(1)%idx, s2(1)%idx
            case("Sigmas")
                allocate(s1(1)%val(nOmN),s2(1)%val(nOmC))
                read(line,*,iostat=ios) key, s1(1)%val, s2(1)%val
        end select
    end do
    close(iuninp)

    write(6,*) "Trot, Tvib =", Trot, Tvib
    write(6,*) "rmin, rmax =", rmin, rmax
    write(6,*) "npoints, nvib, nj =", npoints, nvib, nj
    write(6,*) "fnameN, fnameC =", fnameN, fnameC
    write(6,*) "ZPEN, ZPEC =", ZPEN, ZPEC
    write(6,*) "EN, EC =", EN, EC
    write(6,*) "nLN, nLC =", nLN, nLC
    write(6,*) "nOmN, nOmC =", nOmN, nOmC
    write(6,*) "numvibN, numvibC =", numvibN, numvibC
    write(6,*) "numJN, numJC =", numJN, numJC
    write(6,*) "multN, multC =", multN, multC
    write(6,*) "s1, s2 =", s1(1)%idx, s2(1)%idx
    write(6,*) "s1, s2 =", s1(1)%val, s2(1)%val

    allocate(vibmatN(npoints,nvib),vibmatC(npoints,nvib))
    allocate(rcomp(npoints))

    fNvib = trim(fnameN) // "vibeigenvect_vib.chk"
    fCvib = trim(fnameC) // "vibeigenvect_vib.chk"
    call parse_duo_vib_einfun(fNvib,npoints,nvib,vibmatN)
    call parse_duo_vib_einfun(fCvib,npoints,nvib,vibmatC)

    allocate(JvalsN(nvib,nj,nOmN),JvalsC(nvib,nj,nOmC))
    allocate(indexN(nvib,nj,nOmN,7),indexC(nvib,nj,nOmC,7))
    allocate(keyN(nvib,nj,nOmN,6),keyC(nvib,nj,nOmC,6))
    fNJval = trim(fnameN) // "rovibronic_energies.dat"
    fCJval = trim(fnameC) // "rovibronic_energies.dat"
    call read_J_values(fNJval,nvib,nj,multN,nOmN,JvalsN,indexN,keyN)
    call read_J_values(fCJval,nvib,nj,multC,nOmC,JvalsC,indexC,keyC)
    
    print*, "nvib, nj =", nvib, nj
    

    call read_coefficients(fnameN,nvib,nj,multN,nOmN,coef_dictN,ndictN)
    call read_coefficients(fnameC,nvib,nj,multC,nOmC,coef_dictC,ndictC)
    
    ns1 = size(s1(1)%idx) ; ns2 = size(s2(1)%idx)
    ! s1 = s1(1)%idx ; s2 = s2(1)%idx
    print*, "ns1, ns2 =", ns1, ns2
    print*, "s1, s2 =", s1(1)%idx, s2(1)%idx

    allocate(dysdict(ns1*ns2))
    call read_dys_bySigma(s1(1)%idx,s2(1)%idx,s1(1)%val,s2(1)%val,ns1,ns2,path,dysdict,ndys,rvals)
    print*,ndys

    open(unit=1, file="Dipole_moment_functions.dat", status="old",action="read")
    do i=1,npoints
        read(1,*) rcomp(i), aux, aux
    enddo
    close(1)

    cntmask = 1
    do i=1,npoints
        if (rcomp(i) > rmin .and. rcomp(i) < rmax) then
            cntmask = cntmask + 1
        endif
    enddo
    allocate(mask(cntmask))
    mask = .false.
    mask = (rcomp > rmin).and.(rcomp < rmax)

    allocate(evals(numvibN,numvibC,numJN,numJC,nOmN,nOmC),&
            relInt(numvibN,numvibC,numJN,numJC,nOmN,nOmC))
    ! if (size(DJvals) > 1) then

    print*, indexN(33,81,1,:), indexC(33,81,1,:),keyN(33,81,1,:), keyC(33,81,1,:)
    !     do i=1,size(DJvals)
    !         DJ = DJvals(i)
            call intT_sigma(JvalsN,JvalsC,nvib,rcomp,rvals,npoints,size(rvals),ndys,dysdict, &
                vibmatN,vibmatC,mask,Tvib,Trot,ZPEN,ZPEC,EN,EC,DJ, &
                numvibN,numvibC,numJN,numJC,nOmN,nOmC,indexN,indexC, &
                keyN,keyC,coef_dictN,coef_dictC,ndictN,ndictC,evals,relInt)
            
    !         fileout = trim(fnameN) // trim(fnameC) // trim(character(Trot)) // &
    !                 trim(character(Tvib)) // trim(character(DJvals)) // ".dat"
    !         call dump_spectrum_long(fileout,numvibN,numvibC,numJN,numJC,nOmN,nOmC, &
    !                                 evals,relInt,indexN,indexC,energy_unit="eV",tol_I=0.0d0)
    !     enddo
    ! else
    !     call intT_sigma(JvalsN,JvalsC,nvib,rvals,ndys,npoints,dysdict, &
    !         vibmatN,vibmatC,mask,Tvib,Trot,ZPEN,ZPEC,EN,EC,DJ, &
    !         numvibN,numvibC,numJN,numJC,nOmN,nOmC,indexN,indexC, &
    !         keyN,keyC,coef_dictN,coef_dictC,ndictN,ndictC,evals,relInt)
        
    !     fileout = trim(fnameN) // trim(fnameC) // trim(character(Trot)) // &
    !             trim(character(Tvib)) // trim(character(DJvals)) // ".dat"
    !     call dump_spectrum_long(fileout,numvibN,numvibC,numJN,numJC,nOmN,nOmC, &
    !                             evals,relInt,indexN,indexC,energy_unit="eV",tol_I=0.0d0)
    ! endif

    deallocate(vibmatN,vibmatC)
    deallocate(JvalsN,JvalsC,indexN,indexC,keyN,keyC)

end program calcintsDUO