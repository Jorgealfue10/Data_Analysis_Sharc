module intCalcs
    implicit none

    type coef_key
        real(kind=8) :: J
        real(kind=8) :: Om
        real(kind=8) :: S
        real(kind=8) :: L
        integer      :: parity
        integer      :: indx
    end type

    type coef_dictt
        type(coef_key) :: key
        complex(kind=8), allocatable :: val(:)
        integer :: nread = 0
        logical :: used = .false.
    end type coef_dictt

    type dys_key
        real(kind=8) :: vals(2)
    end type

    type dys_dictt
        type(dys_key) :: key
        real(kind=8), allocatable :: val(:)
        logical :: used = .false.
    end type

    type state_map
        character(len=20) :: name
        integer, allocatable :: idx(:)
        real(8), allocatable :: val(:)
    end type state_map
    
    contains
    
    !---------------------------------------------------------------------------!
    ! Read head of eigenvib vect
    !---------------------------------------------------------------------------!
    logical function is_head_eigenvib(line)
        implicit none

        character(len=256), intent(in) :: line
        integer :: a,c
        real(kind=8) :: b
        
        read(line,*,err=999) a,b,c

        is_head_eigenvib = .true.
        return
        
        999 is_head_eigenvib = .false.
    end function

    !---------------------------------------------------------------------------!
    ! Parsea el archivo de autofunciones de vibración de DUO
    !---------------------------------------------------------------------------!

    subroutine parse_duo_vib_einfun(fname,npoints,nvib,vibmat)
        implicit none

        character(len=256), intent(inout) :: fname
        integer, intent(in) :: npoints,nvib
        real(kind=8), intent(out) :: vibmat(npoints,nvib)

        character(len=256) :: line
        integer :: iunit, ios, ncurrvib, ncurrpoint, auxint
        real(kind=8) :: val, eval

        fname = trim(fname)

        open(newunit=iunit,file=fname,status='old',action='read')

        do
            read(iunit,'(A)', iostat=ios) line
            if (ios /= 0 ) exit

            if (is_head_eigenvib(line)) then
                ncurrpoint = 1
                read(line,*) auxint, eval, auxint, ncurrvib
                ncurrvib = ncurrvib + 1
                if (ncurrvib > nvib) then
                    exit
                endif
            else if (index(line,'End of contracted basis') /= 0) then
                exit
            else
                read(line,*) val
                vibmat(ncurrpoint,ncurrvib) = val
                ncurrpoint = ncurrpoint + 1
            endif
        enddo

        close(iunit)
    end subroutine

    !---------------------------------------------------------------------------!
    ! Bra ket overlap
    !---------------------------------------------------------------------------!

    subroutine bra_ket(npoints, vib_ini, vib_fin, dymat, braket)
        implicit none

        integer, intent(in) :: npoints
        real(kind=8), intent(in) :: vib_ini(npoints), vib_fin(npoints)
        real(kind=8), intent(in) :: dymat(npoints,npoints)
        real(kind=8), intent(out) :: braket

        braket = dot_product(vib_fin,matmul(dymat,vib_ini))
    end subroutine

    !---------------------------------------------------------------------------!
    ! Bra ket overlap with coefficients
    !---------------------------------------------------------------------------!

    subroutine bra_ket_coef(npoints, nvib, vib_ini, vib_fin, dymat, coef_ini, coef_fin, braket)
        implicit none

        integer, intent(in) :: npoints, nvib
        real(kind=8), intent(in) :: vib_ini(npoints), vib_fin(npoints)
        complex(kind=8), intent(in) :: coef_ini(nvib), coef_fin(nvib)
        real(kind=8), intent(in) :: dymat(npoints,npoints)
        real(kind=8), intent(out) :: braket

        integer :: i
        real(kind=8) :: aux

        if (range(coef_ini) == 0 .or. range(coef_fin) == 0) then
            call bra_ket(npoints, vib_ini, vib_fin, dymat, braket)
        else
            braket = dot_product(vib_fin,matmul(dymat,vib_ini))
            braket = braket * real(dot_product(conjg(coef_fin),coef_ini))
        endif
    end subroutine

    !---------------------------------------------------------------------------!
    ! Reading J values E and quantum numbers
    !---------------------------------------------------------------------------!

    subroutine read_J_values(fname,nvib,nj,mult,nOmega,Jvals,index_list,key_dict)
        implicit none

        character(len=256), intent(in) :: fname
        integer, intent(in) :: nvib, nj, nOmega
        real(kind=8), intent(in) :: mult
        real(kind=8), intent(out) :: Jvals(nvib,nj,nOmega), index_list(nvib,nj,nOmega,7),key_dict(nvib,nj,nOmega,6)

        integer :: iunit, ios, indOm
        character(len=256) :: line, aux
        character(len=1) :: parity
        real(kind=8) :: J, eval, stt, vibn, lambda, sval, sigma, omega, pm, ind, Smult


        Smult = (mult-1.d0)/2.d0

        open(newunit=iunit,file=fname,status='old',action='read')        

        Jvals = 0.d0 ; index_list = 0.d0 ; key_dict = 0.d0

        do 
            read(iunit,'(A)', iostat=ios) line
            if (ios /= 0 ) exit

            read(line,*) J, ind, eval, stt, vibn, lambda, sval, sigma, omega, parity, aux, aux
            if (vibn+1 < nvib) then
                if (abs(Smult - nint(Smult)) < 1d-8) then
                    indOm = int(omega+Smult)+1
                else 
                    ! print*,int(vibn)+1, J, int(J)+1, int(omega+Smult)+2
                    indOm = int(omega+Smult)+2
                endif

                ! print*,"joder",int(vibn)+1, int(J)+1, indOm, eval,index
                Jvals(int(vibn)+1, int(J)+1, indOm) = eval

                ! print*,'A'
                if (parity.eq."-") then
                    pm = -1.d0
                else
                    pm = 1.d0
                endif

                if (int(vibn)+1 > nvib .or. int(J)+1 > nj .or. indOm > nOmega) then
                    print*, "OUT OF BOUNDS:"
                    print*, "vibn=", vibn, "J=", J, "omega=", omega
                    print*, "indices=", int(vibn)+1, int(J)+1, indOm
                    print*, "limits=", nvib, nj, nOmega
                    stop
                endif

                index_list(int(vibn)+1, int(J)+1, indOm,1) = vibn
                index_list(int(vibn)+1, int(J)+1, indOm,2) = J
                index_list(int(vibn)+1, int(J)+1, indOm,3) = omega
                index_list(int(vibn)+1, int(J)+1, indOm,4) = sigma
                index_list(int(vibn)+1, int(J)+1, indOm,5) = lambda
                index_list(int(vibn)+1, int(J)+1, indOm,6) = pm
                index_list(int(vibn)+1, int(J)+1, indOm,7) = ind

                key_dict(int(vibn)+1, int(J)+1, indOm,1) = J
                key_dict(int(vibn)+1, int(J)+1, indOm,2) = omega
                key_dict(int(vibn)+1, int(J)+1, indOm,3) = sigma
                key_dict(int(vibn)+1, int(J)+1, indOm,4) = lambda
                key_dict(int(vibn)+1, int(J)+1, indOm,5) = pm
                key_dict(int(vibn)+1, int(J)+1, indOm,6) = ind
            endif
        enddo

        ! close(iunit)
    end subroutine

    !---------------------------------------------------------------------------!
    ! Reading coefficients
    !---------------------------------------------------------------------------!

    logical function same_key(a,b)
        implicit none
        type(coef_key), intent(in) :: a,b

        same_key = .true.

        if (abs(a%J  - b%J ) > 1d-8) same_key = .false.
        if (abs(a%Om - b%Om) > 1d-8) same_key = .false.
        if (abs(a%S  - b%S ) > 1d-8) same_key = .false.
        if (abs(a%L  - b%L ) > 1d-8) same_key = .false.
        if (a%parity /= b%parity)    same_key = .false.
        if (a%indx   /= b%indx)      same_key = .false.
    end function

    ! function find_coef_index(dict, ndict, key) result(idx)
    !     implicit none
    !     integer, intent(in) :: ndict
    !     type(coef_dictt), intent(in) :: dict(ndict)
    !     type(coef_key), intent(in) :: key

    !     integer :: idx, i

    !     idx = -1

    !     do i = 1, ndict
    !         if (dict(i)%used) then
    !             if (same_key(dict(i)%key, key)) then
    !                 idx = i
    !                 return
    !             endif
    !         endif
    !     enddo
    ! end function

    subroutine append_coef(coef_dict,ndict, nvib, key, coeff)
        implicit none
        integer, intent(inout) :: ndict
        integer, intent(in) :: nvib
        type(coef_dictt), allocatable, intent(inout) :: coef_dict(:)
        type(coef_key), intent(in) :: key
        real(kind=8), intent(in) :: coeff

        integer :: i
        logical :: found

        found = .false.

        do i = 1, ndict
            if (coef_dict(i)%used) then
                if (same_key(coef_dict(i)%key, key)) then
                    found = .true.
                    coef_dict(i)%nread = coef_dict(i)%nread + 1
                    coef_dict(i)%val(coef_dict(i)%nread) = cmplx(coeff,0.d0,kind=8)
                    exit
                endif
            endif 
        enddo

        if (.not.found) then
            ndict = ndict + 1
            allocate(coef_dict(ndict)%val(nvib))
            coef_dict(ndict)%used = .true.
            coef_dict(ndict)%key = key
            coef_dict(ndict)%nread = 1
            coef_dict(ndict)%val = cmplx(0.d0,0.d0,kind=8)
            coef_dict(ndict)%val(1) = cmplx(coeff,0.d0,kind=8)
        endif
    end subroutine

    subroutine read_coefficients(fname,nvib,nj,mult,nOmega,coef_dict,ndict)
        implicit none

        character(len=*), intent(in) :: fname
        integer, intent(in) :: nvib, nj, nOmega
        real(kind=8), intent(in) :: mult
        type(coef_dictt), allocatable, intent(out) :: coef_dict(:)
        integer, intent(out) :: ndict

        integer :: iunit, ios, maxdict
        character(len=256) :: line, aux
        real(kind=8) :: J, coeff, stt, vibn, lambda, sval, sigma, omega, pm, index, Smult
        type(coef_key) :: key

        ndict = 0
        Smult = (mult-1.d0)/2.d0

        maxdict = nvib*nj*nOmega
        allocate(coef_dict(maxdict))
        coef_dict(:)%used = .false.
        print*, "maxdict=", maxdict, size(coef_dict)

        open(newunit=iunit,file=fname,status='old',action='read')

        do 
            read(iunit,'(A)', iostat=ios) line
            if (ios /= 0 ) exit
            
            read(line,*,iostat=ios) index, J, pm, coeff, aux, vibn, lambda, sval, sigma, omega, aux
            if (ios /= 0 ) cycle

            if (abs(Smult - nint(Smult)) < 1d-8) then
                    indOm = int(omega+Smult)+1
            else 
                indOm = int(omega+Smult)+2
            endif

            if (lambda < 0.d0) then 
                nl = 1
            coef_mat(int(vibn)+1,int(index),int(J)+1,indOm,int(lambda),int(sigma+Smult)+1,int(pm)) = coeff

            ! key%J = J
            ! key%Om = omega
            ! key%S = sigma
            ! key%L = lambda
            ! key%parity = int(pm)
            ! key%indx = int(index)
                
            ! coef_dict(ndict)%key = key
            call append_coef(coef_dict, ndict, nvib, key, coeff)
            ! print*, vibn, line
        enddo

        print*, "ndict=", ndict
        print*, "nread=", coef_dict(1)%nread

        close(iunit)
    end subroutine

    !---------------------------------------------------------------------------!
    ! Reading Dyson Norms by sigma
    !---------------------------------------------------------------------------!

    subroutine read_dyson_file(fname,r,dyson,n)
        implicit none

        character(len=*), intent(in) :: fname
        integer, intent(inout) :: n
        real(kind=8), allocatable, intent(out) :: r(:), dyson(:)

        integer :: iunit, ios, i
        real(kind=8) :: rtmp, dtmp
        character(len=256) :: line

        n=0

        open(newunit=iunit,file=fname,status='old',action='read')
        do
            read(iunit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            n = n + 1
        enddo

        allocate(r(n),dyson(n))

        rewind(iunit)

        do i = 1, n
            read(iunit,*) rtmp, dtmp
            r(i) = rtmp
            dyson(i) = dtmp
        enddo

        close(iunit)
    end subroutine

    logical function same_dys_key(a,b)
        implicit none
        type(dys_key), intent(in) :: a,b
        same_dys_key = all((a%vals - b%vals) < 1d-8)
    end function

    function find_dys_index(dict, ndict, key) result(idx)
        integer :: idx, i, ndict
        type(dys_dictt), intent(in) :: dict(ndict)
        type(dys_key), intent(in) :: key

        idx = -1
        do i = 1, ndict
            if (dict(i)%used) then
                if (same_dys_key(dict(i)%key, key)) then
                    idx = i
                    return
                endif
            endif
        enddo
    end function

    subroutine append_sigma(sigma_dict, ndict, key, nr, dyson)
        implicit none
        integer, intent(inout) :: ndict, nr
        type(dys_dictt), intent(inout) :: sigma_dict(ndict)
        type(dys_key), intent(in) :: key
        real(kind=8), intent(in) :: dyson(nr)

        integer :: i
        logical :: found

        found = .false.
        do i = 1, ndict
            if (sigma_dict(i)%used) then
                if (same_dys_key(sigma_dict(i)%key, key)) then
                    sigma_dict(i)%val = sigma_dict(i)%val + dyson**2
                    found = .true.
                    exit
                endif 
            endif
        enddo

        if (.not.found) then
            ndict = ndict + 1
            sigma_dict(ndict)%used = .true.
            sigma_dict(ndict)%key = key
            allocate(sigma_dict(ndict)%val(size(dyson)))
            sigma_dict(ndict)%val = dyson**2
        endif
    end subroutine

    subroutine read_dys_bySigma(d1,d2,sd1,sd2,nd1,nd2,path,dict,ndict,r)
        implicit none

        integer, intent(in) :: nd1,nd2
        integer, intent(in) :: d1(nd1),d2(nd2)
        real(kind=8), intent(in) :: sd1(nd1),sd2(nd2)
        character(len=*), intent(in) :: path
        type(dys_dictt), intent(out) :: dict(nd1*nd2)
        integer, intent(out) :: ndict
        real(kind=8), allocatable, intent(out) :: r(:)

        integer :: i,j,nr
        real(kind=8) :: sigma_i, sigma_j
        real(kind=8), allocatable :: dyson(:)
        type(dys_key) :: key

        character(len=256) :: fname

        ndict = 0
        do i = 1, nd1
            do j = 1, nd2
                sigma_i = sd1(i) 
                sigma_j = sd2(j) 
                key%vals = [sigma_i,sigma_j]
                
                write(fname,'(A,"dyson_",I2.2,"_",I2.2,".dat")') trim(path),d1(i),d2(j)
                
                call read_dyson_file(fname,r,dyson,nr)
                ndict = ndict + 1
                dict(ndict)%val = dyson

                ! call append_sigma(dict, ndict, key, nr, dyson)

            enddo
        enddo

    end subroutine

    !---------------------------------------------------------------------------!
    ! Wigner 3j explicit expression
    !---------------------------------------------------------------------------!

    ! Kronecker delta
    integer function kronecker(i,j)
        implicit none
        integer, intent(in) :: i,j
        kronecker = 1
        if (i /= j) kronecker = 0
    end function

    ! Factorial
    integer function fact(i)
        implicit none
        integer, intent(inout) :: i
        fact = 1
        do while (i > 1)
            fact = fact * i
            i = i - 1
        enddo
    end function

    subroutine W3_exp(J1,J2,J3,m1,m2,m3,w3)
        implicit none
        real(kind=8), intent(in) :: J1,J2,J3,m1,m2,m3
        real(kind=8), intent(out) :: w3
        integer :: K,N,i
        integer :: deltaK, m1pf
        integer :: r11, r12, r13, r14
        real(kind=8) :: first_root
        integer :: r21, r22, r23, r24, r25, r26
        real(kind=8) :: second_root
        integer :: r31, r32, r33, r34, r35, m1pk, kf
        real(kind=8) :: third_root, coeff

        K=int(max(0,int(J2-J3-m1),int(J1-J3+m2)))
        N=int(min(int(J1+J2-J3),int(J1-m1),int(J2+m2)))

        deltaK = kronecker(int(m1+m2+m3),0)
        m1pf = int((-1)**(J1-J2-m3))

        r11 = int(J1+J2-J3) ; r11 = fact(r11)
        r12 = int(J1-J2+J3) ; r12 = fact(r12)
        r13 = int(-J1+J2+J3) ; r13 = fact(r13)
        r14 = int(J1+J2+J3+1) ; r14 = fact(r14)
        first_root = sqrt(float((r11*r12*r13)/r14))

        r21 = int(J1-m1) ; r21 = fact(r21)
        r22 = int(J1+m1) ; r22 = fact(r22)
        r23 = int(J2-m2) ; r23 = fact(r23)
        r24 = int(J2+m2) ; r24 = fact(r24)
        r25 = int(J3-m3) ; r25 = fact(r25)
        r26 = int(J3+m3) ; r26 = fact(r26)
        second_root = sqrt(float(r21*r22*r23*r24*r25*r26)) 

        third_root = 0.d0
        do i=K,N+1
            m1pk = (-1)**k
            kf = int(k) ; kf = fact(k)
            r31 = int(J1+J2-J3-k) ; r31 = fact(r31)
            r32 = int(J1-m1-k) ; r32 = fact(r32)
            r33 = int(J2+m2-k) ; r33 = fact(r33)
            r34 = int(J3-J1+m1+k) ; r34 = fact(r34)
            r35 = int(J3-J1-m2+k) ; r35 = fact(r35)
            coeff = m1pk/(kf*r31*r32*r33*r34*r35)
            third_root = third_root + coeff
        enddo

        w3 = deltaK*m1pf*first_root*second_root*third_root
    end subroutine

    !---------------------------------------------------------------------------!
    ! Normalization factors for rotational levels including Projections
    !---------------------------------------------------------------------------!

    subroutine rot_normf(T,Jener,ZPE,Estt,numvib,numJ,Projs,index,normf)
        implicit none
        integer, intent(in)  :: numvib,numJ,Projs
        real(kind=8), intent(in)  :: T,ZPE,Estt
        real(kind=8), intent(in)  :: Jener(numvib,numJ,Projs),index(numvib,numJ,Projs,7)
        real(kind=8), intent(out) :: normf(numvib)

        integer :: i,j,k
        real(kind=8) :: Eh_to_cm, kb, Emin, jval
        Eh_to_cm = 219474.6 ; kb = 3.166811563e-6  ! Eh/K

        Emin = 10000000000000000d0
        do i=1, numvib
            do j=1, numJ
                do k=1, Projs
                    Emin = min(Emin, (Jener(i,j,k) + ZPE)/Eh_to_cm + Estt)
                enddo
            enddo
        enddo

        normf = 0.d0
        do i=1, numvib
            do j=1, numJ
                do k=1, Projs
                    jval = index(i,j,k,2)
                    normf(i) = normf(i) + (2*jval+1)*exp(-kb*(Jener(i,j,k) + ZPE + Emin)/T)
                enddo
            enddo
        enddo
    end subroutine

    !---------------------------------------------------------------------------!
    ! Normalization factors for vibrational J=0 levels including Projections
    !---------------------------------------------------------------------------!

    subroutine vib_normf(T,Jener,ZPE,Estt,numvib,numJ,Projs,normf)
        implicit none
        integer, intent(in)  :: numvib,numJ,Projs
        real(kind=8), intent(in)  :: T,ZPE,Estt
        real(kind=8), intent(in)  :: Jener(numvib,numJ,Projs)
        real(kind=8), intent(out) :: normf

        integer :: i,j,k
        real(kind=8) :: Eh_to_cm, kb, Emin, jvals

        Eh_to_cm = 219474.6 ; kb = 3.166811563e-6  ! Eh/K

        Emin = 10000000000000000d0
        do i=1, numvib
            Emin = min(Emin, (Jener(i,1,2) + ZPE)/Eh_to_cm + Estt)
        enddo

        normf = 0.d0
        do i=1, numvib
            normf = normf + exp(-kb*(Jener(i,1,2) + ZPE + Emin)/T)
        enddo
    end subroutine

    !---------------------------------------------------------------------------!
    !Make intensity matrices including Dyson splines as a function of sigma combs
    !---------------------------------------------------------------------------!

    subroutine intT_sigma(PHener,PHMener,nvib,nj,rcomp,rvals,npoints,nrdys,ndys,dysondict, &
        PHvibs,PHMvibs,mask,Tvib,Trot,ZPEPH,ZPEPHM,PHsttsE,PHMsttsE,DJ, &
        numvibPH,numvibPHM,numJPH,numJPHM,numOmPH,numOmPHM,indexPH,indexPHM, &
        keylistPH,keylistPHM,coefPH,coefPHM,ndictN,ndictC,evals,relInt)

        use spline_module

        implicit none
        integer, intent(in)  :: numvibPH,numvibPHM,numJPH,numJPHM,numOmPH,numOmPHM
        integer, intent(in)  :: ndys,npoints,nvib,nj,nrdys,ndictN,ndictC      
        real(kind=8), intent(in)  :: PHener(nvib,nj,numOmPH),PHMener(nvib,nj,numOmPHM)
        real(kind=8), intent(in)  :: Tvib,Trot,ZPEPH,ZPEPHM,PHsttsE,PHMsttsE,DJ
        real(kind=8), intent(in)  :: PHvibs(npoints,nvib),PHMvibs(npoints,nvib),rcomp(npoints),rvals(nrdys)
        real(kind=8), intent(in)  :: indexPH(nvib,nj,numOmPH,7),indexPHM(nvib,nj,numOmPHM,7)
        real(kind=8), intent(in)  :: keylistPH(nvib,nj,numOmPH,6),keylistPHM(nvib,nj,numOmPHM,6)

        real(kind=8), allocatable :: djlist(:)
        logical, intent(in)  :: mask(npoints)

        complex(kind=8), pointer :: coef_ini(:),coef_fin(:)
        real(kind=8), pointer :: sigma_r(:)

        type(dys_dictt), intent(in), target  :: dysondict(ndys)
        type(dys_key) :: dkey
        type(coef_dictt), intent(in), target  :: coefPH(nvib),coefPHM(nvib)
        type(coef_key)  :: ckeyPH,ckeyPHM
        
        real(kind=8), intent(out)  :: &
            evals(numvibPH,numvibPHM,numJPH,numJPHM,numOmPH,numOmPHM), &
            relInt(numvibPH,numvibPHM,numJPH,numJPHM,numOmPH,numOmPHM)

        integer :: i,j,k,l,m,n,pl,nvalid,djO,djJ
        integer :: idxd,idxPH,idxPHM
        real(kind=8) :: normPHvib,Ehtocm
        real(kind=8) :: bk_val,degJi,degvi,diff
        real(kind=8) :: E_PH, E_PHM, erot, evib
        real(kind=8) :: Jival, Jfval, OmN, OmC, rotcoeff
        real(kind=8), allocatable :: spldys(:),dymat(:,:)
        real(kind=8), allocatable :: vib_PH(:),vib_PHM(:),normPHrot(:),rcomp_masked(:)
        type(spline1d), allocatable :: dys_spline(:)

        nvalid = count(mask)
        allocate(vib_PH(nvalid),vib_PHM(nvalid))
        allocate(dys_spline(ndys),normPHrot(nvib))
        allocate(spldys(nvalid),dymat(nvalid,nvalid),rcomp_masked(nvalid))
        Ehtocm = 219474.6

        n = int(2*DJ) + 1
        allocate(djlist(n))

        do i=0,int(2*DJ)
            djlist(i+1) = real(i) - DJ
        enddo

        do i=1,ndys
            call init_spline(dys_spline(i),rvals*0.52917,dysondict(i)%val(:))
        enddo

        call vib_normf(Tvib,PHener,ZPEPH,PHsttsE,nvib,nj,numOmPH,normPHvib)
        call rot_normf(Trot,PHener,ZPEPH,PHsttsE,nvib,nj,numOmPH,indexPH,normPHrot)

        rcomp_masked = pack(rcomp, mask)
        do i=1,numvibPH

            vib_PH = pack(PHvibs(:,i),mask)

            E_PH = (PHener(i,1,2)+ZPEPH)/Ehtocm
            evib = E_PH - ((minval(PHener(:,1,2)) + ZPEPH)/Ehtocm)
            degvi = exp(-evib/(3.166811563*10**(-6)*Tvib))/normPHvib

            do j=1,numvibPHM
                vib_PHM = pack(PHMvibs(:,j),mask)

                do k=1,numJPH
                    do l=1,numJPHM
                        do m=1,numOmPH

                            E_PH = (PHener(i,k,m)+ZPEPH)/Ehtocm
                            erot = E_PH - ((minval(PHener(i,:,:)) + ZPEPH)/Ehtocm)
                            Jival = indexPH(i,k,m,2) ; OmN = indexPH(i,k,m,3)
                            degJi = ((2*Jival+1)*exp(-erot/(3.166811563*10**(-6)*Trot)))/normPHvib
                            E_PH = E_PH + PHsttsE

                            ! ckeyPH%vals = &
                            ! [keylistPH(i,k,m,1),keylistPH(i,k,m,2),keylistPH(i,k,m,3),&
                            ! keylistPH(i,k,m,4),keylistPH(i,k,m,5),keylistPH(i,k,m,6)]

                            do n=1,numOmPHM

                                ! ckeyPHM%vals = &
                                ! [keylistPHM(j,l,n,1),keylistPHM(j,l,n,2),keylistPHM(j,l,n,3), &
                                ! keylistPHM(j,l,n,4),keylistPHM(j,l,n,5),keylistPHM(j,l,n,6)]

                                dkey%vals = [keylistPH(i,k,m,2),keylistPHM(j,l,n,2)]
                                
                                ! idxPH = find_coef_index(coefPH, nvib, ckeyPH)
                                ! idxPHM = find_coef_index(coefPHM, nvib, ckeyPHM)
                                ! print*,idxPH,idxPHM,ckeyPH%vals,ckeyPHM%vals
                                if (idxPH == 0 .or. idxPHM == 0) cycle

                                coef_ini => coefPH(idxPH)%val
                                coef_fin => coefPHM(idxPHM)%val
                                
                                idxD = find_dys_index(dysondict, ndys, dkey)
                                if (idxD < 0) cycle
                                sigma_r => dysondict(idxD)%val

                                print*,"B"
                                do pl = 1, nvalid
                                    call eval_spline(dys_spline(idxD), rcomp_masked(pl), spldys(pl))
                                end do

                                print*,"C"
                                dymat = 0.d0
                                do pl=1,nvalid
                                    dymat(pl,pl) = spldys(pl)
                                    ! print*,pl,spldys(pl),vib_PH(pl),vib_PHM(pl)
                                enddo
                                
                                print*,size(dymat,1),size(dymat,2)
                                print*,"D"
                                call bra_ket_coef(nvalid, numvibPH, vib_PH, vib_PHM, dymat, &
                                                coef_ini, coef_fin, bk_val)

                                print*,"E"
                                E_PHM = (PHMener(j,l,n)+ZPEPHM)/Ehtocm
                                E_PHM = E_PHM + PHMsttsE

                                diff = abs(E_PH - E_PHM)

                                Jfval = indexPHM(j,l,n,2) ; OmC = indexPHM(j,l,n,3)

                                rotcoeff = 0.d0
                                print*,"D"
                                do djJ=1,n
                                    if ((abs(Jival-Jfval) > djlist(djJ)) .or. (djlist(djJ) > (Jival+Jfval))) cycle
                                    if (abs((Jival + Jfval + (djlist(djJ)) - nint(Jival + Jfval + (djlist(djJ))))) > 1d-8) cycle
                                    do djO=1,n
                                        if (-djlist(djJ) > djlist(djO).and.(djlist(djO) > djlist(djJ))) cycle
                                        if (abs((OmN + OmC + (djlist(djO)) - nint(OmN + OmC + (djlist(djO))))) > 1d-8) cycle
                                        call W3_exp(Jival,Jfval,djlist(djJ),OmN,OmC,-djlist(djO),rotcoeff)
                                        rotcoeff = rotcoeff + abs(degvi*degJi*(rotcoeff*abs(bk_val))**2)
                                    enddo
                                enddo

                                relInt(i,j,k,l,m,n) = rotcoeff
                                evals(i,j,k,l,m,n) = diff

                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine

    !---------------------------------------------------------------------------!
    ! Writing out files
    !---------------------------------------------------------------------------!

    subroutine dump_spectrum_long(fname,numvN,numvC,nJN,nJC,nOmN,nOmC,evals,&
                                relInt,indexN,indexC,energy_unit,tol_I)
        implicit none
        
        character(len=*), intent(in) :: fname
        integer, intent(in) :: numvN,numvC,nJN,nJC,nOmN,nOmC
        real(kind=8), intent(inout) :: evals(numvN,nJN,nOmN,numvC,nJC,nOmC)
        real(kind=8), intent(in) :: relInt(numvN,nJN,nOmN,numvC,nJC,nOmC)
        integer, intent(in) :: indexN(numvN,nJN,nOmN,7), indexC(numvC,nJC,nOmC,7)
        character(len=*), intent(in) :: energy_unit
        real(kind=8), intent(in) :: tol_I

        integer :: i,j,k,l,m,n,iunit
        
        open(newunit=iunit,file=fname,action='write')

        do i=1,numvN
            do j=1,nJN
                do k=1,nOmN
                    do l=1,numvC
                        do m=1,nJC
                            do n=1,nOmC
                                if (abs(relInt(i,j,k,l,m,n)) < tol_I) cycle

                                if (energy_unit == "eV") then
                                    evals(i,j,k,l,m,n) = evals(i,j,k,l,m,n)*27.2114
                                else if (energy_unit == "cm") then
                                    evals(i,j,k,l,m,n) = evals(i,j,k,l,m,n)*219474.63
                                endif
                                    
                                write(iunit,fmt='(6(f2.1,1x),1x,6(f2.1,1x),1x,f12.6,1x,f12.6)') &
                                    indexN(i,j,k,1),indexN(i,j,k,2),indexN(i,j,k,3), &
                                    indexN(i,j,k,4),indexN(i,j,k,5),indexN(i,j,k,6), &
                                    indexC(l,m,n,1),indexC(l,m,n,2),indexC(l,m,n,3),&
                                    indexC(l,m,n,4),indexC(l,m,n,5),indexC(l,m,n,6), &
                                    evals(i,j,k,l,m,n),relInt(i,j,k,l,m,n)

                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        
        close(iunit)
    end subroutine

    !---------------------------------------------------------------------------!
    ! Utility functions
    !---------------------------------------------------------------------------!

    integer function count_values(line)
        implicit none
        character(len=*), intent(in) :: line
        integer :: i
        count_values = 0
        do i = 1, len_trim(line)
            if (line(i:i) /= ' ' .and. (i==1 .or. line(i-1:i-1) /= ' ')) then
                count_values = count_values + 1
            endif
        enddo
    end function

end module intCalcs