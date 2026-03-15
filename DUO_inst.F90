module intCalcs
    implicit none

    type coef_key
        real(kind=8) :: vals(6)
    end type coef_key

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
        complex(kind=8), allocatable :: val(:)
        logical :: used = .false.
    end type
    
    contains
    
    !---------------------------------------------------------------------------!
    ! Read head of eigenvib vect
    !---------------------------------------------------------------------------!
    logical function is_head_eigenvib(line)
        implicit none

        character(len=*), intent(in) :: line
        integer :: a 
        real(kind=8) :: b

        read(line,*,err=999) a,b

        is_head_eigenvib = .true.
        
        999 is_head_eigenvib = .false.
    end function

    !---------------------------------------------------------------------------!
    ! Parsea el archivo de autofunciones de vibración de DUO
    !---------------------------------------------------------------------------!

    subroutine parse_duo_vib_einfun(fname,npoints,nvib,vibmat)
        implicit none

        character(len=*), intent(in) :: fname
        integer, intent(in) :: npoints,nvib
        real(kind=8), intent(out) :: vibmat(npoints,nvib)

        character(len=256) :: line
        integer :: iunit, ios, ncurrvib, ncurrpoint, auxint
        real(kind=8) :: val, eval

        ncurrpoint = 1 ; ncurrvib = 1

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

    subroutine read_J_values(fname,nvib,nj,mult,nLambda,nOmega,Jvals,index_list,key_dict)
        implicit none

        character(len=*), intent(in) :: fname
        integer, intent(in) :: nvib, nj, mult, nLambda, nOmega
        real(kind=8), intent(out) :: Jvals(nvib,nj,nOmega), index_list(nvib,nj,nOmega,7),key_dict(nvib,nj,nOmega,6)

        integer :: iunit, ios, i, Smult
        character(len=256) :: line, aux
        character(len=1) :: parity
        real(kind=8) :: J, eval, stt, vibn, lambda, sval, sigma, omega, pm, index

        Smult = (mult-1)/2

        open(newunit=iunit,file=fname,status='old',action='read')        

        do 
            read(iunit,'(A)', iostat=ios) line
            if (ios /= 0 ) exit

            read(line,*) J, index, eval, stt, vibn, lambda, sval, sigma, omega, parity, aux, aux
            if (vibn < nvib) then
                Jvals(int(vibn)+1, int(J)+1, int(omega)+1) = eval

                if (parity.eq."-") then
                    pm = -1.d0
                else
                    pm = 1.d0
                endif

                index_list(int(vibn)+1, int(J)+1, int(omega+Smult)+1,1) = vibn
                index_list(int(vibn)+1, int(J)+1, int(omega+Smult)+1,2) = J
                index_list(int(vibn)+1, int(J)+1, int(omega+Smult)+1,3) = omega
                index_list(int(vibn)+1, int(J)+1, int(omega+Smult)+1,4) = sigma
                index_list(int(vibn)+1, int(J)+1, int(omega+Smult)+1,5) = lambda
                index_list(int(vibn)+1, int(J)+1, int(omega+Smult)+1,6) = pm
                index_list(int(vibn)+1, int(J)+1, int(omega+Smult)+1,7) = index

                key_dict(int(vibn)+1, int(J)+1, int(omega+Smult)+1,1) = J
                key_dict(int(vibn)+1, int(J)+1, int(omega+Smult)+1,2) = omega
                key_dict(int(vibn)+1, int(J)+1, int(omega+Smult)+1,3) = sigma
                key_dict(int(vibn)+1, int(J)+1, int(omega+Smult)+1,4) = lambda
                key_dict(int(vibn)+1, int(J)+1, int(omega+Smult)+1,5) = pm
                key_dict(int(vibn)+1, int(J)+1, int(omega+Smult)+1,6) = index
            endif
        enddo

        close(iunit)

    end subroutine

    !---------------------------------------------------------------------------!
    ! Reading coefficients
    !---------------------------------------------------------------------------!

    logical function same_key(a,b)
        implicit none
        type(coef_key), intent(in) :: a,b
        same_key = all((a%vals - b%vals) < 1d-8)
    end function

    subroutine append_coef(coef_dict, ndict, nvib, key, coeff)
        implicit none
        integer, intent(inout) :: ndict
        integer, intent(in) :: nvib
        type(coef_dictt), intent(inout) :: coef_dict(ndict)
        type(coef_key), intent(in) :: key
        real(kind=8), intent(in) :: coeff

        integer :: i,n
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
            coef_dict(i)%used = .true.
            coef_dict(i)%key = key
            coef_dict(i)%nread = 1
            allocate(coef_dict(ndict)%val(nvib))
            coef_dict(i)%val = cmplx(0.d0,0.d0,kind=8)
            coef_dict(i)%val(1) = cmplx(coeff,0.d0,kind=8)
        endif

    end subroutine

    subroutine read_coefficients(fname,nvib,nj,mult,nOmega,coef_dict,ndict)
        implicit none

        character(len=*), intent(in) :: fname
        integer, intent(in) :: nvib, nj, mult, nOmega
        type(coef_dictt), allocatable, intent(out) :: coef_dict(:)
        integer, intent(out) :: ndict

        integer :: iunit, ios, maxdict
        character(len=256) :: line, aux
        real(kind=8) :: J, coeff, stt, vibn, lambda, sval, sigma, omega, pm, index
        type(coef_key) :: key

        maxdict = nvib*nj*nOmega
        allocate(coef_dict(maxdict))
        ndict = 0

        open(newunit=iunit,file=fname,status='old',action='read')

        do 
            read(iunit,'(A)', iostat=ios) line
            if (ios /= 0 ) exit

            read(line,*,iostat=ios) index, J, pm, coeff, aux, vibn, lambda, sval, sigma, omega, aux
            if (ios /= 0 ) cycle

            key%vals = [J,omega,sigma,lambda,pm,index]

            call append_coef(coef_dict, ndict, nvib, key, coeff)
        enddo

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

    subroutine read_dys_bySigma(d1,d2,nd1,nd2,path,dict,ndict,r)
        implicit none

        integer, intent(in) :: nd1,nd2
        integer, intent(in) :: d1(nd1),d2(nd2)
        character(len=*), intent(in) :: path
        type(dys_dictt), allocatable, intent(out) :: dict(:)
        integer, intent(out) :: ndict
        real(kind=8), allocatable, intent(out) :: r(:)

        integer :: i,j,nr
        real(kind=8) :: sigma_i, sigma_j
        real(kind=8), allocatable :: dyson(:)
        type(dys_key) :: key

        character(len=256) :: fname

        allocate(dict(nd1*nd2))
        ndict = 0

        do i = 1, nd1
            do j = 1, nd2
                sigma_i = d1(i)
                sigma_j = d2(j)
                key%vals = [sigma_i,sigma_j]
                
                write(fname,'(A,"/dyson_",I2.2,"_",I2.2,".dat")') trim(path),i,j
                
                call read_dyson_file(fname,r,dyson,nr)

                call append_sigma(dict, ndict, key, nr, dyson)

            enddo
        enddo

        do i = 1, ndict
            dict(i)%val = sqrt(dict(i)%val)
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

    subroutine intT_sigma(PHener,PHMener,rvals,ndys,npoints,dysondict,PHvibs,PHMvibs,mask, &
        Tvib,Trot,ZPEPH,ZPEPHM,PHsttsE,PHMsttsE,DJ,numvibPH,numvibPHM, &
        numJPH,numJPHM,numOmPH,numOmPHM,indexPH,indexPHM, &
        keylistPH,keylistPHM,coefPH,coefPHM,evals,relInt)

        implicit none
        integer, intent(in)  :: numvibPH,numvibPHM,numJPH,numJPHM,numOmPH,numOmPHM,ndys,npoints
        real(kind=8), intent(in)  :: PHener(numvibPH,numJPH,numOmPH),PHMener(numvibPHM,numJPHM,numOmPHM)
        real(kind=8), intent(in)  :: Tvib,Trot,ZPEPH,ZPEPHM,PHsttsE,PHMsttsE,DJ
        real(kind=8), intent(in)  :: rvals(npoints),mask
        real(kind=8), intent(in)  :: PHvibs(numvibPH,npoints),PHMvibs(numvibPHM,npoints)
        real(kind=8), intent(in)  :: indexPH(numvibPH,numJPH,numOmPH,7),indexPHM(numvibPHM,numJPHM,numOmPHM,7)
        real(kind=8), intent(in)  :: keylistPH(numvibPH,numJPH,numOmPH,6),keylistPHM(numvibPHM,numJPHM,numOmPHM,6)
        type(dys_dictt), intent(in)  :: dysondict(ndys)
        type(dys_key) :: dkey
        type(coef_dictt), intent(in)  :: coefPH,coefPHM
        type(coef_key)  :: ckeyPH,ckeyPHM

        real(kind=8), intent(out)  :: evals(numvibPH,numvibPHM,numJPH,numJPHM,numOmPH,numOmPHM),relInt(numvibPH,numvibPHM,numJPH,numJPHM,numOmPH,numOmPHM)

        integer :: i,j,k,l,m,n
        real(kind=8) :: normPHvib,normPHrot,Ehtocm

        Ehtocm = 219474.6

        n = int(2*DJ) + 1
        allocate(djlist(n))

        do i=0,int(2*DJ)
            djlist(i+1) = real(i) - DJ
        enddo

        call vib_normf(Tvib,PHener,ZPEPH,PHsttsE,numvibPH,numJPH,numOmPH,normPHvib)
        call rot_normf(Trot,PHener,ZPEPH,PHsttsE,numvibPH,numJPH,numOmPH,indexPH,normPHrot)

        do i=1,numvibPH

            E_PH = (PHener(i,1,2)+ZPEPH)/Ehtocm
            evib = E_PH - ((min(PHener(:,1,2)) + ZPEPH)/Ehtocm)
            degvi = exp(-evib/(3.166811563*10**(-6)*Tvib))/normPHvib

            do j=1,numvibPHM

                ! PHvibsvn = 

                do k=1,numJPH
                    do l=1,numJPHM

                        do m=1,numOmPH

                            E_PH = (PHener(i,k,m)+ZPEPH)/Ehtocm
                            erot = E_PH - ((min(PHener(i,:,:)) + ZPEPH)/Ehtocm)
                            Jival = indexPH(i,k,m,2) ; OmN = indexPH(i,k,m,3)
                            degJi = ((2*Jival+1)*exp(-erot/(3.166811563*10**(-6)*Trot)))/normPHvib
                            E_PH = E_PH + PHsttsE

                            ckeyPH%vals = [keylistPH(i,k,m,1),keylistPH(i,k,m,2),keylistPH(i,k,m,3),keylistPH(i,k,m,4),keylistPH(i,k,m,5),keylistPH(i,k,m,6)]

                            do n=1,numOmPHM

                                ckeyPHM%vals = [keylistPHM(j,l,n,1),keylistPHM(j,l,n,2),keylistPHM(j,l,n,3),keylistPHM(j,l,n,4),keylistPHM(j,l,n,5),keylistPHM(j,l,n,6)]
                                dkey%vals = [indexPH(i,k,m,3),indexPHM(j,l,n,3)]
                                
                                !SPLINES FOR DYSON
                                
                                call bra_ket_coef(npoints, numvibPH, vib_PH, vib_PHM, dymat, coef_PH, coef_PHM, braket)

                                E_PHM = (PHMener(j,l,n)+ZPEPHM)/Ehtocm
                                E_PHM = E_PHM + PHMsttsE

                                diff = abs(E_PH - E_PHM)

                                Jfval = indexPHM(j,l,n,2) ; OmC = indexPHM(j,l,n,3)

                                rotcoeff = 0.d0
                                do dj=1,djlist
                                    if ((abs(Jival-Jfval) > djlist(dj)) .or. (djlist(dj) > (Jival+Jfval))) cycle
                                    if (abs((Jival + Jfval + (djlist(dj)) - nint(Jival + Jfval + (djlist(dj))))) > 1d-8) cycle
                                    do djO=1,djlist
                                        if (-djlist(dj) > djlist(djO).and.(djlist(djO) > djlist(dj))) cycle
                                        if (abs((OmN + OmC + (djlist(djO)) - nint(OmN + OmC + (djlist(djO))))) > 1d-8) cycle
                                        call W3_exp(Jival,Jfval,djlist(dj),OmN,OmC,-djlist(djO),rotcoeff)
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

end module intCalcs