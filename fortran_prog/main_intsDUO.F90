program calcintsDUO
    use intCalcs
    implicit none

    integer :: iuninp,iunout,iunit,ios
    integer :: npoints, nvib
    integer :: numvibN, numvibC, nJN, nJC, nOmN, nOmC, nLN, nLC, multN, multC
    integer :: ndictN, ndictC, ndys, ns1, ns2, nstates

    real(kind=8) :: Tvib, Trot, ZPEN, ZPEC, EN, EC, DJ

    real(kind=8), allocatable :: vibmat(:,:), vibmatN(:,:), vibmatC(:,:)
    real(kind=8), allocatable :: JvalsN(:,:,:), JvalsC(:,:,:)
    real(kind=8), allocatable :: indexN(:,:,:), indexC(:,:,:)
    real(kind=8), allocatable :: keyN(:,:,:), keyC(:,:,:)
    real(kind=8), allocatable :: rvals(:)
    
    type(coef_dictt), allocatable :: coef_dictN(:), coef_dictC(:)
    type(dys_dictt), allocatable :: dysdict(:)
    type(state_map), allocatable :: s1(:), s2(:)
    
    character(len=256) :: fname, fnameN, fnameC, fileout, path
    character(len=2) :: energy_unit
    
    logical :: mask

    open(newunit=iuninp,file='inp_fortran.inp',action='read')

    

    fnameN="/home/jorgebdelafuente/Doctorado/Photoion/DUO/PHPHM/PH/vibeigenvect_vib.chk"
    fnameC="/home/jorgebdelafuente/Doctorado/Photoion/DUO/PHPHM/PHM_GS/vibeigenvect_vib.chk"

    npoints = 10001
    nvib = 35

    allocate(vibmat(npoints,nvib))

    call parse_duo_vib_einfun(fnameN,npoints,nvib,vibmatN)
    call parse_duo_vib_einfun(fnameC,npoints,nvib,vibmatC)

    call read_J_values(fnameN,nvib,nj,multN,nLN,nOmN,JvalsN,indexN,keyN)
    call read_J_values(fnameC,nvib,nj,multC,nLC,nOmC,JvalsC,indexC,keyC)

    call read_coefficients(fnameN,nvib,nj,multN,nOmN,coef_dictN,ndictN)
    call read_coefficients(fnameC,nvib,nj,multC,nOmC,coef_dictC,ndictC)
    
    call read_dys_bySigma(s1,s2,ns1,ns2,path,dysdict,ndys,rvals)

    if size(DJvals) > 1 then

        do i=1,size(DJvals)
            DJ = DJvals(i)
            call intT_sigma(JvalsN,JvalsC,nvib,rvals,ndys,npoints,dysondict, &
                vibmatN,vibmatC,mask,Tvib,Trot,ZPEN,ZPEC,EN,EC,DJ, &
                numvibN,numvibC,numJN,numJC,numOmN,numOmC,indexN,indexC, &
                keyN,keyC,coef_dictN,coef_dictC,ndictN,ndictC,evals,relInt)
            
            call dump_spectrum_long(fileout,numvibN,numbvibC,numJN,numJC,numOmN,numOmC, &
                                    evals,relInt,indexN,indexC,energy_unit="eV",tol_I=0.0d0)
        enddo
    else
        call intT_sigma(JvalsN,JvalsC,nvib,rvals,ndys,npoints,dysondict, &
            vibmatN,vibmatC,mask,Tvib,Trot,ZPEN,ZPEC,EN,EC,DJ, &
            numvibN,numvibC,numJN,numJC,numOmN,numOmC,indexN,indexC, &
            keyN,keyC,coef_dictN,coef_dictC,ndictN,ndictC,evals,relInt)
        
        call dump_spectrum_long(fileout,numvibN,numbvibC,numJN,numJC,numOmN,numOmC, &
                                evals,relInt,indexN,indexC,energy_unit="eV",tol_I=0.0d0)
    endif

    deallocate(vibmat)

end program calcintsDUO