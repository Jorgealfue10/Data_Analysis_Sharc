program calcintsDUO
    use intCalcs
    implicit none

    integer :: npoints, nvib, ios
    real(kind=8), allocatable :: vibmat(:,:)
    character(len=256) :: fname

    fname="/home/jorgebdelafuente/Doctorado/Photoion/DUO/PHPHM/PH/vibeigenvect_vib.chk"
    ! write(*,*) 'Enter file name: '
    ! read(*,'(A)') fname

    ! write(*,*) 'Number of points: '
    ! read(*,*) npoints
    ! write(*,*) 'Number of vibrational states: '
    ! read(*,*) nvib

    npoints = 10001
    nvib = 35

    allocate(vibmat(npoints,nvib))

    call parse_duo_vib_einfun(fname,npoints,nvib,vibmat)

    deallocate(vibmat)

end program calcintsDUO