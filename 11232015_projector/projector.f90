program projector_calc
  
  use, intrinsic :: iso_fortran_env, only: wp => real64

  implicit none
  
  ! Configuration
  integer, parameter :: TITLE_MAX_LENGTH = 32
  integer, parameter :: FILENAME_MAX_LENGTH = 64
  integer, parameter :: BUF_LENGTH = 64
  character(TITLE_MAX_LENGTH) :: title
  integer :: num_type ! number of types of nucleis
  integer :: tot_nuclei
  integer, allocatable :: num_nuclei(:) ! number of nucleis each type
  integer, allocatable :: z_nuclei(:) ! Z of each type
  real(wp), allocatable :: pos_nuclei(:, :) ! position of nucleis
  integer :: num_elec
  real(wp), allocatable :: pos_elec_init(:, :) ! initial position of electrons
  real(wp), allocatable :: pos_elec_final(:, :) ! final position of electrons
  integer :: num_grid
  real(wp) :: grid_xmin, grid_xmax, grid_ymin, grid_ymax
  real(wp) :: tau
  character(BUF_LENGTH) :: tau_str
  integer :: i, j, k
  
  ! Potential Input
  character(FILENAME_MAX_LENGTH) :: file_e, file_Z
  character(BUF_LENGTH) :: buf1, buf2, buf3

  ! Other Variables

  call read_input()
  call read_potential()

  contains

  subroutine read_potential()

    implicit none

    integer :: grid_in
    real(wp) :: r_0, r_n
    integer, parameter :: LINE_SKIP = 2  ! skip lines at the beginning
    real(wp), allocatable :: x_in(:), y_in(:), u_in(:, :)
    real(wp) :: tx, ty, tu  ! temporary variables
    
    ! Read electron potentials

    file_e = get_filename(2, -1)
    
    write (*, '("Reading Data From: ", A)') trim(file_e)

    open (1, file = trim(file_e))

    read (1, *) grid_in, r_0, r_n
    allocate(x_in(grid_in))
    allocate(y_in(grid_in))
    allocate(u_in(grid_in, grid_in))
    
    do i = 1, LINE_SKIP
      read (1, *)
    end do
    
    do i = 1, grid_in
      do j = 1, grid_in
        read (1, *) tx, ty, tu
        u_in(i, j) = tu
        if(i == 1) then
          y_in(j) = ty
        end if 
        if(j == 1) then
          x_in(i) = tx
        end if
      end do
    end do
    
    file_Z = get_filename(2, 2)
    
    write (*, '("Reading Data From: ", A)') trim(file_Z)

  end subroutine read_potential

  function str(num_in)
    integer, intent(in) :: num_in
    character(BUF_LENGTH) :: str
    write (str, *) num_in
    str = adjustl(str)
  end function str

  function get_filename(Z_nuclei, Z_in)
    
    implicit none

    integer :: Z_nuclei, Z_in
    character(FILENAME_MAX_LENGTH) :: get_filename

    ! obtain folder name
    write (buf1, '("Z=", A, "_tau=", A, "_grid=200")'), trim(str(Z_nuclei)), trim(tau_str)

    ! obtain file name
    if(Z_in > 0) then
      write (buf2, '("u_Z=", A, "_tau=", A)') trim(str(Z_in)), trim(tau_str)
    else if(Z_in == -1) then
      write (buf2, '("u_e_tau=", A)') trim(tau_str)
    end if
    
    ! combine folder and file name to filename with path
    write (get_filename, '(A, "/", A)') trim(buf1), trim(buf2)

  end function get_filename

  subroutine read_input()

    implicit none

    write (*, '("==== Start Reading Input ====")')

    ! read info of nucleis
    read (*, '(A)') title
    write (*, '("Title: ", A)') title
    read (*, *) num_type
    write (*, '("Number of Types: ", I5)') num_type
    allocate(num_nuclei(num_type))
    allocate(z_nuclei(num_type))
    
    read (*, *) (num_nuclei(i), i = 1, num_type)
    read (*, *) (z_nuclei(i), i = 1, num_type)
    tot_nuclei = 0
    do i = 1, num_type
      tot_nuclei = tot_nuclei + num_nuclei(i)
      write (*, '("Type, Count, Z: ", 3I5)') i, num_nuclei(i), z_nuclei(i)
    end do
    
    allocate(pos_nuclei(tot_nuclei, 3))
    write (*, '("Nuclei Positions:")') 
    do i = 1, tot_nuclei
      read (*, *) pos_nuclei(i, 1:3)
      write (*, '("X, Y, Z: ", 3F10.6)') pos_nuclei(i, 1:3)
    end do

    ! read info of electrons
    read (*, *) num_elec
    write (*, '("Numbber of Electrons: ", I5)') num_elec
    allocate(pos_elec_init(num_elec, 3))
    allocate(pos_elec_final(num_elec, 3))
    write (*, '("Initial Position of Electrons:")')
    do i = 1, num_elec
      read (*, *) pos_elec_init(i, 1:3)
      write (*, '("X, Y, Z: ", 3F10.6)') pos_elec_init(i, 1:3)
    end do
    write (*, '("Final Position of Electrons:")')
    do i = 1, num_elec - 1
      read (*, *) pos_elec_final(i, 1: 3)
      write (*, '("X, Y, Z: ", 3F10.6)') pos_elec_final(i, 1:3)
    end do

    ! read output grid
    read (*, *) num_grid
    read (*, *) grid_xmin, grid_xmax, grid_ymin, grid_ymax
    write (*, '("Number of Grid Points: ", I5)'), num_grid
    write (*, '("x_min, x_max: ", 2F10.6)') grid_xmin, grid_xmax
    write (*, '("y_min, y_max", 2F10.6)') grid_ymin, grid_ymax
    

    ! read tau
    read (*, '(A)') tau_str
    read (tau_str, *) tau
    write (*, '("Tau: ", F10.6)') tau

    write (*, '("==== Finished Reading Input ====")')  

  end subroutine read_input

end program projector_calc
