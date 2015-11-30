program projector_calc
  
  use bspline_oo_module
  use, intrinsic :: iso_fortran_env, only: wp => real64
  
  implicit none
  
  ! Configuration
  integer, parameter :: TITLE_MAX_LENGTH = 32
  integer, parameter :: FILENAME_MAX_LENGTH = 64
  integer, parameter :: BUF_LENGTH = 64
  real(wp), parameter :: E_T = 0.0_wp
  real(wp), parameter :: grid_z = 0.0_wp
  
  ! Global variables
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
  real(wp) :: grid_dx, grid_dy
  real(wp), allocatable :: proj(:, :)
  real(wp) :: tau
  character(BUF_LENGTH) :: tau_str
  type(bspline_2d) :: pot_ee
  type(bspline_2d), allocatable :: pot_eZ(:)
  
  call read_input()
  call load_potential()
  call evaluate_proj()
  call output_proj()

  contains

  subroutine output_proj()
    ! Output calcualted proj to {title}.proj

    implicit none

    character(FILENAME_MAX_LENGTH) :: filename
    integer :: i, j
    real(wp) :: tx, ty

    write(filename, '(A, ".proj")') trim(title)
    
    open(2, file = trim(filename))

    write (2, '(3A20)') "x", "y", "proj"
    
    do i = 1, num_grid
      tx = grid_dx * (i - 1) + grid_xmin
      do j = 1, num_grid
        ty = grid_dy * (j - 1) + grid_ymin
        write (2, '(3ES20.10)') tx, ty, proj(i, j)
      end do
    end do

    close(2)

    write (*, '("Results Saved To: ", A)') trim(filename)

  end subroutine output_proj
  
  subroutine evaluate_proj()
    ! Calculate position corresponding to each grid point
    ! Change the position of last electron to position calculated
    ! Call evalutate_proj_sub to evaluate proj
    ! Save the result to proj

    implicit none

    integer :: i, j
    real(wp) :: tx, ty
    integer :: elec_id
    
    write (*, '("Evaluate Projector...")')

    elec_id = num_elec

    allocate(proj(num_grid, num_grid))

    grid_dx = (grid_xmax - grid_xmin) / num_grid
    grid_dy = (grid_ymax - grid_ymin) / num_grid

    do i = 1, num_grid
      tx = grid_xmin + grid_dx * (i - 1)
      do j = 1, num_grid
        ty = grid_ymin + grid_dy * (j - 1)
        pos_elec_final(num_elec, 1: 3) = (/tx, ty, grid_z/)
        proj(i, j) = evaluate_proj_sub()
      end do
    end do

  end subroutine evaluate_proj
    
  function evaluate_proj_sub()
    ! Evaluate proj corresponding to the current configuration of particles

    implicit none
    
    real(wp) :: evaluate_proj_sub
    real(wp) :: u_ee
    real(wp), allocatable :: g_eZ(:, :)
    real(wp) :: g_eZ_det
    integer :: i, j
    
    allocate(g_eZ(num_elec, num_elec))
    do i = 1, num_elec
      do j = 1, num_elec
        g_eZ(i, j) = get_g_eZ(i, j)
      end do
    end do
    call matinv(g_eZ, num_elec, g_eZ_det)

    u_ee = get_u_ee()
    
    evaluate_proj_sub = g_eZ_det * exp(E_T - u_ee)
    
    ! write (*, *) evaluate_proj_sub, g_eZ_det, u_ee
    
    
  end function evaluate_proj_sub
  
  function get_g_eZ(elec_to, elec_from)
    
    implicit none
    
    integer :: elec_to, elec_from
    real(wp) :: get_g_eZ
    integer :: i, type_id, type_cnt
    real(wp) :: q, s, t
    real(wp) :: r_jp_a(3), r_i_a(3)
    integer :: iflag

    get_g_eZ = 1.0_wp
    
    type_id = 1
    type_cnt = 0
    
    do i = 1, tot_nuclei
      r_jp_a = pos_elec_final(elec_to, 1:3) - pos_nuclei(i, 1:3)
      r_i_a = pos_elec_init(elec_from, 1:3) - pos_nuclei(i, 1:3)
      q = (norm2(r_jp_a) + norm2(r_i_a)) / 2
      s = norm2(r_jp_a - r_i_a)
      call pot_eZ(type_id)%evaluate(q, s, 0, 0, t, iflag)
      get_g_eZ = get_g_eZ * t

      type_cnt = type_cnt + 1

      if (type_cnt == num_nuclei(type_id)) then
        type_id = type_id + 1
        type_cnt = 0
      end if
    end do
    
  end function get_g_eZ
  
  function get_u_ee()
    ! Evaluate u_ee according to:
    ! Umrigar 2015. "Observation on variational and projector Monte Carlo Method", Eq. 18.

    implicit none

    real(wp) :: get_u_ee
    integer :: i, j
    real(wp) :: q, s
    real(wp) :: t
    real(wp) :: r_ij(3), r_ij_prime(3)
    integer :: iflag
    
    get_u_ee = 0.0_wp

    do i = 1, num_elec
      do j = i + 1, num_elec
        r_ij = pos_elec_init(j, 1:3) - pos_elec_init(i, 1:3)
        r_ij_prime = pos_elec_final(j, 1:3) - pos_elec_final(i, 1:3)
        q = (norm2(r_ij) + norm2(r_ij_prime)) / 2
        s = norm2(r_ij - r_ij_prime)
        call pot_ee%evaluate(q, s, 0, 0, t, iflag)
        get_u_ee = get_u_ee + t
      end do
    end do
    
  end function get_u_ee

  subroutine load_potential()

    implicit none
    
    character(FILENAME_MAX_LENGTH) :: file_e, file_Z
    integer :: i
    
    write (*, '("==== Start Loading Potential ====")')
    file_e = get_filename(z_nuclei(1), -1)
    call load_potential_sub(file_e, pot_ee)
    
    allocate(pot_eZ(num_type))
    do i = 1, num_type
      file_Z = get_filename(z_nuclei(i), z_nuclei(i))
      call load_potential_sub(file_Z, pot_eZ(i))
    end do

    write (*, '("==== Finished Loading Potential ====")')
    
  end subroutine load_potential
  
  subroutine load_potential_sub(filename, bspline_obj)
    implicit none
    
    integer :: grid_in
    real(wp) :: r_0, r_n
    integer, parameter :: LINE_SKIP = 2
    real(wp), allocatable :: x_in(:), y_in(:), u_in(:, :)
    real(wp) :: tx, ty, tu ! temp variables
    character(FILENAME_MAX_LENGTH) :: filename
    integer :: i, j
    type(bspline_2d) :: bspline_obj
    integer :: kx, ky
    integer :: iflag

    write (*, '("Loading Potential From: ", A)') trim(filename)
    
    ! read potential from filename
    open (1, file = trim(filename))
    rewind(1) ! move to the beginning of file
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
        if (i == 1) then
          y_in(j) = ty
        end if
        if (j == 1) then
          x_in(i) = tx
        end if
      end do
    end do
    
    close(1)

    ! spline interpolation
    kx = 3 ! use cubic spline
    ky = 3
    call bspline_obj%initialize(x_in, y_in, u_in, kx, ky, iflag)
    
    if(iflag == 1) then
      write (*, '("Interpolation Succeeded!")')
    else
      write (*, '("Interpolation Failed!")')
      call exit(0)
    end if
  end subroutine load_potential_sub

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
    character(FILENAME_MAX_LENGTH) :: buf1, buf2

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

    integer :: i

    write (*, '("==== Start Reading Input ====")')

    ! read info of nucleis
    read (*, '(A)') title
    write (*, '("Title: ", A)') title
    read (*, *) num_type
    write (*, '("Number of Types: ", A)') trim(str(num_type))
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
    write (*, '("Numbber of Electrons: ", A)') trim(str(num_elec))
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
    write (*, '("Number of Grid Points: ", A)'), trim(str(num_grid))
    write (*, '("x_min, x_max: ", 2F10.6)') grid_xmin, grid_xmax
    write (*, '("y_min, y_max: ", 2F10.6)') grid_ymin, grid_ymax
    

    ! read tau
    read (*, '(A)') tau_str
    read (tau_str, *) tau
    write (*, '("Tau: ", F10.6)') tau

    write (*, '("==== Finished Reading Input ====")')  

  end subroutine read_input

end program projector_calc
