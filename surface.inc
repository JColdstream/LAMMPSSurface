! real kind precision parameters
integer, parameter :: sp = selected_real_kind(6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)

! pi
real, parameter :: pi = 3.14159265

logical :: first_read = .TRUE.

! coordinates, masses
real(dp), dimension(:), allocatable :: x, y, z, mass

! bin widths, histogram normalisation
real(dp) :: binwidth, hdistnrm

! bin number
integer(8) :: nbins

! axis id (x(1), y(2), z(3))
integer(4) :: axisid

! box length
real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax, Lx, Ly, Lz

! index and types for atoms
integer(8), dimension(:), allocatable :: nindex, atomtype

! height distribution
real(dp), dimension(:,:), allocatable :: hdist

! trajectory numbers, timestep, number of total frames, number of frames to skip, number of frames to calculate, 
integer(8) :: timestep, natom, nframe, nskip, nframecalc, nanalyse, ntypes

! system numbers
integer(8) :: ntotal

! file unit numbers
integer(8) :: inputf=10, trajf=11, massf=12, concf=20, densf=21

! input, trajectory, and mass file names
character(LEN=80) :: inputfile, outputfile, trajfile, massfile

! storing strings in the LAMMPS dump file headers
character(LEN=80), dimension(4) :: headertext
