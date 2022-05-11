program surface

implicit none
include 'surface.inc'
integer(8):: frame

call start
do frame = 1, nframecalc
  write(*,*) frame
    call readheader
    call readcoordinates
    call heightdistribution
enddo
call distribution_output

contains

subroutine start

call get_command_argument(1, inputfile)
if (len_trim(inputfile) == 0) then
  write(*,*) 'Provide an input file.'
  call exit
endif

call readinput
! allocate arrays
allocate(nindex(ntotal))
allocate(atomtype(ntotal))
allocate(x(ntotal))
allocate(y(ntotal))
allocate(z(ntotal))

call readmasses

! open lammps trajectory file
open(11, file = trajfile, status='old')

call skipframe

! initialise counter
hdistnrm = 0.0

end subroutine start


! reads analysis.input input file
subroutine readinput
open(10, file = inputfile, status='old')
read(10, *)
read(10, *)
read(10, '(a)') trajfile
read(10, *) outputfile
read(10, *) massfile
read(10, *) nframe
read(10, *) nskip
read(10, *) ntotal
read(10, *) ntypes
read(10, *) binwidth
read(10, *) axisid
write(*,*) trajfile

if (axisid .gt. 3) then
  write(*,*) 'Axis id must take the value 1(x), 2(y) or 3(z).'
  call exit
endif

! calculate relevant quantities for later
! nanalyse = nmol*molsize
nframecalc = nframe-nskip

end subroutine readinput

! skips frames that don't need to be analysed
subroutine skipframe
integer:: i
do i = 1, nskip*(ntotal+9)
  read(11, *)
enddo
write(*,*) "SKIPPED FRAMES : ", nskip
end subroutine skipframe

subroutine readmasses
  integer(sp):: i, j

  allocate(mass(ntypes))

  open(12, file = massfile, status='old')
  read(12, *)
  read(12, *)
  do i = 1, ntypes
    read(12, *) j, mass(i)
    write(*,*) j, mass(i)
  enddo
end subroutine readmasses

! reads in box length, timestep and natoms from the LAMMPS trajectory headers
subroutine readheader
read(11, *) headertext(1)   
read(11, *) timestep
read(11, *) headertext(2)
read(11, *) natom
read(11, *) headertext(3)
read(11, *) xmin, xmax
read(11, *) ymin, ymax
read(11, *) zmin, zmax
read(11, *) headertext(4)

Lx = xmax-xmin
Ly = ymax-ymin
Lz = zmax-zmin

!if(natom  .ne. ntotal) stop 'mismatched number of atoms'
if (first_read) then
  if (axisid .eq. 1) then
    nbins = int((xmax-xmin)/binwidth)+1
  else if (axisid .eq. 2) then
    nbins = int((ymax-ymin)/binwidth)+1
  else if (axisid .eq. 3) then
    nbins = int((zmax-zmin)/binwidth)+1
  endif
  allocate(hdist(ntypes, nbins))
  hdist = 0.0
  first_read = .FALSE.
endif
end subroutine readheader

! reads in the atom types and coordinates from the LAMMPS trajectory file
subroutine readcoordinates
 integer(sp):: i
 do i = 1, ntotal 
   read(11, *) nindex(i), atomtype(i), x(i), y(i), z(i)
 enddo
 ! zero coordinates
 x = x-xmin
 y = y-ymin
 z = z-zmin
end subroutine readcoordinates

subroutine heightdistribution
integer(8):: i, hbin
if (axisid .eq. 1) then
  do i = 1, ntotal
    hbin = int(x(i)/binwidth+1.0)
    hdist(atomtype(i), hbin) = hdist(atomtype(i), hbin) + 1.0
  enddo
else if (axisid .eq. 2) then
  do i = 1, ntotal
    hbin = int(y(i)/binwidth+1.0)
    hdist(atomtype(i), hbin) = hdist(atomtype(i), hbin) + 1.0
  enddo
else if (axisid .eq. 3) then
  do i = 1, ntotal
    hbin = int(z(i)/binwidth+1.0)
    hdist(atomtype(i), hbin) = hdist(atomtype(i), hbin) + 1.0
  enddo
endif
hdistnrm = hdistnrm+1.0
end subroutine heightdistribution

subroutine distribution_output
  integer(8):: i, j
  real(dp):: binvolume
  open(unit = 101, file='conc.dat', status='unknown')
  open(unit = 102, file='density.dat', status='unknown')
  if (axisid .eq. 1) then
    binvolume = Ly*Lz*binwidth
    write(101, *) '#  x', ('    ', j, j = 1, ntypes)
    write(102, *) '#  x', ('    ', j, j = 1, ntypes)
  else if (axisid .eq. 2) then
    binvolume = Lx*Lz*binwidth
    write(101, *) '#  y', ('    ', j, j = 1, ntypes)
    write(102, *) '#  y', ('    ', j, j = 1, ntypes)
  else if (axisid .eq. 3) then
    binvolume = Lx*Ly*binwidth
    write(101, *) '#  z', ('    ', j, j = 1, ntypes)
    write(102, *) '#  z', ('    ', j, j = 1, ntypes)
  endif
  do i = 1, nbins
    write(101, *) (binwidth*(dble(i)-0.5)), (hdist(j, i)/binvolume/hdistnrm, j = 1, ntypes)
    write(102, *) (binwidth*(dble(i)-0.5)), (hdist(j, i)*mass(j)/binvolume/hdistnrm, j = 1, ntypes)
  enddo
end subroutine distribution_output

end program surface
