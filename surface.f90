program surface

implicit none
include 'surface.inc'
integer(8):: frame

call start
do frame = 1, nframecalc
  if (mod(frame, 50) .eq. 0) write(*,*) frame
    call readheader
    call readcoordinates
    call heightdistribution
enddo
close(trajf)
call distribution_output

contains

subroutine start

call get_command_argument(1, inputfile)
if (len_trim(inputfile) == 0) then
  write(*,*) 'Provide an input file.'
  call exit
endif

call readinput

! open lammps trajectory file
open(trajf, file = trajfile, status='old')

! get ntotal
call readheader
! rewind trajectory file
rewind(trajf)
call init_hdist
! allocate arrays
allocate(nindex(ntotal))
allocate(atomtype(ntotal))
allocate(x(ntotal))
allocate(y(ntotal))
allocate(z(ntotal))

call readmasses
call readscatteringlengths
call skipframe

! initialise counter
hdistnrm = 0.0

end subroutine start


! reads analysis.input input file
subroutine readinput
open(inputf, file = inputfile, status='old')
read(inputf, *)
read(inputf, *)
read(inputf, '(a)') trajfile
read(inputf, *) massfile
read(inputf, *) scatteringfile
read(inputf, *) nframe
read(inputf, *) nskip
read(inputf, *) ntypes
read(inputf, *) binwidth
read(inputf, *) axisid
read(inputf, *) fractional_coords
write(*,*) trajfile
close(inputf)

if ((axisid .gt. 3) .or. (axisid .lt. 0)) then
  write(*,*) 'Axis id must take the value 1(x), 2(y) or 3(z).'
  call exit
endif

if (fractional_coords .eq. 0) then
  fractional_bool = .FALSE.
else if (fractional_coords .eq. 1) then
  fractional_bool = .TRUE.
else
  write(*,*) 'Fractional coordinates must take the value of 0 (absolute) or 1 (fractional).'
  call exit
end if  

! calculate relevant quantities for later
! nanalyse = nmol*molsize
nframecalc = nframe-nskip

end subroutine readinput

subroutine readscatteringlengths
  integer(sp):: i, j

  allocate(b(ntypes))
  ! reads scattering lengths input
  if (scatteringfile .eq. '0') then
    do i = 1, ntypes
      b(i) = 1.0_dp
    enddo
  else
    open(12, file = scatteringfile, status='old')
    read(12, *)
    read(12, *)
    do i = 1, ntypes
      read(12, *) j, b(j)
    enddo
    do i = 1, ntypes
      b(i) = b(i) * 0.00001_dp
    enddo
  endif
end subroutine readscatteringlengths

! skips frames that don't need to be analysed
subroutine skipframe
integer:: i
do i = 1, nskip*(ntotal+9)
  read(trajf, *)
enddo
write(*,*) "SKIPPED FRAMES : ", nskip
end subroutine skipframe

subroutine readmasses
  integer(sp):: i, j

  allocate(mass(ntypes))

  open(massf, file = massfile, status='old')
  read(massf, *)
  read(massf, *)
  do i = 1, ntypes
    read(massf, *) j, mass(i)
    ! write(*,*) j, mass(i)
  enddo
  close(massf)
  write(*,*) 'Masses read.'
end subroutine readmasses

! reads in box length, timestep and natoms from the LAMMPS trajectory headers
subroutine readheader
read(trajf, *) headertext(1)   
read(trajf, *) timestep
read(trajf, *) headertext(2)
read(trajf, *) ntotal
read(trajf, *) headertext(3)
read(trajf, *) xmin, xmax
read(trajf, *) ymin, ymax
read(trajf, *) zmin, zmax
read(trajf, *) headertext(4)

Lx = xmax-xmin
Ly = ymax-ymin
Lz = zmax-zmin

end subroutine readheader

! initialize hdist
subroutine init_hdist
  if (fractional_bool) then
    nbins = int(1/binwidth) + 1
  else
    if (axisid .eq. 1) then
      nbins = int((xmax-xmin)/binwidth)+1
    else if (axisid .eq. 2) then
      nbins = int((ymax-ymin)/binwidth)+1
    else if (axisid .eq. 3) then
      nbins = int((zmax-zmin)/binwidth)+1
    endif
  end if
  allocate(hdist(ntypes, nbins))
  hdist = 0.0

end subroutine init_hdist

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
 if (fractional_bool .eqv. .TRUE.) then
   x = x / (xmax-xmin)   
   y = y / (ymax-ymin)   
   z = z / (zmax-zmin)   
 end if  
end subroutine readcoordinates

subroutine heightdistribution
  integer(8):: i, hbin
  ! transform these into a function, and call the function.
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
  open(unit = concf, file='conc.dat', status='unknown')
  open(unit = densf, file='density.dat', status='unknown')
  open(unit = sldf, file='sld.dat', status='unknown')
  ! transform these into a function, and call the function.
  if (axisid .eq. 1) then
    binvolume = Ly*Lz*binwidth
    write(concf, *) '#  x', ('    ', j, j = 1, ntypes)
    write(densf, *) '#  x', ('    ', j, j = 1, ntypes)
  else if (axisid .eq. 2) then
    binvolume = Lx*Lz*binwidth
    write(concf, *) '#  y', ('    ', j, j = 1, ntypes)
    write(densf, *) '#  y', ('    ', j, j = 1, ntypes)
  else if (axisid .eq. 3) then
    binvolume = Lx*Ly*binwidth
    write(concf, *) '#  z', ('    ', j, j = 1, ntypes)
    write(densf, *) '#  z', ('    ', j, j = 1, ntypes)
  endif

  do i = 1, nbins
    write(concf, '(f10.3, 999(e15.6))') (binwidth*(dble(i)-0.5)), &
      (hdist(j, i)/binvolume/hdistnrm, j = 1, ntypes)
    write(densf, '(f10.3, 999(e15.6))') (binwidth*(dble(i)-0.5)), &
      (hdist(j, i)*mass(j)/binvolume/hdistnrm, j = 1, ntypes)
    write(sldf, '(f10.3, 999(e15.6))') (binwidth*(dble(i)-0.5)), &
      (hdist(j, i)*b(j)/binvolume/hdistnrm, j = 1, ntypes)
  enddo

  close(concf)
  close(densf)
  close(sldf)
end subroutine distribution_output

end program surface
