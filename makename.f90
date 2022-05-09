!=============================================================================!
	subroutine makename(base,iver,fname)
!
!.... put a version number on the filename
!
!=============================================================================!
	character*80 base, fname
  
	length = index(base,' ')
	fname = base
	if (iver .lt. 10) then
	  write(fname(length:80),"('.',i1)") iver
	else if (iver .lt. 100) then
	  write(fname(length:80),"('.',i2)") iver
	else if (iver .lt. 1000) then
	  write(fname(length:80),"('.',i3)") iver
	else if (iver .lt. 10000) then
	  write(fname(length:80),"('.',i4)") iver
	else
	  write(*,*) 'Error in MakeName:  iver too large'
	  call exit(1)
	end if
  
	return
	end
