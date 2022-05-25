!=============================================================================!
        program shoot
!=============================================================================!
        use global
        implicit none
        integer :: i, nx
        integer :: imean=10, iparm=11, iout=12
!=============================================================================!

!.... input parameters

        call input

!.... initialize the mean flow

        call initmean(imean)

!.... allocate space for the eigenfunction and adjoint solutions

        allocate( efun(neq,ny), adj(neq,ny) )

!.... read the station information

        open(iparm,file='parm.dat',status='old',err=1000)
        read(iparm,*,err=1000) nx

        if (npcalc.eq.1) goto 10     ! Just compute nonparallel correction

!.... open the output file

        open(iout,file='output.dat',form='unformatted')
        write(iout) nx, ny, ymax, Ma, Re, Pr

!.... Polish the eigenvalues and solve the adjoint problem
!.... Everything you need to compute the nonparallel correction is
!.... written to the output.dat file.

        do i = 1, nx
          read(iparm,*,err=1000) iver, sl, alphar, alphai
          alpha = cmplx(alphar, alphai)
          write(*,"(/,80('-'))")
          write(*,"('Index: ',i3,10x,'s = ',1pe13.6)") iver, sl
          write(*,"(80('-'))")
          call mean(imean, iver)
          call solve
          call adjsolv
          call output(iout)
        end do
        close(imean)
        close(iparm)
        close(iout)

        if (npcalc .eq. 2) call exit(0)  ! don't compute nonparallel correction

  !.... go back through the results and compute the non-parallel growth-rates

  10    call nonpar

        call exit(0)

  1000  write(*,*) 'Error reading parm.dat'
        call exit(1)

        end program shoot
