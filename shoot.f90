!=============================================================================!
        program shoot
!=============================================================================!
        use global
        implicit none
        integer :: i, nx
        integer :: imean=10, iparm=11, iout=12

        character(80) :: pname, tmp
        integer :: is=0, ie=huge(0)
        namelist /parm/ pname, nx, is, ie
        logical :: isParm=.false., isPnml=.false.
!=============================================================================!

!.... input parameters

        call input

!.... initialize the mean flow

        call initmean(imean)

        write(*,'("Completed initmean...")')

!.... allocate space for the eigenfunction and adjoint solutions

        allocate( efun(neq,ny), adj(neq,ny) )

!.... read the station information
        pname = 'parm.dat'
        inquire(file='parm.nml',exist=ispnml)
        if (ispnml) then
          open(iparm,file='parm.nml',status='old',err=1000)
          read(iparm,nml=parm)
          write(*,nml=parm)
          close(iparm)
          open(iparm,file=pname,status='old',err=1000)
          write(*,'("Reading station information from ",a,$)') pname
          nx = 0
  40      continue
          read(iparm,'(a)',end=48) tmp
          if(tmp(1:1).ne.'#') then
            read(tmp,*) iver, sl, alphar, alphai
            if (iver.lt.is .or. iver.gt.ie) goto 40
            nx = nx + 1
            goto 40
          else
            goto 40
          endif
  48      write(*,'("  with Nx = ",i4," stations")') nx
          rewind(iparm)
        else
          pname = 'parm.dat'
          inquire(file=pname,exist=isparm)
          if (isparm) then
            nx = 0
            write(*,'("Reading station information from parm.dat",$)')
            open(iparm,file=pname,status='old',err=1000)
            !read(iparm,*,err=1000) nx
  20        continue
            read(iparm,'(a)',end=30) tmp
            if (tmp(1:1).ne.'#') nx = nx + 1
            goto 20
  30        continue
            write(*,'("  with Nx = ",i4," stations")') nx
            rewind(iparm)
          else
            write(*,*) "ERROR: Either parm.dat or parm.nml must exist"
            call exit(1)
          endif 
        endif
        
        if (npcalc.eq.1) goto 10     ! Just compute nonparallel correction

!.... open the output file

        open(iout,file='output.dat',form='unformatted')
        write(iout) nx, ny, ymax, Ma, Re, Pr

!.... Polish the eigenvalues and solve the adjoint problem
!.... Everything you need to compute the nonparallel correction is
!.... written to the output.dat file.

        i = 1
  45    continue 
          read(iparm,'(a)',end=50) tmp
          !write(*,*) "tmp = ", tmp
          if(tmp(1:1).ne.'#') then
            read(tmp,*) iver, sl, alphar, alphai
            if (iver.lt.is .or. iver.gt.ie) then
              !write(*,*) "skipping : ", i, iver, sl, alphar, alphai
              goto 45
            endif
            i = i + 1
          else
            !write(*,*) "Skipping : ", tmp
            goto 45
          endif
          alpha = cmplx(alphar, alphai)
          write(*,"(/,80('-'))")
          write(*,"('Index: ',i3,10x,'s = ',1pe13.6)") iver, sl
          write(*,"(80('-'))")
#ifdef VERBOSE
          write(*,'("call mean()")')
#endif
          call mean(imean, iver)
#ifdef VERBOSE
          write(*,'("call solve()")')
#endif
          call solve
#ifdef VERBOSE
          write(*,'("call adjsolv()")')
#endif
          call adjsolv
#ifdef VERBOSE
          write(*,'("call output())")')
#endif
          call output(iout)
        goto 45
  50    continue

        close(imean)
        close(iparm)
        close(iout)

        if (npcalc .eq. 2) call exit(0)  ! don't compute nonparallel correction

!.... go back through the results and compute the non-parallel growth-rates

  10    call nonpar

        deallocate( efun, adj )

        call exit(0)
 1000   write(*,'("Error reading ",a)') pname
        call exit(1)
        end program shoot
