 Program Diff

   implicit none


   integer,parameter    :: spmax=5

   integer, allocatable :: OHindex(:)     !Index of the O*

   integer              :: maxstep,    &  !Max Number of production steps (actual)
                           nsp(spmax), &  !Total number of species, for each type nspt
                           mint,       &  !Sampling interval
                           nspt,       &  !Total number of species
                           aindex,     &  !Atom index to be considered displacement, tagged
                           Nat,        &  !Total number of atoms for aindex
                           nsample,    &  !Total Number of samples
                           ndata,      &  !Total Number of datapoints
                           norigin,    &  !Total number of origins to consider in the average
                           nmin,nmax,  &  !Min/Max intervals for each origin's average
                           i,j,k,n,    &  !simple indexes   
                           noOH=0,     &  !THe number of configurations with no OH  
                           dblOH=0,    &  !THe number of configurations with 2  OH  
                           multiOH=0      !The number of configurations with >2 OH   

   real, allocatable    :: tau(:,:),   &  !Atom positions tau(atom-index,dim), for aindex
                                          ! Irrelevant species will be trimed as
                                          ! *.pos file is read 
                           msd(:),     &  !Total mean squared displacement
                           time(:)        !time vector

   real                 :: dx, dy, dz, &  !Displacements 
                           dt,         &  !time step for ONE interval (in fs)
                           dummy,      &  !A dummy variable
                           dum(3),     &  !A dummy array
                           OHbond         !Parameter determining the OH Bond length


   character(len=15)     :: posfile, celfile, outfile

   

   Namelist /input/ maxstep, mint, dt, nspt, nsp, aindex, posfile, &
                     outfile, celfile, OHbond

   read(*,input)

   if( aindex /= 1) then
      write(*,*) ' ERROR: this MSD program will only work with aindex = O* (1) ... '
      stop
   endif

   !Initial Values 
   nsample = maxstep/mint + 1
   !ndata =  nsp(aindex)*nsample
   ndata = nsample

   call print_out(0)


   !Number of Origins is half of the total nummber of samples (configurations)
   norigin = (nsample-1)/2
   !Number of min max samples
   nmin  =  50
   nmax  =  norigin
   !Number of atoms in consideration (tagged)
   Nat = nsp(aindex)

   if(nmin > nmax) then
      write(*,*) ' ERROR: Not enough origins! '
      write(*,*)  ' nmin = ', nmin,  'nmax = ', nmax
      stop
   endif

   allocate( tau(1:ndata,1:3), time(1:norigin), msd(1:norigin) )

   !Open Files
   open(unit=1,file=outfile,form='formatted',status='unknown')
   open(unit=2,file=posfile,form='formatted',status='old')
   open(unit=3,file=celfile,form='formatted',status='old')

   !Find OH indexes
   allocate ( OHindex(nsample) )
   call findOH


   !Read in the correct atoms of aindex, Ostar only
   ! "trim the fat"
   rewind(2)
   n = 0
   do i=1,nsample,1
      read(2,*) dummy, dummy
      do j=1,nspt,1
         do k=1,nsp(j),1
            if( (j==aindex) .and. (OHindex(i) == k) ) then
               n = n + 1
               read(2,*) tau(n,1:3)
               write(98,*) tau(n,1:3)
            else
               read(2,*) dum(1:3)
            endif 
         enddo
      enddo
   enddo
   call print_out(1)

   !calculate the time vector in picosecond
   do i=1,norigin,1
      time(i) = real(i*mint)*dt*0.001
   enddo
   
   msd(1:norigin) = 0.0
   
   call print_out(2)

      !--------------------------------------------------------------
      !Loop over the Origins (Atoms in different configurations)
      !--------------------------------------------------------------
      !Set the Origin to the "Next" atom
      do j=1,norigin,1
         !Loop over the intervals
         do k=nmin,nmax,1
            dx = tau(k,1) - tau(j,1)
            dy = tau(k,2) - tau(j,2)
            dz = tau(k,3) - tau(j,3)
            msd(k) = msd(k) + SQRT(dx**2 + dy**2 + dz**2)
         enddo   
      enddo

   call print_out(3)

   msd(:) = msd/(real(Nat*norigin))

   call print_out(11)   
 
   call print_out(10)

   call print_out(4)

   deallocate(time, msd, tau)

   contains
      
      !---------------------------------------------------------------------------------------
      !Print-out Subroutine
      !io selects from a variety of output functions
      !---------------------------------------------------------------------------------------
      Subroutine print_out(io)
   
         implicit none

         integer, intent(in)  :: io

         select case (io)
            case(0)
                  write(*,*)
                  write(*,*)  '---------------------------------------------------------'
                  write(*,*)  '          Mean Sqaured Displacment Program               '
                  write(*,*)  '                                                         '
                  write(*,*)  '                Hydroixde analysis                       '
                  write(*,*)  '                                                         '
                  write(*,*)  '               Charles W. Swartz VI                      '
                  write(*,*)  '---------------------------------------------------------'
                  write(*,*)  ' nsample = ', nsample, ' ndata = ', ndata
            case(1)
                  write(*,*)  ' Data from ', trim(posfile), ' read in correctly.'
            case(2)
                  write(*,*)
                  write(*,*) ' Begining main loop (this may take a few moments) ...'
            case(3)
                  write(*,*) ' ... Main loop complete!'
            case(4)
                  write(*,*)
                  write(*,*) ' MSD data output-file : ', outfile
                  write(*,*) ' Program Complete!'
                  write(*,*)
            case(10)
                  write(*,*) 
                  write(*,*) 'Total number of configurations with NO  OHminus : ', noOH
                  write(*,*) 'Total number of configurations with TWO OHminus : ', dblOH
                  write(*,*) 'Total number of configurations with >2  OHminus : ', multiOH
                  write(*,*)
            case(11)
                  write(1,*) '#This file contains time and MSD values for ',posfile
                  do i=1,norigin,1
                     write(1,100) time(i), msd(i)
                  enddo

            case default 
         end select

100 FORMAT (2(e14.8,2x))

      end subroutine print_out

      subroutine findOH
   
         use rdf_lib
      
         implicit none


         integer              :: ncount=0, ierror,          &
                                 mic=1                         !From the RDF program that uses multiple supercells
   
         integer              :: na, ns, i, j,              &  !atom/species/general counters
                                 Ostar,                     &  !O* index
                                 Hprime,                    &  !H' index
                                 Hcheck,                    &  !Multiple OH check
                                 prevOstar,                 &  !Last O* index
                                 prevHprime,                &  !Last H' index
                                 H1, H2,                    &  !Hydrogen Counters
                                 total=0,                   &  !Total Atoms
                                 readstep                      !Current step in the read postion loop
                                 


   

         real(DP),allocatable     :: temp(:,:,:), stemp(:,:,:)     !temp coord and scaled coord arrays
                                                                   ! tempx(dim, atoms, species number)


        real(DP)              :: dx, dy, dz,          &  !Coordinate Displacments  
                                 rdist(3), sdist(3),  &  !real/scaled Coordinate Distance
                                 r2,                  &  !squared Distance
                                 aprim(3,3),          &  !primative matrix
                                 apinv(3,3),          &  !inverse of prim matirx
                                 omega                   !Cell volume


         do ns=1,nspt
            total = nat + nsp(ns)
         enddo


         allocate(temp(3,total,nspt), stemp(3,total,nspt) )


         Main_loop: do
            !
            !-----read *.pos file-----
            read(2,*, iostat=ierror) readstep, dummy
            !
            !-----Check for end of *.pos file-----
            if (ierror < 0) then
               write(*,*) ' ERROR: End of File Reached, correct max steps'
               stop
            endif
            !
            !-----Check for stop condition in *.pos file-------
            if(readstep >= maxstep) then
               write(*,*) 'Final step reached : ', maxstep
               write(*,*) 'Total number of steps: ', ncount
               exit
            endif
            !
            !Read in current position values
            do ns=1,nspt,1
               do na=1,nsp(ns),1
                  read(2,*) (temp(j,na,ns),j=1,3)
               enddo
            enddo
            !
            !-----read *cel file------
            read(3, *) dummy, dummy
            do i=1,3,1
               read(3, *) (aprim(i,j),j=1,3)
            enddo
            !
            !calculate the inverse
            call invert(aprim, apinv, omega)
            !
            !Convert to scaled positions 
            do ns=1,nspt,1
               do na=1,nsp(ns),1
                  call r_to_s( temp(1:3,na,ns), stemp(1:3,na,ns), apinv )
               enddo
            enddo
            !
            !---------Find OH--------------
            !  fort.45 is the OH tagging processing
            !  fort.36 tagged OH list
            !  fort.99 OHindex
            write(45,*) 'Configuration :' , (ncount+1)
            Hcheck = 0
            Ostar = 0
            Hprime = 0
            Oloop: do i=1,nsp(1),1
               !Hydrogen counters
               H1 = 0
               H2 = 0
               !
               Hloop: do j=1,nsp(2),1
                  !
                  dx = stemp(1,i,1) - stemp(1,j,2)
                  dy = stemp(2,i,1) - stemp(2,j,2)
                  dz = stemp(3,i,1) - stemp(3,j,2)
                  !
                  !Periodic adjusted distance
                  sdist(1) = dx - NINT( dx/(real(mic)) ) * real(mic)
                  sdist(2) = dy - NINT( dy/(real(mic)) ) * real(mic)
                  sdist(3) = dz - NINT( dz/(real(mic)) ) * real(mic)
                  !
                  !convert from scaled to real
                  call s_to_r( sdist, rdist, aprim )
                  r2 = rdist(1)**2 + rdist(2)**2 + rdist(3)**2
                  !
                  !Check for a O-H
                  if (r2 < (OHbond/0.52917_DP)**2 ) then
                     !How many Hydrogen bonds are already assocated with the O
                     if ( (H1==0) .and. (H2==0) ) then
                        H1 = j
                        write(45,*) 'Hydrogen1 ', i, j, ' ', (sqrt(r2)*0.52917_DP)
                     else
                        H2 = j
                        write(45,*) 'Hydrogen2 ', i, j, ' ', (sqrt(r2)*0.52917_DP)
                     endif
                  endif
                  !
               enddo Hloop
               !
               !Check for a singlar O-H bond
               if (H2 == 0 ) then
                  !
                  write (45, *) 'OH found!!!'
                  Ostar = i
                  OHindex(ncount+1) = i
                  Hprime = H1
                  !
                  !Multiple OH^- check
                  Hcheck = Hcheck + 1
                  !
               endif
               !
            enddo Oloop
            !
            !Check to make sure that there is at least one 
            if ( Hcheck == 0 ) then
               !
               write(36,*) (ncount+1), 'OH not found'
               noOH = noOH + 1
               !
               !Check to be sure that there is a previous OH
               if (prevOstar /= 0 .and. prevHprime /= 0) then
                  Ostar = prevOstar
                  OHindex(ncount+1)= prevOstar
                  Hprime = prevHprime
               else
                  print *, ' ERROR: No previous O* and H'' for configuration ', (ncount+1)
                  stop
               endif
               !
               !Only one OH 
            elseif (Hcheck == 1) then
               !reset the previous OH to this configuration
               prevOstar = Ostar 
               prevHprime = H1
               write (36, *) (ncount+1), 'Ostar : ', Ostar, 'Hprime : ', Hprime
               !
               !Two OH
            elseif ( Hcheck == 2 )  then
               !If there are two OH in the configuration, use the previous 
               !print *, 'WARRINNG: There are two OH- in configuration',
               !(ncount+1)
               write (36, *) (ncount+1), 'Ostar : ', Ostar, 'Hprime : ', Hprime, ' (PT in progress) '
               Ostar = prevOstar
               OHindex(ncount+1) = prevOstar
               Hprime = prevHprime
               dblOH = dblOH + 1
               !
               !More then two OH, NOT GOOD!
            elseif ( Hcheck > 2 ) then
               print *, ' NOTE: There are THREE or more OH- in configuration', (ncount+1), '!'
               Ostar = prevOstar
               OHindex(ncount+1) = prevOstar
               Hprime = prevHprime
               multiOH = multiOH + 1

            endif

            ncount = ncount + 1

            write(99,*) readstep, OHindex(ncount)

         enddo Main_loop


      end subroutine findOH

 End Program Diff
