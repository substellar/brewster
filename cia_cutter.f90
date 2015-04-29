program cia_resample



  implicit none


  integer:: npts, ntemps, i, j, nnew,oldw1,oldw2
  real,allocatable :: temp(:)
  real,allocatable :: oldwavenum(:)
  real,dimension(319188)::newwavenum
  real :: junk, wdiff
  real,allocatable  :: cia(:,:,:)
  real,allocatable  :: newcia(:,:,:)


  
  open(unit=10,file="../LineLists/CH4/ch4_1_1")

  nnew = 319188
  do i = 1, 23
     read(10,*)
  end do
  
  do i = 1, nnew
     read(10,*) newwavenum(i), junk
  end do

  close(10)
  

  open(unit=10, file="final1_abel_CIA.dat")


  read(10,*) npts,ntemps

  allocate(temp(ntemps))
  allocate(cia(ntemps,npts,4))
  allocate(oldwavenum(npts))
  allocate(newcia(ntemps,nnew,4))

  do i = 1, ntemps
     read(10,*)  temp(i)

     do j = 1, npts
        read(10,*) oldwavenum(j),cia(i,j,1), cia(i,j,2), cia(i,j,3), cia(i,j,4)
        

     end do
  end do
  
  close(10)

  do i= 1 , nnew

     wdiff = abs(oldwavenum - newwavenum(i))

     oldw1 = minloc(wdiff,1)

     if (oldwavenum(oldw1) .lt. newwavenum(i)) then
        oldw2 = oldw1 + 1
     else
        oldw2 = oldw1
        oldw1 = oldw2 - 1
     end if


     intfact = (newwavenum(i) - oldwavenum(oldw1)) / (oldwavenum(oldw2) - oldwavenum(oldw1))

     do j = 1, ntemps
        do k = 1, 4
           newcia(j,i,k)  = ((cia(j,oldw2,k) - cia(j,oldw1,k))*intfact) + cia(j,oldw1,k)
        end do
     end do
  end do
  
  open(unit=10,file="high_res_CIA_table"
  write(10,*) nnew, ntemps
  do j = 1, ntemps
     write(10,*) temp(j)

     do i = 1, nnew
        write(10,*) newwavenum(i), newcia(j,i,1),newcia(j,i,2),newcia(j,i,3),newcia(j,i,4)
     end do
     
  end do

  close(10)
  deallocate(cia,newcia,temp,wavenum)


end program cia_resample
  

  
