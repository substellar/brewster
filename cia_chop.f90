program cia_chop



  implicit none


  integer:: nwave, ntemps, i, j
  real,allocatable :: temp(:)
  double precision,dimension(319188)::wavenum
  double precision,allocatable  :: cia(:,:,:)
  character(len=50) :: name
  


  open(unit=10, file="../LineLists/high_res_CIA_table.dat")


  read(10,*) nwave,ntemps

  allocate(temp(ntemps))
  allocate(cia(ntemps,nwave,4))

  
  do j = 1, ntemps
     read(10,*)  temp(j)

     do i = 1, nwave
        read(10,*) wavenum(i),cia(j,i,1), cia(j,i,2), cia(j,i,3), cia(j,i,4)
        

     end do
  end do
  
  close(10)

  open(10,file="cia_temps.dat")
  write(10,*) temp
  close(10)
  
  do j = 1, ntemps
     write(name,"(A,I0)") "../LineLists/CIA/cia_highres_tempindex_",j

     open(10,file=name,status="new")
     write(10,*) temp(j),nwave
     do i = 1, nwave
        write(10,*) wavenum(i), cia(j,i,1),cia(j,i,2),cia(j,i,3),cia(j,i,4)
     end do
     close(10)
     
  end do

  deallocate(cia,temp)


end program cia_chop
  

  
