module global_data_g

  implicit none
  INTEGER*4 i,j, idum
  INTEGER*4 iunit, ios
  INTEGER*4 mrfac, nx4, ny4, nlim, nbdry, nyefit, nxefit 
  INTEGER*4 kxord, kyord, nworkf


  DOUBLE PRECISION xdim, zdim, rcentr, rgrid1, zmid
  DOUBLE PRECISION rmagx, zmagx,simagx,sibdry,bcentr
  DOUBLE PRECISION cpasma, xdum

  DOUBLE PRECISION, dimension(:), allocatable :: fpol, pres, &
       workk1, workk2, qpsi, rbdry, zbdry, xlim, ylim
  DOUBLE PRECISION, dimension(:,:), allocatable :: fold

  character*80 runid
  character*8 label(6)

end module global_data_g



subroutine readg 


  use global_data_g
  INTEGER*4 stat


  !!c *************** read the output of the EFIT code ******************

  print *, 'Reading neqdsk file .......'

  iunit=1 !file is neqdsk
  open (iunit, file='neqdsk', form='formatted', iostat=ios, status='old')

  if (ios .ne. 0) then
     print *, "**** neqdsk file not found, exiting..."
     STOP
  else

     read(iunit,2000) (label(i),i=1,6),idum,nxefit,nyefit
     print *, "   nxefit=", nxefit, ", nyefit=", nyefit

     runid=label(1)//label(2)//label(3)//label(4)//label(5)//label(6)
     print *, "   runid=", runid


     read(iunit,2020) xdim,zdim,rcentr,rgrid1,zmid
     read(iunit,2020) rmagx,zmagx,simagx,sibdry,bcentr
     read(iunit,2020) cpasma,simagx,xdum,rmagx,xdum
     read(iunit,2020) zmagx,xdum,sibdry,xdum,xdum
  
  
     !!-read arrays

     !!-do it in case of a repetitive call to this function
     deallocate (fpol, stat=stat)
     deallocate (pres, stat=stat)
     deallocate (workk1, stat=stat)
     deallocate (workk2, stat=stat)
     deallocate (fold, stat=stat)
     deallocate (qpsi, stat=stat)
     deallocate (rbdry, stat=stat)
     deallocate (zbdry, stat=stat)
     deallocate (xlim, stat=stat)
     deallocate (ylim, stat=stat)


     allocate (fpol(1:nxefit))
     read(iunit,2020) (fpol(i),i=1,nxefit)
     
     allocate (pres(1:nxefit))
     read(iunit,2020) (pres(i),i=1,nxefit)
     
     allocate (workk1(1:nxefit))
     read(iunit,2020) (workk1(i),i=1,nxefit)

     allocate (workk2(1:nxefit))     
     read(iunit,2020) (workk2(i),i=1,nxefit)
     
     allocate (fold(1:nxefit,1:nyefit))
     read(iunit,2020) ((fold(i,j),i=1,nxefit),j=1,nyefit)
     
     allocate (qpsi(1:nxefit))
     read(iunit,2020) (qpsi(i),i=1,nxefit)
     
     read(iunit,2022) nbdry,nlim
     print *, 'nbdry=', nbdry, ', nlim=', nlim
     allocate (rbdry(1:nbdry))
     allocate (zbdry(1:nbdry))
     read(iunit,2020) (rbdry(i),zbdry(i),i=1,nbdry)

     allocate (xlim(1:nlim))
     allocate (ylim(1:nlim))
     read(iunit,2020) (xlim(i),ylim(i),i=1,nlim)

     close (iunit)

2000 format(6a8,3i4)
2020 format(5e16.9)
2022 format(2i5)
     
  endif
  
  print *, '.......finished'

end subroutine readg






subroutine get_nxy (nxefit_, nyefit_, nbdry_, nlim_)
!
! return dimensions so that memory can be allocated at main (i.e. IDL) level
!----------------------------------------------------------------------------!

  use global_data_g
  INTEGER*4, intent (OUT) :: nxefit_, nyefit_, nbdry_, nlim_

  !-can I just return-by-value???
  print *, "***getting g-dims"
  nxefit_=nxefit
  nyefit_=nyefit
  nbdry_=nbdry
  nlim_=nlim

  print *, 'In [f90:get_nxy] nxefit_, nyefit_, nbdry_, nlim_', nxefit_, nyefit_, nbdry_, nlim_

end subroutine get_nxy



subroutine get_psi (&
     nxefit_, nyefit_, nbdry_, nlim_,&
     fold_,& !!-2D
     fpol_, pres_, qpsi_, rbdry_, zbdry_, xlim_, ylim_,& !!-1D
     xdim_, zdim_, rcentr_, rgrid1_, zmid_, rmagx_, zmagx_, simagx_, sibdry_, bcentr_)
!
! return a copy of array fold, the memory is allocated by the caller
!----------------------------------------------------------------------------!

  use global_data_g

  INTEGER*4, intent (IN) :: nxefit_, nyefit_, nbdry_, nlim_
  DOUBLE PRECISION, dimension(nxefit_,nyefit_), intent (OUT)  :: fold_
  DOUBLE PRECISION, dimension(nxefit_), intent (OUT)  :: fpol_, pres_, qpsi_
  DOUBLE PRECISION, dimension(nbdry_), intent (OUT)  :: rbdry_, zbdry_
  DOUBLE PRECISION, dimension(nlim_), intent (OUT)  :: xlim_, ylim_
  DOUBLE PRECISION, intent (OUT)  :: &
       xdim_, zdim_, rcentr_, rgrid1_, zmid_, rmagx_, zmagx_, simagx_, sibdry_, bcentr_

  print *, "***getting g-data"

  fold_=fold !!-copy arrays

  fpol_=fpol
  pres_=pres
  qpsi_=qpsi
  rbdry_=rbdry
  zbdry_=zbdry
  xlim_=xlim
  ylim_=ylim

  xdim_=xdim 
  zdim_=zdim
  rcentr_=rcentr
  rgrid1_=rgrid1
  zmid_=zmid
  rmagx_=rmagx 
  zmagx_=zmagx 
  simagx_=simagx
  sibdry_=sibdry
  bcentr_=bcentr
  
end subroutine get_psi
