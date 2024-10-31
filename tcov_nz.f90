!to compile
!gfortran -fopenmp $finufft/include/finufft_mod.f90 tcov_nz.f90 -o tcov $finufft/lib/libfinufft.so -lfftw3 -lfftw3_omp -lgomp -lstdc++ /usr/lib/liblapack.so.3
      program tcov

      implicit none
      integer*8 tout1,tout2,tout3
      integer*8 out1,out2,out3,out4
      integer*8 i1,i2,i3,i4
      integer*8 ns,ng,ing1,ing2,ing3,ing4
      integer*8 dng,dsq,timep,nc,dngd,ncg
      integer*8 ncol,nrow1,nrow2,nz    
      integer, allocatable :: ia(:),ja(:)

      real*8 k1,h,tfs,c1
      real*8, allocatable :: eig(:),emn(:)


      complex*16, allocatable :: trans(:)


      character (len=30) key01
      character (len=30) tkey1
      character (len=30) tkey2
      character (len=30) tkey3
      character (len=30) key1
      character (len=30) key2
      character (len=30) key3
      character (len=30) key4

      integer allstat


      tout2=88
      write(tkey2,'(a)') "Ueigenvalues_ethLM"
      open(unit=tout2,file=tkey2)
      
 
!**************************************************************************************************
!reading input data
!!!!!!!!!!!!attention!!!!!!!there are few fixed parameters below, which can be changed!!!!!!!!!!!!!

      dng=22176
      dngd=dng+1
      dsq=dng*dng
      nc=50
      ncg=12
      tfs=0.02418880d0
      h=1.0d0/tfs
      timep=1001

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      allocate(eig(dng),stat=allstat)
      if (allstat /= 0) stop "*** not enough memory ***"
     
      do i1=1,dng
        eig(i1)=0.0d0
      end do
      do i1=1,dng
        read(tout2,201) k1
        eig(i1)=k1*h
      end do  
788   format(x,f30.16)
201   format(x,f22.16)
34    format(X,I14)
      
      out1=808
      out2=809
      out3=810
      out4=811
      write(key1,'(a)') "IA_trans"
      write(key2,'(a)') "trans_csr"
      write(key3,'(a)') "JA_trans"
      write(key4,'(a)') "DIMS_trans"
      open(unit=out1,file=key1)
      open(unit=out2,file=key2)
      open(unit=out3,file=key3)
      open(unit=out4,file=key4)

      read(out4,34) nz
      dngd=dng+1
      allocate(trans(nz),emn(nz),ja(nz),ia(dngd),stat=allstat)
      if (allstat /= 0) stop "*** not enough memory ***"

      do i1=1,nz
        trans(i1)=(0.0d0,0.0d0)
        emn(i1)=0.0d0
        ja(i1)=0
      end do

      do i1=1,dngd
        ia(i1)=0
      end do
      do i1=1,dng
        read(out1,34,advance='no') ia(i1)
      end do
      ia(dngd)=nz+1
      do i1=1,nz
        read(out2,788,advance='no') k1
        k1=k1**2
        trans(i1)=dcmplx(k1,0.0d0)
        read(out3,34,advance='no') ja(i1)
      end do


      do i1=1,dng
        nrow1=ia(i1)
        nrow2=ia(i1+1)-1
        do i2=nrow1,nrow2
           ncol=ja(i2)
           emn(i2)=eig(i1)-eig(ncol)
        end do
      end do
!end of reading input data
      call nuft(dngd,nc,ncg,timep,nz,emn,trans,ia,ja)
      
      end
!***************************************END*MAIN*************************************************!
!
!
!
!*************************************************************************************************
!*********************************IMPORTANT NUMERICAL SUBROUTINES********************************* 
!**************************************************************************************************
      subroutine nuft(dngd,nc,ncg,timep,dsq,emn,trans,ia,ja)

!     finufft fortran-header, always needed
      use finufft_mod

      integer*8 i,j,q,nc,nk,ttk
      integer*8 timep,dsq,dngd,ncg
      integer   iflag,ier,out1 
      integer*8 t1,t2,crate
      integer*8 gout1,gout2,gout3,gout4      
      integer*8 k,nz,nd
      integer ia(dngd),ja(dsq)

      real*8 k1,k2,k3,x
      real*8 pi,eps,t,h,tfs

      real*8 emn(dsq),tk(timep),lg(nc,timep)
      real*8 tm(timep,timep),eigs(timep)
      real*8 coef_cos,coef_sin,normc,norms
      real*8 cfcos(timep),cfsin(timep),omega
      
      complex*16 tcov(timep),trans(dsq)
      character (len=30) key1
      character (len=30) key2
      character (len=30) key3
      character (len=30) gkey1
      character (len=30) gkey2
      character (len=30) gkey3
      character (len=30) gkey4

      ! this is how you create the options struct in fortran...
      type(finufft_opts) opts
      ! or this is if you want default opts, make a null pointer...
      type(finufft_opts), pointer :: defopts => null()


786   format(X,F20.8,4(X,F22.16))
788   format(X,F22.16)
789   format(X,F30.16,X,F30.16)
780   format(X,F30.16)
784   format(2(X,F12.8))
781   format(X,F22.16)
761   format(X,F26.8)
722   format(i8,i8,X,F22.16,F22.16)
34    format(X,i8)

!using FFT on the non-uniform grid of frequencies
!c_j=|I_mn(0)|^2 - trans(dsq)
!x_j=(E_m-E_n)*dt- emn(dsq), important: emn should lie within [-pi,pi]                         
!positive sign of the exponent
      pi=dacos(-1.0d0)
      call system_clock(t1)
      iflag=1
      eps=1.0d-15
      do 101 i=1,timep
        tk(i)=dble(i)-1
        tcov(i)=(0.0d0,0.0d0)
101   continue 

      call finufft1d3(dsq,emn,trans,iflag,eps,timep,tk,tcov,defopts,ier)

      call system_clock(t2,crate)
      t = (dble(t2)-dble(t1))/dble(crate)
      if (ier.eq.0) then
         print '("done in ",f10.3," sec, ",e10.2" NU pts/s")',t,dble(dsq)/t
      else
         print *,'failed! ier=',ier
      endif

      do 102 i=1,timep
       write(*,*) tk(i), tcov(i)
       do 103 j=1,timep
         q=abs(i-j)+1         
         tm(i,j)=dble(tcov(q))      
103    continue
102   continue
      out1=245
      write(key1,'(a,i4.4,a)') "svd3d_info_1001fs_nufft"
      open(unit=out1,file=key1)

      call diag(timep,tm,eigs)      
!omega
      do 848 i=1,timep
        k1=dsqrt(dabs(eigs(i)))
        write(out1,761) k1
848   continue
      do 295 i=1,timep
         do 296 j=1,timep
          write(out1,780,advance='no') dble(tm(i,j))
296     continue
        write(out1,'(" ")')
295   continue 
!writing lambdas
      do 888 i=1,nc
        q=timep-(i-1)
        k1=dsqrt(dabs(eigs(q)))
        do 887 j=1,timep
          lg(i,j)=dble(tm(j,q))*k1
887     continue
888   continue
      do 828 i=1,timep
        do 827 j=1,nc
          write(out1,780,advance='no') lg(j,i)
827     continue
        write(out1,'(" ")')
828   continue
!printing Fourier coefs for each lambda
      do 345 i=1,nc
        q=timep-(i-1)
        do 344 j=1,timep
          coef_cos=0.0d0
          coef_sin=0.0d0
          do 343 ttk=1,timep
            x=pi*(dble(ttk)-1.0d0)*(dble(j)-1.0d0)/(dble(timep)-1.0d0)
            coef_cos=coef_cos+dcos(x)*tm(ttk,q)
            coef_sin=coef_sin+dsin(x)*tm(ttk,q)
343       continue
          if (j .eq. 1) then
            coef_cos=coef_cos/2.0d0
          end if
          coef_cos=coef_cos/(dble(timep)-1.0d0)
          coef_sin=coef_sin/(dble(timep)-1.0d0)
          write(out1,781,advance='no') coef_cos          
          write(out1,781,advance='no') coef_sin   
          cfcos(j)=coef_cos
          cfsin(j)=coef_sin       
344     continue
        write(out1,'(" ")')
345   continue
!done with Fourier coefs

      return
      END

!----------------------------------------------------------------------------------------------------
!      
      SUBROUTINE diag(ns,a,d)
      INTEGER*8 ns,info,lwork          
      INTEGER*8 liwork
      real*8 a(ns,ns), d(ns)
      real*8, allocatable :: z(:)
      integer, allocatable :: iwork(:)
      INTRINSIC INT, MIN
      EXTERNAL DSYEVD

!      lwork=-1
!      liwork=-1

!      call DSYEVD('Vectors','Upper',ns,a,ns,d,z,lwork,iwork,liwork,info)

!      lwork=MIN(5000,INT(z(1)))
!      liwork=MIN(5000,iwork(1))
  
      lwork=1+6*ns+2*(ns**2)+200
      liwork=3+5*ns+200
      allocate(z(lwork),iwork(liwork))
      call DSYEVD('Vectors','Upper',ns,a,ns,d,z,lwork,iwork,liwork,info)
      
      return
      end

!----------------------------------------------------------------------------------------------------
