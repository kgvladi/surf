!to compile
!gfortran -fopenmp $finufft/include/finufft_mod.f90 tcov_nz.f90 -o tcov $finufft/lib/libfinufft.so -lfftw3 -lfftw3_omp -lgomp -lstdc++ /usr/lib/liblapack.so.3
      program analytical_constraints

      implicit none
      integer tout1,tout2,tout3
      integer out1,out2,out3,out4,out5
      integer i,j,q,i1,i2,i3,i4
      integer nk,nc
      integer dng,timep,dngd
      integer ncol,nrow1,nrow2,nz    
      integer, allocatable :: ia(:),ja(:)

      real*8 k1,h,tfs,c1,omega
      real*8, allocatable :: eig(:),emn(:)
      real*8, allocatable :: cfcos(:),cfsin(:)
      real*8, allocatable :: trans(:)


      character (len=30) key01
      character (len=30) tkey1
      character (len=30) tkey2
      character (len=30) tkey3
      character (len=30) key1
      character (len=30) key2
      character (len=30) key3
      character (len=30) key4
      character (len=30) key5

      integer allstat


      tout2=88
      write(tkey2,'(a)') "Ueigenvalues_ethLM"
      open(unit=tout2,file=tkey2)

      call getarg(1,key1)
      read(key1,'(i2)') nk
      
 
!**************************************************************************************************
!reading input data
!!!!!!!!!!!!attention!!!!!!!there are few fixed parameters below, which can be changed!!!!!!!!!!!!!

      dng=22176
      nc=50
      dngd=dng+1
      tfs=0.02418880d0
      h=1.0d0/tfs
      timep=1001

      
      allocate(eig(dng),stat=allstat)
      if (allstat /= 0) stop "*** not enough memory ***"
     
      do i1=1,dng
        eig(i1)=0.0d0
      end do
      do i1=1,dng
        read(tout2,201) eig(i1)
      end do  
788   format(x,f30.16)
201   format(x,f22.16)
34    format(X,I14)
786   format(X,F20.8,4(X,F22.16))
789   format(X,F30.16,X,F30.16)
780   format(X,F30.16)
784   format(2(X,F12.8))
781   format(X,F22.16)
761   format(X,F26.8)
722   format(i8,i8,X,F22.16,F22.16)
      
      out1=808
      out2=809
      out3=810
      out4=811
      out5=245
      write(key1,'(a)') "IA_trans"
      write(key2,'(a)') "trans_csr"
      write(key3,'(a)') "JA_trans"
      write(key4,'(a)') "DIMS_trans"
      write(key5,'(a,i4.4,a)') "svd3d_info_1001fs_nufft"
      open(unit=out1,file=key1)
      open(unit=out2,file=key2)
      open(unit=out3,file=key3)
      open(unit=out4,file=key4)
      open(unit=out5,file=key5)

      read(out4,34) nz
      dngd=dng+1
      allocate(trans(nz),emn(nz),ja(nz),ia(dngd),&
      &cfcos(timep),cfsin(timep),&
      &stat=allstat)
      if (allstat /= 0) stop "*** not enough memory ***"

      do i1=1,nz
        trans(i1)=0.0d0
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
        read(out2,788,advance='no') trans(i1)
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

      q=timep-nk
      do i=1,timep
        if ( i .eq. q) then
         read(out5,761) omega
        else
         read(out5,*) 
        end if
      end do
      do i=1,timep*2
         read(out5,*) 
      end do
      do i=1,nc
        q=i-1
        if ( q .eq. nk) then
         do j=1,timep
          read(out5,781,advance='no') cfcos(j)          
          read(out5,781,advance='no') cfsin(j) 
         end do 
        else
          read(out5,*) 
        end if
      end do
!end of reading input data
      call constraints(nk,timep,nz,omega,trans,emn,&
      &cfcos,cfsin)
      
      end
!***************************************END*MAIN*************************************************!
!
!
!
!*************************************************************************************************
!*********************************IMPORTANT NUMERICAL SUBROUTINES********************************* 

!----------------------------------------------------------------------------------------------------
!      
      SUBROUTINE constraints(nk,timep,nz,omega,trans,emn,&
      &cfcos,cfsin)

      integer timep,ng,ns,timep0,k,nk
      integer i,j,n,m,i1,i2,p,p0,q
      integer out1,out2,out3,out4,nd,nz
      
      real*8 t,bigt,h,pi,tfs,omega
      real*8 emn(nz),cfcos(timep),cfsin(timep)
      real*8 trans(nz)
      real*8 k1,k2,k3,k4,thp,thm
      real*8 kn1,kn2,kn3,kn4

      complex*16 im,sumg,sum1,sum2
      complex*16 c1,c2,c3,c4,s3,s4,c8,c10,c11
      complex*16 cn1,cn2,cn3,cn4
      

      character (len=30) key1
      character (len=30) key2
      character (len=30) key3
      character (len=30) key4


      tfs=0.02418880d0
      h=1.0d0/tfs
      out2=923
      write(key2,'(a,i2.2,a)') "G",nk,"_csr"
      open(unit=out2,file=key2)

786   format(X,F20.8,4(X,F22.16))
788   format(X,F22.16)
789   format(X,F30.16,X,F30.16)
    
      timep0=timep-1
      tfs=2.4188843265020d-2
      h=1.0d0/tfs
      bigt=119.0d0/tfs
      pi=dacos(-1.d0)
      im=(0.0d0,1.0d0)
!number of non-zero elements
!      nz=dng+2*nd*dng-nd*(nd+1)
34    format(X,I14)
      do 12 i1=1,nz
         sumg=(0.0d0,0.0d0)
         if ( emn(i1) .eq. 0.0d0) then
          k1=1.0d0/omega
         else
          k1=1.0d0/(2.0d0*omega) 
         end if   
         c1=dcmplx(k1,0.0d0)        
         sum1=(0.0d0,0.0d0)
         sum2=(0.0d0,0.0d0)
         k1=emn(i1)*bigt/2.0d0
         c2=cdexp(-im*dcmplx(k1,0.0d0))
         do 13 p=1,timep
          p0=p-1
          k1=pi*dble(p0)/2.0d0
          c8=cdexp(-im*dcmplx(k1,0.0d0))
          k1=cfcos(p)
          c3=dcmplx(k1,0.0d0)*c8 
          c4=dcmplx(k1,0.0d0)*dconjg(c8) 
          k1=cfsin(p)
          s3=dcmplx(k1,0.0d0)*c8*im 
          s4=dcmplx(k1,0.0d0)*dconjg(c8)*im 
!emn*dt+ *pi/(N-1)
          thp=emn(i1)*h+dble(p0)*pi/dble(timep0)
!emn*dt- *pi/(N-1)
          thm=emn(i1)*h-dble(p0)*pi/dble(timep0)

          kn1=thp/2.0d0
          kn3=dabs(dsin(kn1))
          if (kn3 .gt. 1.0d-14) then
            kn2=dsin(dble(timep)*kn1)/dsin(kn1)
            sum1=sum1+c3*dcmplx(kn2,0.0d0)
            sum2=sum2+s3*dcmplx(kn2,0.0d0)
          end if  

          kn1=thm/2.0d0
          kn3=dabs(dsin(kn1))
          if (kn3 .gt. 1.0d-14) then
           kn2=dsin(dble(timep)*kn1)/dsin(kn1)
           sum1=sum1+c4*dcmplx(kn2,0.0d0)
           sum2=sum2-s4*dcmplx(kn2,0.0d0)
          end if
13       continue
         k1=trans(i1) 
         sumg=-(sum1+sum2)*c1*c2*dcmplx(k1,0.0d0)
         write(out2,789,advance='no') sumg
12    continue
      write(out2,'(" ")')

      return
      END
