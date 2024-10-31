      IMPLICIT NONE
      INTEGER OUT8,OUT9,OUT10,TOUT1,TOUT2
      INTEGER I1,I2,I3,I4,ID3,INP,Q
      INTEGER N1,N2,N3,N4,K
      INTEGER NSTATES,NS,NG,KT,ING1,ING2,ING3,ING4
      INTEGER NSTEP,DNG,ORDER,NGR,NGH
      REAL*8 H,TSTART,TEND,EV
      REAL*8 M1,M2,MU,V0,KR1,KP1,KR2,KP2
      REAL*8 DX,ND,MD,DI,R0,R,ND1,ND2,ND3,ND4
      REAL*8 PI,TFS,WAVENTOAU,EVTOAU
      REAL*8 INFO(11)
      COMMON /INFOS/ INFO
      REAL*8 W(4)
      COMMON /FREQS/ W

      REAL*8, ALLOCATABLE :: EN(:,:),EIG(:)
      CHARACTER (LEN=30) KEY01
      CHARACTER (LEN=30) KEY02
      CHARACTER (LEN=30) KEY03
      INTEGER ALLSTAT

   
    
      WRITE(KEY01,'(A)') "Ueigenvectors_ethLM"
      WRITE(KEY02,'(A)') "Ueigenvalues_ethLM"

      OPEN(UNIT=TOUT1,FILE=KEY01)
      OPEN(UNIT=TOUT2,FILE=KEY02)
      
      
 
!**************************************************************************************************
!READING INPUT DATA
!!!!!!!!!!!!ATTENTION!!!!!!!THERE ARE FEW FIXED PARAMETERS BELOW, WHICH CAN BE CHANGED!!!!!!!!!!!!!
      PI=DACOS(-1.0D0)
      TFS=2.4188843265020D-2
      H=0.50D0/TFS
      TSTART=0.0D0/TFS
      TEND=25.0D0/TFS

      NS=2
      ING1=22
      ING2=12
      ING3=7
      ING4=42
      NGH=10
      NG=ING1*ING2*ING4
      NGR=512
      ORDER=2
      DNG=NG*NS
     
      NSTEP=5000000

      EV=27.21140D0
!TUNING MODES
!6A
      W(2)=0.36/EV
!
      W(1)=0.20540D0/EV
!9A
      W(3)=0.15250D0/EV
!COUPLING MODE  
!10A
      W(4)=0.1100D0/EV
!GAP
      V0=1.90D0/EV

!parameters of the hamiltonian
!k
      INFO(1)=H
!interelectronic coupling linear
      INFO(2)=0.40160D0/EV
      INFO(3)=0.0D0
      INFO(4)=0.0D0
      INFO(5)=0.0D0

!some other useful constants
      INFO(6)=DBLE(ING1)
      INFO(7)=DBLE(ING2)
      INFO(8)=DBLE(ING3)
      INFO(9)=DBLE(ING4)
      INFO(10)=V0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE(EN(DNG,DNG),EIG(DNG),STAT=ALLSTAT)
      IF (ALLSTAT /= 0) STOP "*** Not enough memory ***"
      
      DO I1=1,DNG
        EIG(I1)=0.0D0
        DO I2=1,DNG 
          EN(I2,I1)=0.0D0
        END DO
      END DO
203   FORMAT(X,F30.16,X,F30.16)     
201   FORMAT(X,F22.16)
      DO I1=1,DNG
        READ(TOUT2,201) EIG(I1)
        DO I2=1,DNG
          READ(TOUT1,201,ADVANCE='NO') EN(I1,I2)
        END DO
        READ(TOUT1,*)
      END DO  

!END OF READING INPUT DATA
!*************************************************************************************************!
!                                                                                                 !  
!                                                                                                 !
!                                                                                                 !
!                                                                                                 !
!****************************ELECTRONIC HAMILTONIAN DYN IN DIABATIC BASIS SET*******************************
      CALL INITIAL(K,NG,NS,DNG,EN,EIG,KEY10)
      
      DEALLOCATE(EN,EIG)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END
!***************************************END*MAIN*************************************************!
!
!
!
!****************************************FUNCTIONS*************************************************
!*********************************IMPORTANT NUMERICAL SUBROUTINES*********************************    
!**************************************************************************************************

      SUBROUTINE initial(knum,ng,ns,dng,en,eigs,keyn)

      common /infos/ info
      real*8 info(11)
      common /freqs/ w
      real*8 w(4)
      common /shifts/ shift
      real*8 shift(4,2)

      INTEGER knum,k,l,f,q,n,m,i,j,i1,i2,ng,ns,dng,kickx
      INTEGER out,outg,outgg,i3,i4,nz,nd
      INTEGER factorial,ing1,ing2,ing3,ing4
      INTEGER nstep,order,orderz,ding
      INTEGER filen,imax
      INTEGER out1,out2,out3,out4,out5
      INTEGER out6,out7,out8,out9,out10
      INTEGER tout1,tout2,tout3
      INTEGER out04,out05,out06,out07,out08,out09
      INTEGER allstat,switch,ngh
 
      real*8 t,tstart,tend,h,kcur,TFS,norm,r0
      real*8 k1,k2,k3,k4,k5,kr,kp,c0,alpha,ne,me
      real*8 dx,rmin,ngrd,r,pi
      real*8 ord,nexp,nd1,nd2,nd3,nd4
      real*8 alpha1,alpha2,alpha3,alpha4

      real*8 en(dng,dng),v0
      real*8 eigs(dng)
      
      
      real*8 lam(14,4)
      real*8 c1,c2,c3,c4
      real*8 cn1,cn2,cn3,cn4
      real*8 cn
      real*8 den(dng),trans

      character (len=30) keyn
      character (len=30) key1
      character (len=30) key2
      character (len=30) key3
      character (len=30) key4
      character (len=30) key5
      character (len=30) key6
      character (len=30) key7
      character (len=30) key8
      character (len=30) key9
      character (len=30) key10
      character (len=30) key04
      character (len=30) key05
      character (len=30) key06
      character (len=30) key07
      character (len=30) key08
      character (len=30) key09
      character (len=30) tkey1
      character (len=30) tkey2
      character (len=30) tkey3


      TFS=2.4188843265020D-2

      v0=info(10)
      ing1=floor(info(6))
      ing2=floor(info(7))
      ing3=floor(info(8))
      ing4=floor(info(9))

      pi=dacos(-1.0d0)
      nexp=dexp(1.0d0)
!initial state - for the surprisal of the
!shifted gs wave function v=0 
      do 165 i=1,14
       do 164 j=1,4
         lam(i,j)=0.0d0
164    continue
165   continue

786   format(X,F20.8,4(X,F22.16))
788   format(X,F22.16)
789   format(X,F30.16)
34    format(X,I14)

!initial state
!0.5beta*omega
      alpha1=0.50d0*8000.0d0*w(1)
      alpha2=0.50d0*8000.0d0*w(2)
      alpha4=0.50d0*8000.0d0*w(4)
      k1=alpha1+alpha2+alpha4
!0.5beta*V_0
      c0=-0.50d0*8000.0d0*v0
      lam(1,1)=-c0+k1
      lam(1,2)=c0+k1
!      norm=dcosh(c0)/dsinh(alpha)
!      norm=dlog(norm)+alpha
!a+a
      c0=alpha1*2.0d0 
      lam(2,1)=c0
      lam(2,2)=c0
      c0=alpha2*2.0d0 
      lam(3,1)=c0
      lam(3,2)=c0
!      c0=alpha3*2.0d0 
!      lam(4,1)=dcmplx(c0,0.0d0)
!      lam(4,2)=dcmplx(c0,0.0d0)
      c0=alpha4*2.0d0 
      lam(5,1)=c0
      lam(5,2)=c0
      do 171 i=1,dng
          den(i)=0.0d0
171   continue                         
!coef for electronic                    
!lan rho
!      do 1071 i3=1,ing3
       do 1072 i2=1,ing2
        do 1073 i1=1,ing1
         do 1074 i4=1,ing4
!electronic state 1
         n1=ns*ing4*(i1-1)+ns*(i4-1)+1
         n1=n1+ns*ing1*ing4*(i2-1)
!same vib but electronic state 2
         n2=n1+1
         nd1=dble(i1)-1.0d0
         nd2=dble(i2)-1.0d0
!         nd3=dble(i3)-1.0d0
         nd4=dble(i4)-1.0d0
!diagonal in eig index 
!electronic 
         den(n1)=lam(1,1)+norm
         den(n2)=lam(1,2)+norm

!adag a
         den(n1)=den(n1)+lam(2,1)*nd1
         den(n1)=den(n1)+lam(3,1)*nd2
         den(n1)=den(n1)+lam(5,1)*nd4

         den(n2)=den(n2)+lam(2,2)*nd1
         den(n2)=den(n2)+lam(3,2)*nd2
         den(n2)=den(n2)+lam(5,2)*nd4
1074     continue                           
1073    continue                           
1072   continue
!1071  continue
!checking the norm
      norm=0.0d0
      do 1173 k=1,dng
          norm=norm+dexp(-den(k))
1173  continue                          
      write(*,*) norm
      k2=dlog(norm)
      do 1172 k=1,dng
          den(k)=den(k)+k2
1172  continue                          
      norm=0.0d0
      do 1273 k=1,dng
          norm=norm+dexp(-den(k))
1273  continue                          
      write(*,*) norm
      out1=808
      out2=809
      out3=810
      out4=812
      write(key1,'(a,i4.4,a)') "IA_trans"
      write(key2,'(a,i4.4,a)') "trans_csr"
      write(key3,'(a,i4.4,a)') "JA_trans"
      write(key4,'(a,i4.4,a)') "DIMS_trans"
      open(unit=out1,file=key1)
      open(unit=out2,file=key2)
      open(unit=out3,file=key3)
      open(unit=out4,file=key4)

!transformation to the eigenstates of the evolution operator
!eig - eigenvalues, en - eigenvectors
      k=1
      nd=1000
      nz=0
      do 181 knum=1,dng
       j=0
       write(*,*) knum
       do 184 q=0,nd,1
          m=knum+q
          if (m .ge. 1 .and. m .le. dng) then
            trans=0.0d0
            do 185 l=1,dng
              k1=en(l,knum)*en(l,m)
              trans=trans+den(l)*k1
185         continue
            if (dabs(trans) .gt. 1.0d-12) then 
               if (j .eq. 0) then
                write(out1,34,advance='no') k
               end if

               write(out2,789,advance='no') trans
               nz=nz+1

               write(out3,34,advance='no') m
               k=k+1
               j=1
            end if
          end if
184    continue
181   continue
      write(out4,34) nz

!DYNAMICS
      keyn='true'
      return
      END













 
