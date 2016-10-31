      implicit none
c-----------------------------------------------------------------------
      integer m,Nx,Nz,Lm,Lq,ix,iz,i,iwF,ii,bm,done,jj,kk
      double precision Wx,Wz,pi,wMin,wMax,dw,scr
      double complex delta,lambda,pol,pol0

      double precision x,z,miu,Temp,bt,Zt,tol,fact,accuracy,wF
      double precision accuracyc
c-----------------------------------------------------------------------
c  Things to be allocated later from what we read from myPara:
      double precision, allocatable::  w(:,:)
      double precision, allocatable:: vect(:,:,:,:)
      double precision, allocatable:: pot(:,:)

      double complex, allocatable:: dn(:,:)
      double complex, allocatable:: dnp(:,:)
      double complex, allocatable:: dn0(:,:)
      double complex, allocatable:: vext(:,:)
c-----------------------------------------------------------------------

      pi=4d0*atan(1d0)
      accuracy=0.001d0
      accuracyc = 0.01d0*accuracy
      delta = (0.0,0.01d0)
      bm = 0
      scr = 30d0
      done=0
c-----------------------------------------------------------------------
      open(20,file='myPara.txt') 
      open(21,file='pot.txt')
      open(22,file='e.txt')
      open(23,file='psi.txt')
      open(24,file='testrpaz.txt')
      open(25,file='testlindz.txt')
      open(26,file='dnFinal5-2.txt')
      open(27,file='dnFinalreal.txt')
      open(28,file='dnFinalimag.txt')
      open(29,file='dnFinalabs.txt')
c-----------------------------------------------------------------------
c Determines the number of frequency points:
      iwF=1
      wMax=5.2d0/27.211d0
      wMin=5.2d0/27.211d0
      dw=(wMax-wMin)/(1d0*iwF)
c     ------------------------------------------
c     readin
c     -----------------------------------------
      read(20,*) m
      read(20,*) Nx
      read(20,*) Nz
      read(20,*) Wx
      read(20,*) Wz
      read(20,*) Lq
      read(20,*) miu
      read(20,*) Temp
      read(20,*) bt
      read(20,*) Zt
      read(20,*) tol

      allocate( w(1:m,0:Lq),vect(-Nz-2:Nz+2,-1:Nx+2,1:m,0:Lq),
     &   pot(-Nz-2:Nz+2,-1:Nx+2),dn(-Nz-2:Nz+2,-1:Nx+2),
     &   dnp(-Nz-2:Nz+2,-1:Nx+2),dn0(-Nz-2:Nz+2,-1:Nx+2),
     &   vext(-Nz-2:Nz+2,-1:Nx+2) )

      do Lm = 0,Lq,1
         do i = 1,m,1
            read(22,*) w(i,Lm)
            do ix=-1,Nx+2,1
               do iz = -Nz-2,Nz+2,1
                 read(23,*) vect(iz,ix,i,Lm)
               enddo
            enddo
         enddo
      enddo
      
c      vect = vect/dsqrt(Wx*Wz)

      pot=0d0
      do ix = 1,Nx
         do iz = -Nz,Nz
            read(21,*) pot(iz,ix)
            vext(iz,ix) = -iz*Wz 
         enddo
      enddo

      write(*,*) 'load complete'
c     -----------------------------------------
c     Here we start the frequency loop:
      do ii=1,iwF
         wF=wMin+ii*dw
c        ----------------------
c        computes the un-screen dn
c        -----------------------
         call chi0(m,Lq,Zt,Nx,Nz,Wx,Wz,
     &    miu,bt,tol,accuracyc,delta,
     &    w,vect,pot,
     &    wF,vext,dn0)
c         ------------------------
c         solves the rpa equation
c         -------------------------
          dn = dn0
          dnp= dn0
          call solverpa(bm,m,Lq,Zt,Nx,Nz,Wx,Wz,
     &    miu,bt,tol,accuracy,delta,
     &    w,vect,pot,
     &    wF,dn0,dn,dnp)
c         -------------------------
c         Compute the polarization:
c         -------------------------
          pol = 0d0
          pol0 = 0d0
          do ix=1,Nx,1
             x=ix*wx
             do iz=-Nz,Nz,1
                z = iz*Wz
                pol = pol +2d0*pi* z*dn(iz,ix)*Wx*Wz
                pol0 = pol0 + 2d0*pi*z*dn0(iz,ix)*Wx*Wz
             enddo
         enddo
         write(24,*) wF*27.211d0,2d0*wF/(137d0)*aimag(pol)
         write(25,*) wF*27.211d0,2d0*wF/(137d0)*aimag(pol0)
         write(*,*) wF*27.211d0,2d0*wF/(137d0)*aimag(pol)
      enddo
      
      do ix=1,Nx,1
          x=ix*wx
          do iz=-Nz,Nz,1
              z=iz*wz
              write(26,*) dn(iz,ix)
              write(27,*) real(dn(iz,ix))
              write(28,*) aimag(dn(iz,ix))
              write(29,*) abs(dn(iz,ix))
          end do
      end do
      
      
c  We're done writing RPA results.     
c------------------------------------------
      close(20)
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)
      
c ------------------------------------------------------
      stop
      end
