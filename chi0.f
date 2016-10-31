      subroutine chi0(m,Lq,Zt,Nx,Nz,Wx,Wz,
     & miu,bt,tol,accuracy,delta,
     & w,vect,pot,
     & freq,vext,dn)

      implicit none
      integer m,Lq,Nx,Nz
      double precision Wx,Wz,Zt
      double precision miu,bt,tol,accuracy,freq
      double complex delta
      double precision w(1:m,0:Lq),vect(-Nz-2:Nz+2,-1:Nx+2,1:m,0:Lq)
      double precision pot(-Nz-2:Nz+2,-1:Nx+2)
      double complex vext(-Nz-2:Nz+2,-1:Nx+2),dn(-Nz-2:Nz+2,-1:Nx+2)

c-----------------------------------------------------------------------
      integer ix,iz,Lm,i
      double precision pi,fact,x,z
      double complex lambda

      double complex, allocatable:: F(:,:),start(:,:),
     &   T_P(:,:),T_M(:,:)
      allocate(  F(-Nz-2:Nz+2,-1:Nx+2),
     &                T_P(-Nz-2:Nz+2,-1:Nx+2),     
     &                T_M(-Nz-2:Nz+2,-1:Nx+2),     
     &                start(-Nz-2:Nz+2,-1:Nx+2) )


      pi=4d0*atan(1d0)
c-----------------------------------------------------------------------


      start = 0d0

c     --------------------------------
c     start building dn
c     -------------------------------
      dn = 0d0    
      do Lm=0,Lq
         do i=1,m
            if(Lm.eq.0) then
              fact = 2d0/(1d0+dexp((w(i,Lm)-miu)*bt)) 
            else
              fact = 4d0/(1d0+dexp((w(i,Lm)-miu)*bt))
            endif
c            write(*,*) fact/Zt
            if(fact/Zt.lt.tol) cycle 
            do ix=1,Nx
               do iz=-Nz,Nz
                   F(iz,ix)= vext(iz,ix)*vect(iz,ix,i,Lm)
               enddo
            enddo

c           ------------------------------------
c           solve (H-[Enm+wF+delta)T_P=F
c           -- ---------------------------------
            lambda = w(i,Lm) + freq + delta
      call green(Lm,Nx,Nz,Wx,Wz,T_P,F,start,accuracy,pot,lambda)
c           ------------------------------------
c           solve (H-[Enm-wF-delta)T_M=F
c           -----------------------------------
            lambda = w(i,Lm) - freq - delta
      call green(Lm,Nx,Nz,Wx,Wz,T_M,F,start,accuracy,pot,lambda)

c           ---------------------------------------------
c           compute dn
c           ----------------------------------------------
            do ix=1,Nx
               x=ix*wx
               do iz=-Nz,Nz
                    dn(iz,ix)= dn(iz,ix)
     &            -  1d0/(2d0*pi)*fact*vect(iz,ix,i,Lm)
     &            *  (T_P(iz,ix)+T_M(iz,ix))
                enddo
            enddo 
         enddo
      enddo

c      dn = 1.2234d0 * dn

      return

      end 
