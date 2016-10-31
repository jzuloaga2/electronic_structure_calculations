      subroutine green(Lm,Nx,Nz,Wx,Wz,var,b,in,accuracy,pot,lambda)
      implicit none
      integer Lm,Nr,Nx,Nz
      double precision Wx,Wz,accuracy
      double complex lambda
      double complex var(-Nz-2:Nz+2,-1:Nx+2)
      double complex in(-Nz-2:Nz+2,-1:Nx+2)
      double complex b(-Nz-2:Nz+2,-1:Nx+2)
      double precision pot(-Nz-2:Nz+2,-1:Nx+2)

c     **********************************************
      integer ix,iz,it,maxit
      double complex beta, z, rr
      double precision eps,norm,gx,gz
      double precision x,check,bc,bc1,bc2

      double complex, allocatable:: p(:,:),pp(:,:),
     &  q(:,:),qp(:,:),r(:,:),rp(:,:)
      allocate( p(-Nz-2:Nz+2,-1:Nx+2),pp(-Nz-2:Nz+2,-1:Nx+2),
     & q(-Nz-2:Nz+2,-1:Nx+2),qp(-Nz-2:Nz+2,-1:Nx+2),
     & r(-Nz-2:Nz+2,-1:Nx+2),rp(-Nz-2:Nz+2,-1:Nx+2) )

c     **********************************************
c     define some constants
c     **********************************************

      Nr = Nx*(2*Nz+1)
c      maxit = 5000
      maxit = 5000
      gx = 1d0/(24d0*Wx*Wx)
      gz = 1d0/(24d0*Wz*Wz)
      if(Lm.eq.0) then
        bc = 2d0**(Lm+0.5d0)+0.25d0-Lm*Lm
        bc1 = (3d0*Lm*Lm+30d0-0.75d0)*2d0**(Lm+0.5d0)+2d0**(2d0*Lm+1d0)
     & -16d0*3d0**(Lm+0.5d0)
        bc2 = 12d0*(Lm*Lm-0.25d0)+3d0**(Lm+0.5d0)-bc1*2d0**(Lm+0.5d0)
      else
        bc=2d0
        bc1=16d0
        bc2=-29d0
      end if
c     *********************************************************
c     calculates the norm of b
c     *********************************************************

      norm=0d0
      do ix = 1,Nx
         do iz =-Nz,Nz
            norm=norm+abs(b(iz,ix))
         enddo
      enddo
c      write(*,*) 'green norm',norm
      if(norm.lt.0.00000000000000001) then
        var = 0d0
        goto 111
      end if
c     ********************************************************
c     this takes care of the boundary conditions
c     ********************************************************
      in(:,0) = (16d0-bc1)*in(:,1)
      in(:,-1) = 16d0*in(:,0)-(30d0+bc2)*in(:,1)
     & +(16d0-bc1)*in(:,2)

      in(:,Nx+1) = 0d0
      in(:,Nx+2) = -in(:,Nx)

      in(-Nz-1,:) = 0d0
      in(-Nz-2,:) = -in(-Nz,:)

      in(Nz+1,:) = 0d0
      in(Nz+2,:) = -in(Nz,:)
c     *********************************************************
c     calculates the rest
c     *********************************************************
      r = 0d0
      do ix = 1,Nx
         x = ix*Wx
         do iz = -Nz,Nz
            r(iz,ix)= ( 30d0*(gx+gz)  
     &     +(Lm*Lm-0.25d0)/(2d0*x*x)+pot(iz,ix)-lambda)*in(iz,ix)
     &     - 16d0*gx*(in(iz,ix+1)+in(iz,ix-1))
     &     + gx*(in(iz,ix+2)+in(iz,ix-2))
     &     - 16d0*gz*(in(iz+1,ix)+in(iz-1,ix))
     &     + gz*(in(iz+2,ix)+in(iz-2,ix))
         enddo
      enddo
      r=b-r
      rp=r
c     ********************************************************
c     starts the iterative process
c     ********************************************************
      p=r
      pp=p
      var=in

      do it = 1,maxit
c        ******************
c        computes q=A*p
c        *************************************************
c        this takes care of the boundary conditions
c        *************************************************
         p(:,0) = (16d0-bc1)*p(:,1)
         p(:,-1) = 16d0*p(:,0)-(30d0+bc2)*p(:,1)
     &   +(16d0-bc1)*p(:,2)

         p(:,Nx+1) = 0d0
         p(:,Nx+2) = -p(:,Nx)

         p(-Nz-1,:) = 0d0
         p(-Nz-2,:) = -p(-Nz,:)
      
         p(Nz+1,:) = 0d0
         p(Nz+2,:) = -p(Nz,:)
c        *************************************************
c        calculates q=A*p
c        *************************************************
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
            q(iz,ix)= ( 30d0*(gx+gz)  
     &     +(Lm*Lm-0.25d0)/(2d0*x*x)+pot(iz,ix)-lambda )*p(iz,ix)
     &     - 16d0*gx *(p(iz,ix+1)+p(iz,ix-1))
     &     + gx*(p(iz,ix+2)+p(iz,ix-2))
     &     - 16d0*gz*(p(iz+1,ix)+p(iz-1,ix))
     &     + gz*(p(iz+2,ix)+p(iz-2,ix))
            enddo
         enddo
c        ******************
c        computes qp=A*pp
c        *************************************************
c        this takes care of the boundary conditions
c        *************************************************
         pp(:,0) = (16d0-bc1)*pp(:,1)
         pp(:,-1) = 16d0*pp(:,0)-(30d0+bc2)*pp(:,1)
     &   +(16d0-bc1)*pp(:,2)

         pp(:,Nx+1) = 0d0
         pp(:,Nx+2) = -pp(:,Nx)

         pp(-Nz-1,:) = 0d0
         pp(-Nz-2,:) = -pp(-Nz,:)

         pp(Nz+1,:) = 0d0
         pp(Nz+2,:) = -pp(Nz,:)
c        *************************************************
c        calculates qp=A^**pp
c        *************************************************
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
            qp(iz,ix)= ( 30d0*(gx+gz)
     &     +(Lm*Lm-0.25d0)/(2*x*x)+pot(iz,ix)-dconjg(lambda))*pp(iz,ix)
     &     - 16d0*gx*(pp(iz,ix+1)+pp(iz,ix-1))
     &     + gx*(pp(iz,ix+2)+pp(iz,ix-2))
     &     - 16d0*gz*(pp(iz+1,ix)+pp(iz-1,ix))
     &     + gz*(pp(iz+2,ix)+pp(iz-2,ix))
            enddo
         enddo
c        *************************************************
c        makes the updates
c        *************************************************
         rr=0d0
         z=0d0  
         do ix = 1,Nx
            do iz = -Nz,Nz
               rr=rr+dconjg(rp(iz,ix))*r(iz,ix)
               z=z+dconjg(pp(iz,ix))*q(iz,ix)
            enddo
         enddo
         z=rr/z
         var=var+z*p

         norm = 0d0
         check=0d0
         do ix = 1,Nx
            do iz = -Nz,Nz
               check=check+abs(z*p(iz,ix))
               norm = norm + abs(var(iz,ix))
            enddo
         enddo

         r=r-z*q
         rp=rp-dconjg(z)*qp

         beta=0d0
         do ix = 1,Nx
            do iz = -Nz,Nz
               beta=beta+ dconjg(rp(iz,ix))*r(iz,ix)
            enddo
         enddo
         beta=beta/rr
         p=r+beta*p
         pp=rp+dconjg(beta)*pp         

c        ************************************************
c        checks accuracy
c        ************************************************
         if (check/norm.lt.accuracy) then
            goto 111
         end if
c         write(*,*) it,check/norm    
      enddo

      write(*,*) 'accuracy not achieved in green.f'

111   return     

      end

