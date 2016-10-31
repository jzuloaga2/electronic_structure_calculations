      subroutine solverpa(bm,m,Lq,Zt,Nx,Nz,Wx,Wz,
     & miu,bt,tol,accuracy,delta,
     & w,vect,pot,
     & freq,dn0,dn,dnp)
      integer bm,m,Lq,Nx,Nz,Nr
      double precision Wx,Wz,Zt,scr,pi
      double precision miu,bt,tol,accuracy,freq
      double complex delta,deltap
      double precision w(1:m,0:Lq),vect(-Nz-2:Nz+2,-1:Nx+2,1:m,0:Lq)
      double precision pot(-Nz-2:Nz+2,-1:Nx+2)
      double complex dn0(-Nz-2:Nz+2,-1:Nx+2)
      double complex dn(-Nz-2:Nz+2,-1:Nx+2)
      double complex dnp(-Nz-2:Nz+2,-1:Nx+2)
      double complex cden(-Nz-2:Nz+2,-1:Nx+2)
      double complex vh(-Nz-2:Nz+2,-1:Nx+2)
      double complex start(-Nz-2:Nz+2,-1:Nx+2)
c     **********************************************
      double precision accuracyh,accuracyc,accuracyPedo
      integer ix,iz,it,maxit
      double complex beta, z, rr
      double precision eps,norm,gx,gz
      double precision x,check,bc,bc1,bc2

      double complex, allocatable:: p(:,:),pp(:,:),
     &  q(:,:),qp(:,:),r(:,:),rp(:,:)
      allocate( p(-Nz-2:Nz+2,-1:Nx+2),pp(-Nz-2:Nz+2,-1:Nx+2),
     & q(-Nz-2:Nz+2,-1:Nx+2),qp(-Nz-2:Nz+2,-1:Nx+2),
     & r(-Nz-2:Nz+2,-1:Nx+2),rp(-Nz-2:Nz+2,-1:Nx+2) )

      maxit = 100
      scr = 200d0
      Nr = Nx*(2*Nz+1)
      pi=4d0*atan(1d0)
      start = 0d0
      accuracyh = 0.001d0*accuracy
      accuracyc = 0.01d0*accuracy
      accuracyPedo = 1d-9
c      write(*,*) 'delta',delta
c     *********************************************************
c     calculates the norm of dn0
c     *********************************************************
      deltap = -delta
      norm=0D0
      do ix = 1,Nx
         do iz =-Nz,Nz
            norm=norm+abs(dn0(iz,ix))
         enddo
      enddo
      if(norm.lt.0.00000000000000001) then
        dn = 0d0
        goto 111
      end if
c     ********************************************************
c     calculates the rests, r and rp  assuming initials = dn, dnp
c     ********************************************************
c     ---------------------------------
c     compute the cden and then the hartree potential
c     ---------------------------------
      do ix = 1,Nx
         x = ix*Wx
         do iz = -Nz,Nz
               cden(iz,ix)=-4d0*pi*dn(iz,ix)/dsqrt(x)
         enddo
      enddo
      call chartree(bm,Nx,Nz,Wx,Wz,scr,vh,cden,start,accuracyh) 
      do ix = 1,Nx
         x = ix*Wx
         do iz = -Nz,Nz
            vh(iz,ix)=vh(iz,ix)/dsqrt(x)
         enddo
      enddo
      call chi0(m,Lq,Zt,Nx,Nz,Wx,Wz,
     & miu,bt,tol,accuracyc,delta,
     & w,vect,pot,
     & freq,vh,r)
      r = dn0-dn+r

      call chi0(m,Lq,Zt,Nx,Nz,Wx,Wz,
     & miu,bt,tol,accuracyc,deltap,
     & w,vect,pot,
     & freq,dnp,rp)
      do ix = 1,Nx
         x = ix*Wx
         do iz = -Nz,Nz
            cden(iz,ix)=-4d0*pi*rp(iz,ix)/dsqrt(x)
         enddo
      enddo
      call chartree(bm,Nx,Nz,Wx,Wz,scr,vh,cden,start,accuracyh)
      do ix = 1,Nx
         x = ix*Wx
         do iz = -Nz,Nz
            vh(iz,ix)=vh(iz,ix)/dsqrt(x)
         enddo
      enddo
      rp=dn0-dnp+vh
c     ********************************************************
c     starts the iterative process
c     ********************************************************
      p=r
      pp=rp
      do it = 1,maxit
c        ******************
c        computes q=(1-M)*p
c        *************************************************
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
               cden(iz,ix)=-4d0*pi*p(iz,ix)/dsqrt(x)
            enddo
         enddo
         call chartree(bm,Nx,Nz,Wx,Wz,scr,vh,cden,start,accuracyh)
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
               vh(iz,ix)=vh(iz,ix)/dsqrt(x)
            enddo
         enddo
         call chi0(m,Lq,Zt,Nx,Nz,Wx,Wz,
     &   miu,bt,tol,accuracyc,delta,
     &   w,vect,pot,
     &   freq,vh,q)
         q=p-q
c        ******************
c        computes qp=(1-M*)pp
c        *************************************************
         call chi0(m,Lq,Zt,Nx,Nz,Wx,Wz,
     &   miu,bt,tol,accuracyc,deltap,
     &   w,vect,pot,
     &   freq,pp,qp)
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
               cden(iz,ix)=-4d0*pi*qp(iz,ix)/dsqrt(x)
            enddo
         enddo
         call chartree(bm,Nx,Nz,Wx,Wz,scr,vh,cden,start,accuracyh)  
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
               vh(iz,ix)=vh(iz,ix)/dsqrt(x)
            enddo
         enddo
         qp=pp-vh
c        *************************************************
c        makes the updates
c        *************************************************
         rr=0D0
         z=0D0
         do ix = 1,Nx
            do iz = -Nz,Nz
               rr=rr+dconjg(rp(iz,ix))*r(iz,ix)
               z=z+dconjg(pp(iz,ix))*q(iz,ix)
            enddo
         enddo
c         write(*,*) 'z',z
         z=rr/z
         dn=dn+z*p
         dnp=dnp+dconjg(z)*pp
         norm = 0d0
         check=0d0
         do ix = 1,Nx
            do iz = -Nz,Nz
               check=check+abs(z*p(iz,ix))
               norm = norm + abs(dn(iz,ix))
            enddo
         enddo
c         write(*,*) 'norm', norm
         r=r-z*q
         rp=rp-dconjg(z)*qp
c         write(*,*) 'absr',sum(abs(r)),'absrp',sum(abs(rp))
         beta=0D0
         do ix = 1,Nx
            do iz = -Nz,Nz
               beta=beta + dconjg(rp(iz,ix))*r(iz,ix)
            enddo
         enddo
         beta=beta/rr
         p=r+beta*p
         pp=rp+dconjg(beta)*pp     
c         write(*,*) 'absp',sum(abs(p)),'abspp',sum(abs(pp))    
c        ************************************************
c        checks accuracy
c        ************************************************
c         write(*,*) it,check/norm
         if (check/norm.lt.accuracy) then
            goto 111
         end if

      enddo
      write(*,*) 'accuracy not achieved in solverpa'

111   continue

      return
      end




