c     --------------------------------------------------
      program eqdavidson                                              

      implicit none
      integer Nx,Hx,Nz,Hz,Nr,Natoms
      integer dftit,m,Lq,Lmax
      integer charge,dosp
      double precision Lx,Lz,Wx,Wz,Zt
      double precision Temp,bt,prec,scr,tol
      double precision nj,vj,cut,r1,r2,def1,def2
      double precision cent
c     *******************************************************
c     Nx = x grid points; Nz = z grid points(odd nr)
c     Nr = total grid points
c     Lx, Lz = radius and height of the super-cell
c     Wx, Wz = grid constants
c     b = cell of the chain
c     Natoms = atoms in the chain
c     miuiterations = chemical pot. iterations
c     dftit = dft iterations
c     Natoms = Nr of atoms
c     Temp = temperature
c     prec = precision in the density
c     tol = dNel/Nel
c     scr = screening length in the coulomb kernel
c     Ll = leads length
c     charge = atomic valance charge
c     dosp = nr. of points in dos plot
c     *******************************************************

c     *********************************************************
c     Define size of computational box and other parameters:
c     *********************************************************
      parameter(Nx=99,Hx=2*Nx,Nz=99,Hz=Nz,Nr=Nx*(2*Nz+1),
     &  Lx=50D0,Lz=50D0,Wx=Lx/(Nx+1D0),Wz=Lz/(Nz+1D0),Lq=15,
     &  tol=0.00001d0,dftit=90,Temp=1000d0,
     &  scr=5d0,nj=0.0088d0,vj=-0.1d0,dosp=100)
      parameter(r1=12d0,r2=24d0,cent=0d0)
c    ***********************************************************
      integer Lm,miuit,Nel,flag,count
      integer i,j,ix,iz,kx,kz,jj,iteration
      double precision iTemp,miu,miup,vol,voleads
      double precision pi,dN,fact,ek,lim,dlt
      double precision x,z,mixing,xm,checkp
      double precision Nelectrons,NelectronsP,norm,check,n0
      double precision para,rans,delta,delta0,dos,proj,en
      integer nv(-1:Lq)
      double precision Ze(100),z0(100)
      double precision den(-Nz-2:Nz+2,-1:Nx+2)
      double precision nden(-Nz-2:Nz+2,-1:Nx+2)
      double precision newden(-Nz-2:Nz+2,-1:Nx+2)
      double precision cden(-Hz-2:Hz+2,-1:Hx+2)
      double precision potD(-Hz-2:Hz+2,-1:Hx+2)
      double precision potN(-Hz-2:Hz+2,-1:Hx+2)
      double precision varh(-Hz-2:Hz+2,-1:Hx+2)
      double precision var(-Nz-2:Nz+2,-1:Nx+2)
      double precision pot(-Nz-2:Nz+2,-1:Nx+2)
      double precision tpot(-Nz-2:Nz+2,-1:Nx+2)
      double precision pseudov(-Nz-2:Nz+2,-1:Nx+2)
      double precision, allocatable:: vect(:,:,:,:),w(:,:),
     & work(:,:,:)      
      double precision vv(-Nz-2:Nz+2,-1:Nx+2,100)
c      double precision vect(-Nz-2:Nz+2,-1:Nx+2,m,0:10)
c      double precision w(m,0:10)
      external spectrum,dhartree,nhartree
      double precision vxc,jell,shell,ran
      double precision Sodium3S,PseudoPotential


      pi=4d0*atan(1d0)
c     *****************************************************
c     initialize
c     *****************************************************
      den = 0d0
      nden = 0d0
      cden = 0d0
      newden = 0d0
      potD = 0d0
      potN = 0d0
      var = 0d0
      varh = 0d0
      pot = 0d0
      tpot = 0d0
      pseudov = 0d0
c ***********************************************************
c sets the jellium
c ***********************************************************
      vol=0d0
      voleads=0d0
      Zt=0d0
      do ix = 1,Nx
         x = ix*Wx
         do iz = -Nz,Nz
            z=iz*Wz
c           ***************************************************
c           sets the nuclear, electron and net charge densities
c           ***************************************************
            nden(iz,ix) = nden(iz,ix)
     &     +nj*shell(x,z,r1,r2,cent)
            den(iz,ix) = nden(iz,ix)
            cden(iz,ix) = nden(iz,ix)-den(iz,ix)
c           ***************************************************
c           sets the jellium pseudopotential
c           ***************************************************
            pseudov(iz,ix)=pseudov(iz,ix)
     &     -0.16*shell(x,z,r1,r2,cent)
c           ***************************************************
c           calculates the total nr of electrons (rs = 3 au)
c           ***************************************************
            Zt=Zt+nj*shell(x,z,r1,r2,cent)
     &     *2d0*pi*x*Wx*Wz
         enddo
      enddo
      write(*,*) 'nj=',nj,'Lx=',Lx,'Lz=',Lz,'Zt =',Zt
c     *********************************************************
c     rough guess for miu;
c     *********************************************************
      miu=-0.1d0
      delta0=0.01d0
c     *********************************************************
c     Starts dft iteration
c     *********************************************************
c     generates the starting vect from plane waves
c     *********************************************************
      count = 0
      lim = 4.5d0
      do kx = 1,100
         do kz = 0,100
            ek=dsqrt(kx*kx+kz*kz+0d0)
            if(ek.lt.lim) then
              if(kz.eq.0) then
                 count = count + 1
                 do ix = 1,Nx
                    x = ix*Nx
                    do iz = -Nz,Nz
                       z = iz*Nz
                       vv(iz,ix,count)=dsin(x*kx*pi/Lx)
                    enddo
                 enddo
              else
                count = count + 1
                do ix = 1,Nx
                   x = ix*Nx
                   do iz = -Nz,Nz
                      z = iz*Nz
                      vv(iz,ix,count)=dsin(x*kx*pi/Lx)
     &               *dsin(z*kz*pi/Lz)
                   enddo
                enddo
                count = count + 1
                do ix = 1,Nx
                   x = ix*Nx
                   do iz = -Nz,Nz
                      z = iz*Nz
                      vv(iz,ix,count)=dsin(x*kx*pi/Lx)
     &               *dcos(z*kz*pi/Lz)
                   enddo
                enddo
              end if
            end if
         enddo
      enddo
      write(*,*) '***********************************'
      write(*,*) 'we start with m =',count,'vectors'
      write(*,*) '***********************************'
      m = count
      nv = count
      allocate(vect(-Nz-2:Nz+2,-1:Nx+2,2*m,0:Lq),w(2*m,0:Lq),
     & work(-Nz-2:Nz+2,-1:Nx+2,2*m))
      do Lm = 0,Lq
         do i = 1,m
            vect(:,:,i,Lm) = vv(:,:,i)
         enddo
       enddo
      potD = 0d0
      potN = 0d0
      pot = 0d0
      flag=0
      do iteration = 1,dftit
         write(*,*) '******************************'
         write(*,*) 'iteration',iteration
         bt=11604.8d0*27.211d0/Temp
         prec = 10d0**(-4d0-3d0*(iteration-1d0)/(dftit-1d0))
         write(*,*) 'precision =',prec
c        ***************************************************************
c        go to the cylindrical cdensity
c        ***************************************************************
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
               z = iz*Wz
               cden(iz,ix)=4d0*pi*cden(iz,ix)*sqrt(x)
            enddo
         enddo
c        ***********************************************************
c        constructs the Dirichlet solution of the poisson eq
c        ***********************************************************
         call dhartree(Hx,Hz,Wx,Wz,scr,varh,cden,potD,prec)
         potD=varh
c        ***********************************************************
c        constructs the Neumann solution of the poisson eq
c        ***********************************************************
         call nhartree(Hx,Hz,Wx,Wz,scr,varh,cden,potN,prec)
         potN=varh
c        **********************************************************
c        constructs the effective potential: 
c        (VD+VN)/2+vxc+pseudov
c        **********************************************************    
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
               z = iz*Wz
               pot(iz,ix)=0.5d0*(potD(iz,ix)+potN(iz,ix))/dsqrt(x)
     &         +vxc(den(iz,ix))+pseudov(iz,ix)
            enddo
         enddo
c        ***************************************************************
c        the Veff is updated now
c        calculates the new den
c        ***************************************************************    
         do Lm = 0,Lq
            write(*,*) '     ***************************'
            write(*,*) '     Lm =',Lm
c           ***********************************************************
c           Updates the potential for Lm (adds the centripetal part)
c           ***********************************************************
            do ix = 1,Nx
               x = ix*Wx
               do iz = -Nz,Nz
                  z = iz*Wz
                  tpot(iz,ix)=pot(iz,ix)+(Lm*Lm-0.25d0)/(2d0*x*x)
               enddo
            enddo
            if(iteration.eq.1) then
              if(Lm.eq.0) then
                work = vect(:,:,:,Lm)
              else
                work = vect(:,:,:,Lm-1)
                nv(Lm)=nv(Lm-1)
              end if
            else
              work = vect(:,:,:,Lm)
            end if
            write(*,*) '     nv =',nv(Lm)
c           ***********************************************************
c           calls the Davidson subroutine
c           ***********************************************************
	    call spectrum(Lm,Nx,Nz,Wx,Wz,m,nv(Lm),tpot,w(:,Lm),
     &      work,prec)
            vect(:,:,:,Lm) = work
c           write(*,*) (w(i,Lm), i=1,m)
c           write(*,*) (1d0/(1d0+dexp((w(i,Lm)-miu)*bt)), i=1,m)
c           **********************************************************
c           checks if Lm contribution to the density is in the noise
c           **********************************************************
            check=0d0
            count = 0
            do i = 1,m
               if(Lm.eq.0) then
                 fact = 2d0/(1d0+dexp((w(i,Lm)-miu)*bt))/Zt
               else
                 fact = 4d0/(1d0+dexp((w(i,Lm)-miu)*bt))/Zt
               end if
               check = check + fact
               if(fact.gt.tol) count = count +1
            enddo
            if(iteration.eq.1) nv(Lm) = count + 1
            write(*,*) '     nv =',nv(Lm),'Lm=',Lm,'check=',check
            Lmax = Lm
            if(check.lt.tol) goto 23
         enddo
c        **************************************************************
c        the cycle over Lm is closed
c        **************************************************************
 23      continue 
c        **************************************************************
c         delta=delta0
          do miuit = 1,10
c           *************************************************************
c           checks the neutrality and search for the new miu
c           **************************************************************
            Nelectrons=0d0
            do Lm =0,Lmax
               do i = 1,nv(Lm)
                  if(Lm.eq.0) then
                    fact=2d0/(1d0+dexp((w(i,Lm)-miu)*bt))
                  else
                    fact=4d0/(1d0+dexp((w(i,Lm)-miu)*bt))
                  endif
                  Nelectrons=Nelectrons+fact
               enddo
            enddo
            if (abs(Nelectrons/Zt-1d0).lt.tol) then
               goto 100
            end if
c           ****************************************************
c           if this is the first time we change miu, add a delta
c           ****************************************************
            if(flag.eq.0) then
              if (Nelectrons.gt.Zt) then
                  miup=miu
                  miu=miu-delta0
              else
                  miup=miu
                  miu=miu+delta0
              end if
              NelectronsP=Nelectrons
              flag=1
              goto 99
            end if
c           ************************************************
c           this is all that is done for the first time
c           ************************************************
c           for the rest of the iterations
c           if miuit=1 use dN from the previous iteration
c           if miuit>1 calculate dN
c           ************************************************
            if(miuit.gt.1) then
              dN=abs((Nelectrons-NelectronsP)/(miu-miup))
            end if
            miup=miu
            miu=miu+(Zt-Nelectrons)/dN
            NelectronsP=Nelectrons
 99         continue
         enddo
 100     continue
         write(*,*) '     Nelectrons/Zt =',Nelectrons/Zt,'miu =',miu
c        *******************************************************
c        the correct miu was found
c        *******************************************************
         newden=0d0
         do Lm = 0,Lmax
            do i = 1,nv(Lm)
               if(Lm.eq.0) then
                 fact=2d0/(1d0+dexp((w(i,Lm)-miu)*bt))
               else
                 fact=4d0/(1d0+dexp((w(i,Lm)-miu)*bt))
               endif
               fact = fact/(2d0*pi*Wx*Wz)
               newden=newden+fact*vect(:,:,i,Lm)*vect(:,:,i,Lm)
            enddo
         enddo
c        ******************************************************
c        final touch and check
c        ******************************************************
         check = 0d0
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
               newden(iz,ix) = newden(iz,ix)/x
               check = check + 2d0*pi*newden(iz,ix)*x*Wx*Wz
            enddo
          enddo

c            do iz = 1,Nz
c               newden(iz,ix) = (newden(iz,ix)+newden(-iz,ix))/(2d0*x)
c               newden(-iz,ix) = newden(iz,ix)
c               check = check + 4d0*pi*newden(iz,ix)*x*Wx*Wz
c            enddo
c         enddo
         checkp=0d0
         do Lm =0,Lmax
            do i = 1,nv(Lm)
               if(Lm.eq.0) then
                 fact=2d0/(1d0+dexp((w(i,Lm)-miu)*bt))
               else
                 fact=4d0/(1d0+dexp((w(i,Lm)-miu)*bt))
               endif
               checkp=checkp+fact
            enddo
         enddo
c        ******************************************************
c        the new density was found
c        ******************************************************
         write(*,*) 'check newden',check/Zt,checkp/Zt
c        ******************************************************
c        puts the feedback
c        *******************************************************
         para=0d0
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
               z = iz*Wz
               para=para+abs(newden(iz,ix)-den(iz,ix))/Nr
            enddo
         enddo
         mixing=0.05d0
         write(*,*) 'new-old=',para,'mixing =',mixing
         den = mixing*newden+(1d0-mixing)*den        
         do ix = 1,Nx
            do iz = -Nz,Nz
               cden(iz,ix)=nden(iz,ix)-den(iz,ix)
            enddo
         enddo
         if(para.lt.0.00001) goto 113
      enddo
 113  continue
c     ***********************************************************
c     equilibrium is completed
c     ***********************************************************
c     calculates the projected density of states
c     ***********************************************************
      dlt = 0.01d0
      do iz = -Nz,Nz
         z = iz*Wz
         do j = 1,dosp
            en = w(1,0)-0.05d0+(w(nv(Lm),Lm)-w(1,0)+0.05d0)*j/dosp
            dos = 0d0
            do Lm = 0,Lmax
               do i = 1,nv(Lm)
                  proj = 0d0
                  do ix = 1,Nx
                     x = ix*Wx
                     proj = proj + 2d0*pi*Wx*vect(iz,ix,i,Lm)**2
                  enddo
                  fact = (pi/dlt)/((w(i,Lm)-en)**2+dlt*dlt)
                  dos = dos + fact*proj
               enddo
            enddo
c            write(28,*) z,en,dos
         enddo
      enddo
c     ***********************************************************


c-----------------------------------------------------------------------
c  In this section we write everything out into text files:

      open(19,file='myPara.txt')
      open(20,FILE='den.txt')
      open(21,FILE='pot.txt')
      open(22,FILE='e.txt')
      open(23,FILE='psi.txt')

c  Parameters:
      write(19,*) m
      write(19,*) Nx
      write(19,*) Nz
      write(19,*) Wx
      write(19,*) Wz
      write(19,*) Lq
      write(19,*) miu
      write(19,*) Temp
      write(19,*) bt
      write(19,*) Zt
      write(19,*) tol
c  I don't think I use these to calculate screening, but I
c  write them anyway.
      write(19,*) Lx
      write(19,*) Lz
      write(19,*) nj
      write(19,*) vj
      write(19,*) Lmax       
      write(19,*) prec
      write(19,*) scr
      write(19,*) dosp


c  den,pot  ------ Things that do NOT depend on Lm:
      do ix = 1,Nx
         do iz = -Nz,Nz
            write(20,*) den(iz,ix)
            write(21,*) pot(iz,ix)
         enddo
      enddo

c  e,psi  ------ Things that DO depend on Lm:
         do Lm = 0,Lq,1
            do i = 1,m,1
               write(22,*) w(i,Lm)
               do ix=-1,Nx+2,1
                  do iz = -Nz-2,Nz+2,1
                     write(23,*) vect(iz,ix,i,Lm)
     &                /dsqrt(Wx*Wz)                     
                  enddo
               enddo
            enddo
         enddo


      close(19)
      close(20)
      close(21)
      close(22)
      close(23)
c-----------------------------------------------------------------------      
c End of the program:     
      stop

      end
c-----------------------------------------------------------------------
c  Here we start writing the functions we will use:

      function shell(x,z,r1,r2,cent)
      double precision shell,x,z,r1,r2,cent
      double precision rIn,rOut
      rOut  = dsqrt(x*x+z*z)
      rIn   = dsqrt(x*x+(z-cent)*(z-cent))
      if(rIn.ge.r1.and.rOut.le.r2) then
         shell = 1d0
      else
         shell = 0d0
      endif
      return
      end


      function jellium(xm,zm,x,z)
      double precision jellium,x,z
      double precision xm,zm
c      jellium=1d0/(1+exp((x-xm)/0.5d0))*
c     &  1d0/(1+exp((z-zm)/0.5d0)) 
      jellium=0d0
      if (abs(z).lt.zm) then
         if (x.lt.xm) then
               jellium=1d0
         end if
      end if
      return 
      end


      function vxc(den)
      double precision vxc,den,rs
      if(den.gt.0.0000001) then
        rs=(3d0/(16d0*atan(1d0)*den))**(1d0/3d0)
      else
        rs=(3d0/(16d0*atan(1d0)*0.0000001))**(1d0/3d0)
      end if
       
      if (rs.lt.1d0) then
        vxc=-0.611d0/rs+0.03109d0*log(rs)-0.0584d0+
     &  0.0013d0*rs*log(rs)-0.0071d0*rs
      else
        vxc=-0.611d0/rs-0.1423D0*(1d0+1.23d0*dsqrt(rs)+0.445d0*rs)/
     &   (1d0+1.0592d0*dsqrt(rs)+0.3334d0*rs)**2d0
      end if
      return
      end




      subroutine dhartree(Nx,Nz,Wx,Wz,scr,var,b,in,accuracy)
      integer Nr,Nx,Nz
      double precision Wx,Wz,scr,accuracy
      double precision var(-Nz-2:Nz+2,-1:Nx+2)
      double precision in(-Nz-2:Nz+2,-1:Nx+2)
      double precision b(-Nz-2:Nz+2,-1:Nx+2)
c     **********************************************
      integer ix,iz,it,maxit
      double precision eps,beta,z,norm,rr,gx,gz
      double precision x,check,bc,bc1,bc2
      double precision, allocatable:: p(:,:),q(:,:),r(:,:)
      allocate(p(-Nz-2:Nz+2,-1:Nx+2),q(-Nz-2:Nz+2,-1:Nx+2),
     & r(-Nz-2:Nz+2,-1:Nx+2))
c     **********************************************
c     define some constants
c     **********************************************
      Nr = Nx*(2*Nz+1)
      maxit = 500
      gx = 1d0/(12d0*Wx*Wx)
      gz = 1d0/(12d0*Wz*Wz)
      bc1 = (30d0-0.75d0)*dsqrt(2d0)-16d0*dsqrt(3d0)+2d0
      bc2 = -3d0+dsqrt(3d0)-bc1*dsqrt(2d0)
c     *********************************************************
c     calculates the norm of b
c     *********************************************************
      norm=0D0
      do ix = 1,Nx
         do iz =-Nz,Nz
            norm=norm+b(iz,ix)**2d0
         enddo
      enddo
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
            r(iz,ix)=(-30d0*(gx+gz)+0.25d0/(x*x)
     &     -1d0/(scr*scr))*in(iz,ix)
     &     +16d0*gz*(in(iz+1,ix)+in(iz-1,ix))
     &     -gz*(in(iz+2,ix)+in(iz-2,ix))
     &     +16d0*gx*(in(iz,ix+1)+in(iz,ix-1))
     &     -gx*(in(iz,ix+2)+in(iz,ix-2))
         enddo
      enddo
      r=b-r
c     ********************************************************
c     starts the iterative process
c     ********************************************************
      p=r
      var=in
      do it = 1,maxit
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
c        calculates q=s*p
c        *************************************************
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
               q(iz,ix)=(-30d0*(gx+gz)+0.25d0/(x*x)
     &        -1d0/(scr*scr))*p(iz,ix)
     &        +16d0*gz*(p(iz+1,ix)+p(iz-1,ix))
     &        -gz*(p(iz+2,ix)+p(iz-2,ix))
     &        +16d0*gx*(p(iz,ix+1)+p(iz,ix-1))
     &        -gx*(p(iz,ix+2)+p(iz,ix-2))
            enddo
         enddo
c        *************************************************
c        makes the updates
c        *************************************************
         rr=0D0
         z=0D0
         do ix = 1,Nx
            do iz = -Nz,Nz
               rr=rr+r(iz,ix)*r(iz,ix)
               z=z+p(iz,ix)*q(iz,ix)
            enddo
         enddo
         z=rr/z
         var=var+z*p
         check=0d0
         do ix = 1,Nx
            do iz = -Nz,Nz
               check=check+abs(z*p(iz,ix))
            enddo
         enddo
         r=r-z*q
         beta=0D0
         do ix = 1,Nx
            do iz = -Nz,Nz
               beta=beta+r(iz,ix)*r(iz,ix)
            enddo
         enddo
         beta=beta/rr
         p=r+beta*p
c        ************************************************
c        checks accuracy
c        ************************************************
         if (check/Nr.lt.accuracy) then
            goto 111
         end if
      enddo
      write(*,*) 'accuracy not achieved in dhartree'
111   return
      end


      subroutine nhartree(Nx,Nz,Wx,Wz,scr,var,b,in,accuracy)
      integer Nr,Nx,Nz
      double precision Wx,Wz,scr,accuracy
      double precision var(-Nz-2:Nz+2,-1:Nx+2)
      double precision in(-Nz-2:Nz+2,-1:Nx+2)
      double precision b(-Nz-2:Nz+2,-1:Nx+2)
c     **********************************************
      integer ix,iz,it,maxit
      double precision eps,beta,z,norm,rr,gx,gz
      double precision x,check,bc,bc1,bc2
      double precision, allocatable:: p(:,:),q(:,:),r(:,:)
      allocate(p(-Nz-2:Nz+2,-1:Nx+2),q(-Nz-2:Nz+2,-1:Nx+2),
     & r(-Nz-2:Nz+2,-1:Nx+2))
c     **********************************************
c     define some constants
c     **********************************************
      Nr = Nx*(2*Nz+1)
      maxit = 500
      gx = 1d0/(12d0*Wx*Wx)
      gz = 1d0/(12d0*Wz*Wz)
      bc1 = (30d0-0.75d0)*dsqrt(2d0)-16d0*dsqrt(3d0)+2d0
      bc2 = -3d0+dsqrt(3d0)-bc1*dsqrt(2d0)
c     *********************************************************
c     calculates the norm of b
c     *********************************************************
      norm=0D0
      do ix = 1,Nx
         do iz =-Nz,Nz
            norm=norm+b(iz,ix)**2d0
         enddo
      enddo
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

      in(:,Nx+1) = in(:,Nx)
      in(:,Nx+2) = in(:,Nx-1)

      in(-Nz-1,:) = in(-Nz,:)
      in(-Nz-2,:) = in(-Nz+1,:)

      in(Nz+1,:) = in(Nz,:)
      in(Nz+2,:) = in(Nz-1,:)
c     *********************************************************
c     calculates the rest
c     *********************************************************
      r = 0d0
      do ix = 1,Nx
         x = ix*Wx
         do iz = -Nz,Nz
            r(iz,ix)=(-30d0*(gx+gz)+0.25d0/(x*x)
     &     -1d0/(scr*scr))*in(iz,ix)
     &     +16d0*gz*(in(iz+1,ix)+in(iz-1,ix))
     &     -gz*(in(iz+2,ix)+in(iz-2,ix))
     &     +16d0*gx*(in(iz,ix+1)+in(iz,ix-1))
     &     -gx*(in(iz,ix+2)+in(iz,ix-2))
         enddo
      enddo
      r=b-r
c     ********************************************************
c     starts the iterative process
c     ********************************************************
      p=r
      var=in
      do it = 1,maxit
c        *************************************************
c        this takes care of the boundary conditions
c        *************************************************
         p(:,0) = (16d0-bc1)*p(:,1)
         p(:,-1) = 16d0*p(:,0)-(30d0+bc2)*p(:,1)
     &   +(16d0-bc1)*p(:,2)

         p(:,Nx+1) = p(:,Nx)
         p(:,Nx+2) = p(:,Nx-1)

         p(-Nz-1,:) = p(-Nz,:)
         p(-Nz-2,:) = p(-Nz+1,:)
      
         p(Nz+1,:) = p(Nz,:)
         p(Nz+2,:) = p(Nz-1,:)
c        *************************************************
c        calculates q=s*p
c        *************************************************
         do ix = 1,Nx
            x = ix*Wx
            do iz = -Nz,Nz
               q(iz,ix)=(-30d0*(gx+gz)+0.25d0/(x*x)
     &        -1d0/(scr*scr))*p(iz,ix)
     &        +16d0*gz*(p(iz+1,ix)+p(iz-1,ix))
     &        -gz*(p(iz+2,ix)+p(iz-2,ix))
     &        +16d0*gx*(p(iz,ix+1)+p(iz,ix-1))
     &        -gx*(p(iz,ix+2)+p(iz,ix-2))
            enddo
         enddo
c        *************************************************
c        makes the updates
c        *************************************************
         rr=0D0
         z=0D0
         do ix = 1,Nx
            do iz = -Nz,Nz
               rr=rr+r(iz,ix)*r(iz,ix)
               z=z+p(iz,ix)*q(iz,ix)
            enddo
         enddo
         z=rr/z
         var=var+z*p
         check=0d0
         do ix = 1,Nx
            do iz = -Nz,Nz
               check=check+abs(z*p(iz,ix))
            enddo
         enddo
         r=r-z*q
         beta=0D0
         do ix = 1,Nx
            do iz = -Nz,Nz
               beta=beta+r(iz,ix)*r(iz,ix)
            enddo
         enddo
         beta=beta/rr
         p=r+beta*p
c        ************************************************
c        checks accuracy
c        ************************************************
         if (check/Nr.lt.accuracy) then
            goto 111
         end if
      enddo
      write(*,*) 'accuracy not achieved in nhartree'
111   return
      end

      subroutine spectrum(Lm,Nx,Nz,Wx,Wz,m,nv,pot,eig,vec,prec)
      integer Nr,Nx,Nz,m,nv,Lm
      double precision Wx,Wz,bc,bc1,bc2,T,miu,prec
      double precision pot(-Nz-2:Nz+2,-1:Nx+2)
      double precision eig(2*m)
      double precision vec(-Nz-2:Nz+2,-1:Nx+2,2*m)
c     ****************************************************************
      integer maxit,ix,iz,it,i,j,mdown,mup
      double precision test,pi,check,gx,gz
      integer n,kx,kz,info
      double precision, allocatable:: s(:,:),z(:,:),work(:),
     & x(:,:),y(:,:),ortho(:,:,:),w(:),wg(:),oldeig(:)
      double precision element

      allocate(
     & s(2*nv,2*nv),
     & x(-Nz-2:Nz+2,-1:Nx+2),
     & y(-Nz-2:Nz+2,-1:Nx+2),
     & work(6*nv-2),
     & w(2*nv),
     & wg(2*nv),
     & oldeig(2*nv),
     & z(2*nv,2*nv),
     & ortho(-Nz-2:Nz+2,-1:Nx+2,2*nv))
      pi=4d0*atan(1d0)
c     ****************************************************************
c     sets the max iterations and other constants
c     ****************************************************************
      Nr = Nx*(2*Nz+1)
      maxit=5000
      accuracy=0.0001d0
      gx=1d0/(24d0*Wx*Wx)
      gz=1d0/(24d0*Wz*Wz)
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
      oldeig=eig(1:2*nv)
      do it = 1,maxit
c     *****************************************************************
c     calculate vec nv+1,...,mup
c     *****************************************************************
      do i=1,nv
         x = vec(:,:,i)
         x(:,0) = (16d0-bc1)*x(:,1)
         x(:,-1) = 16d0*x(:,0)-(30d0+bc2)*x(:,1)+(16d0-bc1)*x(:,2)

         x(:,Nx+1) = 0d0
         x(:,Nx+2) = -x(:,Nx)

         x(-Nz-1,:) = x(Nz,:)
         x(-Nz-2,:) = x(Nz-1,:)

         x(Nz+1,:) = x(-Nz,:)
         x(Nz+2,:) = x(-Nz+1,:)

         do ix = 1,Nx
            do iz = -Nz,Nz
               vec(iz,ix,2*nv-i+1)=(30d0*(gx+gz)+pot(iz,ix))*x(iz,ix)
     &        -16d0*gz*(x(iz+1,ix)+x(iz-1,ix))
     &        +gz*(x(iz+2,ix)+x(iz-2,ix))
     &        -16d0*gx*(x(iz,ix+1)+x(iz,ix-1))
     &        +gx*(x(iz,ix+2)+x(iz,ix-2))
            enddo
         enddo
      enddo
c     *****************************************************************
c     calculate the overlap matrix for the mup vectors
c     *****************************************************************
      do j = 1,2*nv
         do i = j,2*nv
            element=0d0
            do ix = 1,Nx
               do iz = -Nz,Nz
                  element=element
     &           +vec(iz,ix,i)*vec(iz,ix,j)
               enddo
            enddo
            s(1+i-j,j)=element
         enddo
      enddo
c     *****************************************************************
c     diagonalze the overlap matrix
c     *****************************************************************
      call dsbev('v','l',2*nv,2*nv-1,s,2*nv,w,z,2*nv,work,info)
c      write(*,*) 'eig of overlap'
c      write(*,*) (w(i), i=1,2*nv)
c     *****************************************************************
c     orthonormalize vec
c     *****************************************************************
      mup=0
      do i = 1,2*nv
        if(w(i).gt.prec*prec) then
         mup = mup+1
         ortho(:,:,mup) = 0d0
         do j = 1,2*nv
           ortho(:,:,mup)=ortho(:,:,mup)+vec(:,:,j)*z(j,i)/dsqrt(w(i))
         enddo
        end if
      enddo
c     ****************************************************************
c     test
c     ****************************************************************
c      test = 0d0
c      do ix = 1,Nx
c         do iz = -Nz,Nz
c            test = test + ortho(iz,ix,1)*ortho(iz,ix,1)
c         enddo
c      enddo
c      write(*,*) 'orthotest =',test
c     *****************************************************************
c     compute the matrix elements
c     *****************************************************************
      do j = 1,mup
         x = ortho(:,:,j)

         x(:,0) = (16d0-bc1)*x(:,1)
         x(:,-1) = 16d0*x(:,0)-(30d0+bc2)*x(:,1)+(16d0-bc1)*x(:,2)

         x(:,Nx+1) = 0d0
         x(:,Nx+2) = -x(:,Nx)

         x(-Nz-1,:) = x(Nz,:)
         x(-Nz-2,:) = x(Nz-1,:)

         x(Nz+1,:) = x(-Nz,:)
         x(Nz+2,:) = x(-Nz+1,:)

         do i = j,mup
            element = 0d0
            do ix = 1,Nx
               do iz = -Nz,Nz
                  element=element+ortho(iz,ix,i)
     &           *((30d0*(gx+gz)+pot(iz,ix))*x(iz,ix)
     &           -16d0*gz*(x(iz+1,ix)+x(iz-1,ix))
     &           +gz*(x(iz+2,ix)+x(iz-2,ix))
     &           -16d0*gx*(x(iz,ix+1)+x(iz,ix-1))
     &           +gx*(x(iz,ix+2)+x(iz,ix-2)))
               enddo
            enddo
            s(1+i-j,j)=element
         enddo
      enddo
c     ******************************************************************
c     diagonalize
c     ******************************************************************
      call dsbev('v','l',mup,mup-1,s,2*nv,wg,z,2*nv,work,info)
c      write(*,*) 'eig of H'
c      write(*,*) (eig(i), i=1,mup)
c     ******************************************************************
c     recalculate the eigenvectors
c     ******************************************************************
      do i = 1,mup
         vec(:,:,i) = 0d0
         do j = 1,2*nv
            vec(:,:,i) = vec(:,:,i) + ortho(:,:,j)*z(j,i)
         enddo
      enddo
c     *****************************************************************
c     test
c     *****************************************************************
c      test = 0d0
c      do ix = 1,Nx
c         do iz = -Nz,Nz
c            test = test + vec(iz,ix,1)*vec(iz,ix,1)
c         enddo
c      enddo
c      write(*,*) 'test =',test
c     *****************************************************************
c     convergence
c     *****************************************************************
      check = 0d0
      do i = 1,nv
         check = check + abs(oldeig(i)-wg(i))
      enddo
      write(*,*) 'it =',it,'check =',check/nv
      if(check/nv.lt.prec) goto 101
      oldeig = wg
      enddo
      write(*,*) 'accuracy not achieved in spectrum'
 101  continue
      write(*,*) '     precision achieved in',it,'steps'
      eig(1:2*nv) = wg
      return
      end
