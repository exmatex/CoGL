c Copyright (c) 2013, Los Alamos National Security, LLC
c All rights reserved.
c Copyright 2013. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL),
c which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.

c NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

c If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

c Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
c ·         Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
c ·         Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other
c           materials provided with the distribution.
c ·         Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used
c           to endorse or promote products derived from this software without specific prior written permission.

c THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
c WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
c DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
c OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
c ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

        implicit none
        integer i,j,k,l,nip,nop,nmp,ne1,ncount,i1,j1,i2,j2
     *     ,ntime,itime,interval,iseed,np1,nr1,i3,id
	parameter(k=64)

	real e1(k,k,k),e2(k,k,k),e3(k,k,k)
	real e4(k,k,k),e5(k,k,k),e6(k,k,k)
	real ux(k,k,k),uy(k,k,k),uz(k,k,k)
	real uxdt(k,k,k),uydt(k,k,k),uzdt(k,k,k)
	real gxx(k,k,k),gyy(k,k,k),gzz(k,k,k)
	real gxy(k,k,k),gyz(k,k,k),gxz(k,k,k)
	real gxx1(k,k,k),gyy1(k,k,k),gzz1(k,k,k)
	real gxy1(k,k,k),gyz1(k,k,k),gxz1(k,k,k)
        real G1(k,k,k),G2(k,k,k),G3(k,k,k)
	real lpe2(k,k,k),lpe3(k,k,k),lpvx(k,k,k),lpvy(k,k,k)
        real lpvz(k,k,k)
        real sxxx(k,k,k),sxyy(k,k,k),sxzz(k,k,k)
        real sxyx(k,k,k),syyy(k,k,k),syzz(k,k,k)
        real sxzx(k,k,k),syzy(k,k,k),szzz(k,k,k)
        real etaxx(k,k,k),etayy(k,k,k),etazz(k,k,k)
        real etaxy(k,k,k),etaxz(k,k,k),etayz(k,k,k) 
        real uxx(k,k,k),uyy(k,k,k),uzz(k,k,k)
        real uxy(k,k,k),uxz(k,k,k),uzy(k,k,k)
        real uyx(k,k,k),uyz(k,k,k),uzx(k,k,k)
        real uz_defect(k,k,k)
        real uxxdt(k,k,k),uyydt(k,k,k),uzzdt(k,k,k)
        real uxydt(k,k,k),uxzdt(k,k,k),uzydt(k,k,k)
        real uyxdt(k,k,k),uyzdt(k,k,k),uzxdt(k,k,k)
        real uxx_applied(k,k,k)
        real uxxdt_applied(k,k,k)
        real uxydt_applied(k,k,k),uxzdt_applied(k,k,k)
        real e1dt(k,k,k),e2dt(k,k,k),e3dt(k,k,k)
        real e4dt(k,k,k),e5dt(k,k,k),e6dt(k,k,k)
        real etaxxdt(k,k,k),etayydt(k,k,k),etazzdt(k,k,k)
        real etaxydt(k,k,k),etaxzdt(k,k,k),etayzdt(k,k,k)
        real eta_xx_zero(k,k,k),eta_yy_zero(k,k,k),eta_xy_zero(k,k,k)
        real eta_xx_prime(k,k,k),eta_yy_prime(k,k,k)
     *,eta_xy_prime(k,k,k)
        real eta_xx_ave,eta_yy_ave,eta_xy_ave,eta_xx_ave_theta
        real eta0_xx,eta0_yy,eta0_zz,e1_0,e3_0
        real sigma_xx, sigma_xx_prime
        real A1_prime, A2_prime, A3_prime
        real G11(k,k,k),G21(k,k,k),G31(k,k,k)
        real G41(k,k,k),G51(k,k,k),G61(k,k,k)
        real temp,theta,k_cube,stress_level,stress_orient
        real uapp,udtapp 
        real float_time,eps_0
        real gxx_ave,gyy_ave,gzz_ave,gxy_ave,gxz_ave,gyz_ave
        real gxx_ave_prime,gyy_ave_prime,gzz_ave_prime
        real gxy_ave_prime,gxz_ave_prime,gyz_ave_prime
        real sqrt2, sqrt3, sqrt6


        real an1,delta,h,AS,AC,Ao,D2,Eo,vs,tau,eta,pi,
     *     A_shear,A_bulk,ran,c_cprime,alpha,u_o
        real ran1
        character*20 command

c        READ INPUT
        ! open (11,file='input1',status='unknown')
           ntime= 200000
           interval=20000
           temp=255.
           D2=1.0
           eta=10.
           delta=0.01 
           iseed=3173 
           h=1.
           eps_0=0.0 
           ! READ(11,*) h    
           ! READ(11,*) eps_0
           ! READ(11,'(a3)') command
         !close(11)

         
c      Define A1_prime,A2_prime,A3_prime
        A1_prime=eta
        A2_prime=eta
        A3_prime=eta

c      Define pi
         pi=2.0*asin(1.0)
c        theta=pi/stress_orient
         theta=0.0

         k_cube=float(k*k*k)

c       PARAMETERS
c      F=A(e2^2+e3^2)+Be3(e3^2-3e2^2)+C(e2^2+e3^2)^2
c        + D[(Del e2)^2 + (Del e3)^2]
c        + A_bulk/2 e1^2 + A_comp(e4^2+e5^2+e6^2)

c       grid spacing is h, delta is time step 
c       AS is scaled (A_shear/Ao) modulus
        A_shear=28.0*(10.0**(10.0))
c       AC is scaled (A_bulk/Ao)  modulus  (see notes)
        A_bulk=14.0*(10.0**(10.0))
c       Ao=2.32*(10.0**(9.0))
        Ao=1.97*(10.0**(10.0))
c       D2 is scaled (2*D/Ao) gradient coeffn.
c       Eo=E/B where the term Ee1(e2^2+e3^2) is added to
c       the energy for arresting the transformation.
c       the volume change e1 = -E(e2^2+e3^2)/A_bulk
c        c_cprime=c/c' where c'=c-(E^2/A_bulk), where c' is the
c       renormalized coefficient of c for the single variant case
c       (e2=0,e3=eps_0)

c       Evaluate AS from A_shear,Ao
c       Evaluate AC from A_bulk,Ao

        AS = A_shear/Ao
        AC = A_bulk/Ao
        write(*,*) 'AS=',AS,'AC=',AC

c       TAU, LATTICE PARAMETERS AND BULK MODULI

        tau=(temp-270.)/(295.-270.)
         !tau=-1.

        eta0_xx=(3.795-3.756)/3.756
        eta0_yy=eta0_xx
        eta0_zz=(3.725-3.756)/3.756

        e1_0=(eta0_xx+eta0_yy+eta0_zz)/sqrt(3.)
        e3_0=(eta0_xx+eta0_yy-2.0*eta0_zz)/sqrt(6.)


        A_bulk=19.23*(10.0**(10.0))
c       Ao=2.32*(10.0**(9.0))
        Ao=1.97*(10.0**(10.0))

        write(*,*)'temp',temp,'tau',tau
        write(*,*)'strains',eta0_xx,eta0_yy,eta0_zz
        write(*,*)'e1_0',e1_0,'e3_0',e3_0


c       calculate alpha_0, Eo, c_cprime


c       Eo=-0.5000*AC*e1_0/e3_0
        Eo=0.0
c       c_cprime=1.0+ (alpha_0**2)*(e3_0**2)*AC
        c_cprime=1.0+ ((2.0*Eo**2)/AC)


        sigma_xx=0.0

c------------------------------------------------

c	HERE WE SET THE INITIAL CONDITIONS

          !if (command.eq.'ini') then
 
c         initialize displacements ux,uy,uz and
c         velocities uxdt,uydt,uzdt

c           uz_defect
            uz_defect=0.0000
            u_o=0.4
            alpha=0.5
            write(*,*) 'k',k
             do i=1,k
              do j=1,k
               do l=1,k
 	      ux(i,j,l)=(2.0*ran1(iseed)-1.0)*.1
 	      uy(i,j,l)=(2.0*ran1(iseed)-1.0)*.1
 	      uz(i,j,l)=(2.0*ran1(iseed)-1.0)*.1
	      uxdt(i,j,l)=0.0
	      uydt(i,j,l)=0.0
	      uzdt(i,j,l)=0.0
         write(31,*) ux(i,j,l)
         write(32,*) uy(i,j,l)
         write(33,*) uz(i,j,l)
	      enddo
	     enddo
            enddo

         close(31)
         close(32)
         close(33)

         write(*,*) ux(1,1,1)
         write(*,*) ux(2+1, 1+1, 3+1)
         write(*,*) uy(2+1, 1+1, 3+1)
         write(*,*) uz(2+1, 1+1, 3+1)

          !else
            open (12,file='ux',status='unknown')
            open (13,file='uy',status='unknown')
            open (14,file='uz',status='unknown')
            open (15,file='uxdt',status='unknown')
            open (16,file='uydt',status='unknown')
            open (17,file='uzdt',status='unknown')

c	HERE WE START THE ITERATIONS.


c       EVOLVE
        ncount=0

      sqrt2 = 1.41421356237
      sqrt3 = 1.73205080757
      sqrt6 = 2.44948974278

	DO itime=1,ntime

               write(*,*) 'itime ', itime

c             if (itime .eq. ntime) then
               id = 128871
               i3 = MODULO(id,k)
               i2 = MODULO(id/k,k)
               i1 = id/(k*k)
                write(*,*) 't=', itime, 'e2=', e2(i1+1,i2+1,i3+1)

               id = 0*k*k;
               i3 = MODULO(id,k)
               i2 = MODULO(id/k,k)
               i1 = id/(k*k)
                write(*,*) 't=', itime, 'e2=', e2(i1+1,i2+1,i3+1)

               id = 1*k*k;
               i3 = MODULO(id,k)
               i2 = MODULO(id/k,k)
               i1 = id/(k*k)
                write(*,*) 't=', itime, 'e2=', e2(i1+1,i2+1,i3+1)

               id = 30*k*k;
               i3 = MODULO(id,k)
               i2 = MODULO(id/k,k)
               i1 = id/(k*k)
                write(*,*) 't=', itime, 'e2=', e2(i1+1,i2+1,i3+1)

               id = 63*k*k;
               i3 = MODULO(id,k)
               i2 = MODULO(id/k,k)
               i1 = id/(k*k)
                write(*,*) 't=', itime, 'e2=', e2(i1+1,i2+1,i3+1)

c            Apply strain.
c           2*u/L=eps_0*(t/ntime). Thus eps_0 is max strain after t=ntime
c            u= L*eps_0*(t/ntime)/2. We are assuming that the displacement u
c           at right edge is +u and left edge is -u
                  float_time=float(itime)/float(ntime)
                
               uxx_applied=0.

c       Evaluate gradients of ux,uy,uz


c               write(*,*)'ux',ux(1,1,1)
c               write(*,*)'h',h

               call GRADIENT(ux,k,1,uxx,h)

c               write(*,*)'uxx',uxx(1,1,1)

               call GRADIENT(ux,k,2,uxy,h)
               call GRADIENT(ux,k,3,uxz,h)
               call GRADIENT(uy,k,1,uyx,h)
               call GRADIENT(uy,k,2,uyy,h)
               call GRADIENT(uy,k,3,uyz,h)
               call GRADIENT(uz,k,1,uzx,h)
               call GRADIENT(uz,k,2,uzy,h)
               call GRADIENT(uz,k,3,uzz,h)


c      Evaluate strain components eta_ij (or etaij)
               etaxx=uxx+uxx_applied
               etayy=uyy
               etazz=uzz
               etaxy=0.5*(uxy+uyx)/h
               etaxz=0.5*(uxz+uzx)/h
               etayz=0.5*(uyz+uzy)/h

c      Evaluate symmetry adapted strains
               e1=(etaxx+etayy+etazz)/(sqrt3)
               e2=(etaxx-etayy)/(sqrt2)
               e3=(etaxx+etayy-2.*etazz)/(sqrt6)
               e6=etaxy
               e5=etaxz
               e4=etayz
               
c      Evaluate velocity gradients (ie uxdt, ...)



c               write(*,*)'3'
!c      Evaluate strain components eta_ij (or etaij)

c      Evaluate 
               call LAPLACIAN(e2,k,lpe2,h)
               call LAPLACIAN(e3,k,lpe3,h)
               call LAPLACIAN(uxdt,k,lpvx,h)
               call LAPLACIAN(uydt,k,lpvy,h)
               call LAPLACIAN(uzdt,k,lpvz,h)
c              write(*,*)'5'
c              write(*,*)'6'
                 G1=AC*e1 

                 G2=2.*tau*e2-D2*lpe2-12.*e2*e3+
     .           4.*e2*(e2*e2+e3*e3)
                 G3=2.*tau*e3-D2*lpe3+6.*(e3*e3-e2*e2)+
     .           4.*e3*(e2*e2+e3*e3)
            gxx=(G1/sqrt3)+(G2/sqrt2)+(G3/sqrt6)
            gyy=(G1/sqrt3)-(G2/sqrt2)+(G3/sqrt6)
            gzz=(G1/sqrt3) - (2.*G3/sqrt6)
            gxy=AS*e6
            gxz=AS*e5
            gyz=AS*e4

c------------------------------------------------------------

c           Evaluate (Div. sigma) in component form 
             call GRADIENT(gxx,k,1,sxxx,h)
             call GRADIENT(gyy,k,2,syyy,h)
             call GRADIENT(gzz,k,3,szzz,h)
             call GRADIENT(gxy,k,1,sxyx,h)
             call GRADIENT(gxy,k,2,sxyy,h)
             call GRADIENT(gxz,k,1,sxzx,h)
             call GRADIENT(gxz,k,3,sxzz,h)
             call GRADIENT(gyz,k,3,syzz,h)
             call GRADIENT(gyz,k,2,syzy,h)

         do i=1,k
          do j=1,k
           do l=1,k 
c            UPDATE VELOCITY FIELD u_idt
          uxdt(i,j,l)=uxdt(i,j,l)+delta*(sxxx(i,j,l)+sxyy(i,j,l)
     *    +sxzz(i,j,l)+eta*lpvx(i,j,l))
          uydt(i,j,l)=uydt(i,j,l)+delta*(sxyx(i,j,l)+syyy(i,j,l)
     *    +syzz(i,j,l)+eta*lpvy(i,j,l))
          uzdt(i,j,l)=uzdt(i,j,l)+delta*(sxzx(i,j,l)+syzy(i,j,l)
     *    +szzz(i,j,l)+eta*lpvz(i,j,l))
c            UPDATE DISPLACEMENTS u_i   
          ux(i,j,l)=ux(i,j,l)+delta*uxdt(i,j,l)
          uy(i,j,l)=uy(i,j,l)+delta*uydt(i,j,l)
          uz(i,j,l)=uz(i,j,l)+delta*uzdt(i,j,l)

 
          end do
         end do
        end do


c         DIAGNOSTICS and OUTPUT
                  
                  if(mod(itime,interval) .eq. 0.) then
                    ncount=ncount+1
                    nip=100+ ncount
                       do l=1,k
                        do j=1,k
                         do i=1,k
                          write(nip,*) e3_0*e2(i,j,l), e3_0*e3(i,j,l)
                         enddo
                        enddo
                       enddo

                       do l=1,k
                        do j=1,k
                         do i=1,k
                           write(22,*)  e3_0*e2(i,j,l)
                           write(23,*)  e3_0*e3(i,j,l)
                         enddo
                        enddo
                       enddo

                       
                       close(nip)
                       close(22)
                       close(23)
                 endif
                 if(itime .eq. ntime) then
                       do l=1,k
                        do j=1,k
                         do i=1,k
                          write(12,*)  ux(i,j,l)
                          write(13,*)  uy(i,j,l)
                          write(14,*)  uz(i,j,l)
                          write(15,*)  uxdt(i,j,l)
                          write(16,*)  uydt(i,j,l)
                          write(17,*)  uzdt(i,j,l)
                          enddo
                          enddo
                          enddo
                     

                      close(12)
                      close(13)
                      close(14)
                      close(15)
                      close(16)
                      close(17)
                endif

           
c          END OF TIME ITERATION
             END DO


	stop
        end



      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END


      SUBROUTINE LAPLACIAN(phipp,k,lap,r_mesh)
      integer i,j,i1,j1,i2,j2,k,dim
      integer l,l1,l2,l3,l4,i3,i4,j3,j4
      real r_mesh,r_mesh2

      real phipp(k,k,k),lap(k,k,k)

      r_mesh2=r_mesh**2


       do i=1,k
        do j=1,k
        do l=1,k
              i1=i+1
              i2=i-1
              j1=j+1
              j2=j-1
              l1=l+1
              l2=l-1
             if(i.eq.1)  i2=k
             if(i.eq.k)  i1=1
             if(j.eq.1)  j2=k
             if(j.eq.k)  j1=1
             if(l.eq.1)  l2=k
             if(l.eq.k)  l1=1

             if(i.eq.1) then
             lap(i,j,l)=(phipp(i+1,j,l)+phipp(i,j,l)-2.*phipp(i,j,l)
     .        +phipp(i,j1,l)+phipp(i,j2,l)-2.*phipp(i,j,l)
     .        +phipp(i,j,l1)+phipp(i,j,l2)-2.*phipp(i,j,l))/(r_mesh2)
             end if

             if(i.eq.k) then
             lap(i,j,l)=(phipp(i-1,j,l)+phipp(i,j,l)-2.*phipp(i,j,l)
     .        +phipp(i,j1,l)+phipp(i,j2,l)-2.*phipp(i,j,l)
     .        +phipp(i,j,l1)+phipp(i,j,l2)-2.*phipp(i,j,l))/(r_mesh2)
             end if

             if(i.ne.1.OR.i.ne.k) then
             lap(i,j,l)=(phipp(i1,j,l)+phipp(i2,j,l)-2.*phipp(i,j,l)
     .        +phipp(i,j1,l)+phipp(i,j2,l)-2.*phipp(i,j,l)
     .        +phipp(i,j,l1)+phipp(i,j,l2)-2.*phipp(i,j,l))/(r_mesh2)
             end if
           
         enddo
        enddo
       enddo



       return
       end


      SUBROUTINE GRADIENT(phipp,k,dim,gra,r_mesh)

      integer i,j,l,k,j1,j2,j3,j4, dim
      real r_mesh

      real phipp(k,k,k),gra(k,k,k) 



c      1 is x, 3  is z, 2 is y

      IF (dim .eq. 1) then

         do i=1,k
         do j=1,k
         do l=1,k
              i1=i+1
              i2=i-1
             if(i.eq.1)  i2=k
             if(i.eq.k)  i1=1
c        this is second order central difference
        if (i .eq. 1) then 
         gra(i,j,l)=(phipp(i+1,j,l)-phipp(i,j,l))/(r_mesh)
        endif
        if (i .eq. k) then 
         gra(i,j,l)=(phipp(i,j,l)-phipp(i-1,j,l))/(r_mesh)
        endif
         if (i.ne.1.OR.i.ne.k) then
         gra(i,j,l)=(phipp(i1,j,l)-phipp(i2,j,l))/(2.*r_mesh)
         end if

         end do
         end do
         end do



      ELSE IF (dim .eq. 2) then

         do i=1,k
         do j=1,k
         do l=1,k
              j1=j+1
              j2=j-1
             if(j.eq.1)  j2=k
             if(j.eq.k)  j1=1
         gra(i,j,l)=(phipp(i,j1,l)-phipp(i,j2,l))/(2.*r_mesh)
         end do
         end do
         end do



      ELSE IF (dim .eq. 3) then


         do j=1,k
         do i=1,k
         do l=1,k
              l1=l+1
              l2=l-1
             if(l.eq.1)  l2=k
             if(l.eq.k)  l1=1
         gra(i,j,l)=(phipp(i,j,l1)-phipp(i,j,l2))/(2.*r_mesh)
         end do
         end do
         end do



      ELSE
         print *, 'Wrong Input of "dim"'
         STOP
      ENDIF


      RETURN
      end

      SUBROUTINE SGRADIENT(phipp,k,dim,gra,r_mesh)

      integer i,j,l,k,j1,j2,j3,j4, dim
      real r_mesh

      real phipp(k,k,k),gra(k,k,k)



c      1 is x, 3  is z, 2 is y

      IF (dim .eq. 1) then

         do i=2,k-1
         do j=1,k
         do l=1,k
              i1=i+1
              i2=i-1
c        this is second order central difference
         gra(i,j,l)=(phipp(i1,j,l)-phipp(i2,j,l))/(2.*r_mesh)
         end do
         end do
         end do


      ELSE IF (dim .eq. 2) then

         do i=2,k-1
         do j=1,k
         do l=1,k
              j1=j+1
              j2=j-1
             if(j.eq.1)  j2=k
             if(j.eq.k)  j1=1
         gra(i,j,l)=(phipp(i,j1,l)-phipp(i,j2,l))/(2.*r_mesh)
         end do
         end do
         end do

      ELSE IF (dim .eq. 3) then

         do i=2,k-1
         do j=1,k
         do l=1,k
              l1=l+1
              l2=l-1
             if(l.eq.1)  j2=k
             if(l.eq.k)  j1=1
         gra(i,j,l)=(phipp(i,j,l1)-phipp(i,j,l2))/(2.*r_mesh)
         end do
         end do
         end do


      ELSE
         print *, 'Wrong Input of "dim"'
         STOP
      ENDIF


      RETURN
      end


       SUBROUTINE FOURN(DATA,NN,NDIM,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION NN(NDIM),DATA(4*128*128)
      NTOT=1
      DO 11 IDIM=1,NDIM
        NTOT=NTOT*NN(IDIM)
11    CONTINUE
      NPREV=1
      DO 18 IDIM=1,NDIM
        N=NN(IDIM)
        NREM=NTOT/(N*NPREV)
        IP1=2*NPREV
        IP2=IP1*N
        IP3=IP2*NREM
        I2REV=1
        DO 14 I2=1,IP2,IP1
          IF(I2.LT.I2REV)THEN
            DO 13 I1=I2,I2+IP1-2,2
              DO 12 I3=I1,IP3,IP2
                I3REV=I2REV+I3-I2
                TEMPR=DATA(I3)
                TEMPI=DATA(I3+1)
                DATA(I3)=DATA(I3REV)
                DATA(I3+1)=DATA(I3REV+1)
                DATA(I3REV)=TEMPR
                DATA(I3REV+1)=TEMPI
12            CONTINUE
13          CONTINUE
          ENDIF
          IBIT=IP2/2
1         IF ((IBIT.GE.IP1).AND.(I2REV.GT.IBIT)) THEN
            I2REV=I2REV-IBIT
            IBIT=IBIT/2
          GO TO 1
            ENDIF
          I2REV=I2REV+IBIT
14      CONTINUE
        IFP1=IP1
2       IF(IFP1.LT.IP2)THEN
          IFP2=2*IFP1
          THETA=ISIGN*6.28318530717959D0/(IFP2/IP1)
          WPR=-2.D0*DSIN(0.5D0*THETA)**2
          WPI=SIN(THETA)
          WR=1.D0
          WI=0.D0
          DO 17 I3=1,IFP1,IP1
            DO 16 I1=I3,I3+IP1-2,2
              DO 15 I2=I1,IP3,IFP2
                K1=I2
                K2=K1+IFP1
                TEMPR=SNGL(WR)*DATA(K2)-SNGL(WI)*DATA(K2+1)
                TEMPI=SNGL(WR)*DATA(K2+1)+SNGL(WI)*DATA(K2)
                DATA(K2)=DATA(K1)-TEMPR
                DATA(K2+1)=DATA(K1+1)-TEMPI
                DATA(K1)=DATA(K1)+TEMPR
                DATA(K1+1)=DATA(K1+1)+TEMPI
15            CONTINUE
16          CONTINUE
            WTEMP=WR
            WR=WR*WPR-WI*WPI+WR
            WI=WI*WPR+WTEMP*WPI+WI
17        CONTINUE
          IFP1=IFP2
        GO TO 2
        ENDIF
        NPREV=N*NPREV
18    CONTINUE
      RETURN
      END



