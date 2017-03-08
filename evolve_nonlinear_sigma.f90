!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden, CITA, University of Toronto
!>
!>@brief
!> A simple program to evolve homogeneous trajectories of for nonlinear sigma models
!> of multi-scalar field dynamics.
!> Uses a Gauss-Legendre integrator (whose order can be selected at compile time)
!> Also evolves a series of nearby trajectories in order to estimate the local value of the
!> Hessian for the purposes of density calculations and caustic finding.
!> In a future iteration, this will be supplemented with solving the Geodesic Deviation Equations, with the neighbouring trajectory codes converted to bred vector solutions
!>
!>
!>@todo
!>@arg Debug better (in particular, off-diagonal metrics)
!>@arg Modularise so it's more readable (but probably slower).
!>     Also makes dropping new equations of motion in simpler.
!>@arg Improve file output to use binary (preferably using a higher level library such as FITS, NetCDF or HDF5) instead of ASCII output
!
! Compile with:
!   gfortran -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -O3 -cpp -o trajectory evolve_nonlinear_sigma.f90
! Then run with
!   ./trajectory
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define CONFORMAL 1
program scan_hom
  implicit none

  integer :: i,j, ii, jj
  real*8 :: chi
  integer, parameter :: nstep = 2**12
  integer, parameter :: stepsize=2**3
  real*8, parameter :: tstep = 0.01 !7.5/640./4.

! Pack the phase space variables
  integer, parameter :: dl = kind(1.d0)
  integer, parameter :: nfld=3, ntraj=2
  integer, parameter :: nvar = 2*(nfld+1)

  real*8, dimension(nvar, ntraj) :: yvec
  real*8, dimension(nvar) :: delvec
! Preprocessor macros for identifying individual field components
#define FLD yvec(1:nfld,l)
#define FLDP yvec(nfld+1:2*nfld,l)
#define AVAR yvec(2*nfld+1,l)
#define PIA yvec(2*(nfld+1),l)

  real*8, parameter :: lam=1., g2=1.875, g2_2=2., zeta=2.
  real*8, parameter :: zetafac = zeta*(1.+6.*zeta)

! Parameters controlling the rescaling, etc.
  integer :: k
  real*8 :: len
  real*8, parameter :: rescale = 0.1
  real*8, parameter :: eps = 1.e-3
  real*8, dimension(1:ntraj) :: norm=1.
  real*8 :: rho
  real*8 :: chi0, lnchi

! Parameters for computing floquet exponent
  real*8, parameter :: period = 7.4162987092054875
  real*8, parameter :: floquet = 0.2376 !0.235895932659  ! this is a function of g2 above
  real*8, parameter :: growfactor = 1./period/floquet
  integer :: pstep

  open(unit=99,file='bg_traj.out')
  open(unit=98, file='pert_traj.out')
  chi0 = 1.e-7  
  print*, chi0

  do i=0,500
     ! log sample the chi direction
     lnchi = log(chi0) + dble(i)*0.002*log(5.)
     chi = exp(lnchi)
     ! normalize lnchi for output
!     lnchi = (lnchi - log(2.3393837654714997732962993666073) ) * growfactor
     call init(chi)
     lnchi = (log(yvec(2,1)) - log(yvec(1,1)) ) * growfactor
     norm = 1.

     write(99,*) lnchi, 0., yvec(:,1), evalh(yvec(:,1)), modelv(yvec(1:nfld,1))
     delvec = yvec(:,2)-yvec(:,1)
     len = sum(delvec(:)**2)**0.5
     write(98,*) lnchi, 0., delvec(:)/len

     do j=1,nstep
        do ii=1, ntraj
        do jj=1, stepsize
           call gl10( yvec(:,ii), tstep )
        enddo
        enddo
           
        write(99,*) lnchi, j*stepsize*tstep, yvec(:,1), evalh(yvec(:,1)), modelv(yvec(1:nfld,1))
        delvec = yvec(:,2)-yvec(:,1)
        len = sum(delvec(:)**2)**0.5
        write(98,*) lnchi, j*stepsize*tstep, delvec(:)/len

        ! check vectors for rescaling
        do k=2,ntraj
           delvec = yvec(:,2)-yvec(:,1)
           len = sum( abs(delvec(:)) )

           if (len .gt. eps) then
              delvec = rescale*delvec
              norm(k) = norm(k)*rescale
              yvec(:,k) = yvec(:,1)+delvec(:)
           endif
        enddo
     enddo

     write(99,*)
     write(98,*)

  enddo

  contains
#define SCL yscl(l)
#define SCLP ysclp(l)
!
! 10th order Gauss-Legendre integrator
! Used for pieces of Hamiltonian without exact solutions
! To do : implement entire evolution using this
! To do : pack into smaller vector as required (just chop a, so probably not worth it)
!
    subroutine gl10( y, dt )
      real*8 :: y(nvar)
      real*8 :: dt

      integer, parameter :: s = 5
      integer, parameter :: niter = 8
      real*8 :: g(nvar,s)

      ! Butcher tableau for 10th order Gauss-Legendre method
      real*8, parameter :: a(s,s) = reshape( (/ &
           0.5923172126404727187856601017997934066D-1, -1.9570364359076037492643214050884060018D-2, &
           1.1254400818642955552716244215090748773D-2, -0.5593793660812184876817721964475928216D-2, &
           1.5881129678659985393652424705934162371D-3,  1.2815100567004528349616684832951382219D-1, &
           1.1965716762484161701032287870890954823D-1, -2.4592114619642200389318251686004016630D-2, &
           1.0318280670683357408953945056355839486D-2, -2.7689943987696030442826307588795957613D-3, &
           1.1377628800422460252874127381536557686D-1,  2.6000465168064151859240589518757397939D-1, &
           1.4222222222222222222222222222222222222D-1, -2.0690316430958284571760137769754882933D-2, &
           4.6871545238699412283907465445931044619D-3,  1.2123243692686414680141465111883827708D-1, &
           2.2899605457899987661169181236146325697D-1,  3.0903655906408664483376269613044846112D-1, &
           1.1965716762484161701032287870890954823D-1, -0.9687563141950739739034827969555140871D-2, &
           1.1687532956022854521776677788936526508D-1,  2.4490812891049541889746347938229502468D-1, &
           2.7319004362580148889172820022935369566D-1,  2.5888469960875927151328897146870315648D-1, &
           0.5923172126404727187856601017997934066D-1 /) , [s,s])
      real, parameter :: b(s) = (/ &
           1.1846344252809454375713202035995868132D-1,  2.3931433524968323402064575741781909646D-1, &
           2.8444444444444444444444444444444444444D-1,  2.3931433524968323402064575741781909646D-1, &
           1.1846344252809454375713202035995868132D-1 /)

      integer :: i,k

      g = 0.
      do k=1,niter
         g = matmul(g,a)
         do i=1,s
            call derivs( y+g(:,i)*dt, g(:,i) )
         enddo
      enddo
      y = y + matmul(g,b)*dt

    end subroutine gl10
!
! Evolution vector in phase space for GL integration
!
    subroutine derivs(yc, yp)
      real*8, intent(in) :: yc(nvar)
      real*8, intent(out) :: yp(nvar)

      real*8 :: dgi(nfld,nfld,nfld)
      real*8 :: gi(nfld,nfld), metric(nfld,nfld), vtmp(nfld,nfld)
      real*8 :: vprime(nfld)
      real*8 :: KE2, acur
      integer :: i

      gi = ginv(yc(1:nfld))
      acur = yc(2*nfld+1)
      dgi = gprime(yc(1:nfld))
      vprime = modeldv(yc(1:nfld))

      do i=1,nfld
         metric = dgi(1:nfld,1:nfld,i)
         vtmp(:,i) = matmul(metric,yc(nfld+1:2*nfld))
      enddo
      KE2 = sum( yc(nfld+1:2*nfld)* matmul(gi,yc(nfld+1:2*nfld)) )
      
      yp(1:nfld) = matmul( gi,yc(nfld+1:2*nfld) ) / acur**2
      do i=1,nfld
         yp(nfld+i) = -0.5*sum( yc(nfld+1:2*nfld)*vtmp(1:nfld,i)  )/acur**2 - acur**4*vprime(i)
      enddo
      yp(2*nfld+1) = -yc(2*nfld+2) / 6.
      yp(2*nfld+2) = KE2 / acur**3 - 4.*acur**3*modelv(yc(1:nfld))
#ifdef COSMIC

#endif
#ifdef MINKOWSKI

#endif
    end subroutine derivs

    function ginv(f)
      real*8, dimension(nfld, nfld) :: ginv
      real*8, dimension(nfld) :: f
      real*8 :: fac, ftmp
      
      fac = 1. + zeta*f(1)**2
      ftmp = 1. + zetafac*f(1)**2

      ginv = 0.
      ginv(1,1) = fac**2/ftmp
      ginv(2,2) = fac
      ginv(3,3) = fac
!      ginv(1,1) = 1.
!      ginv(2,2) = 1.
    end function ginv

    function gprime(f)
      real*8, dimension(nfld, nfld, nfld) :: gprime
      real*8, dimension(nfld) :: f
      real*8 :: fac, fac2

      fac = 1.+zeta*f(1)**2
      fac2 = 1.+zetafac*f(1)**2

      gprime = 0.
      gprime(1,1,1) = 4.*zeta*f(1)*fac / fac2 - 2.*zetafac*f(1)*fac**2/fac2**2
      gprime(2,2,1) = 2.*zeta*f(1)
      gprime(3,3,1) = 2.*zeta*f(1)
    end function gprime

    function modelv(f)
      real*8 :: modelv
      real*8, dimension(1:nfld) :: f
      real*8 :: denom
      
!      denom = (1.+zeta*f(1)**2)**2
!      modelv = 0.5*lam*f(1)**4 + g2*f(1)**2*f(2)**2
!      modelv = modelv / denom
      modelv = 0.25*lam*f(1)**4 + 0.5*g2*f(1)**2*f(2)**2 + 0.5*g2_2*f(1)**2*f(3)**2
      return
    end function modelv

    function modeldv(f)
      real*8, dimension(1:nfld) :: modeldv
      real*8, dimension(1:nfld) :: f
      real*8 :: denom
      
!      denom = (1.+zeta*f(1)**2)
!      modeldv(1) = ( lam*f(1)**3 + g2*f(1)*f(2)**2 )/denom**2   &
!                  - zeta*f(1)*( lam*f(1)**4 + 2.*g2*f(1)**2*f(2)**2   ) /denom**3
!      modeldv(2) =  g2*f(1)**2*f(2) / denom**2

      modeldv(1) = lam*f(1)**3+g2*f(1)*f(2)**2+g2_2*f(1)*f(3)**2
      modeldv(2) = g2*f(1)**2*f(2)
      modeldv(3) = g2_2*f(1)**2*f(3)
      return
    end function modeldv

    subroutine init(chirat)
      real*8 :: chirat
      real*8 :: rhotmp, ftmp(1:nfld)
      real*8 :: H0=1.9348974397391251388968698880012

      real*8, dimension(nfld,nfld) :: gi
      integer :: l,m

      if (abs(zeta).lt.1.e-3) then
         yvec(1,:) = 2.3393837654714997732962993666073
         yvec(nfld+1,:) = -2.7363582010758065274616992909302
      else
         yvec(1,:) = (0.5*(sqrt((1+24*zeta)*(1.+8.*zeta))-1.)/(1.+6.*zeta)/abs(zeta))**0.5
         yvec(nfld+1,:) = 0.
      endif
      yvec(2,:) = chirat * yvec(1,:)
      yvec(3,:) = chirat * yvec(1,:)
      ftmp = modeldv( yvec(1:nfld,1) )
      yvec(nfld+2,:) = 0.

      gi = ginv(yvec(1:nfld,1))

! Fix rhotmp in here (only works since a=1 initially)
      do l=1,10
         rhotmp = scalar_rho(yvec(1:nfld,1),yvec(nfld+1:2*nfld,1),1.)
!0.5*yvec(3,1)**2*gi(1,1) + 0.5*yvec(4,1)**2*gi(2,2) + modelv( yvec(1:nfld,1) )
         H0 = rhotmp**0.5/3.**0.5
         do m=1,nfld
            yvec(nfld+m,:) = -ftmp(m) / 3. / H0
         enddo
      enddo

!fix rhotmp here
      rhotmp = scalar_rho(yvec(1:nfld,1),yvec(nfld+1:2*nfld,1),1.)
!0.5*yvec(nfld+1,1)**2*gi(1,1) + 0.5*yvec(nfld+2,1)**2*gi(2,2) + modelv( yvec(1:nfld,1) )
      H0 = rhotmp**0.5 / 3.**0.5
 
      yvec(2,2) = chi0 + 1.e-9

      yvec(2*nfld+1,:) = 1.
#ifdef CONFORMAL
      yvec(2*nfld+2,:) = -6.*H0*yvec(2*nfld+1,:)  ! H0 is a^{-1}da/dtau, tau conformal time?
#else
      yvec(2*nfld+2,:) = -4.*H0*yvec(2*nfld+1,:)
#endif
    end subroutine init

    function scalar_rho(fld, fldp, acur)
      real*8 :: scalar_rho
      real*8, intent(in), dimension(1:nfld) :: fld, fldp
      real*8, intent(in) :: acur

      real*8, dimension(nfld,nfld) :: gi
      real*8, dimension(nfld) :: vtmp
      real*8 :: KE

      gi = ginv(fld)
      vtmp = matmul(gi,fldp)
      KE = sum(fldp*vtmp)
#ifdef CONFORMAL
      scalar_rho = 0.5*KE/acur**6 + modelv(fld)
#endif
    end function scalar_rho

! Check Hamiltonian constraint
    function evalh(fv)
      real*8 :: evalh
      real*8, dimension(2*(nfld+1)) :: fv

      real*8, dimension(nfld,nfld) :: gi
      gi = ginv(fv(1:nfld))
      
#ifdef CONFORMAL
      evalh = fv(2*nfld+2)**2/12.   &
              - fv(2*nfld+1)**4*scalar_rho( fv(1:nfld),fv(nfld+1:2*nfld),fv(2*nfld+1) )
#else
      evalh = 
#endif


    end function evalh
    
  end program scan_hom
