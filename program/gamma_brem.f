!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++ calculate the bremsstrahlung-gamma from the e+/e- in ISM field
!++ Eelec(GeV), Felec(cm^-3 GeV^-1), Eg in GeV, Fg in cm^-3 s^-1 GeV^-1, 
!++ line of sight integral divided by 4pi gives the flux with unit 
!++ cm^-2 s^-1 sr^-1 GeV^-1
!++ use bremss_spec.f, inter.f
      subroutine GAMMA_BREM(Eelec,Felec,nelec,NH,Eg,Fg,ng)
      implicit real*8 (a-z)
      integer nelec,ng,ielec,ig
      real*8 Eelec(nelec),Felec(nelec),Eg(ng),Fg(ng)
      ME=0.511d-3
      if (abs(Eg(ng)-Eg(1))<1d-50) then
         Eg_min=1.d-3               ! GeV
         Eg_max=1.d6                ! GeV
         gfactor=exp(log(Eg_max/Eg_min)/(ng-1))
         do ig=1,ng
            Eg(ig)=Eg_min*gfactor**(ig-1) ! GeV
         enddo
      endif
      elecfactor=Eelec(nelec)/Eelec(nelec-1)
      do ig=1,ng
         Eg_MeV=Eg(ig)*1d3
         HI=0.d0
         He=0.d0
         do ielec=1,nelec
            Ee_MeV=(Eelec(ielec)+ME)*1d3
            ve=sqrt(1-(0.511/Ee_MeV)**2)*3d10
            call bremss_spec(Eg_MeV,Ee_MeV,1,1,dSdK) ! HI
!++ dSdK in barn/MeV, times 1e-24*1e3 gives cm^2/GeV
            HI=HI+dSdK*Felec(ielec)*Eelec(ielec)*ve
            call bremss_spec(Eg_MeV,Ee_MeV,2,2,dSdK) ! He
            He=He+dSdK*Felec(ielec)*Eelec(ielec)*ve
         enddo
         Fg(ig)=NH*(HI+0.11*He)*log(elecfactor)*1d-21
      enddo
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * bremss_spec.f *                               galprop package * 4/14/2000 
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

      subroutine bremss_spec(Egam,E0,IZ1,Ne1,dSdK) 
c***********************************************************************
c          *** I.Moskalenko (MPE,Garching) *** version of 14.05.1998 ***
c PURPOSE:
c Calculation of the spectrum of e-bremsstrahlung in the fully ionized,
c hydrogen-like, or helium-like gas as well as in neutral H, He, and 
c heavier in Born approximation; valid for E0,Egam >~ 0.01 MeV and
c beta0, beta1 >> 2*pi*Af*Z, corresponds to gamma=1.001
c REFERENCES:
c [1] Blumenthal & Gould 1970,Rev.Mod.Phys.42,237 (BG70); 
c [2] Gould 1969, Phys.Rev.185,72 (G69);
c [3] Koch & Motz 1959, Rev.Mod.Phys.31,920 (KM59).
c [4] Strong A.W., Moskalenko I.V., Reimer O. 2000, ApJ 537, 763
c INPUT/OUTPUT parameters:
c Egam, E0 [MeV] are the energy of a gamma and total electron energy,
c respectively;
c IZ1 is the nucleus charge; 
c Ne1 = 0 for unshielded charge, = 1 for 1-e atoms, = 2 for 2-e atoms.
c dSdK [barn/MeV] is the differential cross section (output).
c Note:
c #1 a one-parameter Hylleraas function is used for helium-like ions,
c    while Hartree-Fock functions are used for heutral He atoms;
c #2 a contribution of free electrons could be estimated by calling
c    the subroutine with Z=1, Ne=0 and multiplying dSdk by ionization level.
c++ using inter.f, pp_meson.f  yuanqiang 20090312
c***********************************************************************
      implicit real*8 (a-h,m,o-z)
      real*8 DD1(11)
      common/delta/ delta, IZ, Ne
     2      /Hartree/ DD(11),PH1(11),PH2(11)  ! DD = delta/(2.d0*Af)
      external FPHI1, FPHI2

      dSdK = 0.d0
      mc2 = 0.511   ! MeV, electron rest mass
      E1 = E0-Egam  ! total energy of an electron after scattering
      if(E0 .le. 1.02*mc2 .or. E1 .le. mc2) return
         
      IZ = IZ1
      Ne = Ne1
      Pi = 3.14159265359d0
      Af = 7.2974d-3      ! fine structure constant
      R02 = 2.81794d-1**2 ! classical electron radius square (barn)
      EB1 = 0.07 ! MeV,kinetic energy; the upper bound for NonRelat.approx.
      EB2 = 2.  ! MeV, -"-; the lower bound for HE approx. (nonscreened)
      delta = Egam*mc2/(2.d0*E0*E1)
      Ef = 0.
      FS = 0.d0
      FS90 = 0.90
      T0 = E0-mc2   ! Ekin
      P0 = mc2*dsqrt((E0/mc2)**2-1.d0) ! initial electron momentum
      P1 = mc2*dsqrt((E1/mc2)**2-1.d0) ! final -"-
      beta0 = P0/E0
      beta1 = P1/E1

c### Fano-Sauter high-frequency limit (Egam ->E0-mc2, KM59,p.933) ###
      if(Egam/T0 .ge. FS90)
     2   FS =4.*Pi*IZ**3*R02*Af**2/T0*beta0*E0/mc2/(E0/mc2-1.d0)**2 
     3      *(4./3. +E0/mc2*(E0/mc2-2.d0)/(E0/mc2+1.d0) *(1.d0
     4      -dlog(1.d0+2.d0*beta0/(1.d0-beta0))*(mc2/E0)**2/2.d0/beta0))

c## Elwert factor (Ef) is valid if   Z*Af*(1./beta1-1./beta0)<<1 ###
      if(T0 .lt. EB2 .and. IZ*Af*(1./beta1-1./beta0) .le. 0.5)
     2   Ef = beta0/beta1 *(1.d0-dexp(-2.d0*Pi*Af*IZ/beta0))
     3      /(1.d0-dexp(-2.d0*Pi*Af*IZ/beta1))

c### NONRELATIVISTIC ENERGIES ( T0 < 0.07 MeV) ###
c# Born approx.(KM59,p.925,3BNa) with Elwert factor (p.931)
c# Born approximation: 2*Pi*Z*Af/beta{0,1}<<1; 
c# nonscreened: 2*delta/[Af*Z^(1/3)] >> 1;
      if(T0 .lt. EB1) then
         A = IZ*IZ                    ! no electron contribution
         dSdK = Ef*R02*Af*A*16./3. /Egam 
     2      *(mc2/P0)**2*dlog(1.d0+2.d0*P1/(P0-P1))
         if(dSdK .lt. FS) dSdK = (FS+dSdK)/2.   !# 1/2 Fano-Sauter limit
         RETURN
      endif

c### LOW ENEGRY (T0 < 2. MeV) OR delta/(2.*Af*Z) > 4. ###
c# Born approximation (KM59,p.925,3BN): 2*Pi*Z*Af/beta{0,1}<<1;
c# nonscreened: 2*delta/Af >> 1
      if(T0 .lt. EB2 .or. delta/(2.*Af*IZ) .ge. 4.) then
         if (E0>mc2*1d3) then
            eps0 = dlog(1+4*P0*P0/mc2/mc2)
         else
            eps0 = dlog(1.d0+2.d0*P0/(E0-P0))  !=dlog((E0+P0)/(E0-P0))
         endif
         if (E1>mc2*1d3) then
            eps1 = dlog(1+4*P1*P1/mc2/mc2)
         else
            eps1 = dlog(1.d0+2.d0*P1/(E1-P1))  !=dlog((E1+P1)/(E1-P1))
         endif
         AL = 2.d0*dlog((E0*E1+P0*P1-mc2**2)/(Egam*mc2))
         T1 = 2.d0*E0*E1*((P0/P1)**2+1.d0)/P0**2
         T2 = mc2**2*(eps0*E1/P0**3+eps1*E0/P1**3-eps0*eps1/P0/P1)
         T10 = 8.d0/3.*(E0*E1)/(P0*P1)
         T11 = Egam**2*((E0/P0*E1/P1)**2+1.d0)/(P0*P1)
         T12 = mc2**2 *Egam/(2.d0*P0*P1) *(
     2      eps0*(E0/P0*E1/P0+1.d0)/P0 -eps1*(E0/P1*E1/P1+1.d0)/P1
     3      +2.d0*Egam*E0/P0*E1/P1/(P0*P1) )
         A = IZ*IZ +Ne            ! correction for atomic electrons
         if(T0 .lt. EB2)
     2      A =( IZ*IZ +Ne*(1.d0-dexp(-(T0-EB1)/9./EB1)) ) ! screening
     3         *(1.-0.3*dexp(-Egam/0.033)) *Ef  ! left & right connection
         dSdK = R02*Af*A/Egam*P1/P0*(4.d0/3.-T1+T2+AL*(T10+T11+T12))
         if(dSdK .lt. FS) dSdK = (FS+dSdK)/2.   !# 1/2 Fano-Sauter limit
         RETURN                    
      endif
      
c### HIGH  ENERGY (T0 > 2. MeV) AND delta/(2.*Af*Z) < 4. ###

c# Schiff approx. for neutral atoms heavier then He (KM59,p.925,3BSe) 
      if(IZ .gt. 2 .and. IZ .eq. Ne) then
         b = IZ**(1.d0/3.d0)/111.d0/delta
         OM = 1.d0/delta**2/(1.d0+b*b)
         dSdK = 2.d0*IZ*IZ*R02*Af/Egam *( 
     2      (1.d0+(E1/E0)**2-2.d0/3.*E1/E0)
     3         *(dlog(OM)+1.d0-2.d0/b*datan(b))
     4      +E1/E0*(2.d0/b/b*dlog(1.d0+b*b)
     5         +4./3.*(2.d0-b*b)/b**3*datan(b)-8.d0/3./b/b+2.d0/9.) )
         RETURN
      endif

c# arbitrary screening (G69, BG70; KM59,p.925,3BSb)
      A = IZ*IZ +Ne            !correction for atomic electrons
      phiu = -dlog(delta)-1.d0/2.d0
c# Hartree-Fock approximation for neutral He atoms
      if(IZ .eq. 2 .and. IZ .eq. Ne) then
         DD1(1) = dlog(1.d-3)
         do i=2,11
            DD1(i) = dlog(DD(i))
         enddo
         if(delta/(2.d0*Af) .ge. DD(2)) then
            if(delta/(2.d0*Af) .le. DD(11)) then
               call DINTER(DD1,PH1,11,dlog(delta/(2.d0*Af)),phi1)
               call DINTER(DD1,PH2,11,dlog(delta/(2.d0*Af)),phi2)
            else
               phi1 = 4.*A *phiu  ! asymptotics
               phi2 = 4.*A *phiu
            endif
         else
            phi1 = PH1(1)
            phi2 = PH2(1)
         endif
         phi1 = phi1/4.
         phi2 = phi2/4.
      else
         if(Ne .eq. 0) then  !# UNSHIELDED charge
            phi1 = A*phiu
            phi2 = A*phiu
         else
c# H-like atoms & Hylleraas-1 approximation for He-like atoms
            if(delta/(2.d0*Af) .le. 1.d-4) delta = 1.d-4*(2.d0*Af)
            call SIM1(1.d0,delta,1.d-3,1.d-5,1.d-4,FPHI1,p1)
            call SIM1(1.d0,delta,1.d-3,1.d-5,1.d-4,FPHI2,p2)
            phi1 = 2.*IZ*(-p1       +1.+(Ne-1.)/IZ)  +(IZ-Ne)**2*phiu
            phi2 = 2.*IZ*(-p2+5./6.*(1.+(Ne-1.)/IZ) )+(IZ-Ne)**2*phiu
         endif         
      endif
      if(phi1 .gt. A*phiu) phi1 = A *phiu  ! asymptotics
      if(phi2 .gt. A*phiu) phi2 = A *phiu

      dSdK =4.*R02*Af/Egam *((1.d0+(E1/E0)**2)*phi1 -2./3.*E1/E0*phi2)
      if(dSdK .lt. FS) dSdK = (FS+dSdK)/2.      !# 1/2 Fano-Sauter limit
      return
      end

      real*8 function FPHI1(Q) 
c***********************************************************************
c used for calculation of PHI1 function (G69,p.74,76)
c          *** I.Moskalenko (MPE,Garching) *** version of 23.09.1997 ***
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common/delta/ delta, IZ, Ne

      Af = 7.2974d-3          ! fine structure constant
      if(Ne .eq. 1) then      ! one-electron atoms
         AZ = 1.d0/(2.d0*Af*IZ)
         FZ = 1.d0/(1.d0+(AZ*Q)**2)**2 ! form-factor for one-electron atom
         FPHI1 = (Q-delta)**2/Q**3 *(1.d0-FZ)
         return
      endif                   ! two-electron atoms (Hylleraas-1 function)
      AZ = 1.d0/(2.d0*Af*(IZ-5.d0/16.d0))
      FZ = 1.d0/(1.d0+(AZ*Q)**2)**2 ! form-factor for two-electron atom
      FPHI1 = (Q-delta)**2/Q**3 *(2.*(1.d0-FZ) -(1.d0-FZ**2)/IZ)
      return
      end

      real*8 function FPHI2(Q) 
c***********************************************************************
c used for calculation of PHI2 function (G69,p.74,76)
c          *** I.Moskalenko (MPE,Garching) *** version of 23.09.1997 ***
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common/delta/ delta, IZ, Ne

      Af = 7.2974d-3          ! fine structure constant
      if(Ne .eq. 1) then      ! one-electron atoms
         AZ = 1.d0/(2.d0*Af*IZ)
         FZ = 1.d0/(1.d0+(AZ*Q)**2)**2 ! atomic form-factor
         FPHI2=(Q**3+(3.d0-6.d0*dlog(Q/delta))*delta**2*Q 
     2      -4.d0*delta**3)/Q**4 *(1.d0-FZ)
         return
      endif                   ! two-electron atoms (Hylleraas-1 function)
      AZ = 1.d0/(2.d0*Af*(IZ-5.d0/16.d0))
      FZ = 1.d0/(1.d0+(AZ*Q)**2)**2 ! form-factor for two-electron atom
      FPHI2=(Q**3+(3.d0-6.d0*dlog(Q/delta))*delta**2*Q 
     2   -4.d0*delta**3)/Q**4 *(2.*(1.d0-FZ) -(1.d0-FZ**2)/IZ)
      return
      end
c***********************************************************************

!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * inter.f *                                     galprop package * 4/14/2000 
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

      SUBROUTINE INTER(XX,FF,NN,X0,F1)
C********************************************************************
C* INTERPOLATION BY THE CUBIC POLINOM, REAL*4 version               *
C* XX,FF  ARE THE ARGUMENT AND FUNCTION ARRAYS (j=1..NN);           *
C* X0  IS THE ARGUMENT AT WHICH THE FUNCTION IS TO BE EVALUATED;    *
C* F1  ON EXIT, IS THE EVALUATED VALUE OF THE FUNCTION.             *
C********************************************************************
      IMPLICIT real (A-H,O-Z)
      IMPLICIT integer (I-N)
      DIMENSION XX(NN),FF(NN)
      IN=0
      N1=NN-1
      DO 1 M=1,N1
         IF(X0.GE.XX(M).AND.X0.LE.XX(M+1)) GOTO 2  
1           CONTINUE
2     IN=M
      IF(IN.EQ.1.OR.X0.LT.XX(1)) IN=2
      IF(IN.EQ.N1.OR.X0.GT.XX(NN)) IN=NN-2
      X1=XX(IN-1)
      X2=XX(IN)
      X3=XX(IN+1)
      X4=XX(IN+2)
      Y1=FF(IN-1)
      Y2=FF(IN)
      Y3=FF(IN+1)
      Y4=FF(IN+2)
      X01=X0-X1             
      X02=X0-X2       
      X03=X0-X3       
      X04=X0-X4       
      X12=X1-X2       
      X13=X1-X3       
      X14=X1-X4       
      X23=X2-X3       
      X24=X2-X4       
      X34=X3-X4
      X1= X02*X03*X04/X12/X13/X14
      X2=-X01*X03*X04/X12/X23/X24
      X3= X01*X02*X04/X13/X23/X34
      X4=-X01*X02*X03/X14/X24/X34
      F1= Y1*X1+Y2*X2+Y3*X3+Y4*X4
      RETURN
      END

      SUBROUTINE DINTER(XX,FF,NN,X0,F1)
C********************************************************************
C* INTERPOLATION BY CUBIC POLINOM, REAL*8 version                   *
C* XX,FF  ARE THE ARGUMENT AND FUNCTION ARRAYS (j=1..NN);           *
C* X0  IS THE ARGUMENT AT WHICH THE FUNCTION IS TO BE EVALUATED;    *
C* F1  ON EXIT, IS THE EVALUATED VALUE OF THE FUNCTION.             *
C********************************************************************
      IMPLICIT real*8  (A-H,O-Z)
      IMPLICIT integer (I-N)
      DIMENSION XX(NN),FF(NN)
      IN=0
      N1=NN-1
      DO 1 M=1,N1
         IF(X0.GE.XX(M).AND.X0.LE.XX(M+1)) GOTO 2  
1           CONTINUE
2     IN=M
      IF(IN.EQ.1.OR.X0.LT.XX(1)) IN=2
      IF(IN.EQ.N1.OR.X0.GT.XX(NN)) IN=NN-2
      X1=XX(IN-1)
      X2=XX(IN)
      X3=XX(IN+1)
      X4=XX(IN+2)
      Y1=FF(IN-1)
      Y2=FF(IN)
      Y3=FF(IN+1)
      Y4=FF(IN+2)
      X01=X0-X1             
      X02=X0-X2       
      X03=X0-X3       
      X04=X0-X4       
      X12=X1-X2       
      X13=X1-X3       
      X14=X1-X4       
      X23=X2-X3       
      X24=X2-X4       
      X34=X3-X4
      X1= X02*X03*X04/X12/X13/X14
      X2=-X01*X03*X04/X12/X23/X24
      X3= X01*X02*X04/X13/X23/X34
      X4=-X01*X02*X03/X14/X24/X34
      F1= Y1*X1+Y2*X2+Y3*X3+Y4*X4
      RETURN
      END
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
