!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++ calculate the IC-gamma from the e+/e- in photon field
!++ Eelec(GeV), Felec(cm^-3 GeV^-1), Ebk(eV), Fbk(cm^-3 eV^-1)
!++ Eg in GeV, Fg in cm^-3 s^-1 GeV^-1, line of sight integral
!++ divided by 4pi gives the flux with unit cm^-2 s^-1 sr^-1 GeV^-1
      subroutine GAMMA_IC(Eelec,Felec,nelec,Ebk,Fbk,nbk,Eg,Fg,ng)
      implicit real*8 (a-z)
      integer nelec,ng,nbk,ielec,ig,ibk
      real*8 Eelec(nelec),Felec(nelec),Eg(ng),Fg(ng),Ebk(nbk),Fbk(nbk)
      ME=0.511d-3
      if (abs(Eg(ng)-Eg(1))<1d-50) then
         Eg_min=1.d-3               ! GeV
         Eg_max=1.d6                ! GeV
         gfactor=exp(log(Eg_max/Eg_min)/(ng-1))
         do ig=1,ng
            Eg(ig)=Eg_min*gfactor**(ig-1) ! GeV
         enddo
      endif
      if (nbk==1) then
         bkfactor=exp(1/Ebk(nbk))
      else
         bkfactor=Ebk(nbk)/Ebk(nbk-1)
      endif
      if (nelec==1) then
         elecfactor=exp(1/Eelec(nelec))
      else
         elecfactor=Eelec(nelec)/Eelec(nelec-1)
      endif
      do ig=1,ng
         gamma_emiss=0.d0
         do ibk=1,nbk
            E_tar=Ebk(ibk)*1d-9/ME
            do ielec=1,nelec
               e_lorentz=Eelec(ielec)/ME+1
               gamma_emiss=gamma_emiss+Fbk(ibk)*Felec(ielec)*
!     &FJONES(e_lorentz,E_tar,Eg(ig)/ME)*Ebk(ibk)*Eelec(ielec)
     &F_IC(e_lorentz,E_tar,Eg(ig)/ME)*Ebk(ibk)*Eelec(ielec)
!++ FJONES in cm^3 s^-1 GeV^-1
            enddo
         enddo
         Fg(ig)=gamma_emiss*log(bkfactor)*log(elecfactor)
      enddo
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++ calculate the blackbody background radiation field with temperature
!++ T [K] and energy density u [eV cm^-3]
!++ Ebk(nbk) in eV, Fbk(nbk) in cm^-3 eV-1
      subroutine bkg_rad(Ebk,Fbk,nbk,T,u)
      implicit real*8 (a-z)
      integer nbk,ibk
      real*8 Ebk(nbk),Fbk(nbk)
      if (abs(Ebk(2)-Ebk(1))<1d-50) then
         Ebk_min=2.d-6*T    ! eV
         Ebk_max=2.d-2*T    ! eV
      else
         Ebk_min=Ebk(1)
         Ebk_max=Ebk(nbk)
      endif
      bkfactor=exp(log(Ebk_max/Ebk_min)/(nbk-1))
      do ibk=1,nbk
         eps=Ebk_min*bkfactor**(ibk-1) ! eV
         Fbk(ibk)=eps*eps/(exp(eps/(8.62e-5*T))-1)*u*2.7891e15/T**4
         Ebk(ibk)=eps
      enddo
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++ calculate the integral flux (cm^-2 s^-1)
      subroutine inte_flux(Eg,Fg,ng,ptflag)
      implicit real*8 (a-z)
      integer ng,ig,ptflag
      real*8 Eg(ng),Fg(ng)
      s=0
      temp=Fg(ng)
      do ig=ng-1,1,-1
         deltaE=Eg(ig+1)-Eg(ig)
         s=s+(Fg(ig)+temp)/2*deltaE
         temp=Fg(ig)
         Fg(ig)=s
         if (ptflag.ne.0) write(*,*) Eg(ig),Fg(ig)
      enddo
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 function FJONES(gam,E1,E4)
c***********************************************************************
c                            *** I.Moskalenko, version of 26.03.1998 ***
c calculation of the spectrum from INVERSE COMPTON scattering
c of ISOTROPIC background photons off isotropic electrons
c    INPUT:
c gam is the electron Lorentz factor;
c E1 is the energy of a background photon (in electron rest mass units mc^2);
c E4 is the energy of scattered photon (in electron rest mass units mc^2);
c    OUTPUT:
c FJONES - the differential energy spectrum of gamma-rays
c+++ photon cm^3 s^-1 GeV^-1, yuanqiang 20081119 +++
c    (phot ccm/sec/energy) per one electron in ccm divided by Pi*r0^2*c;
c    REFERENCES:
c F.C.Jones 1968, Phys.Rev. 167, 1159
c***********************************************************************
      implicit real*8 (a-h,o-z)

      FJONES = 0.d0
      if(E4/gam .gt. 4.d0*E1*gam/(1.d0+4.d0*E1*gam)) return
      q = E4/(4.d0*E1*gam**2)*(1.d0+4.d0*E1*gam)
      if(4.d0*E1*gam .lt. 1.d10) q = E4/(4.d0*E1*gam*(gam-E4))
      FJONES = 2.d0*q*dlog(q)+(1.d0+2.d0*q)*(1.d0-q)
     #   +(1.d0-q)/2.d0*(4.d0*E1*gam*q)**2/(1.d0+4.d0*E1*gam*q)
      FJONES = FJONES*2.d0/E1/gam**2  ! in units Pi*r0^2*c  eq.(12)
      FJONES = FJONES*1.46458184E-11  ! 1.46458184E-11=pi*r0*r0*c/m_e

      return
      end

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++ IC radiation spectrum for single electron-photon interaction
!++ Eq.(18) of Fang & Zhang, arXiv:0711.4173
!++ e_lorentz: lorentz factor of electron
!++ E_tar: photon energy in unit of electron mass
!++ Eg: scattering gamma energy in unit of electron mass
!++ F_IC in cm^3 s^-1 GeV^-1
      function F_IC(e_lorentz,E_tar,Eg)
      implicit real*8 (a-z)
      SIGMAT=6.652d-25
      Gam=4*E_tar*e_lorentz
      E1=Eg/e_lorentz
      q=E1/(Gam*(1-E1))
      F_IC=0.
      if (q<1.and.q>1/(4*e_lorentz*e_lorentz)) then
         F_IC=3*SIGMAT*3d10/(4*e_lorentz*e_lorentz*E_tar*0.511d-3)*
     &(2*q*log(q)+(1+2*q)*(1-q)+(Gam*q)**2*(1-q)/(2*(1+Gam*q)))
      endif
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
