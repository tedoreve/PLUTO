!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++ calculate the gamma, e+ or pbar flux from pp collision
!++ key=0 (gamma), +1 (e+), -1 (pbar)
!++ input: Ep(np) in GeV, Fp(np) in cm^-3 GeV^-1, den_ISM in cm^-3
!++ output: Esec(nsec) in GeV, Fsec(nsec) in cm^-3 s^-1 GeV^-1
!++ line of sight integral (or times the propagator for charged particle) 
!++ divided by 4pi gives the flux with unit cm^-2 s^-1 sr^-1 GeV^-1
!++ use antiproton.f, aa_kamae.f, pp_meson.f, nucleon_cs.f, crn6.f
      subroutine SEC_PION_PP(Ep,Fp,np,den_ISM,key,Esec,Fsec,nsec)
      implicit real*8 (a-z)
      integer np,nsec,ip,isec,key
      real*8 Ep(np),Fp(np),Esec(nsec),Fsec(nsec)
      MP=0.9382                  ! GeV
      if (abs(Esec(nsec)-Esec(1))<1d-50) then
         Esec_min=1.d-3          ! GeV
         Esec_max=1.d6           ! GeV
         secfactor=exp(log(Esec_max/Esec_min)/(nsec-1))
         do isec=1,nsec
            Esec(isec)=Esec_min*secfactor**(isec-1) ! [GeV]
         enddo
      endif
      pfactor=Ep(np)/Ep(np-1)
      do isec=1,nsec
         p_HI=0.d0
         do ip=1,np
            pp=sqrt(Ep(ip)*(2*MP+Ep(ip))) ! momentum, GeV
            vp=sqrt(1-(MP/(MP+Ep(ip)))**2)*3d10
            if (key==0) then !++ gamma-rays
!               p_HI=p_HI+PP_MESON(Esec(isec),pp,1,1,0) ! old GALPROP
               p_HI=p_HI+AA_KAMAE(Esec(isec),pp,1,1,0) ! Kamae et al.
     &*Fp(ip)*Ep(ip)*vp*1d-24 ! [GeV^-1 s^-1]
            else if (key==1) then !++ positrons
!               p_HI=p_HI+PP_MESON(Esec(isec),pp,1,1,4)
               p_HI=p_HI+AA_KAMAE(Esec(isec),pp,1,1,1)
     &*Fp(ip)*Ep(ip)*vp*1d-24 ! [GeV^-1 s^-1]
            else if (key==-1) then !++ antiprotons
               ppbar=sqrt(Esec(isec)*(2*MP+Esec(isec)))  ! momentum
               p_HI=p_HI+ANTIPROTON(2,ppbar,pp,1,1,1,1)
     &*Fp(ip)*Ep(ip)*vp*(Esec(isec)+MP)/ppbar*1d-24 ! [GeV^-1 s^-1]
!++ factor (Esec(isec)+MP)/ppbar convert dN/dp to dN/dE_tot (or dN/dE_k)
            endif
         enddo
         Fsec(isec)=den_ISM*p_HI*log(pfactor)
      enddo
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++ secondary production of AA collision
!++ NZ1 & NA1 are charge number and mass number for beam particle
!++ NZ2 & NA2 are charge number and mass number for target particle
!++ key=0 (gamma), +1 (e+), -1 (pbar)
!++ input: Ecr(ncr) in GeV, Fcr(ncr) in cm^-3 GeV^-1, den2 in cm^-3
!++ output: Esec(nsec) in GeV, Fsec(nsec) in cm^-3 s^-1 GeV^-1
      subroutine SEC_PION(Ecr,Fcr,ncr,den2,key,Esec,Fsec,nsec
     &,NZ1,NA1,NZ2,NA2)
      implicit real*8 (a-z)
      integer ncr,nsec,icr,isec,key,NZ1,NA1,NZ2,NA2
      real*8 Ecr(ncr),Fcr(ncr),Esec(nsec),Fsec(nsec)
      MP=0.9382 ! GeV
      MCR=MP*NA1
      if (abs(Esec(nsec)-Esec(1))<1d-50) then
         Esec_min=1.d-3          ! GeV
         Esec_max=1.d6           ! GeV
         secfactor=exp(log(Esec_max/Esec_min)/(nsec-1))
         do isec=1,nsec
            Esec(isec)=Esec_min*secfactor**(isec-1) ! [GeV]
         enddo
      endif
      crfactor=Ecr(ncr)/Ecr(ncr-1)
      do isec=1,nsec
         tmp=0.d0
         do icr=1,ncr
            pcr=sqrt(Ecr(icr)*(2*MCR+Ecr(icr)))
            vcr=sqrt(1-(MCR/(MCR+Ecr(icr)))**2)*3d10
            con=Fcr(icr)*Ecr(icr)*vcr*1d-24
            if (key==0) then !++ gamma-rays
               tmp=tmp+AA_KAMAE(Esec(isec),pcr,NA1,NA2,0)*con 
            else if (key==1) then !++ positrons
               tmp=tmp+AA_KAMAE(Esec(isec),pcr,NA1,NA2,1)*con 
            else if (key==-1) then !++ antiprotons
               ppbar=sqrt(Esec(isec)*(2*MP+Esec(isec)))
               tmp=tmp+ANTIPROTON(2,ppbar,pcr,NZ1,NA1,NZ2,NA2)*con
     &*(Esec(isec)+MP)/ppbar 
!++ factor (Esec(isec)+MP)/ppbar convert dN/dp to dN/dE_tot (or dN/dE_k)
            endif
         enddo
         Fsec(isec)=den2*tmp*log(crfactor)
      enddo
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function AA_KAMAE(Esec, Pp1, NA1, NA2, key1)
c INPUT parameters:
c Esec [GeV] - (kinetic) energy of secondary particle or photon;
c Pp1 [GeV/c] - beam momentum per nucleus;
c NA1 & NA2 are the atomic numbers of beam and target nuclei,
c correspondingly (NA1=NA2=1 for pp-collisions);
c key1 - choice of which kind of secondary particle
c key1= 0:  gamma rays
c key1=-1:  electrons
c key1= 1:  positrons
c OUTPUT: differential cross section in [barn/GeV]
      implicit real*8 (a-z)
      integer NA1,NA2,key1
      mp=0.9382
      Pp=Pp1/NA1
      Tp=sqrt(Pp*Pp+mp*mp)-mp
      AA_KAMAE=PP_KAMAE(Esec,Tp,key1)*(NA1**(3./8.)+NA2**(3./8.)-1.)**2
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++ The secondary production cross section from pp collisions
!++ key= 0:  gamma rays
!++ key=-1:  electrons
!++ key= 1:  positrons
!++ key=-2:  electron antineutrinos
!++ key= 2:  electron neutrinos
!++ key=-3:  muon antineutrinos
!++ key= 3:  muon neutrinos
!++ inputs: the kinetic energy of secondary particle Eksec (GeV)
!++ energy of proton Tp (GeV), flag key
!++ output: dsigma/dE (barn/GeV)
!++ Reference: Kamae et al., 2006, ApJ, 647, 692
!++ 2008-10-27
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function PP_KAMAE(Eksec,Tp,key)
      implicit real*8 (a-z)
      integer key
      PP_KAMAE=0.
      if (Eksec<1d-4.or.Tp>1d6) return
      sig_ND=F_ND(Eksec,Tp,key)*F_NDKL(Eksec,Tp,key)
      if (sig_ND<0) sig_ND=0
      sig_DIFF=F_DIFF(Eksec,Tp,key)*F_KL(Eksec,Tp)
      if (sig_DIFF<0) sig_DIFF=0
      sig_RES=F_RES(Eksec,Tp,key)*F_KL(Eksec,Tp)
      if (sig_RES<0) sig_RES=0
      PP_KAMAE=(sig_ND+sig_DIFF+sig_RES)*1d-3/Eksec
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function F_ND(Eksec,Tp,key) ! Eq.(6) of Kamae et al.
      implicit real*8 (a-z)
      integer key
      x=log10(Eksec)
      p=sqrt(Tp*Tp+2*0.9382*Tp)
      F_ND=0.0
      if (p<1.0) return ! Eq.(1)
      yr=log10(Tp)
      y=yr-3.
      if (key==0) then 
         if (Tp>0.4) then
         a0=-0.51187*(y+3.3)+7.6179*(y+3.3)**2-2.1332*(y+3.3)**3
     &+0.22184*(y+3.3)**4
         a1=-1.2592d-5+1.4439d-5*exp(-0.2936*(y+3.4))+5.9363d-5/
     &(y+4.1485)+2.2640d-6*y-3.3723d-7*y*y
         a2=-174.83+152.78*log10(1.5682*(y+3.4))-808.74/(y+4.6157)
         a3=0.81177+0.56385*y+0.0040031*y*y-0.0057658*y**3+
     &0.00012057*y**4
         a4=0.68631*(y+3.32)+10.145*(y+3.32)**2-4.6176*(y+3.32)**3+
     &0.86824*(y+3.32)**4-0.053741*(y+3.32)**5
         a5=9.0466d-7+1.4539d-6*log10(0.015204*(y+3.4))+1.3253d-4/
     &(y+4.7171)**2-4.1228d-7*y+2.2036d-7*y*y
         a6=-339.45+618.73*log10(0.31595*(y+3.9))+
     &250.2/(y+4.4395)**2
         a7=-35.105+36.167*y-9.3575*y*y+0.33717*y*y*y
         a8=0.17554+0.373*y-0.014938*y*y+0.0032314*y**3+0.0025579*y**4
         else
         a0=0
         a1=0
         a2=0
         a3=0
         a4=0
         a5=0
         a6=0
         a7=0
         a8=0
         endif
         r=1.01
         if (Tp<=1.95) r=3.05*exp(-107*((yr+3.25)/
     &(1+8.08*(yr+3.25)))**2)
      else if (key==-1) then
         if (Tp>0.4) then
         a0=-0.018639*(y+3.3)+2.4315*(y+3.3)**2-0.57719*(y+3.3)**3
     &+0.063435*(y+3.3)**4
         a1=7.1827d-6-3.5067d-6*y+1.3264d-6*y**2-3.3481d-7*y**3
     &+2.3551d-8*y**4+3.4297d-9*y**5
         a2=563.91-362.18*log10(2.7187*(y+3.4))-2.8924d4/(y+7.9031)**2
         a3=0.52684+0.57717*y+0.0045336*y*y-0.0089066*y**3
         a4=0.36108*(y+3.32)+1.6963*(y+3.32)**2-0.074456*(y+3.32)**3-
     &0.071455*(y+3.32)**4+0.010473*(y+3.32)**5
         a5=9.7387d-5+7.8573d-5*log10(0.0036055*(y+4.3))+0.0002466/
     &(y+4.939)-3.8097d-7*y*y
         a6=-273-106.22*log10(0.341*(y+3.4))+89.037*y-12.546*y*y
         a7=432.53-883.99*log10(0.19737*(y+3.9))-4.1938d4/(y+8.5518)**2
         a8=-0.12756+0.43478*y-0.0027797*y*y-0.0083074*y**3
         else
         a0=0
         a1=0
         a2=0
         a3=0
         a4=0
         a5=0
         a6=0
         a7=0
         a8=0
         endif
         r=1.01
         if (Tp<=15.6) r=3.63*exp(-106*((yr+3.26)/
     &(1+9.21*(yr+3.26)))**2)-0.182*yr-0.175*yr*yr
      else if (key==1) then
         if (Tp>0.4) then
         a0=-0.79606*(y+3.3)+7.7496*(y+3.3)**2-3.9326*(y+3.3)**3
     &+0.80202*(y+3.3)**4-0.054994*(y+3.3)**5
         a1=6.7943d-6-3.5345d-6*y+6.0927d-7*y**2+2.0219d-7*y**3
     &+5.1005d-8*y**4-4.2622d-8*y**5
         a2=44.827+81.378*log10(0.027733*(y+3.5))-1.3886d4/(y+8.4417)**2
         a3=0.52010+0.59336*y+0.012032*y*y-0.0064242*y**3
         a4=2.1361*(y+3.32)+1.8514*(y+3.32)**2-0.47572*(y+3.32)**3+
     &0.0032034*(y+3.32)**4+0.0082955*(y+3.32)**5
         a5=1.0845d-6+1.4336d-6*log10(0.0077255*(y+4.3))+0.00013018/
     &(y+4.8188)**2+9.3601d-8*y
         a6=-267.74+14.175*log10(0.35391*(y+3.4))+64.669*y-7.7036*y*y
         a7=138.26-529.84*log10(0.12467*(y+3.9))-1.9869d4/(y+7.6884)**2
     &+1.0675*y*y
         a8=-0.14707+0.40135*y+0.0039899*y*y-0.0016602*y**3
         else
         a0=0
         a1=0
         a2=0
         a3=0
         a4=0
         a5=0
         a6=0
         a7=0
         a8=0
         endif
         r=1
         if (Tp<=5.52) r=2.22*exp(-98.9*((yr+3.25)/
     &(1+1.04*(yr+3.25)))**2)
      else if (key==-2) then
         if (Tp>0.4) then
         a0=0.0013113+0.36538*(y+3.31)+1.5178*(y+3.31)**2-
     &0.20668*(y+3.31)**3+0.024255*(y+3.31)**4
         a1=-4.7833e-6+4.5837e-5*exp(-0.42980*(y+3.4))+6.1559e-6/
     &(y+4.1731)+1.1928e-6*y
         a2=-245.22+73.223*y-19.652*y*y+0.83138*y*y*y+0.71564*y*y*y*y
         a3=0.45232+0.52934*y+0.010078*y*y-0.0017092*y*y*y
         a4=-0.0025734+0.38424*(y+3.32)+1.5517*(y+3.32)**2+
     &0.17336*(y+3.32)**3-0.17160*(y+3.32)**4+0.021059*(y+3.32)**5
         a5=4.7673e-5+5.4936e-5*log10(0.0067905*(y+4.3))+0.00020740/
     &(y+4.9772)
         a6=-270.30-114.47*log10(0.34352*(y+3.4))+80.085*y-7.9240*y*y
         a7=3271.9-2.9161e5/(y+87.847)-6.2330*y*y
         a8=-0.17787+0.36771*y-0.025397*y*y+0.0019238*y*y*y+
     &0.0032725*y*y*y*y
         else
         a0=0
         a1=0
         a2=0
         a3=0
         a4=0
         a5=0
         a6=0
         a7=0
         a8=0
         endif
         r=1
         if (Tp<=15.6) r=2.67*exp(-45.7*((yr+3.27)/
     &(1+6.59*(yr+3.27)))**2)-0.301*yr-0.208*yr*yr
      else if (key==2) then
         if (Tp>0.4) then
         a0=0.0074087+2.9161*(y+3.31)+0.99061*(y+3.31)**2-
     &0.28694*(y+3.31)**3+0.038799*(y+3.31)**4
         a1=-3.2480e-5+7.1944e-5*exp(-0.21814*(y+3.4))+2.0467e-5/
     &(y+4.1640)+5.6954e-6*y-3.4105e-7*y*y
         a2=-230.50+58.802*y-9.9393*y*y+1.2473*y*y*y-0.26322*y*y*y*y
         a3=0.45064+0.56930*y+0.012428*y*y-0.0070889*y*y*y
         a4=-0.011883+1.7992*(y+3.32)+3.5264*(y+3.32)**2-
     &1.7478*(y+3.32)**3+0.32077*(y+3.32)**4-0.017667*(y+3.32)**5
         a5=-1.6238e-7+1.8116e-6*exp(-0.30111*(y+3.4))+9.6112e-5/
     *(y + 4.8229)**2
         a6=-261.30-43.351*log10(0.35298*(y+3.4))+70.925*y-8.7147*y*y
         a7=184.45-1473.6/(y+6.8788)-4.0536*y*y
         a8=-0.24019+0.38504*y+0.0096869*y*y-0.0015046*y*y*y
         else
         a0=0
         a1=0
         a2=0
         a3=0
         a4=0
         a5=0
         a6=0
         a7=0
         a8=0
         endif
         r=1
         if (Tp<=7.81) r=0.329*exp(-249*((yr+3.26)/
     &(1+6.56*(yr+3.26)))**2)-0.957*yr-0.229*yr*yr
      else if (key==-3) then
         if (Tp>0.4) then
         a0=-1.5243*(y+3.3)+10.107*(y+3.3)**2-4.3126*(y+3.3)**3+
     &0.80081*(y+3.3)**4-0.048724*(y+3.3)**5
         a1=-2.6297e-5+9.3858e-5*exp(-0.32384*(y+3.4))+7.7821e-6/
     &(y+4.0560)+7.6149e-6*y-8.4091e-7*y*y
         a2=-243.62+59.374*y-5.7356*y*y+1.9815*y*y*y-1.0478*y*y*y*y
         a3=0.50807+0.60221*y+0.0034120*y*y-0.011139*y*y*y
         a4=2.6483*(y+3.32)+4.4585*(y+3.32)**2-1.2744*(y+3.32)**3+
     &0.11659*(y+3.32)**4+0.0030477*(y+3.32)**5
         a5=9.1101e-7+1.3880e-6*log10(0.016998*(y+4.3))+0.00012583/
     &(y+4.7707)**2
         a6=-272.11+53.477*log10(0.35531*(y+3.4))+56.041*y-6.0876*y*y
         a7=6431.8+893.92*log10(5.7013e-9*(y+3.9))+2103.6/(y+5.6740)-
     &6.1125*y*y
         a8=-0.1112+0.38144*y-0.040128*y*y+0.0047484*y**3+0.0054707*y**4
         else
         a0=0
         a1=0
         a2=0
         a3=0
         a4=0
         a5=0
         a6=0
         a7=0
         a8=0
         endif
         r=1
         if (Tp<=15.6) r=2.56*exp(-107*((yr+3.25)/
     &(1+8.34*(yr+3.25)))**2)-0.385*yr-0.125*yr*yr
      else if (key==3) then
         if (Tp>0.4) then
         a0=-0.63611*(y+3.3)+9.9015*(y+3.3)**2-4.5897*(y+3.3)**3+
     &0.91778*(y+3.3)**4-0.060724*(y+3.3)**5
         a1=6.8700e-6-2.8245e-6*y+7.6032e-7*y*y-3.2953e-7*y*y*y+
     &7.4292e-8*y*y*y*y
         a2=-240.46+58.405*y-9.8556*y*y+3.1401*y*y*y-0.88932*y*y*y*y
         a3=0.49935+0.60919*y+0.0024963*y*y-0.0099910*y*y*y
         a4=2.5094*(y+3.32)+4.1350*(y+3.32)**2-0.89534*(y+3.32)**3-
     &2.7577e-3*(y+3.32)**4+0.014511*(y+3.32)**5
         a5=8.2046e-7+1.4085e-6*log10(0.016793*(y+4.3))+0.00013340/
     &(y+4.7136)**2
         a6=-267.55-0.21018*log10(0.35217*(y+3.4))+69.586*y-9.9930*y*y
         a7=2741.8+222.01*log10(9.7401e-13*(y+3.9))-4772.5/(y+19.773)-
     &6.1001*y*y
         a8=-0.11857+0.39072*y-0.037813*y*y+0.0022265*y**3+
     &0.0046931*y**4
         else
         a0=0
         a1=0
         a2=0
         a3=0
         a4=0
         a5=0
         a6=0
         a7=0
         a8=0
         endif
         r=1
         if (Tp<=15.6) r=2.23*exp(-93.4*((yr+3.25)/
     &(1+8.38*(yr+3.25)))**2)-0.376*yr-0.121*yr*yr
      endif
      F_ND=(a0*exp(-a1*(x-a3+a2*(x-a3)**2)**2)+a4*exp(-a5*(x-a8+
     &a6*(x-a8)**2+a7*(x-a8)**3)**2))*r
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function F_NDKL(Eksec,Tp,key) ! Eq.(7)
      implicit real*8 (a-z)
      integer key
      x=log10(Eksec)
      Lmin=-2.6
      if (key==0) then
         Lmax=0.96*log10(Tp)
         W_NDL=15
         W_NDH=44
      else if (key==-1) then
         Lmax=0.96*log10(Tp)
         W_NDL=20
         W_NDH=45
      else if (key==1) then
         Lmax=0.94*log10(Tp)
         W_NDL=15
         W_NDH=47
      else if (key==-2) then
         Lmax=0.94*log10(Tp)
         W_NDL=20
         W_NDH=45
      else if (key==2) then
         Lmax=0.98*log10(Tp)
         W_NDL=15
         W_NDH=42
      else if (key==-3) then
         Lmax=0.98*log10(Tp)
         W_NDL=15
         W_NDH=40
      else if (key==3) then
         Lmax=0.98*log10(Tp)
         W_NDL=15
         W_NDH=40
      endif
      F_NDKL=1/((exp(W_NDL*(Lmin-x))+1)*(exp(W_NDH*(x-Lmax))+1))
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function F_DIFF(Eksec,Tp,key) ! Eq.(9)
      implicit real*8 (a-z)
      integer key
      x=log10(Eksec)
      p=sqrt(Tp*Tp+2*0.9382*Tp)
      F_DIFF=0.0
      if (p<2.25) return ! Eq.(2)
      y=log10(Tp*1d-3)
      if (key==0) then
         if (Tp>1.38) then
         if (Tp<5.52) then
            b0=0
            b1=0
            b2=0
            b3=0
         else
            b0=60.142*tanh(-0.37555*(y+2.2))-5.9564*(y+0.59913)**2+
     &6.0162d-3*(y+9.4773)**4
            b1=35.322+3.8026*tanh(-2.5979*(y+1.9))-
     &2.1870d-4*(y+369.13)**2
            b2=-15.732-0.082064*tanh(-1.9621*(y+2.1))+
     &2.3355d-4*(y+252.43)**2
            b3=-0.086827+0.37646*exp(-0.53053*((y+1.0444)/(1+
     &0.33933*(y-0.83562)))**2)
         endif
         b4=2.5982+0.39131*(y+2.95)**2-0.0049693*(y+2.95)**4+
     &0.94131*exp(-24.347*(y+2.45-0.19717*(y+2.45)**2)**2)
         b5=0.11198-0.64582*y+0.16114*y*y+2.2853*exp(-0.0032432*(
     &(y-0.83562)/(1.0+0.33933*(y-0.83562)))**2)
         b6=1.7843+0.91914*y+0.050118*y**2+0.038096*y**3-
     &0.027334*y**4-0.0035556*y**5+0.0025742*y**6
         b7=-0.19870-0.071003*y+0.019328*y**2-
     &0.28321*exp(-6.0516*(y+1.8441)**2)
         else
         b0=0
         b1=0
         b2=0
         b3=0
         b4=0
         b5=0
         b6=0
         b7=0
         endif
      else if (key==-1) then
         if (Tp>1.38) then
         if (Tp<5.52) then
            b0=0
            b1=0
            b2=0
            b3=0
         else
            b0=0.20463*tanh(-6.237*(y+2.2))-0.16362*(y+1.6878)**2+
     &3.5183d-4*(y+9.64)**4
            b1=1.6537+3.8530*exp(-3.2027*((y+2.0154)/
     &(1+0.62779*(y+2.0154)))**2)
            b2=-10.722+0.082672*tanh(1.8879*(y+2.1))+
     &1.4895d-4*(y+256.63)**2
            b3=-0.023752-0.51734*exp(-3.3087*((y+1.9877)/(1+
     &0.403*(y+1.988)))**2)
         endif
         b4=0.94921+0.1228*(y+2.9)**2-7.1585d-4*(y+2.9)**4+
     &0.5213*log10(y+2.9)
         b5=-4.2295-1.0025*tanh(9.0733*(y+1.9))-0.11452*(y-62.382)
         b6=1.4862+0.99544*y-0.042763*y**2-0.0040065*y**3+
     &0.0057987*y**4
         b7=6.2629+6.9517*tanh(-0.3648*(y+2.1))-0.26033*(y-2.8542)**2
         else
         b0=0
         b1=0
         b2=0
         b3=0
         b4=0
         b5=0
         b6=0
         b7=0
         endif
      else if (key==1) then
         if (Tp>1.38) then
         if (Tp<11.05) then
            b0=0
            b1=0
            b2=0
            b3=0
         else
            b0=29.192*tanh(-0.37879*(y+2.2))-3.2196*(y+0.675)**2+
     &3.6687d-3*(y+9.0824)**4
            b1=-142.97+147.86*exp(-0.37194*((y+1.8781)/
     &(1+3.8389*(y+1.8781)))**2)
            b2=-14.487-4.2223*tanh(-13.546*(y+2.2))+
     &0.00016988*(y+234.65)**2
            b3=-0.0036974-0.41976*exp(-6.1527*((y+1.8194)/(1+
     &0.99946*(y+1.8194)))**2)
         endif
         b4=1.8108+0.18545*(y+2.95)**2-2.0049d-3*(y+2.95)**4+
     &0.85084*exp(-14.987*(y+2.29-0.18967*(y+2.29))**2)
         b5=2.0404-0.51548*tanh(2.2758*(y+1.9))-0.035009*(y-6.6555)
         b6=1.5258+1.0132*y-0.064388*y**2-0.0040209*y**3+
     &0.0082772*y**4
         b7=3.0551+3.5240*tanh(-0.36739*(y+2.1))-0.13382*(y-2.7718)**2
         else
         b0=0
         b1=0
         b2=0
         b3=0
         b4=0
         b5=0
         b6=0
         b7=0
         endif
      else if (key==-2) then
         if (Tp>1.38) then
         if (Tp<11.0) then
            b0=0
            b1=0
            b2=0
            b3=0
         else
            b0=41.307*tanh(-0.37411*(y+2.2))-4.1223*(y+0.55505)**2+
     &0.0042652*(y+9.2685)**4
            b1=-132.50+142.12*exp(-8.0289*((y+1.9196)/(1.0+11.530*(y+
     &1.9196)))**2)
            b2=-17.223+0.011285*tanh(69.746*(y+1.9))-0.048233*y+
     &0.00025881*(y+250.77)**2
            b3=8.1991-9.6437*exp(-45.261*((y+1.9292)/(1.0+16.682*(y+
     &1.9292)))**2)
         endif
         b4=0.55919+0.36647*(y+2.95)+0.056194*(y+2.95)**2+
     &0.49957*exp(-5.5317*(y+2.2+0.43867*(y+2.2)**2)**2)
         b5=1.2544-0.52362*tanh(2.7638*(y+1.9))-0.055837*(y-17.638)
         b6=1.4788+1.0278*y-0.092852*y*y-0.0062734*y**3+0.011920*y**4
         b7=5.1651+5.7398*tanh(-0.37356*(y+2.1))-0.22234*(y-2.7889)**2
         else
         b0=0
         b1=0
         b2=0
         b3=0
         b4=0
         b5=0
         b6=0
         b7=0
         endif
      else if (key==2) then
         if (Tp>1.38) then
         if (Tp<11.0) then
            b0=0
            b1=0
            b2=0
            b3=0
         else
            b0=53.809*tanh(-0.41421*(y+2.2))-6.7538*(y+0.76010)**2+
     &0.0088080*(y+8.5075)**4
            b1=-50.211+55.131*exp(-1.3651*((y+1.8901)/(1.0+4.4440*(y
     &+1.8901)))**2)
            b2=-17.231+0.041100*tanh(7.9638*(y+1.9))-0.055449*y+
     &0.00025866*(y+250.68)**2
            b3=12.335-12.893*exp(-1.4412*((y+1.8998)/(1.0+5.5969*(y+
     &1.8998)))**2)
         endif
         b4=1.3558+0.46601*(y+2.95)+0.052978*(y+2.95)**2+
     &0.79575*exp(-5.4007*(y+2.2+4.6121*(y+2.2)**2)**2)
         b5=1.8756-0.42169*tanh(1.6100*(y+1.9))-0.051026*(y-3.9573)
         b6=1.5016+1.0118*y-0.072787*y*y-0.0038858*y**3+0.0093650*y**4
         b7=4.9735+5.5674*tanh(-0.36249*(y+2.1))-0.20660*(y-2.8604)**2
         else
         b0=0
         b1=0
         b2=0
         b3=0
         b4=0
         b5=0
         b6=0
         b7=0
         endif
      else if (key==-3) then
         if (Tp>1.38) then
         if (Tp<11.0) then
            b0=0
            b1=0
            b2=0
            b3=0
         else
            b0=70.430*tanh(-0.35816*(y+2.2))-6.6796*(y+0.52273)**2+
     &0.0065659*(y+9.5266)**4
            b1=-8.1145+7686.0*exp(-44046*((y+2.2190)/(1.0+
     &81.105*(y+2.2190)))**2)
            b2=-1.3095+0.071270*tanh(0.0075463*(y+1.9))+
     &0.067759*(y+5.3433)-0.0044205*(y-1.8683)**2
            b3=0.082149-2190.1*exp(-533.75*((y+2.8363)/(1.0+
     &7.0976*(y+2.8363)))**2)
         endif
         b4=2.7540+0.33859*(y+2.95)**2-0.0034274*(y+2.95)**4+
     &1.1679*exp(-10.408*(y+2.28-0.18922*(y+2.2)**2)**2)
         b5=2.1817-0.59584*tanh(2.7054*(y+1.9))-0.010909*(y-14.880)
         b6=1.4591+1.0275*y-0.074949*y*y-0.0060396*y**3+0.0097568*y**4
         b7=3.7609+4.2843*tanh(-0.37148*(y+2.1))-0.16479*(y-2.7653)**2
         else
         b0=0
         b1=0
         b2=0
         b3=0
         b4=0
         b5=0
         b6=0
         b7=0
         endif
      else if (key==3) then
         if (Tp>1.38) then
         if (Tp<11.0) then
            b0=0
            b1=0
            b2=0
            b3=0
         else
            b0=64.682*tanh(-0.34313*(y+2.2))-5.5955*(y+0.44754)**2+
     &0.0050117*(y+9.9165)**4
            b1=-7.6016+3042.7*exp(-1.0134e4*((y+2.3066)/(1.0+
     &41.612*(y+2.3066)))**2)
            b2=-1.4978-0.58163*tanh(-0.36488*(y+1.9))+
     &0.031825*(y+2.8097)+0.022796*(y-1.8861)**2
            b3=-0.0061483-65.799*exp(-4.8239*((y+3.8835)/(1.0+
     &0.53343*(y+3.8835)))**2)
         endif
         b4=2.8009+0.35351*(y+2.95)**2-0.0039779*(y+2.95)**4+
     &1.3012*exp(-10.592*(y+2.28-0.19149*(y+2.28)**2)**2)
         b5=1.8016-0.69847*tanh(2.8627*(y+1.9))-0.015722*(y-45.410)
         b6=1.4617+1.0167*y-0.078617*y*y-0.0038336*y*y*y+0.010141*y**4
         b7=3.5599+4.0041*tanh(-0.41889*(y+2.1))-0.18182*(y-2.4209)**2
         else
         b0=0
         b1=0
         b2=0
         b3=0
         b4=0
         b5=0
         b6=0
         b7=0
         endif
      endif
      F_DIFF=b0*exp(-b1*((x-b2)/(1+b3*(x-b2)))**2)+
     &b4*exp(-b5*((x-b6)/(1+b7*(x-b6)))**2)
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function F_KL(Eksec,Tp) ! Eq.(7)
      implicit real*8 (a-z)
      x=log10(Eksec)
      Lmax=log10(Tp)
      W_DIFF=75
      F_KL=1/(exp(W_DIFF*(x-Lmax))+1)
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function F_RES(Eksec,Tp,key) ! Eq.(12)
      implicit real*8 (a-z)
      integer key
      x=log10(Eksec)
      p=sqrt(Tp*Tp+2*0.9382*Tp) 
      y=log10(Tp*1d-3)
      if (key==0) then
         if (Tp>0.4.and.Tp<2) then
         c0=2.4316*exp(-69.484*((y+3.1301)/(1.0+1.24921*(y+
     &3.1301)))**2)-6.3003-9.5349/y+0.38121*y*y
         c1=56.872+40.627*y+7.7528*y*y
         c2=-5.4918-6.7872*tanh(4.7128*(y+2.1))+0.68048*y
         c3=-0.36414+0.039777*y
         c4=-0.72807-0.48828*y-0.092876*y**2
         else
         c0=0
         c1=0
         c2=0
         c3=0
         c4=0
         endif
         if (Tp>0.6.and.Tp<2.9) then
         d0=3.2433*exp(-57.133*((y+2.9507)/(1.0+1.2912*(y+2.9507)))**2)
     &-1.0640-0.43925*y
         d1=16.901+5.9539*y-2.1257*y**2-0.92057*y**3
         d2=-6.6638-7.5010*tanh(30.322*(y+2.1))+0.54662*y
         d3=-1.50648-0.87211*y-0.17097*y**2
         d4=0.42795+0.55136*y+0.20707*y**2+0.027552*y**3
         else
         d0=0
         d1=0
         d2=0
         d3=0
         d4=0
         endif
      else if (key==-1) then
         c0=0
         c1=0
         c2=0
         c3=0
         c4=0
         if (Tp>0.6.and.Tp<2.9) then
         d0=0.3779*exp(-56.836*((y+2.9537)/(1+1.5221*(y+2.9537)))**2)
     &-0.059458+0.0096583*y*y
         d1=-5.5135-3.3988*y
         d2=-7.1209-7.1850*tanh(30.801*(y+2.1))+0.35108*y
         d3=-6.7841-4.8358*y-0.91523*y**2
         d4=-134.03-139.63*y-48.316*y**2-5.5526*y**3
         else
         d0=0
         d1=0
         d2=0
         d3=0
         d4=0
         endif
      else if (key==1) then
         if (Tp>0.4.and.Tp<2) then
         c0=2.9841*exp(-67.857*((y+3.1272)/(1.0+0.22831*(y+3.1272)))**2)
     &-(6.5855+9.6984/y-0.41256*y*y)
         c1=6.8276+5.2236*y+1.4630*y*y
         c2=-6.0291-6.4581*tanh(5.0830*(y+2.1))+0.46352*y
         c3=0.59300+0.36093*y
         c4=0.77368+0.44776*y+0.056409*y*y
         else
         c0=0
         c1=0
         c2=0
         c3=0
         c4=0
         endif
         if (Tp>0.6.and.Tp<2.9) then
         d0=1.9186*exp(-56.544*((y+2.9485)/(1+1.2892*(y+2.9485)))**2)
     &-0.2372+0.041315*y*y
         d1=-4.9866-3.1435*y
         d2=-7.0550-7.2165*tanh(31.033*(y+2.1))+0.38541*y
         d3=-2.8915-2.1495*y-0.45006*y**2
         d4=-1.297-0.13947*y+0.41197*y**2+0.10641*y**3
         else
         d0=0
         d1=0
         d2=0
         d3=0
         d4=0
         endif
      else if (key==-2) then
         c0=0
         c1=0
         c2=0
         c3=0
         c4=0
         if (Tp>0.6.and.Tp<2.9) then
         d0=0.36459*exp(-58.210*((y+2.9537)/(1.0+1.4320*(y+2.9537)))**2)
     &-(0.11283+0.046244*y)
         d1=-9.5066-5.4655*y-0.31769*y*y
         d2=-7.1831-7.1551*tanh(30.354*(y+2.1))+0.33757*y
         d3=2.7938+1.6992*y+0.20161*y*y
         d4=0.61878+0.62371*y+0.18913*y*y+0.019118*y*y*y
         else
         d0=0
         d1=0
         d2=0
         d3=0
         d4=0
         endif
      else if (key==2) then
         if (Tp>0.4.and.Tp<2) then
         c0=2.8290*exp(-71.339*((y+3.1282)/(1.0+0.48420*(y+3.1282)))**2)
     &-(9.6339+15.733/y-0.52413*y*y)
         c1=-24.571-15.831*y-2.1200*y*y
         c2=-5.9593-6.4695*tanh(4.7225*(y+2.1))+0.50003*y
         c3=0.26022+0.24545*y
         c4=0.076498+0.061678*y+0.0040028*y*y
         else
         c0=0
         c1=0
         c2=0
         c3=0
         c4=0
         endif
         if (Tp>0.6.and.Tp<2.9) then
         d0=1.7951*exp(-57.260*((y+2.9509)/(1.0+1.4101*(y+2.9509)))**2)
     &-(0.58604+0.23868*y)
         d1=-2.6395-1.5105*y+0.22174*y*y
         d2=-7.0512-7.1970*tanh(31.074*(y+2.1))+0.39007*y
         d3=-1.4271-1.0399*y-0.24179*y*y
         d4=0.74875+0.63616*y+0.17396*y*y+0.017636*y*y*y
         else
         d0=0
         d1=0
         d2=0
         d3=0
         d4=0
         endif
      else if (key==-3) then
         if (Tp>0.4.and.Tp<2) then
         c0=2.8262*exp(-62.894*((y+3.1250)/(1.0-0.47567*(y+3.1250)))**2)
     &+5.6845+13.409/y-0.097296*y*y
         c1=16.721+11.750*y+2.4637*y*y
         c2=-6.0557-6.3378*tanh(21.984*(y+2.1))+0.43173*y
         c3=0.37009+0.27706*y
         c4=0.047507+0.061570*y+0.0070117*y*y
         else
         c0=0
         c1=0
         c2=0
         c3=0
         c4=0
         endif
         if (Tp>0.6.and.Tp<2.9) then
         d0=2.2400*exp(-57.159*((y+2.9492)/(1.0+1.2994*(y+2.9492)))**2)
     &-(0.66521+0.27554*y)
         d1=-7.0650-4.2773*y-0.17648*y*y
         d2=-7.0410-7.1977*tanh(31.095*(y+2.1))+0.40238*y
         d3=-1.2354-0.87581*y-0.20829*y*y
         d4=-0.11395+0.34418*y+0.27103*y*y+0.050248*y*y*y
         else
         d0=0
         d1=0
         d2=0
         d3=0
         d4=0
         endif
      else if (key==3) then
         if (Tp>0.4.and.Tp<2) then
         c0=3.6052*exp(-60.914*((y+3.1278)/(1.0-0.19497*(y+3.1278)))**2)
     &-0.92514+2.1315/y+0.23548*y*y
         c1=95.310+70.497*y+13.636*y*y
         c2=-6.2158-6.2939*tanh(21.592*(y+2.1))+0.37440*y
         c3=2.7485+1.1692*y
         c4=-2.7568-1.8461*y-0.31376*y*y
         else
         c0=0
         c1=0
         c2=0
         c3=0
         c4=0
         endif
         if (Tp>0.6.and.Tp<2.9) then
         d0=2.5489*exp(-58.488*((y+2.9509)/(1.0+1.3154*(y+2.9509)))**2)
     &-(0.83039+0.34412*y)
         d1=88.173+65.148*y+12.585*y*y
         d2=-7.0962-7.1690*tanh(30.890*(y+2.1))+0.38032*y
         d3=-4.1440-3.2717*y-0.70537*y*y
         d4=2.2624+1.1806*y+0.0043450*y*y-0.043020*y*y*y
         else
         d0=0
         d1=0
         d2=0
         d3=0
         d4=0
         endif
      endif
      F_1232=c0*exp(-c1*((x-c2)/(1+c3*(x-c2)+c4*(x-c2)**2))**2)
      F_1600=d0*exp(-d1*((x-d2)/(1+d3*(x-d2)+d4*(x-d2)**2))**2)
      if (p<1.4.or.p>10.0) F_1232=0.0 ! Eq.(3)
      if (p<1.6.or.p>20.0) F_1600=0.0 ! Eq.(4)
      F_RES=F_1232+F_1600
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * pp_meson.f *                                  galprop package * 9/09/2003 
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

      function PP_MESON(Esec, Pp1, NA1, NA2, key1)
c***********************************************************************
c### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c### PURPOSE: MAIN subprogram for calculation of the LS spectrum of
c### secondary electrons, positrons, and gamma-rays [barn/GeV]
c### produced in AA-, AP-, and PP-collisions (1 target nucleus per cm^3).
c### REFERENCES:
c### [1] Badhwar, Stephens, & Golden 1977, Phys. Rev. D 15, 820
c### [2] Dermer 1986, A&A 157, 223;  ApJ 307, 47
c### [3] Mori 1997, ApJ 478, 225
c### [4] Moskalenko & Strong 1998, ApJ 493, 694
c### [5] Stecker 1970, Astrophys. Spa. Sci. 6, 377
c### [6] Stephens & Badhwar 1981, Astrophys. Spa. Sci. 76, 213
c INPUT/OUTPUT parameters:
c Esec [GeV] - secondary particle or photon energy;
c Pp1 [GeV/c] -beam momentum per nucleUS;
c NA1 & NA2 are the atomic numbers of beam and target nuclei,
c correspondingly (NA1=NA2=1 for pp-collisions);
c key = 0 for pp->pi^0 +X;
c     < 0 for pp->pi^- +X react.;
c     > 0 for pp->pi^+ +X react.,
c if |key| = 3, 4 then KAON decay included,
c if |key| >= 4 then an isotropic distribution for electrons/positrons
c    in the muon rest system is assumed.
c The values |key| = 1, 3 are recommended.
c NB Don't use |key|=2, this value is reserved for the internal purposes !
c Delta-isobar mass has the Breit-Wigner distr. (M0 -mean mass, G -width);
c History:
c !! uses common "common"-block with ANTIPROTON: /en/...,Pap !! =5.6.98=
c common/key/ added "muonkey" to avoid inf at gamma_pi=1 
c muonkey=1 to use in STECKER & D_PION, muonkey=0 in SB_MESON IMOS20030909
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     1      /thres/Pth0,Pth1,Pth2,Pth3,Pth4,Pth5
     1      /branch/BR1,BR2
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c
     3      /key/key, kaonkey, muonkey
      external STECKER, SB_MESON, D_PION
      
      PP_MESON = 0.
      Pp  = Pp1/NA1
      if(Pp .le. Pth0) return ! the lowest threshold momentum for pp->pi+X
      AI1 = 0.d0
      AI2 = 0.d0
      AI3 = 0.d0
      AI4 = 0.d0

      P1 = 3.d0   ! GeV/c, lower boundary of the interpolation region
      P2 = 7.d0   ! GeV/c, upper boundary of the interpolation region
      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) !Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      key = key1
      Egam= Esec
      Pe  = Esec  ! massless electron/positron
      if(key .ne. 0 .and. Pe .le. 0.d0) return
      
      if(key .eq. 0) then  !--- GAMMA-RAYS ---!
         Mmes= Mpi0        ! Neutral pion rest mass
         MX = 2.d0*Mp      ! The mass in the channel X (=2Mp)
         gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s) ! max pion Lf in CMS
         if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
         betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)         ! max beta*gamma -"-
         EL = dmax1(Mmes,Egam+Mmes**2/4.d0/Egam) !Lower limit for pion LS Energy
         EU = Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max) !Upper limit
      else                 !--- ELECTRONS/POSITRONS ---!
         Mmes= Mpi1        ! Charged pion rest mass
         if(Pp .le. Pth3) return   !the lowest thres. momentum for pi^+ produc.
         MX = Mp+Mn                ! The mass in the channel X (=Mp+Mn)
         if(key .lt. 0) then          !--- electrons ---!
            if(Pp .le. Pth1) return   ! the threshold momentum for pp->pi^- X
            MX = 2.*Mp+Mmes           ! The mass in the channel X (=2Mp+Mmes)
         endif
         gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s) ! max pion Lf in CMS
         if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
         betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)         ! max beta*gamma -"-
         Emu1 = (Mmes**2 +Mmu**2)/2.d0/Mmes ! Muon energy in the pion rest mass
         Pmu1 = dsqrt(Emu1**2 -Mmu**2)      ! Muon momentum -"-  ~30 MeV/c
         EL=Mmes*dmax1(1.d0,Pe/(Emu1+Pmu1)+(Emu1+Pmu1)/4.d0/Pe) !Lo PION Lf
         if(Pe .lt. (Emu1+Pmu1)/2.d0) EL = Mmes
c^         EL=Mmes*dmax1(1.d0,(Emu1-Pmu1)*Pe/Mmu**2+(Emu1+Pmu1)/4.d0/Pe)!Lo PION Lf
         EU=Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max)  !Up PION Lf
      endif

      if(EL .lt. EU) then
c STECKER: integral over the Breit-Wigner mass distribution for Delta-isobara
         call CS_TOT(Pp,CS_S,CS_SB,CS_D,CS_K)
         if(Pp .lt. P2) then
c#         if(dsqrt(s)-Mp .ge. M0) AI1 = -STECKER(M0) ! delta-function approx.
            if(key.le.0 .or. Pp.gt.Pth2) ! the threshold momentum for pp->pi^+ X
     1      call SIM2(dsqrt(s)-Mp,Mp+Mmes,G/9.,1.d-3,1.d-10,STECKER,AI1)
            AI1 = -AI1
         endif
c STEPHENS & BADHWAR: integral over the pion momentum in LS
         if(Pp .gt. P1) then
            call SIM2(EL, EU, EL/9.,1.d-3,1.d-10,SB_MESON,AI2)
            AI2=AI2*2.*3.1415926  /CS_SB
         endif
c positron production in pp->pi^+ +d reaction (at >4 GeV it's negligible)
         if(key .gt. 0 .and. Pe .le. 4.d0) then
c#         if(Pp .le. Pth3) return !the threshold momentum for pp->pi^+ d
            MX = Md                ! The mass in the channel X (=Md)
            gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s)!max pion Lf in CMS
            if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
            betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)        !max beta*gamma -"-
            EL=Mmes*dmax1(1.d0,Pe/(Emu1+Pmu1)+(Emu1+Pmu1)/4.d0/Pe) !Lo PION Lf
            if(Pe .lt. (Emu1+Pmu1)/2.d0) EL = Mmes
            if(EL .eq. Mmes) EL = EL * (1.d0+1.d-8) ! IMOS20030909
c^         EL=Mmes*dmax1(1.d0,(Emu1-Pmu1)*Pe/Mmu**2+(Emu1+Pmu1)/4.d0/Pe)!Lo PION Lf
            EU=Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max)  !Up PION Lf
            if(EL .lt. EU)
     #         call SIM2(EU,EL,EL/9.,1.d-4,1.d-11,D_PION,AI3)
         endif
      endif

c Interpolation between P1 & P2 GeV/c; the total cross section is corrected
      A12 = (AI2-AI1)*(Pp-P1)/(P2-P1) +AI1  !(CS_S/CS_SB) ~ 1.09
      if(Pp .le. P1 .or. Pp .ge. P2) A12 = AI1+AI2

      PP_MESON = (A12*CS_S-AI3*CS_D) *1.d-3  !*1.d-24 ! 1 barn=10^{-24}cm^2
     1     *(NA1**(3./8.)+NA2**(3./8.)-1.)**2

c KAON: electron/positron spectrum from pp -> K +X reaction
      if(iabs(key) .ge. 3 .and. iabs(key1) .le. 4) then !KAON decay mode K->mu +nu
         Mmes= MK    ! KAON rest mass
         if(Pp .le. Pth5) return    !the lowest thres. momentum for K produc.
         MX = Mp+Mn                 ! The mass in the channel X (=Mp+Mn)
         if(key .lt. 0) then        !--- electrons ---!
            if(Pp .le. Pth4) return ! the threshold momentum for pp->K^- X
            MX = 2.*Mp+Mmes         ! The mass in the channel X (=2Mp+Mmes)
         endif
         gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s) ! max pion Lf in CMS
         if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
         betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)         ! max beta*gamma -"-
         Emu1 = (Mmes**2 +Mmu**2)/2.d0/Mmes ! Muon energy in the KAON rest syst.
         Pmu1 = dsqrt(Emu1**2 -Mmu**2)      ! Muon momentum -"- 
         EL=Mmes*dmax1(1.d0,Pe/(Emu1+Pmu1)+(Emu1+Pmu1)/4.d0/Pe) !Lo PION Lf
         if(Pe .lt. (Emu1+Pmu1)/2.d0) EL = Mmes
         EU=Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max)  !Up KAON Lf
         if(EL .ge. EU) return

c STEPHENS & BADHWAR: integral over KAON momentum in LS
         kaonkey = 1
         call SIM2(EL, EU, EL/9.,1.d-3,1.d-10,SB_MESON,AI4)
         kaonkey = 0
         AI4=AI4*2.*3.1415926

         PP_MESON=(A12*CS_S-AI3*CS_D+AI4*BR1) *1.d-3 !barn, 1 barn=10^{-24}cm^2
     1     *(NA1**(3./8.)+NA2**(3./8.)-1.)**2
      endif
ccc      print *, "PP_MESON: ",AI1,AI2,AI3,AI4
      return
      end

      function MUON(gam_mes)
c***********************************************************************
c PURPOSE: calculation of the spectrum of electrons/positrons [1/GeV] from
c charged pion/kaon decay (MESON ->muon +neutrino, for the fixed meson Lf).
c key > 0 then the positron spectrum; < 0 for the electron spectrum.
c if |key| >3 then isotropic distribution for electrons/positrons,
c if |key| =2 then KAON decay.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c
     3      /key/key, kaonkey, muonkey
       Gmu1(gam,bet) = (                     !# mu^- decay
     1  (Pe/Mmu)**3*(-32.d0*(1.d0-bet)*gam**3+(24.d0-32.d0*bet)*gam)
     2 +(Pe/Mmu)**2*((27.d0-9.d0*bet)*gam**2+9.d0*dlog((1.d0+bet)*gam))
     3                 )*4.d0/9.d0      
       Gmu2(gam,bet) = (                    !# mu^+ decay
     1  (Pe/Mmu)**3*(-32.d0*(1.d0+bet)*gam**3+(24.d0+32.d0*bet)*gam)
     2 +(Pe/Mmu)**2*((27.d0+9.d0*bet)*gam**2-9.d0*dlog((1.d0+bet)*gam))
     3                 )*4.d0/9.d0
       Gmu0(gam) = (Pe/Mmu)**3*(-128.d0/9.d0*gam**3+32.d0/3.d0*gam)
     1            +(Pe/Mmu)**2*12.d0*gam**2

      Mmes= Mpi1    ! Charged PION rest mass
      if(kaonkey .eq. 1) Mmes= MK ! KAON rest mass

      MUON = 0.d0
      Emu1 = (Mmes**2 +Mmu**2)/2.d0/Mmes ! Muon energy in the MESON rest syst.
      Pmu1 = dsqrt(Emu1**2 -Mmu**2)      ! Muon momentum in the MESON rest syst.
      if(gam_mes .le. 1.d0) gam_mes = 1.d0
      betgam_mes = dsqrt(gam_mes**2 -1.d0) ! Pion betgam in the LS

      EmuL = gam_mes*Emu1 -betgam_mes*Pmu1 !Lower limit for the muon LS energy
      if(EmuL .le. Mmu) EmuL = Mmu
      PmuL = dsqrt(EmuL**2-Mmu**2)         ! -"- muon momentum
      EmuU = gam_mes*Emu1 +betgam_mes*Pmu1 !Upper limit for the muon LS energy
      if(EmuU .le. Mmu) EmuU = Mmu
      PmuU = dsqrt(EmuU**2-Mmu**2)         ! -"- muon momentum

      if (key*key .gt. 10) then            ! TEST with isotropic distribution
         if(Pe .le. Mmu**2/(2.d0*(EmuU+PmuU))) then           ! (I)    
            MUON =Gmu0(EmuU/Mmu) -Gmu0(EmuL/Mmu)
         else
            EmuB = Pe +Mmu**2/4.d0/Pe          ! A boundary from electron energy
            if(EmuB .le. Mmu) EmuB = Mmu
            PmuB = dsqrt(EmuB**2-Mmu**2)       ! -"- muon momentum
            if(Pe .le. Mmu**2/(2.d0*(EmuL+PmuL))) then        ! (II)
               MUON =FMU(EmuU,key) -FMU(EmuB,key)
     2              +Gmu0(EmuB/Mmu)-Gmu0(EmuL/Mmu)
            else
               if(Pe .le. Mmu**2/(2.d0*(EmuL-PmuL))) then     ! (III)
                  MUON =FMU(EmuU,key)-FMU(EmuL,key)
               else
                  if(Pe .le. Mmu**2/(2.d0*(EmuU-PmuU))) then  ! (IV)
                     MUON =FMU(EmuU,key)-FMU(EmuB,key)
                  endif
               endif
            endif
         endif
         MUON = MUON /Pmu1                         ! IMOS20030909
         if(muonkey.eq.1) MUON = MUON /betgam_mes  ! IMOS20030909
         return
      endif

c~~      if(key      .ne.       0) then  
      if (key .lt. 0) then               ! electron spectrum
         if(Pe .le. Mmu**2/(2.d0*(EmuU+PmuU))) then           ! (I)    
            MUON =Gmu1(EmuU/Mmu,PmuU/EmuU) -Gmu1(EmuL/Mmu,PmuL/EmuL)
         else
            EmuB = Pe +Mmu**2/4.d0/Pe          ! A boundary from electron energy
            if(EmuB .le. Mmu) EmuB = Mmu
            PmuB = dsqrt(EmuB**2-Mmu**2)       ! -"- muon momentum
            if(Pe .le. Mmu**2/(2.d0*(EmuL+PmuL))) then        ! (II)
               MUON =FMU(EmuU,key) -FMU(EmuB,key)
     2              +Gmu1(EmuB/Mmu,PmuB/EmuB)-Gmu1(EmuL/Mmu,PmuL/EmuL)
            else
               if(Pe .le. Mmu**2/(2.d0*(EmuL-PmuL))) then     ! (III)
                  MUON =FMU(EmuU,key)-FMU(EmuL,key)
               else
                  if(Pe .le. Mmu**2/(2.d0*(EmuU-PmuU))) then  ! (IV)
                     MUON =FMU(EmuU,key)-FMU(EmuB,key)
                  endif
               endif
            endif
         endif
      else                               ! positron spectrum
         if(Pe .le. Mmu**2/(2.d0*(EmuU+PmuU))) then           ! (I)    
            MUON =Gmu2(EmuU/Mmu, PmuU/EmuU) -Gmu2(EmuL/Mmu,PmuL/EmuL)
         else
            EmuB = Pe +Mmu**2/4.d0/Pe          ! A boundary from positron energy
            PmuB = dsqrt(EmuB**2-Mmu**2)       ! -"- muon momentum
            if(Pe .le. Mmu**2/(2.d0*(EmuL+PmuL))) then        ! (II)
               MUON =FMU(EmuU,key) -FMU(EmuB,key)
     2              +Gmu2(EmuB/Mmu,PmuB/EmuB) -Gmu2(EmuL/Mmu,PmuL/EmuL)
            else
               if(Pe .le. Mmu**2/(2.d0*(EmuL-PmuL))) then     ! (III)
                  MUON =FMU(EmuU,key)-FMU(EmuL,key)
               else
                  if(Pe .le. Mmu**2/(2.d0*(EmuU-PmuU))) then  ! (IV)
                     MUON =FMU(EmuU,key)-FMU(EmuB,key)
                  endif
               endif
            endif
         endif
      endif
      MUON = MUON /Pmu1                        ! IMOS20030909
      if(muonkey.eq.1) MUON = MUON /betgam_mes ! IMOS20030909
      if(MUON .lt. 0.d0) MUON = 0.d0
      return
      end

      function FMU(Emu, key)
c***********************************************************************
c PURPOSE: a subroutine used for calculation of the electron/positron
c spectrum from charged pion/kaon decay.
c key > 0 then the positron spectrum; < 0 for the electron spectrum.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c

      gam = Emu/Mmu
      bet = dsqrt(gam**2-1.d0)/gam

      if(key*key .ge. 10) then       !# TEST with uniform distribution
         if(gam .lt. 2.d1) then
            FMU = ( 
     1 (Pe/Mmu)**3*(-128.d0*(1.d0-bet)*gam**3 +(96.d0-32.d0*bet)*gam)
     2+(Pe/Mmu)**2*108.d0*(1.d0-bet)*gam**2 +15.d0*dlog(gam*(1.d0+bet))
     3                  )/18.d0
         else
            FMU = 
     1-(Pe/Mmu)**3*(2.d0/9.d0 +(1.d0/6.d0 +(1.d0/8.d0 +(7.d0/72.d0 
     2         +5.d0/64.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**3
     3+(Pe/Mmu)**2*(3.d0+(3.d0/4.d0 +(3.d0/8.d0 +(15.d0/64.d0 
     4         +(21.d0/128.d0 +(63.d0/512.d0 +99.d0/1024.d0/gam**2)
     5         /gam**2)/gam**2)/gam**2)/gam**2)/gam**2)
     6  -(5.d0/24.d0 +(5.d0/64.d0 +(25.d0/576.d0 +(175.d0/6144.d0 
     7  +(21.d0/1024.d0 +385.d0/24576.d0/gam**2)/gam**2)/gam**2)
     8  /gam**2)/gam**2)/gam**2  +5.d0/6.d0*dlog(2.d0*gam)  
     
         endif
         return
      endif

c~~      if(key      .ne.       0) then  
      if(key .lt. 0) then            !# mu^- decay
c         if(gam .lt. 2.d4) then
         if(gam .lt. 2.d1) then
            FMU = ( 
     1 (Pe/Mmu)**3*( -512.d0*(1.d0-bet)*gam**3 +(576.d0-320.d0*bet)*gam
     2                +48.d0*dlog((gam-1.d0)/(gam+1.d0))  )
     3+(Pe/Mmu)**2*(288.d0*(1.d0-bet)*gam**2+72.d0*dlog((1.d0+bet)/bet))
     4  +6.d0*dlog(gam*bet) +30.d0*dlog((1.d0+bet)*gam) )/36.d0
         else
            FMU =
     1-(Pe/Mmu)**3*(2.d0/3.d0 +(8.d0/15.d0 +(71.d0/168.d0+(149.d0/432.d0
     2         +611.d0/2112.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**3 
     3+(Pe/Mmu)**2*(4.d0+dlog(4.d0) +(3.d0/2.d0 +(13.d0/16.d0 
     4         +(13.d0/24.d0 +(205.d0/512.d0 +(403.d0/1280.d0 
     5         +1585.d0/6144.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**2)
     6         /gam**2)
     7  +dlog(gam)+5.d0/6.d0*dlog(2.d0) -(7.d0/24.d0 +(23.d0/192.d0
     8  +(41.d0/576.d0 +(101.d0/2048.d0 +(571.d0/15360.d0
     9  +2179.d0/73728.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**2)/gam**2
         endif      
      else                           !# mu^+ decay
c         if(gam .lt. 2.d4) then
         if(gam .lt. 2.d1) then
            FMU = (                   
     1 (Pe/Mmu)**3*(16.d0*dlog((gam+1.d0)/(gam-1.d0))
     2             -64.d0*(1.d0-bet)*gam)
     3+(Pe/Mmu)**2*(48.d0*(1.d0-bet)*gam**2 +24.d0*dlog(bet/(1.d0+bet)))
     4  -2.d0*dlog(gam*bet) +10.d0*dlog((1.d0+bet)*gam) )/12.d0
     5                 
         else
            FMU = 
     1 (Pe/Mmu)**3*(2.d0/9.d0 +(1.d0/5.d0 +(29.d0/168.d0 +(65.d0/432.d0
     2         +281.d0/2112.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**3 
     3+(Pe/Mmu)**2*(2.d0-dlog(4.d0) -(1.d0/16.d0 +(7.d0/96.d0 
     4         +(37.d0/512.d0 +(11.d0/160.d0 +397.d0/6144.d0/gam**2) 
     5         /gam**2)/gam**2)/gam**2)/gam**4)
     6  +2.d0/3.d0*dlog(gam) +5.d0/6.d0*dlog(2.d0) -(1.d0/8.d0 
     7  +(7.d0/192.d0 +(1.d0/64.d0 +(47.d0/6144.d0 +(59.d0/15360.d0
     8  +131.d0/73728.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**2)/gam**2
         endif
      endif
      return
      end

      subroutine CS_TOT(Pp1,CS_S,CS_SB,CS_D,CS_K)
c***********************************************************************
c PURPOSE: calculation of the total inclusive CROSS SECTION [mbarn]
c of the meson production in the reaction p+p->(pi,K) +X
c (Dermer 1986, A&A 157,223; ApJ 307, 47).
c INPUT/OUTPUT parameters:
c E_gam [GeV] - energy of secondary particle/photon;
c P_p [GeV/c] - proton beam momentum
c key = 0 for pi^0;
c     < 0 for pp->pi^- +X react.;
c     > 0 for pp->pi^+ +X react.
c 2< |key| <4 also calculates KAON production cross section.
c CS_S [mbarn] - PION production cross section (the Dermer's fit);
c CS_SB [mbarn] - PION cross section in the Stephens & Badhwar formalism;
c CS_D [mbarn] - PION cross section in the deutron channel (p+p->pi+d);
c CS_K [mbarn] - KAON cross section in the Stephens & Badhwar formalism.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     1      /thres/Pth0,Pth1,Pth2,Pth3,Pth4,Pth5
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c
     3      /key/key, kaonkey, muonkey
      external SB_EPI
      
      CS_S = 0.d0  ! PION product. cross sect. [mbarn], DERMER's fit
      CS_SB = 0.d0 ! PION product. cross sect. [mbarn], SB's approximation
      CS_D = 0.d0  ! PION product. cross sect. [mbarn], deutron channel (DERMER)
      CS_K = 0.d0  ! KAON product. cross sect. [mbarn], SB's approximation
      P1 = 3.d0
      if(Pp .le. Pth0) return    ! the lowest threshold momentum for pp->pi+X
      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) !Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      Pp = Pp1

c ### GAMMA-RAYS ###
      if(key .eq. 0) then
         Mmes= Mpi0    ! NEUTRAL PION rest mass
         MX = 2.d0*Mp  ! The mass in the channel X (=2Mp)
c Dermer's fit
         if(Pp .le. 0.96) then
            eta = dsqrt((s-Mmes**2-MX**2)**2-4.d0*(Mmes*MX)**2)
     #           /dsqrt(s)/2.d0/Mmes
            CS_S = eta**2*( 0.032 +eta**4*(0.040 +eta**2*0.047) )
         else
            if(Pp .le. 1.27) then
               CS_S = 32.6*(Pp-0.8)**3.21
            else
               if(Pp .le. 8.) then
                  CS_S = 5.40*(Pp-0.8)**0.81
               else
                  if(Pp .le. 1.d3) then
                     CS_S = 32.0*dlog(Pp)+48.5/dsqrt(Pp)-59.5
                  else
                     CS_S = 163.*(s/1.876d3)**0.21 ! M.Mori (1997, ApJ)
                  endif
               endif
            endif
         endif
      endif

c ### ELECTRONS ###
      if(key .lt. 0) then
         if(Pp .le. Pth1) return   ! the threshold momentum for pp->pi^- X
         Mmes= Mpi1                ! CHARGED PION rest mass
         MX = 2.*Mp+Mmes           ! The mass in the channel X (=2Mp+Mmes)
c Dermer's fit
         if(Pp .le. 2.81) then
            CS_S = 2.33*(Pp-1.65)**1.2
         else
            if(Pp .le. 5.52) then
               CS_S = 0.32*Pp**2.1
            else
               CS_S = 28.2*dlog(Pp)+74.2/dsqrt(Pp)-69.3
            endif
         endif
      endif

c ### POSITRONS ###
      if(key .gt. 0) then
         if(Pp .le. Pth3) return   ! the lowest thres. momentum for pi^+ produc.
         Mmes= Mpi1                ! CHARGED PION rest mass
c DERMER; charged pions from pp->pi^+ +d
         MX = Md                   ! The mass in the channel X (=Md)
         Tp = dsqrt(Pp**2+Mp**2) - Mp
         if(Tp .le. 0.65) then
            eta = dsqrt((s-Mmes**2-MX**2)**2-4.d0*(Mmes*MX)**2)
     #             /dsqrt(s)/2.d0/Mmes
            CS_D = eta*( 0.18 +eta**2*(0.95 -eta**6*0.016) )
         else
            if(Tp .le. 1.43) then
               CS_D = 0.56/Tp**3.9
            else
               CS_D = 0.34/Tp**2.5
            endif
         endif
c Dermer's fit
         if(Pp .le. Pth2) return   ! the threshold momentum for pp->pi^+ X
         MX = Mp+Mn                ! The mass in the channel X (=Mp+Mn)
         if(Pp .le. 1.29) then
            eta = dsqrt((s-Mmes**2-MX**2)**2-4.d0*(Mmes*MX)**2)
     #            /dsqrt(s)/2.d0/Mmes
            CS_S = eta**4*( 0.95 +eta**2*(0.099 +eta**2*0.204) )
            if(Pp .gt. 0.95) CS_S = 0.67*eta**4.7 +0.3
         else
            if(Pp .le. 4.0) then
               CS_S = 22.0*(Pp-1.27)**0.15
            else
               CS_S = 27.0*dlog(Pp)+57.9/dsqrt(Pp)-40.9
            endif
         endif
      endif

c STEPHENS & BADHWAR's formalism for PIONs
      if(Pp .ge. P1) then
         gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s) ! max pion Lf in CMS
         if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
         betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)         ! max beta*gamma -"-
         EU = Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max) ! Upper limit
c integral over pion momentum in LS
         call SIM2(EU,Mmes,Mmes/9.,1.d-3,1.d-3,SB_EPI,B1)
         CS_SB = -B1 *2.*3.1415926
      endif

      if(iabs(key) .lt. 2 .or. iabs(key) .gt. 4) return
c STEPHENS & BADHWAR's formalism for KAONs
      Mmes= MK                     ! KAON rest mass
      if(key .lt. 0) then
         if(Pp .le. Pth4) return   ! the threshold momentum for pp->K^- X
         MX = 2.*Mp+Mmes           ! The mass in the channel X (=2Mp+Mmes)
      else
         if(Pp .le. Pth5) return   ! the threshold momentum for pp->K^+ X
         MX = Mp+Mn                ! The mass in the channel X (=2Mp+Mmes)
      endif
      gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s) ! max pion Lf in CMS
      if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
      betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)         ! max beta*gamma -"-
      EU = Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max) ! Upper limit
c integral over pion momentum in LS
      call SIM2(EU,Mmes,Mmes/9.,1.d-3,1.d-3,SB_EPI,B1)
      CS_K = -B1 *2.*3.1415926
      return
      end

      function D_PION(Epi)
c***********************************************************************
c PURPOSE: calculation of the positron spectrum [1/GeV] from p+p->pi^+ +d
c (for Pp, Ppi fixed).
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c
     3      /key/key, kaonkey, muonkey

      D_PION = 0.d0
      if(key .le. 0) return
      Mpi= Mpi1    ! Charged pion rest mass
      MX = Md      ! The mass in the channel X (=Md)
      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) !Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      Ppi = dsqrt(Epi**2-Mpi**2)

      gam_pi_cms = (s-MX**2+Mpi**2)/2.d0/dsqrt(s)/Mpi !Pion Lf in CMS
      if(gam_pi_cms .le. 1.d0) gam_pi_cms = 1.d0
      betgam_pi_cms= dsqrt(gam_pi_cms**2-1.d0)        !gamma*beta -"-
      Emu1 = (Mpi**2 +Mmu**2)/2.d0/Mpi   ! Muon energy in the pion rest mass
      Pmu1 = dsqrt(Emu1**2 -Mmu**2)      ! Muon momentum -"-
      gamL = dmax1((Emu1-Pmu1)*Pe/Mmu**2+(Emu1+Pmu1)/4.d0/Pe, ! L-limit PION Lf
     1       gam_cms*gam_pi_cms-betgam_cms*betgam_pi_cms)
      if(Pe .lt. (Emu1+Pmu1)/2.d0) gamL = 1.d0
      gamU = gam_cms*gam_pi_cms+betgam_cms*betgam_pi_cms      ! U-limit PION Lf
      if(gamL .ge. gamU) return
      muonkey = 1                                                ! IMOS20030909
      D_PION = MUON(Epi/Mpi)/(2.d0*Mpi*betgam_cms*betgam_pi_cms) ! *Mpi/Mpi
      return
      end

      function SB_EPI(Epi)
c***********************************************************************
c PURPOSE: calculation of the differential cross section of the pion/kaon
c production [mbarn/GeV] vs. LS meson energy from p+p->meson+X (for Pp fixed).
c Used for calculation of the total inclusive cross section (CS_SB & CS_K).
c Stephens & Badhwar formalism (1981, Ap.Spa.Sci. 76,213)
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c
     3      /key/key, kaonkey, muonkey
      external SB

      Mmes= Mpi0                  ! NEUTRAL PION rest mass
      if(key .ne. 0) Mmes= Mpi1   ! CHARGED PION rest mass
      if(kaonkey .eq. 1) Mmes= MK ! KAON rest mass
      Ppi = dsqrt(Epi**2-Mmes**2)
      call SIM1(0.d0,+1.d0,1.d-3,5.d-4,1.d-15,SB,AI1)
      call SIM1(0.d0,-1.d0,1.d-3,5.d-4,1.d-15,SB,AI2)
      SB_EPI = (AI1-AI2)  *Ppi                                !<==(*)
      return
      end

      function SB_MESON(Epi)
c***********************************************************************
c PURPOSE: calculates the product of the differential (vs. LS meson energy,
c Pp fixed) cross section of the pion/kaon production [mbarn/GeV], and the
c spectrum of photons or electrons/positrons [1/GeV] from pion/kaon decay.
c Stephens & Badhwar formalism (1981,Ap.Spa.Sci.76,213)
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c
     3      /key/key, kaonkey, muonkey
      external SB

      Mmes= Mpi0                       ! Neutral pion rest mass ## GAMMA-RAYS ##
      MX = 2.d0*Mp                     ! Mass in the channel X (=2Mp)
      if(key .ne. 0) then 
         MX = Mp+Mn                    !## POSITRONS ##
         Mmes= Mpi1                    ! Charged PION rest mass
         if(kaonkey .eq. 1) Mmes= MK   ! KAON rest mass
         if(key .lt. 0) MX= 2.*Mp+Mmes !## ELECTRONS ##
      endif

      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) ! Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          ! CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) ! CMS beta*gamma
      Ppi = dsqrt(Epi**2-Mmes**2)         ! meson momentum
      
      Epi_cms = (s-MX**2+Mmes**2)/2.d0/dsqrt(s) ! max pion E in CMS
      COSmax =-1.d0
      if(betgam_cms .gt. 0.d0 .and. Ppi .gt. 0.d0)
     1   COSmax =(gam_cms*Epi-Epi_cms)/betgam_cms/Ppi
      if(COSmax .lt. -1.d0) COSmax = -1.d0      
      AI = 0.d0
      if(COSmax .lt. 1.d0)
     #   call SIM1(1.d0,COSmax,1.d-3,5.d-4,1.d-15,SB,AI)
      SB_MESON= -AI  *2.d0                     ! gamma-rays (2 photons per pion)
c      if(key .ne. 0) SB_MESON= -AI*Ppi*MUON(Epi/Mmes) ! positrons/electrons
      muonkey = 0                                      ! IMOS20030909
      if(key .ne. 0) SB_MESON= -AI*Mmes*MUON(Epi/Mmes) ! IMOS20030909
      return
      end

      function SB(COSX)
c***********************************************************************
c Invariant cross section [Epi*(d^3 sigma/d^3 Ppi)]=[mbarn/...] for the
c inclusive PION/KAON production in pp-collisions at energies ~6 - 24 GeV
c Stephens & Badhwar (1981, Ap.Spa.Sci. 76,213)
c COSX = cos(theta) is the LS polar angle.
c key = 0 for gamma-rays; > 0 for positrons; < 0 for electrons.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c
     3      /key/key, kaonkey, muonkey
      
      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) !Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      SB = 0.d0
c inclusive cross section parameters
      if(key .eq. 0) then        ! GAMMA-RAYS
         Mmes= Mpi0              ! Neutral pion rest mass
         MX = 2.d0*Mp            ! Mass in the channel X (=2Mp)
         A  = 140.
         B  = 5.43
         R  = 2.
         C1 = 6.1
         C2 = 3.3
         C3 = 0.6
         Fp=(1.d0+23.d0/dsqrt(Pp**2+Mp**2)**2.6)*(1.d0-4.d0*Mp**2/s)**R
      else
         if(kaonkey .eq. 1) then ! KAON production
            Mmes= MK             ! KAON rest mass
            if(key .gt. 0) then  ! K+ (POSITRONS)
               MX = Mp+Mn        ! The mass in the channel X (=Mp+Mn)
               A = 8.85
               B = 4.05
               C = 2.5
            else                 ! K- (ELECTRONS)
               MX = 2.*Mp+Mmes   ! The mass in the channel X (=2Mp+Mmes)
               A = 9.3
               B = 3.8
               C = 8.3
            endif
         else                    ! charged PION production
            Mmes= Mpi1           ! Charged PION rest mass
            if(key .gt. 0) then  ! pi+ (POSITRONS)
               MX = Mp+Mn        ! The mass in the channel X (=Mp+Mn)
               A  = 153.
               B  = 5.55
               R  = 1.
               C1 = 5.3667
               C2 = 3.5
               C3 = 0.8334
            else                 ! pi- (ELECTRONS)
               MX = 2.*Mp+Mmes   ! The mass in the channel X (=2Mp+Mmes)
               A  = 127.
               B  = 5.3
               R  = 3.
               C1 = 7.0334
               C2 = 4.5
               C3 = 1.667
            endif
            Fp=(1.d0+4.d0*Mp**2/s)**(-R)
         endif
      endif

      PT= Ppi*dsqrt(1.d0-COSX**2)   !Transverse momentum of PION/KAON (inv) 
c Parallel CMS pion momentum in LS variables 
      PZ=gam_cms*Ppi*COSX -betgam_cms*dsqrt(Ppi**2+Mmes**2)
      X1=PZ*2.d0*dsqrt(s) /dsqrt((s-Mmes**2-MX**2)**2-4.d0*(Mmes*MX)**2)
      X = dsqrt( X1**2+4.d0/s*(PT**2+Mmes**2) )
      if(X .lt. 1.d0) then
         if(kaonkey .eq. 1) then  ! KAON production
            if(C*dlog10(1.d0-X)-0.44*B*PT.gt.-200.)            !IMOS20001124.1
     #         SB = A *(1.d0-X)**C *dexp(-B*PT)
         else                     ! PION production
            Q = (C1 -C2*PT +C3*PT**2)/dsqrt(1.d0+4.d0*Mp**2/s)
            if(Q*dlog10(1.d0-X)-0.44*B*PT/(1.d0+4.d0*Mp**2/s).gt.-200.) !IMOS20001124.2
     #         SB = A*Fp *(1.d0-X)**Q *dexp(-B*PT/(1.d0+4.d0*Mp**2/s))
         endif
      endif
      return
      end
                                
      function STECKER(MD)
c***********************************************************************
c spectrum of gamma-rays and electrons/positrons [1/GeV] from p+p->pi+X
c (for Pp, MD fixed). Stecker's formalism (1970, Ap.Spa.Sci. 6,377)
c MD [GeV/c^2]- Delta-isobar mass,Breit-Wigner distr.(M0 -mean mass,G -wigth)
c Egam [GeV] - gamma-ray energy; Pe [GeV] - electron/positron energy (massless);
c Pp [GeV/c] - proton beam momentum.
c key = 0 for gamma-rays; > 0 for positrons; < 0 for electrons.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md1,Mpi0,Mpi1,Mmu,MK
     1      /thres/Pth0,Pth1,Pth2,Pth3,Pth4,Pth5
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c
     3      /key/key, kaonkey, muonkey
      external MUON
      FF(x,y)=dlog( (x + dsqrt(x*x - 1.d0))/(y + dsqrt(y*y - 1.d0)) )
      
      STECKER = 0.
      if(Pp .le. Pth0) return ! the lowest threshold momentum for pp->pi+X

      A1 = 0.d0
      A2 = 0.d0
      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) !Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma

      gam_d_cms = (s+MD**2-Mp**2)/2./dsqrt(s)/MD   !Delta CMS Lf
      if(gam_d_cms .le. 1.d0) gam_d_cms = 1.d0
      betgam_d_cms = dsqrt(gam_d_cms**2-1.d0)      !Delta beta*gamma in CMS
      gam_d_ls1=gam_cms*gam_d_cms-betgam_cms*betgam_d_cms !Delta LS Lf back
      if(gam_d_ls1 .le. 1.d0) gam_d_ls1 = 1.d0
      gam_d_ls2=gam_cms*gam_d_cms+betgam_cms*betgam_d_cms !Delta LS Lf fwrd
      if(gam_d_ls2 .le. 1.d0) gam_d_ls2 = 1.d0
      betgam_d_ls1 = dsqrt(gam_d_ls1**2-1.d0)             !gamma*beta -"-
      betgam_d_ls2 = dsqrt(gam_d_ls2**2-1.d0)             !gamma*beta -"-

      if(key .eq. 0) then  ! --- GAMMA-RAYS ---
         Mpi= Mpi0    ! Neutral pion rest mass
         gam_pi_drs = (MD**2-Mp**2+Mpi**2)/2.d0/MD/Mpi !Pion Lf in Delta Rest Sys.
         if(gam_pi_drs .le. 1.d0) gam_pi_drs = 1.d0
         betgam_pi_drs= dsqrt(gam_pi_drs**2-1.d0)      !gamma*beta -"-
         gamL1 = dmax1(Egam/Mpi+Mpi/4.d0/Egam,    !Lower limit for PION Lf
     #           gam_d_ls1*gam_pi_drs-betgam_d_ls1*betgam_pi_drs)
         gamU1 = gam_d_ls1*gam_pi_drs+betgam_d_ls1*betgam_pi_drs !Upper lim.
         gamL2 = dmax1(Egam/Mpi+Mpi/4.d0/Egam,    !Lower limit for PION Lf
     #           gam_d_ls2*gam_pi_drs-betgam_d_ls2*betgam_pi_drs)
         gamU2 = gam_d_ls2*gam_pi_drs+betgam_d_ls2*betgam_pi_drs !Upper lim.
         if(gamL1 .lt. gamU1)   !Factor 2 is from 2 photons per pion
     #      A1 = 2.d0/(betgam_d_ls1*betgam_pi_drs) *FF(gamU1,gamL1)
         if(gamL2 .lt. gamU2)   !Factor 2 is from 2 photons per pion
     #      A2 = 2.d0/(betgam_d_ls2*betgam_pi_drs) *FF(gamU2,gamL2)
      else                 ! --- POSITRONS/ELECTRONS ---
         if(key.lt.0 .and. Pp.le.Pth1) return !the threshold momentum, pp->pi^- X
         if(key.gt.0 .and. Pp.le.Pth2) return !the threshold momentum, pp->pi^+ X
         Mpi= Mpi1    ! Charged pion rest mass
         Emu1 = (Mpi**2 +Mmu**2)/2.d0/Mpi   ! Muon energy in the pion rest mass
         Pmu1 = dsqrt(Emu1**2 -Mmu**2)      ! Muon momentum in the pion rest mass
         gam_pi_drs = (MD**2-Mp**2+Mpi**2)/2.d0/MD/Mpi !Pion Lf in Delta Rest Sys.
         if(gam_pi_drs .le. 1.d0) gam_pi_drs = 1.d0
         betgam_pi_drs= dsqrt(gam_pi_drs**2-1.d0)      !gamma*beta -"-

         gamL1=dmax1(Pe/(Emu1+Pmu1)+(Emu1+Pmu1)/4.d0/Pe, !Lo-limit PION Lf
     1               gam_d_ls1*gam_pi_drs-betgam_d_ls1*betgam_pi_drs)
         if(Pe .lt. (Emu1+Pmu1)/2.d0) 
     1        gamL1 = gam_d_ls1*gam_pi_drs-betgam_d_ls1*betgam_pi_drs
         gamU1=gam_d_ls1*gam_pi_drs+betgam_d_ls1*betgam_pi_drs !Up-limit PION Lf
         gamL2=dmax1(Pe/(Emu1+Pmu1)+(Emu1+Pmu1)/4.d0/Pe, !Lo-limit PION Lf
     1               gam_d_ls2*gam_pi_drs-betgam_d_ls2*betgam_pi_drs)
         if(Pe .lt. (Emu1+Pmu1)/2.d0) 
     1        gamL2 = gam_d_ls2*gam_pi_drs-betgam_d_ls2*betgam_pi_drs
         gamU2=gam_d_ls2*gam_pi_drs+betgam_d_ls2*betgam_pi_drs !Up-limit PION Lf
         muonkey = 1                                   ! IMOS20030909
         if(gamL1 .lt. gamU1) then
            call SIM1(gamU1,gamL1,1.d-3,5.d-4,1.d-15,MUON,AI)
            A1 = -AI*Mpi/(betgam_d_ls1*betgam_pi_drs)          
         endif
         if(gamL2 .lt. gamU2) then
            call SIM1(gamU2,gamL2,1.d-3,5.d-4,1.d-15,MUON,AI)
            A2 = -AI*Mpi/(betgam_d_ls2*betgam_pi_drs)
         endif
      endif
      BW = G/((MD-M0)**2+G*G)       !Delta mass spectrum (Breit-Wigner)
     #     /( atan( (dsqrt(s)-Mp-M0)/G )-atan( (Mp+Mpi-M0)/G ) )
c      STECKER = (A1+A2)/4.d0/Mpi        
      STECKER = (A1+A2)*BW/4.d0/Mpi        
      return
      end

      SUBROUTINE SIM1(A1,B1,H1,REPS1,AEPS1,FU,AI)
c***********************************************************************
c calculation the definite integral by Simpson's method with the automatic
c choice of the step of integration 
C INPUT: A1,B1 - the limits of integration; H1 - the initial step;
C REPS1,AEPS1 - the relative and absolute precision; FU - the name of the 
C user function f(x); OUTPUT: AI - the value of the integral;
C AIH - the value of integral with one more step of integration;
C AIABS - the value of the integral for module of the integrand;
C # NOTES # the subprogram returns the value of integral as one of the
C precise conditions (AEPS1,EPS1) are reached; when AEPS1=EPS1=0, 
c then it is calculated with the constant step H1.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      IMPLICIT real*8 (A-H,O-Z), integer (I-N)
      DIMENSION F(7),P(5)
      H=dSIGN(H1,B1-A1)
      S=dSIGN(1.d0,H) 
      A=A1
      B=B1
      AI=0.d0
      AIH=0.d0
      AIABS=0.d0
      P(2)=4.d0
      P(4)=4.d0
      P(3)=2.d0
      P(5)=1.d0
      IF(B-A) 1,2,1
    1 REPS=dABS(REPS1)
      AEPS=dABS(AEPS1)
      DO 3 K=1,7
    3    F(K)=1.d20
      X=A
      C=0.d0
      F(1)=FU(X)/3.d0
    4 X0=X
      IF((X0+4.d0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.d0
      IF(H) 7,2,7
    7 DO 8 K=2,7
    8    F(K)=1.d20
      C=1.d0
    5 DI2=F(1)
      DI3=dABS(F(1))
      DO 9 K=2,5
         X=X+H
         IF((X-B)*S) 23,24,24
   24    X=B
   23    IF(F(K)-1.d20) 10,11,10
   11    F(K)=FU(X)/3.d0
   10    DI2=DI2+P(K)*F(K)
    9    DI3=DI3+P(K)*dABS(F(K))
      DI1=(F(1)+4.d0*F(3)+F(5))*2.d0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=dABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=dABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.d0) 17,14,14
   17 H=H*2.d0
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
   19    F(K)=1.d20
      GOTO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=1.d20
      F(4)=1.d20
      F(6)=1.d20
      F(7)=1.d20
   18 DI1=DI2+(DI2-DI1)/15.d0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GOTO 22
   21 H=H/2.d0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=1.d20
      F(4)=1.d20
      X=X0
      C=0.d0
      GOTO 5
   22 IF(C) 2,4,2
    2 RETURN
      END

      SUBROUTINE SIM2(A1,B1,H1,REPS1,AEPS1,FU,AI)
c***********************************************************************
c calculation the definite integral by Simpson's method with the automatic
c choice of the step of integration 
C INPUT: A1,B1 - the limits of integration; H1 - the initial step;
C REPS1,AEPS1 - the relative and absolute precision; FU - the name of the 
C user function f(x); OUTPUT: AI - the value of the integral;
C AIH - the value of integral with one more step of integration;
C AIABS - the value of the integral for module of the integrand;
C # NOTES # the subprogram returns the value of integral as one of the
C precise conditions (AEPS1,EPS1) are reached; when AEPS1=EPS1=0, 
c then it is calculated with the constant step H1.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      IMPLICIT real*8 (A-H,O-Z), integer (I-N)
      DIMENSION F(7),P(5)
      H=dSIGN(H1,B1-A1)
      S=dSIGN(1.d0,H)
      A=A1
      B=B1
      AI=0.d0
      AIH=0.d0
      AIABS=0.d0
      P(2)=4.d0
      P(4)=4.d0
      P(3)=2.d0
      P(5)=1.d0
      IF(B-A) 1,2,1
    1 REPS=dABS(REPS1)
      AEPS=dABS(AEPS1)
      DO 3 K=1,7
    3    F(K)=1.d20
      X=A
      C=0.d0
      F(1)=FU(X)/3.d0
    4 X0=X
      IF((X0+4.d0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.d0
      IF(H) 7,2,7
    7 DO 8 K=2,7
    8    F(K)=1.d20
      C=1.d0
    5 DI2=F(1)
      DI3=dABS(F(1))
      DO 9 K=2,5
         X=X+H
         IF((X-B)*S) 23,24,24
   24    X=B
   23    IF(F(K)-1.d20) 10,11,10
   11    F(K)=FU(X)/3.d0
   10    DI2=DI2+P(K)*F(K)
    9    DI3=DI3+P(K)*dABS(F(K))
      DI1=(F(1)+4.d0*F(3)+F(5))*2.d0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=dABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=dABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.d0) 17,14,14
   17 H=H*2.d0
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
   19    F(K)=1.d20
      GOTO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=1.d20
      F(4)=1.d20
      F(6)=1.d20
      F(7)=1.d20
   18 DI1=DI2+(DI2-DI1)/15.d0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GOTO 22
   21 H=H/2.d0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=1.d20
      F(4)=1.d20
      X=X0
      C=0.d0
      GOTO 5
   22 IF(C) 2,4,2
    2 RETURN
      END

      BLOCK DATA
c***********************************************************************
c a common BLOCK DATA for PION production and decay code and for BREMSS. code
c MASSes, THRESHOLD momenta, and BRANCHING ratio; phi1 & phi2 for He-like atoms
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     1      /thres/Pth0,Pth1,Pth2,Pth3,Pth4,Pth5
     1      /branch/BR1,BR2
     2      /Hartree/ DD(11),PH1(11),PH2(11) !! DD = delta/(2.d0*Af)
      data
     #   M0   /1.232d0/,   ! GeV, Delta(1232)-isobar mean mass 
     #   G    /0.0575d0/,  ! GeV, Delta-isobar wight (0.115d0/2.)
     #   Mp   /0.938d0/,   ! GeV, Proton rest mass
     #   Mn   /0.9396d0/,  ! GeV, Neutron rest mass
     #   Md   /1.8756d0/,  ! GeV, Deutron rest mass
     #   Mpi0 /0.135d0/,   ! GeV, Neutral pion rest mass
     #   Mpi1 /0.1396d0/,  ! GeV, Charged pion rest mass
     #   Mmu  /0.10566d0/, ! GeV, Muon rest mass
     #   MK   /0.49365d0/, ! GeV, Kaon rest mass
c
     #   Pth0 /0.78/,      ! GeV/c, the threshold momentum for pp->pi^o (+2p)
     #   Pth1 /1.65/,      !(1.22)GeV/c, -"- for pp->pi^-(+2p+pi^+)
     #   Pth2 /0.80/,      ! GeV/c, -"- for pp->pi^+ (+p+n)
     #   Pth3 /0.791/,     ! GeV/c, -"- for pp->pi^+ (+d)
     #   Pth4 /3.302/,     ! GeV/c, -"- for pp->K^- (+2p+K^+)
     #   Pth5 /1.8332/,    ! GeV/c, -"- for pp->K^+ (+p+n)
c
     #   BR1  /0.635/,     ! branching ratio for channel K -> mu +nu
     #   BR2  /0.212/      ! branching ratio for channel K -> pi^(+/-) +pi^0

c He formfactors
      data DD/0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10./,
     #     PH1/134.60, 133.85, 133.11, 130.86, 127.17, 
     #         120.35, 104.60,  89.94,  74.19,  54.26,  40.94/,
     #     PH2/131.40, 130.51, 130.33, 129.26, 126.76,
     #         120.80, 105.21,  89.46,  73.03,  51.84,  37.24/ 
      end

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * antiproton.f *                                galprop package * 2001/05/11
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

      real*8 function ANTIPROTON(key,Pap1,Pp1,NZ1,NA1,NZ2,NA2)        ! IMOS20010511
c***********************************************************************
c             ###   I.Moskalenko   ###    version of 22.06.1998     ###
c Antiproton (+antineitron) production spectrum vs. momentum [barn c/GeV] for 
c pp-, pA-, Ap-, and AA-collisions (Pp1, Pap1 fixed) per 1 target nucleus/cm^3. 
c Refs: Moskalenko I.V. et al. 2002, ApJ 565, 280
c (Tan & Ng 1983, J.Phys.G:Nucl.Phys.9,227; ibid.,p.1289;
c Letaw et al.1983,ApJS,51,271; Gaisser & Schaeffer 1992,ApJ,394,174;
c Westfall et al.1979,PRC,19,1309)
c
c Pap1 [GeV/c] - secondary anti-p momentum; Pp1 [GeV/c] -beam momentum/nucleus
c NA1 & NA2 are the atomic numbers of beam and target nuclei, correspondingly
c (NA1=NA2=1 for pp-collisions)
c++ using nucleon_cs.f, pp_meson.f  yuanqiang 20090312
c***********************************************************************
      implicit real*8 (a-h,m,o-z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c
      external TAN_NG

      Pi = 3.1415926
      ANTIPROTON = 0.d0
      Pap = Pap1
      Pp = Pp1/NA1
c      Pp = Pp1  ! use if Pp1 is the beam momentum per nucleon
      s2 = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp)!square of the total CMS energy
      s1 = dsqrt(s2)                      !Total pp-energy in CMS
      if(s1 .le. 4.d0*Mp) return          !anti-p production threshold

      MX = 3.d0*Mp                        !pp->ppp(anti-p) channel
      gam_cms = s1/2.d0/Mp                !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      Eap = dsqrt(Pap**2+Mp*Mp)           !anti-p LS energy
c###
      if(NA2 .gt. 1) then                 ! IMOS20010215
         Eap = Eap+0.06                   ! gives more pbars at low energies
         Pap = dsqrt(Eap**2-Mp*Mp)        ! in case of nucleus target
      endif 
c###
      Eap_max = (s2-MX*MX+Mp*Mp)/2.d0/s1  !max anti-p CMS energy
      Ek = dsqrt(Pp*Pp+Mp*Mp)-Mp          !kinetic energy per nucleon
      COSmax =-1.d0
      if(betgam_cms*Pap .gt. 0.d0)
     1   COSmax =(gam_cms*Eap-Eap_max)/betgam_cms/Pap
      if(COSmax .lt. -1.d0) COSmax = -1.d0      
      AI = 0.d0
      if(COSmax .lt. 1.d0)
     #   call SIM1(1.d0,COSmax,1.d-3,1.d-4,1.d-40,TAN_NG,AI)
c dF/dEap - production spectrum vs. energy; factor 2 accounts for antineutrons
      ANTIPROTON = -AI*Pap  *2.d0  *2.d0*Pi*1.d-3  
     #   *Pap/Eap                         !tranformation to dF/dPap
      if(NA1*NA2 .eq. 1) return           !return if pp-collisions

c pA-, Ap-, and AA-collisions: 
      if(Ek .le. 0.3) return              !kinetic energy must be>0.3GeV
      call NUCLEON_CS(key,Ek,1,NZ2,NA2
     #       ,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann)  ! IMOS20010511 IMOS20000606
c nuclei: pA total inelastic cross section, mb (Letaw et al.1983,ApJS,51,271)
      CS2 = PP_inel
      if(NA2 .gt. 1) CS2 = PA_inel
      CS1 = PP_inel
      if(NA1 .gt. 1) then 
         call NUCLEON_CS(key,Ek,1,NZ1,NA1
     #       ,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann)  ! IMOS20010511 IMOS20000606
         CS1 = PA_inel
      endif
c multiplicity of anti-p in pA-,Ap-,AA-coll.; Gaisser&Schaeffer 1992,ApJ,394,174
c###
      MULT=1.2*  (NA1*CS2+NA2*CS1)/2./PP_inel   ! IMOS20010223 put factor 1.2
c###                                            ! from comparison w. Simon et al. 
      ANTIPROTON = ANTIPROTON*MULT

      return
      end

      real*8 function TAN_NG(COSX)
c***********************************************************************
c The invar. cross section of the inclusive anti-proton production 
c [mbarn c^3/GeV^2] from p+p->(anti-p)+X reaction (for Pp,(Pap,COSX) fixed);
c (Tan & Ng 1983,J.Phys.G:Nucl.Phys.9,227; ibid., p.1289)
c COSX = cos(theta) is the LS polar angle.
c***********************************************************************
      implicit real*8 (a-h,m,o-z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c

      TAN_NG = 0.d0
      s2 = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp)!square of the total CMS energy
      s1 = dsqrt(s2)                      !Total pp-energy in CMS
      MX = 3.d0*Mp                        !pp->ppp(anti-p) channel
      if(s1 .le. 4.d0*Mp) return          !anti-p production threshold

      gam_cms = s1/2.d0/Mp                !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      Eap = dsqrt(Pap**2+Mp*Mp)           !anti-p LS energy
      SINX = 0.d0
      if(1.d0-COSX**2 .gt. 0.d0) SINX = dsqrt(1.d0-COSX**2)
      PT = Pap*SINX                       !anti-p transverse momentum
c radial variable: XR=(E*)/(Emax*), where E*,Emax* - anti-p CMS energies
      XR=(gam_cms*Eap-betgam_cms*Pap*COSX)*2.d0*s1/(s2-MX*MX+Mp*Mp)
      if(XR .gt. 1.d0) return

c      goto 999
c >> parametrization used in Tan & Ng 1983,J.Phys.G:Nucl.Phys.9,1289 <<
c inclusive cross section at s^(1/2) >~ 10 GeV; pp.1296-8
888   A = 0.465*dexp(-3.70d-2*XR)+2.31*dexp(1.40d-2*XR)   ! c/GeV
      B = 3.02d-2*dexp(-3.19*(XR+0.399))*(XR+0.399)**8.39 !(c/GeV)**2
      F = (3.15-1.05d-4)*(1.d0-XR)**7.90                  !mbarn c^3/GeV^2
      if(XR .le. 0.5) F = F+1.05d-4*dexp(-10.1*XR)
      TAN_NG = F*dexp(-(A*PT+B*PT*PT))
      if(s1 .ge. 1.d1) return

c inclusive cross section at s^(1/2) <~ 10 GeV; pp.1303-5
      XT =((gam_cms*Eap-betgam_cms*Pap*COSX)-Mp)
     #   /((s2-MX*MX+Mp*Mp)/2.d0/s1-Mp)   !XT=(T*)/(Tmax*)
      Q = s1-4.d0*Mp                      !4Mp = anti-p production threshold
      A = 0.306*dexp(-0.120*XT)
      B = 0.0552*dexp(2.72*XT)
      C = 0.758-0.680*XT+1.54*XT*XT
      D = 0.594*dexp(2.87*XT)
      delta = 1.-dexp(-(1.-dexp(-A*Q**B))*dexp(C*Q-D))
      TAN_NG = TAN_NG/delta
      return

c ----------------------------------------------------------------------
c >> parametrization used in Tan & Ng 1983,J.Phys.G:Nucl.Phys.9,227 <<
c inclusive cross section at s^(1/2) >~ 10 GeV
999   A = 3.95*dexp(-2.76*XR)             ! c/GeV
      B = 40.5*dexp(-3.21*XR)*XR**2.13    !(c/GeV)**2
      F = 2.10*(1.d0-XR)**7.80            !mbarn c^3/GeV^2
      if(XR .le. 0.5) F = F+3.34*dexp(-17.6*XR)
      TAN_NG = F*dexp(-(A*PT+B*PT*PT))
      if(s1 .ge. 1.d1) return

c inclusive cross section at s^(1/2) <~ 10 GeV
      DXR =XR -dsqrt(PT**2+Mp**2) *2.d0*s1/(s2-MX*MX+Mp*Mp) !DXR=XR-XRmin
      Q = s1-4.d0*Mp                      !4Mp = anti-p production threshold
      delta = 6.25d-3*(dexp(-0.592*Q)+493.*dexp(-5.40*Q))
     #   *(dexp(6.08+2.57*DXR+7.95*DXR**2)-1.d0)*dexp(3.00*DXR*(3.09-Q))
      TAN_NG = (delta+1.d0)*TAN_NG
      return
      end

!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
! * nucleon_cs.f *                   adapted from galprop package * 2008/10/06 
!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

! parametrization of the pp-, pA-, AA-, (anti_p)p-, and (anti_p)A total
! inelastic cross sections, as well as (anti_p) annihilation cross sect.
! A sample driver program is in the bottom of this file.
!    INPUT:
! option -pA total inelastic cross section 
!       =0 -[L83]; 
!       =1 -[WA96] for Zt>5 and [BP01] for Zt<=5;
!       =2 -[BP01];
! Ek -kinetic energy per nucleon of the beam momentum particles (GeV);
! Zp=+/-1 is the charge of beam and At is the atomic number of target nuclei
! (Zp=At=1 for pp-collisions, Zp= -1 for antiprotons);
!    OUTPUT: 
! PP_inel  [mbarn] is the proton-proton inelastic cross sect.;
! PA_inel  [mbarn] is the proton-nucleus inelastic cross sect. (p+A2);
! aPP_non  [mbarn] is the antiproton-proton inelastic cross sect.(nonannihil.);
! aPA_non  [mbarn] put =PA_inel, is the antiproton-nucleus inelastic cross sect.(nonan.);
! aPP_ann  [mbarn] is the antiproton-proton annihilation cross sect.,
! aPA_ann  [mbarn] is the antiproton-nucleus annihilation cross sect.,
!    REFERENCES:
! [W79] Westfall et al. 1979, PRC, 19, 1309
! [L83] Letaw et al. 1983, ApJS, 51, 271
! [TN83] Tan & Ng 1983, J.Phys.G:Nucl.Phys., 9, 227
! [PDG00] D.E.Groom et al. 2000, Europ. Phys. J. C15, 1 
! [MO97] Moiseev & Ormes 1997, Astropart. Phys.6, 379
! [WA96] Wellisch H.P., Axen D. 1996, PRC 54, 1329; Wellish 1999, private comm. 
!    (typo corrections & code)
! [BP01] V.S.Barashenkov, A.Polanski code used for calc. of the pA total inelastic cross sections
!++ using SIGHAD in crn6.f     yuanqiang 20090312
!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
      subroutine nucleon_cs(option,Ek,Zp,Zt,At,PP_inel,PA_inel,aPP_non,
     &aPA_non,aPP_ann,aPA_ann)
      implicit real*8 (a-z)
      integer option,Zp,Zt,At,ISS

      PP_inel=0
      PA_inel=0
      aPP_non=0
      aPA_non=0
      aPP_ann=0
      aPA_ann=0
      if( Ek<= 0.) return
      aPP_inel=0
      aPA_inel=0
      Emev=Ek
      PZ=abs(1.*Zp)
      Em=1000.*Ek
      TZ=Zt
      TA=At
      ISS=2

!## Proton-Proton INELASTIC cross section, mb [TN83,p.234]
      if(Ek > 0.3) then
         U = log((Ek+Mp)/200.)
         PP_inel = 32.2*(1.+0.0273*U)
         if(U >= 0.) PP_inel = PP_inel+32.2*0.01*U*U
         if(Ek < 3.) PP_inel = PP_inel/(1.+0.00262/Ek**(17.9+
     &13.8*log(Ek)+4.41*(log(Ek))**2))
      endif
      if(Zp*At.eq.1) return

!## Proton-Nucleus INELASTIC cross section, mb
      if (option.eq.0) then
         if (At.eq.4) then ! a correction for He
            C1=0.8
         else
            C1=1
         endif
         if (At.eq.9) C1 = 1.+0.75*exp(-Ek/0.075) ! a correction for Be
         PA_inel = C1 *45. *TA**0.7*(1.+0.016*sin(5.3-2.63*log(TA)))
         if(Ek < 5.) PA_inel = PA_inel*(1.-0.62*exp(-Ek/0.2)*
     &sin(10.9/(Ek*1.e3)**0.28))
         if(At.eq.4) then ! pHe, Moskalenko
            if (Ek>0.01) then
               PA_inel=111.*(1.-(1.-sin(9.72*(log10(Ek*1000.))**0.319
     &-4.14))*exp(-3.84*(Ek-0.1)))
            else
               PA_inel=0.
            endif
         endif
      else if (option.eq.1) then ! Zt>5 [WA96], Zt<=5 [BP01]
      if(Zt>5) then
         b0 = 2.247-0.915*(1.-TA**(-1./3.));
         Fcorr = (1.-0.15*exp(-Emev))/(1.-0.0007*At); ! high-energy correction
         if (At-Zt>1.5) then
            rN = log(TA-Zt)
         else
            rN = 1.
         endif
         s0 = Pi*10.*1.36*1.36*Fcorr*rN*(1.+TA**(1./3.)-
     &b0*(1.-TA**(-1./3.)))
         p1 = 8.-8./At-0.008*At
         p2 = 2.*(1.17-2.7/At-0.0014*At)
         p3 = 0.8+18./At-0.002*At
         p4 = 5.6-0.016*At
         p5 = 1.37*(1.+1./At)
         x = log10(Emev)
         f1 = 1./(1.+exp(-p1*(x+p2))) ! low-energy return to zero
         f2 = 1. +p3 *( 1. -1./(1.+exp(-p4*(x+p5)))) ! low-energy rise
         PA_inel = f1*f2*s0
      endif
      else if (option.eq.2) then ! [BP01]
         if (Em<14.)  Em=14.
         if (Em>1.e6) Em=1.e6
         PA_inel = sighad(ISS,PZ,PZ,TA,TZ,Em) ! IMOS20020502
      endif
      if(Zp*At >= 1) return

!## AntiProton-Proton ANNIHILATION cross section [TN83]
      if(Ek < 10.) then
         aPP_ann=661.*(1.+0.0115/Ek**0.774-0.948*Ek**0.0151) ! 0.1GeV<Ek<12GeV
      else
! assuming aPP_ann = aPP_tot -PP_tot (i.e., aPP_elast = PP_elast); (aPP_tot-PP_tot) from [PDG00]
         s2 = 2.*Mp*(Ek+2*Mp) ! square of the total CMS energy
         aPP_ann = 2*35.43/s2**0.560
      endif

!## AntiProton-Proton TOTAL INELASTIC cross section
      aPP_inel = PP_inel + aPP_ann;
      if(Ek <= 14.) then
         aPP_inel = 24.7*(1.+0.584/Ek**0.115+0.856/Ek**0.566)
         if(aPP_ann > aPP_inel) aPP_ann = aPP_inel
      endif

!## AntiProton-Proton TOTAL INELASTIC NON-ANNIHILATION cross section
      aPP_non = aPP_inel - aPP_ann
      if(aPP_non < 0.) aPP_non = 0.

!## AntiProton-NUCLEUS cross sections
      if(At > 1) then
!## AntiProton-NUCLEUS TOTAL INELASTIC NON-ANNIHILATION cross section
         aPA_non = PA_inel
!## AntiProton-NUCLEUS ANNIHILATION cross section on 12C-nucleus [mb] (0.4<Pp<300) [MO97]
         A = At
         Z = Zt ! Z = 0.59*pow(A,.927);  for Z > 2 nuclei
         if(At.eq.4) then ! Modified to agree with HE p-He cs / imos
            Z = 2.
            A = 3.30
         endif
         aPA_ann = A**(2./3.)*(48.2 +19./Ek**0.55 ! Scaling to other nuclei
     &+(0.1-0.18/Ek**1.2)*Z+0.0012/Ek**1.5*Z*Z)-aPA_non
         if(aPA_ann < 0.) aPA_ann = 0.
         if(aPA_ann < aPP_ann) aPA_ann = aPP_ann
         if(At.eq.4.and.Ek>5.) aPA_ann = aPP_ann
      endif
      end
!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!/**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!/ * crn6.f *                                      galprop package * 2001/05/11
!/**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

!/**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
*              Joint Institute for Nuclear Research
*        Laboratory of Computing Techniques and Automation
*	            V.S.Barashenkov,A.Polanski
*
*                        Electronic Guide
*             	  for Nuclear Cross Sections
*
*   A short description of a fortran CROSEC code providing the integral cross
*-sections for pion-nucleus,nucleon-nucleus and nucleus-nucleus interactions
*(total,nonelastic,elastic).The hadron-nucleus cross-sections are obtained,by
*means of interpolation between evaluated experimental data at target mass 
*numbers A>=4 and energies from 14(20)MeV up to 1 TeV.The nucleus-nucleus cross
*sections are calculated with the help of approximation formula with fitted
*coefficient at energies above several Mev/nucleon.The CROSEC code uses 168 KB
*memory and can be used in interactive way to generate separate values or ta-
*bles of cross sections on a display screen or file.Fractions of the CROSEC
*code can be used as subroutines employed by other codes.A set of all available
*experimental integral hadron-nucleus cross-sections at energies exceeding 14
*MeV and plots of evaluated total,nonelastic and elastic cross sections for pio
*and nucleon interactions are presented in/3-5/.Two methods have been employed
*to calculate dependence of cross-sections vs energy.At high energies,where
*the projectile de Brogle wave length is significantly smaller then the size
*of the target nucleus quasioptical model is used.The parameters of the
*model have been fitted to obtain best agreement of calculated and experimental
*data.The high-energy region has been divided into separate intervals with
*characteristic behavior of cross sections.(For example,the region near the
*minimum of nucleon cross-sections at energy about 200 MeV - the resonance
*region.In case of pion-nucleus cross sections,the interval of smooth cross
*section alterations at energy above 1 GeV).A set of parameters has been
*defined for each interval.Phenomenological approximation of cross sections
*was used at lower energies/3-5/.
*  The known experimental information on nucleus-nucleus cross-sections is
*insufficient to compile a detailed plots curves,especially if one considers
*the great number of particular interesting nuclei pairs.In this case approxi-
*mation relations can be used with coefficients fitted by means of comparison
*with the known experimental data/5/.Of course,the accuracy of the results is
*lower then that for hadron-nucleus interactions.
*   The CROSEC code includes following modules:
** main module-reading input and printing separate values or tables of hadron
*  -nucleus and nucleus-nucleus cross sections(mbarns)
** FUNCTION SIGHAD-calculation of pion-nucleus and nucleon-nucleus cross
*  -sections.
** FUNCTION SIGION-calculation of nucleus-nucleus cross-sections
** FUNCTION FHS-   calculation of high-energy nucleus-nucleus cross
*  -sections
** FUNCTION FC- calculation of low-energy parameters
** BLOCK DATA A-target nucleus mass numbers
** BLOCK DATA B-projectile kinetic energies
** BLOCK DATA C-nucleus-nucleus cross section parameters
** BLOCK DATA S1-neutron-nucleus cross-sections
** BLOCK DATA S2-proton-nucleus cross-sections
** BLOCK DATA S3-pi-meson-nucleus cross sections
** BLOCK DATA S4-pi+meson-nucleus cross-sections
** Integer BLOCKS B,S1,S2,S3,S4 there are on the BARPOL.FOR or on the BARPOL.DAT 
*After starting the program you must insert the eight numbers in free format:
*ITYPE
*PA_PZ_TA_TZ_E1_E2_ES
*representing respectively:
*ITYPE=1 OR ITYPE=2 - Calculation of total or nonelastic cross-sections,
*ITYPE=3            - Calculation at once of total,nonelastic and elastic cross
*                     sections,
* ITYPE=-1            -Read cross sections from BARPOL.DAT
* ITYPE=0             -Exit
*PA,PZ-projectile mass and charge numbers(for pions PA<0.2)
*TA,TZ-the same for target nucleus,
*E1,E2-the lowest and the highest energies in the considered interval,
*ES   -energy step,
*   For example if one needs to obtain cross sections for pion(pi-)+nucleus
*Pb-207.19 in the energy range from 20 MeV to 200 MeV with a step 10 Mev one
*must enter(in free format):
*3 
*0.15 -1 207.19 82 20 200 10
*
*The SIGHAD and SIGION functions can be used separately with others codes
*to calculate a current value of hadron-nucleus or nucleus-nucleus cross
*-section:
*CS=SIGHAD(IS,PA,PZ,TA,TZ,E)
*CS=SIGION(IS,PA,PZ,TA,TZ,E)
*where:E-projectile energy in MeV (for hadrons) or in MeV/nucleon(for nucleus)
*IS=1 or IS=2 -calculation of total or nonelastic cross-sections,
*PA,PZ,TA,TZ-explained above.
*   The CROSEC code is written in FORTRAN.
*   Further information can be requested at e-mail addresses:
*   barashen@lcta30.jinr.dubna.su
*   polanski@cyf.gov.pl
*   olek@neutron.kth.se
*References:
*1.Barashenkov V.S.,Polanski A.,Comm. JINR E2-94-417,Dubna, 1994.
*2.Barashenkov V.S.Interaction cross-sections of particles and nuclei
*  JINR.Dubna,1993.(In Russian)
*3.Barashenkov V.S.Gareeva G.F.,Polanski A.,Comm. JINR 10-92-214,Dubna 1992.(In Russ.)
*4.Barashenkov V.S.,Polanski A.,Comm. JINR B2-90-489,Dubna, 1990.(In Russ.)
*5.V.S.Barashenkov, A. Polanski, A.N. Sosnin. Comm. JINR P2- 90-159,Dubna 1990.(In Russ.)
*6.V.S.Barashenkov, A. Polanski, A.N. Sosnin.Comm. JINR P2-89-753, Dubna 1989.(In Russ.)

C     PROGRAM CROSEC(INPUT,OUTPUT)
      subroutine CROSEC

      IMPLICIT REAL*8 (A-H,O-Z)
 5000 FORMAT(' ****************************************************' 
     *  /,   ' CODE FOR CALCULATION OF NUCLEON-NUCLEUS,PION-NUCLEUS'
     *  /,   ' AND NUCLEUS-NUCLEUS TOTAL,NONELASTIC AND ELASTIC    '
     *  /,   ' ++++++++++++++++ CROSS-SECTIONS(MBARNS)+++++++++++++'
     *  /,   ' ****************************************************'
     *  /,   ' WRITTEN BY VLADILEN S.BARASHENKOV-JINR-DUBNA(RUSSIA)' 
     *  /,   ' AND ALEKSANDER POLANSKI-SINS-SWIERK(POLAND)         '
     *  /,   ' ****************************************************'
     *  /,   ' *********** e mail:polanski@ipj.gov.pl ************')

C     Input parameters:
* ITYPE
* PA_PZ_TA_TZ_E1_E2_ES
* representing respectively:
* ITYPE=1 OR ITYPE=2 - Calculation of total or inelastic cross-sections,
* ITYPE=3            - Calculation at once of total,inelastic and elastic cross
*                      sections,
* ITYPE=-1            -Read cross sections from BARPOL.DAT
* ITYPE=0             -Exit
* PA,PZ-projectile mass and charge numbers(for pions PA<0.2)
* TA,TZ-the same for target nucleus,
* E1,E2-the lowest and the highest energies in the considered interval,
* ES   -energy step,
* To break the procedure one must insert ITYPE=0 or PA=0 
*    For example if one needs to obtain cross sections for proton+nucleus
* Pb-207.19 in the energy range from 20 MeV to 40 MeV with a step 2.5 MeV one
* must enter(in free format):
* 3 
* 1 1 207.19 82 20 40 2.5 
C     REGION OF APPLICABILITY OF THIS CODE:
C     FROM 14 MEV UP TO 1 TEV FOR NUCLEON-NUCLEUS COLISIONS
C     FROM 20 MEV UP TO 1 TEV FOR PION-NUCLEUS COLISIONS
C     FROM 1.0 MEV/NUCLEON  UP TO 1 TEV/NUCLEON FOR NUCLEUS-NUCLEUS
C     COLISIONS 
C     A1,A2<240
c     A2>=4
 6000 FORMAT(' ENTER:'/' ITYPE (ITYPE=0 exit,ITYPE>0 continue)')
     
 6001 FORMAT(' ENTER:'/' PA_PZ_TA_TZ_E1_E2_ES'/)
      ITYPE=-1
      CALL SIGTAP2(ITYPE)
      PRINT 5000
 100  PRINT 6000
      IN=5
      READ(IN,*) ITYPE
      IF(ITYPE.LT.0) CALL SIGTAP2(ITYPE)
      IF(ITYPE.EQ.0) GO TO 200
      IF(ITYPE.LT.0) GO TO 100
      PRINT 6001
      READ(IN,*) PA,PZ,TA,TZ,E1,E2,ES
c      WRITE(6,*) PA,PZ,TA,TZ,E1,E2,ES
      NE1=0
      ND=1
      ISS=ITYPE
      IF(ES.GT.0.0.AND.E2.GT.0.0) NE1=(E2-E1)/ES
      NE=1+ABS(NE1)
      IF(E1.GT.E2) ES=-ES
      IF(ITYPE.GT.2) ND=2
      IF(PA.EQ.1.) WRITE(6,1000)
      IF(PA.LT.1.) WRITE(6,2000)
      IF(PA.GT.1.) WRITE(6,3000)
      DO 10 IE=1,NE
      CR1=0.
      CR2=0.
      CR3=0.
      DO 20 I=1,ND
      IF(ITYPE.GT.2) ISS=I
      E=E1+(IE-1)*ES
      IF(PA-1.0) 11,11,12
   11 CS=SIGHAD(ISS,PA,PZ,TA,TZ,E)
      GO TO 13
   12 CS=SIGION(ISS,PA,PZ,TA,TZ,E)
   13 CONTINUE
      IF(ISS.EQ.1) CR1=CS
      IF(ISS.EQ.2) CR2=CS
   20 CONTINUE
      IF(ITYPE.GT.2) CR3=CR1-CR2
      WRITE(6,7000) E,CR1,CR2,CR3
   10 CONTINUE
      GO TO 100 
 1000 FORMAT('     NUCLEON NUCLEUS CROSS-SECTIONS (MBARNS) '/
     *        ' ENERGY(MEV)  TOTAL   NONELASTIC     ELASTIC ')
 2000 FORMAT('         PION-NUCLEUS CROSS-SECTIONS(MBARNS) '/
     *       ' ENERGY(MEV)   TOTAL   NONELASTIC     ELASTIC ')
 3000 FORMAT('         NUCLEUS-NUCLEUS CROSS-SECTIONS(MBARNS)'/
     *       ' ENERGY(MEV/NUC) TOTAL   NONELASTIC     ELASTIC ')
 7000 FORMAT(4(F8.1,4X))
 200  CONTINUE
      END


      REAL*8 FUNCTION BINT(U,E,F,N,IS)
      IMPLICIT REAL*8 (A-H,O-Z)
C     LINEAR INTERPOLATION  IS=1
C     QUADRATIC INTERPOLATION  IS=2
      DIMENSION E(N),F(N)
      IF(IS.LE.0.or.n.eq.1) BINT=U
      IF(IS.LE.0.or.n.eq.1) RETURN
      IF(N.GT.2) GO TO 10
      X1=E(1)
      Y1=F(1)
      X2=E(2)
      Y2=F(2)
      GO TO 8
 10   CONTINUE     
      IF(U-E(1))1,1,2
  1   X1=E(1)
      Y1=F(1)
      X2=E(2)
      Y2=F(2)
      X3=E(3)      
      Y3=F(3)  
      GO TO 7
  2   IF(U-E(N-1)) 3,4,4
 4    X2=E(N-1)
      X3=E(N)      
      Y2=F(N-1)
      Y3=F(N)
      X1=E(N-2)
      Y1=F(N-2)  
      GO TO 7
 3    CONTINUE
      IF(N.LE.2) GO TO 7
      N1=N-1
      DO 5 J=2,N1
      IF(U-E(J)) 6,5,5
  6   X1=E(J-1)
      X2=E(J)
      X3=E(J+1)
      Y1=F(J-1)
      Y2=F(J)
      Y3= F(J+1)
      GO TO 7
 5    CONTINUE
 7    CONTINUE
      IF(IS.NE.2.OR.N.EQ.2)  GOTO 8
       BINT=Y1*(((U-X2)*(U-X3))/((X1-X2)*(X1-X3)))+
     *      Y2*(((U-X1)*(U-X3))/((X2-X1)*(X2-X3)))+
     *      Y3*(((U-X1)*(U-X2))/((X3-X1)*(X3-X2)))
  8   CONTINUE
      IF(IS.EQ.1.OR.N.EQ.2) BINT=Y1+(U-X1)*(Y2-Y1)/(X2-X1)
      RETURN
      END


      REAL*8 FUNCTION SIGION(ISS,A1,Z1,A2,Z2,T)
      IMPLICIT REAL*8 (A-H,O-Z)
C     FOR CALCULATION OF NUCLEUS-NUCLEUS TOTAL (ISS=1)
C     AND INELASTIC (ISS=2) CROSS SECTIONS
C     A1,Z1 - PROJECTILE MASS AND CHARGE NUMBERS (A1>1)
C     A2,Z2 - THE SAME FOR TARGET NUCLEUS (3<A2<240)
C     T - LAB. KINETIC ENERGY OF PROGECTALE (1 MEV/NUCLEON< T <1 TEV/NUCLEON)

      COMMON/CX/CX(38)
      COMMON /FH/AMP,AMT,AP,AT,B0,R0
C      WRITE(*,*)' SIGION',ISS,A1,Z1,A2,Z2,T
      IS=3-ISS
      IF(A1.LT.1.0D0.OR.A1.GT.240.0D0.OR.A2.LT.3.0D0.OR.A2.GT.240.0D0)
     *GO TO 101
      IF(DABS(Z1).LT.1.0D0) GO TO 101
      IF(T.LT.1.0D0) GO TO 101
      SIGION=0.0D0
      TP=T/A1
      AP=A1**0.333333
      AT=A2**0.333333
      AMP=A1*930.63D0
      AMT=A2*930.63D0
C     PARAMETER FOR CALCULATION OF NUCLEAR RADIUS
      R0=1.4D0
      IF(DABS(A1-4.0D0) .LT. 0.1D0) R0=1.3D0
      B0=1.44D0*Z1*Z2
      I=1
       IF(IS.EQ.2) I=20
C     SELECTION OF PROGECTALES
C     HEVY ION
      IF(A1.GT.4.1D0) N=I
C     ALFA,HELION,TRITON
      IF(A1.GT.2.1D0 .AND. A1.LT.4.1D0) N=I+6
C     DEUTRON
      IF(A1.LT.2.1D0) N=I+12
C     HIGH-ENERGY CROSS-SECTION
C     SELECTION OF PROGECTALE ENERGY
      IF(TP.LT.CX(N+1)) K=2
      IF(TP.LT.CX(N+4)) K=5
C     CROSS-SACTION PARAMRTERS
      C=CX(I)
      IF(TP.LT.CX(N+1)) C=CX(N+K)+CX(N+K+1)*DLOG10(TP)
      CP=CX(N+5)+CX(N+6)
      IF(TP.LT.10.) GO TO 1
C     HIGH-ENERGY CROSS-SECTION
      SIGION=FHS(IS,T,C)
      RETURN
C     CALCULATION OF LOW-ENERGY CROSS-SECTION
C     NORMALUSED HIGH-ENERGY CROSS-SECTION
 1     SH10=FHS(IS,10.*A1,CP)
      R0=1.45D0
      IF(DABS(A1-4.0D0).LT. 0.1D0) R0=1.4D0
C     RENORMALUSED COULOMB BARRIER
      B=B0/R0/(AP+AT)
C     LOW-ENERGY CROSS-SECTION
      SIGION=SH10*FC(T,B)/FC(10.*A1,B)
      IF(SIGION) 101,100,100
  101 CONTINUE
C      WRITE(*,1001)
C      WRITE(*,*)' ISS',ISS,' A1',A1,' Z1',Z1,' A2',A2,' Z2',Z2,' T',T
C      PAUSE
      SIGION=1.0D-07
C 1001 FORMAT(' ERROR IN INPUT OF PARAMETERS OF FUNCTION SIGION')
  100 CONTINUE
      RETURN
      END
      REAL*8 FUNCTION FC(T,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /FH/AMP,AMT,AP,AT,B0,R0
C     CMS ENERGY
      TC=T*AMT/(AMP+AMT)
      X=(TC-B)/1.2
      IF(X.GT.5) GO TO 1
      D=1.+EXP(X)
      FC=DLOG(D)/TC
      RETURN
  1   FC=X/TC
      RETURN
      END


      REAL*8 FUNCTION FHS(IS,E,C)
      IMPLICIT REAL*8 (A-H,O-Z)
C     CALCULATION OF HIGH-ENERGY TOTAL (IS=2) AND
C     INELASTIC (IS=1) CROSS-SECTIONS

C      E - LAB. KINETIC ENERGY OF PROGECTALE(MEV)
      COMMON/FH/AMP,AMT,AP,AT,B0,R0
C     SQUED PROGECTALE CMS MOMENTUM
      PPC=AMT*AMT*E*(E+2.*AMP)/((AMP+AMT)**2+2.*AMT*E)
C     DE BROGLE WAVE LANGTH
      AL=1.41*140./SQRT(PPC)
      EC=SQRT(PPC+AMP*AMP)-AMP
C     COULOMB BARRIER
      B=B0/R0/(AP+AT+AL)
      FHS=31.416*1.21*(1.-B/EC)*(AP+AT+1.85*AP*AT/(AP+AT)
     *+AL-C)**2*IS
      RETURN
      END


      BLOCK DATA C
      IMPLICIT REAL*8 (A-H,O-Z)
C     NUCLEUS-NUCLEUS CROSS-SECTION PARAMETERS

      COMMON/CX/CX(38)
      DATA  CX/2.07,560.,0.8,0.426,
     *           100.,-2.05,1.9,
     *           200.,0.07,0.87,
     *           20.,-1.55,2.1,
     *           700.,-1.01,1.08,
     *           400.,-0.59,0.94,
     *        2.45,225.,-2.25,2.,
     *           100.,-4.61,3.18,
     *           185.,-3.,2.4,
     *           185.,-3.,2.4,
     *           185.,-4.77,3.18,
     *           185.,-4.77,3.18/
      END


      REAL*8 FUNCTION SIGHAD(IS,A1,Z1,A2,Z2,T)
      IMPLICIT REAL*8 (A-H,O-Z)
C     CODE FOR CALCULATION OF NUCLEON-NUCLEUS AND PION-NUCLEUS
C     TOTAL(IS=1) AND NONELASTIC(IS=2)CROSS-SECTIONS(MBARNS) 
C     IS=0 EXIT
C     A1=1.0 FOR NUCLEON OR 0<A1<0.2 FOR PIONS
C     Z1-PROJECTILE CHARGE NUMBER
C     A2,Z2- TARGET NUCLEUS MASS AND CHARGE NUMBERS (4.0<=A2<=239.0)
C     T- PROJECTILE PARTICLE KINETIC ENERGY(MEV;14(20)MEV<T<1TEV)
C     AAI1(K),AAI2(K)- TARGET MASS NUMBERS for nucleons and pions
C     IENER(L,JS)-PROJECTILE ENERGY(MEV)
C     CROSS-SECTIONS ARE STORED IN ISIG1(K,L),ISIG2(K,L),ISIG3(K,L),ISIG4(K,L) 
C      ISIG(JS,K,L)- JS=1 NEUTRON CROSS SECTIONS  
C                    JS=2 PROTON CROSS SECTIONS 
C                    JS=3PI- MESON CROSS SECTIONS  
C                    JS=4PI+ MESON CROSS SECTIONS 
C      K- VS. TARGET MASS NUMBERS
C      L- VS. ENERGY
      COMMON /BARPO1/ AAI1(24),AAI2(20)
      COMMON/BARPO/ NEL(2),NE(4),IENER(62,4),ISIG(4,24,62)
      DIMENSION SS(24),EE(3),FF(3),SIG(2),JSIG(2)
      SIGHAD=0.0
      IF(IS.LE.0) RETURN
      IF(A1.LE.0.OR.A1.GT.1.0) GO TO 101
      IF(A1.LT.1.0.AND.A1.GT.0.2) GO TO 101  
      IF(A2.LT.1.0.OR.A2.GT.250.) GO TO 103 
      IF(ABS(Z1).GT.1.0) GO TO 101
      IF(Z1.EQ.0.0.and.T.LT.10.0) GO TO 102
      IF(T.LT.1.0) go to 102
      NS=4
C     NE(JS) - NUMBER OF ENERGY POINTS (62 for protons,53 for neutrons, 43 FOR PIONS-+)
C     NEL(IN) - NUMBER OF TARGETS (23 FOR NUCLEONS AND 19 FOR PIONS)
      IF(Z1.EQ.0.0.AND.A1.EQ.1.) JS=1
      IF(Z1.EQ.1.0.AND.A1.EQ.1.) JS=2
      IF(Z1.LT.0.0.AND.A1.LE.0.2) JS=3
      IF(Z1.EQ.1.0.AND.A1.LE.0.2) JS=4 
      IN=1
      IF(JS.GT.2)IN=2
C     FLAG JS=1  CROSS SECTIONS FOR NEUTRONS 
C     FLAG JS=2  CROSS SECTIONS FOR PROTONS 
C     FLAG JS=3  CROSS SECTIONS FOR MESONS PI- 
C     FLAG JS=4  CROSS SECTIONS FOR MESONS PI+
      IS1=1
      IS2=2
c      write(6,*) NEL(1),NEL(2),NE(1),NE(2),NE(3),NE(4),JS
      IF(T.GE.1.E6) T=1.E6
C     SELECTION OF ENERGY IENER(LP,JS)<T<IENER(LP+1,JS)
      DO 80 L=1,NE(JS)
      TI=IENER(L,JS)
      IF(TI-T) 80,90,90
  80  CONTINUE
  90  LP=L-1
      IF(L.LE.2) LP=1
      IF(L.GE.NE(JS))LP=NE(JS)-1
      NPAR=1
      NEE=3
C     NPAR=1 FOR NEUTRON,PROTON,PI-,P+.
C     NPAR=2 FOR MESON PI-0
      IF(Z1.EQ.0.0.AND.A1.LE.0.2) NPAR=2
      DO 55 IPAR=1,NPAR
      ZP1=Z1
      IF(IPAR.EQ.1.AND.NPAR.EQ.2) ZP1=-1.
      IF(IPAR.EQ.2.AND.NPAR.EQ.2) ZP1=1.
      DO 50 LI=1,NEE
      L1=LP-1+LI 
      EE(LI)=IENER(L1,JS)
      DO 30 K=1,NEL(IN)
C     JSIG(1)-TOTAL CROSS SECTIONS
C     JSIG(2)-NONELASTIC CROSS SECTIONS
      JSIG(1)=ISIG(JS,K,L1)/10000
      JSIG(2)=ISIG(JS,K,L1)-10000*JSIG(1)
 30   SS(K)=JSIG(IS)      
C     QUADRATIC INTERPOLATION OF CROSS SECTONS VS TARGET MASS
      IF(IN.EQ.1) FF(LI)=BINT(A2,AAI1,SS,NEL(IN),IS2)
      IF(IN.EQ.2) FF(LI)=BINT(A2,AAI2,SS,NEL(IN),IS2)
c      IF(IN.EQ.1) write(6,*)A2,AAI1,SS,NEL(IN),IS2 
c      IF(IN.EQ.2) write(6,*)A2,AAI2,SS,NEL(IN),IS2
  50  CONTINUE  
C     LINEAR INTERPOLATION OF CROSS SECTIONS VS ENERGY
 55   SIG(IPAR)=BINT(T,EE,FF,NEE,IS1)
c      write(6,*) T,EE,FF,NEE,IS1 
      IF(NPAR.EQ.1) SIGHAD=SIG(1)
      IF(NPAR.EQ.2) SIGHAD=(SIG(1)+SIG(2))/2
c      write(6,*) ' sig(1),sig(2) ',sig(1),sig(2)
      IF(SIGHAD) 101,100,100
  101 SIGHAD=0.0
c      WRITE(6,1001)A1,Z1,A2,T,sig(1),sig(2)
      GO TO 100
  102 WRITE(6,1002)A1,Z1,A2,T,sig(1),sig(2)
      GO TO 100
  103 WRITE(6,1003)A1,Z1,A2,T,sig(1),sig(2)
 1001 FORMAT(3F5.1,3E10.4,' ERROR IN INPUT (A1,Z1) OF FUNCTION SIGHAD')
 1002 FORMAT(3F5.1,3E10.4,' ERROR IN INPUT(Proj.ene)OF FUNCTION SIGHAD')
 1003 FORMAT(3F5.1,3E10.4,' ERROR IN INPUT (A2,Z2) OF FUNCTION SIGHAD')
  100 CONTINUE
      RETURN
      END


      BLOCK DATA ASIG
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /BARPO1/ AAI1(24),AAI2(20)     
C     AAI1- for nucleus, AAI2- for pions TARGET NUCLEUS MASS NUMBERS

C    1   H, D, He,  Li,  Be,  C,   N,   O,     Na,  Al,   S,   Ca, 
c       Ti,   Fe,   Cu,   Br,   Mo,  Cd,   Sn,   Ba,     W,     Pb, 
c        U,  Cf      
c
c    3   He,  Be,  C,   N,    O,   Na,  Al,  Ca,  Ti,  Fe,  Cu,  Br, 
c         Mo,  Cd,    Sn,    Ba,     W,    Pb,     U,     Cf          
      DATA AAI1/
     *1.0,2.0,4.0,6.90,9.01,12.00,14.00,16.00,23.00,26.98,32.08,40.08,  
     *47.90,55.85,63.50,79.90,95.94,112.40,118.7,137.34,183.85, 207.19,
     *238.03,250.0/
      DATA AAI2/
     *4.0,9.01,12.0,14.0,16.0,23.0,26.98,40.08,47.90,55.85,63.55,79.90,
     *95.94, 112.40, 118.69, 137.34,183.85,207.19,238.03,250.0/
      END


      SUBROUTINE SIGTAP2 (ITYPE)
      IMPLICIT REAL*8 (A-H,O-Z)
C     READ CROSS SECTIONS FROM BARPOL.DAT
C 2  4 55 62 43 43 24 20 Barashenkov and Polanski 
C (1)Neutrons,(2)Protons (3)Pi- (4)Pi+ Cross Sections
      COMMON /BARPO1/ AAI1(24),AAI2(20)
      COMMON/BARPO/ NEL(2),NE(4),IENER(62,4),ISIG(4,24,62)
      CHARACTER TITLE(60)
      DIMENSION ISIGN(24),ISIGE(24)
      IF(ITYPE.GE.0) RETURN
C     N=2 - LOGICAL UNIT OF TAPE
C     NS=4
C     NE(J) - NUMBER OF ENERGY POINTS 
C     NEL(IN) - NUMBER OF TARGETS
C     IN=1 FOR NUCLEONS
C     IN=2 FOR PIONS 
C     IENER(L,J) THE ENERGY BEAMS IN MEV
C     ISIG(J,K,L)- CROSS SECTIONS
C     (J=1) CROSS SECTIONS FOR NEUTRONS 
C     (J=2) CROSS SECTIONS FOR PROTONS
C     (J=3) CROSS SECTIONS FOR PI- 
C     (J=4) CROSS SECTIONS FOR PI+
C     K- VS. TARGET MASS NUMBERS
C     L- VS. ENERGY
C
      
      OPEN(2,FILE='barpol.dat',STATUS='OLD',FORM='FORMATTED')  
c      OPEN(8,FILE='BARPOL8.DAT',STATUS='NEW',FORM='FORMATTED')
C------------------------------------------------------------------------------------
      INP1=2
      INP2=8
      REWIND INP1       
      IX0=0
      READ(INP1,3001) N,NS,(NE(JS),JS=1,4),(NEL(K),K=1,2),
     1 (TITLE(I),I=1,60)
c      WRITE(INP2,3001)N,NS,(NE(JS),JS=1,4),(NEL(K),K=1,2),
c     1 (TITLE(I),I=1,60) 
      DO 170 JS=1,NS
      IN=1
      IF(JS.GT.2)IN=2
      READ(INP1,3000) IDUM 
c      WRITE(INP2,3000) JS
      DO 170 L=1,NE(JS)
      IF(IN.EQ.1) READ(INP1,5100)IENE,(ISIGE(K),K=1,NEL(IN)),IDUM,IEN,
     1 (ISIGN(K),K=1,NEL(IN)) 
      IF(IN.EQ.2) READ(INP1,4300)IENE,(ISIGE(K),K=1,NEL(IN)),IDUM,IEN,
     1 (ISIGN(K),K=1,NEL(IN))
c      WRITE(INP2,3030)IENE,(ISIGE(K),K=1,NEL(IN)),JS,IEN,
c     1 (ISIGN(K),K=1,NEL(IN))
      IENER(L,JS)=IENE
      IF(L.GT.(NE(JS)-3)) IENER(L,JS)=IENE*1000
      DO 170 K=1,NEL(IN)
      ISIG(JS,K,L)=(ISIGN(K)+ ISIGE(K))*10000+ISIGN(K)
  170 CONTINUE
C---------------------------------------------------------------------
      IF(ITYPE.NE.-2) GO TO 200
      WRITE(INP2,3001)N,NS,(NE(JS),JS=1,4),(NEL(K),K=1,2),
     1 (TITLE(I),I=1,60)     
      DO 180 JS=1,NS
      WRITE(INP2,1015) (IENER(L,JS),L=1,NE(JS))
      IN=1
      IF(JS.GT.2)IN=2
      DO 180 K=1,NEL(IN)
      WRITE(INP2,1015)(ISIG(JS,K,L),L=1,NE(JS)),K,JS 
  180 CONTINUE 
  200 CONTINUE
 3000 FORMAT(20I4)
 5100 FORMAT(51I5)  
 4300 FORMAT(43I5) 
 3030 FORMAT(52I5) 
 3020 FORMAT(A1,2I8)
 3001 FORMAT(8I3,60A1)
 1015 FORMAT(8I9)
      RETURN
      END 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
