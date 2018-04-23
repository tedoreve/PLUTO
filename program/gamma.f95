program main
    real*8 Eelec1(2000),Felec1(2000),Eg1(100),Fg1(100),Ebk1(100),Fbk1(100),Ep1(100),Fp1(100)
    character(len = 512) :: cfile
    character :: xx
    integer nn,mm,ii,jj,kk,key1
    real*8 T1,u1,B1,nISM
    !real, allocatable, dimension(:) :: ISM
    !real, allocatable, dimension(:) :: B
    !real, dimension(2) :: ISM
    
	!!!!!!!!!!!!!!!!!!!!!!!!!set the density and the magnetic field in every pixel!!!!!!!!!!
    !open(12,file="rho.txt")
    !    read(12,*)xx,mm
    !    allocate(ISM(mm))
    !    read(12,*)(ISM(ii),ii=1,mm)
    !close(12)

    !open(12,file="B.txt")
    !    read(12,*)xx,nn
    !    allocate(B(nn))
    !    read(12,*)(B(ii),ii=1,nn)
    !close(12)
    !!!!!!!!!!!!!!!!!!!!!!!!!!set the constants for various mechanism!!!!!!!!!!!!!!!!!!!!!!!
    T1=2.7d0
    u1=2.5d-1 
    key1=1
    !write(*,*)ISM
!do  kk=1,mm
!do  jj=1,nn
    
    nISM=1.0    
    B1=1.0e-6

    nelec1=2000    
    do ii=1,nelec1
       Eelec1(ii)=(1.d0)*1.d1**(ii/200.d0)
       Felec1(ii)=1.d2*nISM*Eelec1(ii)**(-2)*exp(-Eelec1(ii)/30000.d0)
    end do
    
    ng1=100
    do ii=1,ng1
       Eg1(ii)=(1.d-15)*1.d1**(ii/4.d0)
       !Eg1(ii)=1
       Fg1(ii)=0
    end do

    np1=100
    do ii=1,np1
       Ep1(ii)=(1.d0)*(10d0**(ii/10.d0))
       Fp1(ii)=1.d3*nISM*(Ep1(ii)**(-2))*exp(-Ep1(ii)/30000.d0)
    end do

    nbk1=100
    do ii=1,nbk1
       Ebk1(ii)=1
       Fbk1(ii)=1
    end do
    call bkg_rad(Ebk1,Fbk1,nbk1,T1,u1)

!write(cfile,*)kk,jj


Open(12,File="Eelec.txt")
Do ii = 1,2000
   write(12,*)Eelec1(ii)
End Do
Close(12)

Open(12,File="Felec.txt")
Do ii = 1,2000
   write(12,*)Felec1(ii)
End Do
Close(12)

Open(12,File="Ep.txt")
Do ii = 1,100
   write(12,*)Ep1(ii)
End Do
Close(12)

Open(12,File="Fp.txt")
Do ii = 1,100
   write(12,*)Fp1(ii)
End Do
Close(12)

Open(12,File="Ebk.txt")
Do ii = 1,100
   write(12,*)Ebk1(ii)
End Do
Close(12)

Open(12,File="Fbk.txt")
Do ii = 1,100
   write(12,*)Fbk1(ii)
End Do
Close(12)

Open(12,File="Eg.txt")
Do ii = 1,100
   write(12,*)Eg1(ii)
End Do
Close(12)

call RADIATION_SYN(Eelec1,Felec1,nelec1,B1,Eg1,Fg1,ng1)
!Open(12,File="Fsyn"//Trim(AdjustL(cfile))//".txt")
Open(12,File="Fsyn.txt")
Do ii = 1,100
   write(12,*)Fg1(ii)
End Do
!write(*,*)Fg1
Close(12)

call GAMMA_BREM(Eelec1,Felec1,nelec1,nISM,Eg1,Fg1,ng1)
Open(12,File="Fbrem.txt")
Do ii = 1,100
   write(12,*)Fg1(ii)
End Do
Close(12)

call SEC_PION_PP(Ep1,Fp1,np1,nISM,key1,Eg1,Fg1,ng1)
Open(12,File="Fpp.txt")
Do ii = 1,100
   write(12,*)Fg1(ii)
End Do
Close(12)

call GAMMA_IC(Eelec1,Felec1,nelec1,Ebk1,Fbk1,nbk1,Eg1,Fg1,ng1)
Open(12,File="Fic.txt")
Do ii = 1,100
   write(12,*)Fg1(ii)
End Do
Close(12)


!end do
!end do
end
