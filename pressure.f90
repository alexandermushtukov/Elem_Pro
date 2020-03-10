!=========================================================================================
! The function calculates the energy of photon after the scattering in non-magnetized case.
!=========================================================================================
real*8 function ComptonEf(Ei,beta,theta_i,theta_f,ksi_f)
implicit none
real*8,intent(in)::Ei,beta,theta_i,theta_f,ksi_f
real*8::cos_a,gamma
  gamma=1.d0/sqrt(1.d0-beta**2)
  cos_a=sin(theta_i)*sin(theta_f)*cos(ksi_f)+cos(theta_i)*cos(theta_f)
  ComptonEf=Ei*(1.d0-beta*cos(theta_i)) / (Ei/gamma*(1.d0-cos_a)+(1.d0-beta*cos(theta_f)))
return
end function ComptonEf


real*8 function ComptPressureM(Ei,beta,theta_i)
implicit none
real*8,intent(in)::Ei,beta,theta_i
real*8::res,gamma
  gamma=1.d0/sqrt(1.d0-beta**2)
  res=0.5d0*( gamma*beta*log(2.*gamma*Ei*(1.d0-cos(theta_i))+1.d0) +gamma*(cos(theta_i)-1.d0)/(1.d0-beta*cos(theta_i))*&
              ((gamma*Ei*(1.d0-cos(theta_i))+1.d0)/gamma/Ei/(1.-cos(theta_i))*log(2.*gamma*Ei*(1.-cos(theta_i))+1.) -2.))
  ComptPressureM=res
return
end function ComptPressureM




!================================================================================================
! The function calculates the avarage momentum of the photon after the scattering by electron gas
! of temperature TkeV and velocity beta_0 (beta_0<0 - gas moves towards the NS surface).
! 1d relativistic Maxwellian distribution is taken to describe distribution of electrons.
! Ei - photon energy in units mc^2
! SBD: limits of integration should be reconsidered.
!================================================================================================
real*8 function ComptPressure_Max(Ei,theta_i,beta_0,TkeV,s_i)
implicit none
real*8,intent(in)::Ei,theta_i,beta_0,TkeV
integer,intent(in)::s_i
real*8::res,mas
dimension mas(5)
real*8::int_simpson_SD
  mas(1)=Ei
  mas(2)=theta_i
  mas(3)=beta_0
  mas(4)=TkeV/511.d0
  mas(5)=s_i*1.d0
  res=int_simpson_SD(f_int,-40.d0,40.d0,1.d-2,mas,5)
  ComptPressure_Max=res
return
contains
  real*8 function f_int(Z,mas,n_mas)
  implicit none
  real*8,intent(in)::Z,mas
  integer,intent(in)::n_mas
  dimension mas(n_mas)
  real*8::beta,Ei,theta_i,beta_0,T
  real*8::MaxwellRel1d_gen,ComptPressure01,ComptPressure_mag00,MaxwellRel_1d_BF  !==functions==!
  integer::s_i
    beta=sign(1.d0,Z)*sqrt(1.d0-1.d0/(1.d0+Z**2))
    Ei=mas(1)
    theta_i=mas(2)
    beta_0=mas(3)
    T=mas(4)
    s_i=int(mas(5))
    !f_int=MaxwellRel1d_gen(Z,T,beta_0)*ComptPressure01(Ei,beta,theta_i)
    !f_int=MaxwellRel1d_gen(Z,T,beta_0)*ComptPressure_mag00(Ei,beta,theta_i,s_i)
    f_int=MaxwellRel_1d_BF(Z,T,beta_0)*ComptPressure01(Ei,beta,theta_i)
    !f_int=MaxwellRel_1d_BF(Z,T,beta_0)*ComptPressure_mag00(Ei,beta,theta_i,s_i)
  return
  end function f_int
end function ComptPressure_Max


subroutine Test_ComptPressure_Max()
implicit none
real*8::Ei,theta_i,beta_0,TkeV,res
integer::s_i
real*8::ComptPressure_Max,ComptPressure01,ComptPressure_mag00  !==functions==!
  beta_0=-0.d0
  Ei=1.d0/511
  theta_i=0.d0
  TkeV=5.d0
  do while(Ei.le.(200./511))
    res=ComptPressure01(Ei,beta_0,theta_i)
    !res=ComptPressure_mag00(Ei,beta_0,theta_i,1)
    !res=ComptPressure_Max(Ei,theta_i,beta_0,TkeV,1)
    write(*,*)Ei,res
    Ei=Ei*1.05d0
  end do
return
end subroutine Test_ComptPressure_Max
!================================================================================================



!================================================================================================
! Non-magnetic case.
! The function calculates the average photon momentum after the scattering by electron moving with 
! velocity beta.
! Ei - photon energy in units mc^2
! beta - electron velosity (beta<0 - electron moves towards NS surface)
! SBD: the number of inervals in integation process should be choosen authomatically.
!================================================================================================
real*8 function ComptPressure01(Ei,beta,theta_i)
implicit none
real*8,intent(in)::Ei,beta,theta_i
real*8::gamma
real*8::pi=3.141592653589793d0
integer::i,j,ni,nj
real*8::dtheta,dksi,res,ksi_f,theta_f_,nu_i_,nu_f_
  nu_i_=(cos(theta_i)-beta)/(1.d0-beta*cos(theta_i))
  ni=2000
  nj=100
  dtheta=pi/(ni-1)
  dksi=2*pi/(nj-1)
  gamma=1.d0/sqrt(1.d0-beta**2)

  !res=CP00(Ei,0.8d0,nu_i_,cos(3.d0),0.d0)
  !write(*,*)res
  !read(*,*)

!write(*,*)0.01d0,0.d0,0.d0,pi/4,0.d0,CP00(0.01d0,0.d0,1.d0,cos(pi/4),0.d0)
!write(*,*)0.01d0,0.d0,0.d0,pi/3,0.d0,CP00(0.01d0,0.d0,1.d0,cos(pi/3),0.d0)
!write(*,*)0.01d0,-0.5d0,0.d0,pi/4,0.d0,CP00(0.01d0,-0.5d0,1.d0,cos(pi/4),0.d0)
!write(*,*)0.01d0,-0.5d0,0.d0,pi/3,0.d0,CP00(0.01d0,-0.5d0,1.d0,cos(pi/3),0.d0)
!read(*,*)

  !==integration over the solid angle==!
  res=0.d0
  i=1
  do while(i.lt.ni)
    theta_f_=(i-1)*dtheta
    nu_f_=cos(theta_f_)
    j=1
    do while(j.lt.nj)
      ksi_f=(j-1)*dksi
      res=res+dtheta*dksi*CP00(Ei,beta,nu_i_,nu_f_,ksi_f)   !==isotropic scattering in electron RF==!
      !res=res+dtheta*dksi*CP00_(Ei,beta,nu_i_,nu_f_,ksi_f)    !==non-isotropic scattering in electron RF==!
      write(*,*)theta_f_,ksi_f,CP00(Ei,beta,nu_i_,nu_f_,ksi_f)
      j=j+1
    end do
    i=i+1
  end do
  ComptPressure01=res*Ei/4/pi
return
contains
  !== The scattering in the elecctron reference frame is assumed to be isotropic ==!
  real*8 function CP00(Ei,beta,nu_i_,nu_f_,ksi_f)
  implicit none
  real*8,intent(in)::Ei,beta,ksi_f,nu_i_,nu_f_
  real*8::gamma,sin_i_,sin_f_,nu_i,sin_i,sin_f,nu_f,cos_a
    gamma=1.d0/sqrt(1.d0-beta**2)
    sin_i_=sqrt(1.d0-nu_i_**2)
    sin_f_=sqrt(1.d0-nu_f_**2)
    CP00=sin_f_*(beta+nu_f_)
    CP00=CP00/( (1.d0+beta*nu_f_)*(1.d0+gamma*Ei*(1.d0+beta*nu_f_))&
                -Ei/gamma*sin_i_*sin_f_*cos(ksi_f)-gamma*Ei*(nu_f_+beta)*(nu_i_+beta) )
  return
  end function CP00

  !== The scattering in the elecctron reference frame is assumed to be isotropic ==!
  real*8 function CP00_(Ei,beta,nu_i_,nu_f_,ksi_f)
  implicit none
  real*8,intent(in)::Ei,beta,ksi_f,nu_i_,nu_f_
  real*8::gamma,sin_i_,sin_f_,nu_i,sin_i,sin_f,nu_f,cos_a,cos_gamma
    gamma=1.d0/sqrt(1.d0-beta**2)
    sin_i_=sqrt(1.d0-nu_i_**2)
    sin_f_=sqrt(1.d0-nu_f_**2)
    CP00_=sin_f_*(beta+nu_f_)
    CP00_=CP00_/( (1.d0+beta*nu_f_)*(1.d0+gamma*Ei*(1.d0+beta*nu_f_))&
                -Ei/gamma*sin_i_*sin_f_*cos(ksi_f)-gamma*Ei*(nu_f_+beta)*(nu_i_+beta) )
    cos_gamma=sin_i_*sin_f_*(cos(ksi_f)+0.d0*sin(ksi_f))+nu_i_*nu_f_
    cos_gamma=max(-1.0,min(1.d0,cos_gamma))
    CP00_=CP00_*0.75d0*(1.d0+cos_gamma**2)
  return
  end function CP00_

end function ComptPressure01


subroutine Test_ComptPressure01()
implicit none
real*8::ComptPressure01  !==function==!
real*8::x,Ei,beta,theta_i
  Ei=0.01d0
  beta=0.d0
  theta_i=0.d0
  x=ComptPressure01(Ei,beta,theta_i)
  write(*,*)x
return
end subroutine Test_ComptPressure01





!================================================================================================
! Magnetic case.
! The function calculates the average photon momentum after the scattering by electron moving with
! velocity beta.
! B-field is taken into account in the conservation laws and in angular distrubution of photons
! after the scattering. We used Herols approximation of differential cross sections.
! Ei - photon energy in units mc^2
! beta - electron velosity (beta<0 - electron moves towards NS surface)
! theta_i - the angle between the direction of B-fiels and photon momentum (theta_i=0 - photon propagates away from the NS surface).
! s_i - polarization state of a photon
! SBD: Ecyc should be a parameter in this problem
!================================================================================================
real*8 function ComptPressure_mag00(Ei,beta,theta_i,s_i)
implicit none
real*8,intent(in)::Ei,beta,theta_i
integer,intent(in)::s_i
real*8::gamma
real*8::pi=3.141592653589793d0
integer::i,j,ni,nj
real*8::dtheta,dksi,res,ksi_f,theta_f_,nu_i_,nu_f_,nu_f,nu_i,Zi,Ecyc,theta_i_,SigmaCS,Ei_
real*8::SigmaTot,dSigmadOmega,dSigmadOmega_pl2  !==functions==!
  !open(unit = 23, file = './res/res1', status = 'old')
  nu_i=cos(theta_i)
  nu_i_=(nu_i-beta)/(1.d0-beta*nu_i)  ! directin of photon momentum in electron's reference frame
  theta_i_=acos(nu_i_)
  ni=100
  nj=50
  dtheta=pi/(ni-1)
  dksi=2*pi/(nj-1)
  gamma=1.d0/sqrt(1.d0-beta**2)
  Ei_=Ei*gamma*(1.d0-beta*nu_i)
  Zi=sqrt(gamma-1.d0)*sign(1.d0,beta)

  Ecyc=.1d0

  SigmaCS=SigmaTot(s_i,1,theta_i_,0.d0,Ei_,Ecyc)+SigmaTot(s_i,2,theta_i_,0.d0,Ei_,Ecyc)
  !==integration over the solid angle in the electron's reference frame==!
  res=0.d0
  i=1
  do while(i.lt.ni)
    theta_f_=(i-1)*dtheta                 !== the angle in the electron's reference frame==!
    nu_f_=cos(theta_f_)                   !== cosine in electron's reference frame==!
    nu_f =(nu_f_+beta)/(1.d0+beta*nu_f_)  !== cosine in observer's reference frame==!
    j=1
    do while(j.lt.nj)
      !write(*,*)"sss"
      ksi_f=(j-1)*dksi
      !==non-isotropic case==!
      res=res+dtheta*dksi*sin(theta_f_)&
          *(dSigmadOmega_pl2(s_i,1,Ei_,Ecyc,theta_i_,theta_f_,0.d0,ksi_f)&
           +dSigmadOmega_pl2(s_i,2,Ei_,Ecyc,theta_i_,theta_f_,0.d0,ksi_f))&
          /SigmaCS&
          *kf_00(Ei,theta_i,acos(nu_f),Zi)*nu_f
      !==isotropic case==!
      !res=res+dtheta*dksi*sin(theta_f_)&
      !        *kf_00(Ei,theta_i,acos(nu_f),Zi)*nu_f
      j=j+1
    end do
    i=i+1
  end do
  ComptPressure_mag00=res/4/pi
return
contains
  !==the photon momentum after the ground-to-ground state scattering in a strong magnetic field==!
  real*8 function kf_00(ki,theta_i,theta_f,Zi)
  implicit none
  real*8,intent(in)::ki,theta_i,theta_f,Zi
  real*8::Ei,E_T,Z_T,Z_f,Ef,Y_f,dkfdzi
  integer::det
    call calc(2,0,Zi,0.d0,ki,theta_i,0.d0,0,Z_f,Y_f,kf_00,theta_f,0.d0,dkfdzi,1.d0,det)
    Ei=sqrt(1.d0+Zi**2)
    E_T=ki+Ei
    Z_T=ki*cos(theta_i)+Zi
    Z_f=Z_T-kf_00*cos(theta_f)
    Ef=sqrt(1.d0+Z_f**2)
  return
  end function kf_00
end function ComptPressure_mag00



!================================================================================================
! The subroutine calculates the average photon momentum and averaged photon energy
! after the scattering by electron moving with velocity beta.
! B-field is taken into account in the conservation laws and in angular distrubution of photons
! after the scattering. We used Herols approximation of differential cross sections.
! Ei - photon energy in units mc^2
! beta - electron velosity (beta<0 - electron moves towards NS surface)
! theta_i - the angle between the direction of B-fiels and photon momentum
! (theta_i=0 - photon propagates away from the NS surface).
! s_i - polarization state of a photon.
! SBD: add non-zero temperature; add dispersion measure to the problem.
!!! B-field should be a parameter in the problem.
!================================================================================================
subroutine Compt_Ave_mag00(Ecyc,Ei,beta,theta_i,s_i,p_ave,E_ave)
implicit none
real*8,intent(in)::Ecyc,Ei,beta,theta_i
integer,intent(in)::s_i
real*8::gamma
real*8::pi=3.141592653589793d0
integer::i,j,ni,nj
real*8::dtheta,dksi,res,ksi_f,theta_f_,nu_i_,nu_f_,nu_f,nu_i,Zi,theta_i_,SigmaCS,Ei_
real*8::p_ave,E_ave
real*8::SigmaTot,dSigmadOmega,dSigmadOmega_pl2  !==functions==!
  !open(unit = 23, file = './res/res1', status = 'old')
  nu_i=cos(theta_i)
  nu_i_=(nu_i-beta)/(1.d0-beta*nu_i)  !==directin of photon momentum in electron's reference frame==!
  theta_i_=acos(nu_i_)
  ni=100
  nj=50
  dtheta=pi/(ni-1)
  dksi=2*pi/(nj-1)
  gamma=1.d0/sqrt(1.d0-beta**2)
  Ei_=Ei*gamma*(1.d0-beta*nu_i)
  Zi=sqrt(gamma-1.d0)*sign(1.d0,beta)

  !Ecyc=1.d0

  SigmaCS=SigmaTot(s_i,1,theta_i_,0.d0,Ei_,Ecyc)+SigmaTot(s_i,2,theta_i_,0.d0,Ei_,Ecyc)
  !==integration over the solid angle in the electron's reference frame==!
  p_ave=0.d0
  E_ave=0.d0
  i=1
  do while(i.lt.ni)
    theta_f_=(i-1)*dtheta                 !== the angle in the electron's reference frame==!
    nu_f_=cos(theta_f_)                   !== cosine in electron's reference frame==!
    nu_f =(nu_f_+beta)/(1.d0+beta*nu_f_)  !== cosine in observer's reference frame==!
    j=1
    do while(j.lt.nj)
      !write(*,*)"sss"
      ksi_f=(j-1)*dksi
      p_ave=p_ave+dtheta*dksi*sin(theta_f_)&
                  *(dSigmadOmega_pl2(s_i,1,Ei_,Ecyc,theta_i_,theta_f_,0.d0,ksi_f)&
                   +dSigmadOmega_pl2(s_i,2,Ei_,Ecyc,theta_i_,theta_f_,0.d0,ksi_f))&
                 /SigmaCS&
                 *kf_00(Ei,theta_i,acos(nu_f),Zi)*nu_f
      E_ave=E_ave+dtheta*dksi*sin(theta_f_)&
                  *(dSigmadOmega_pl2(s_i,1,Ei_,Ecyc,theta_i_,theta_f_,0.d0,ksi_f)&
                   +dSigmadOmega_pl2(s_i,2,Ei_,Ecyc,theta_i_,theta_f_,0.d0,ksi_f))&
                  /SigmaCS&
                  *( kf_00(Ei,theta_i,acos(nu_f),Zi))
      j=j+1
    end do
    i=i+1
  end do
  p_ave=p_ave/4/pi
return
contains
  !==the photon momentum after the ground-to-ground state scattering in a strong magnetic field==!
  real*8 function kf_00(ki,theta_i,theta_f,Zi)
  implicit none
  real*8,intent(in)::ki,theta_i,theta_f,Zi
  real*8::Ei,E_T,Z_T,Z_f,Ef,Y_f,dkfdzi
  integer::det
    call calc(2,0,Zi,0.d0,ki,theta_i,0.d0,0,Z_f,Y_f,kf_00,theta_f,0.d0,dkfdzi,1.d0,det)
    Ei=sqrt(1.d0+Zi**2)
    E_T=ki+Ei
    Z_T=ki*cos(theta_i)+Zi
    Z_f=Z_T-kf_00*cos(theta_f)
    Ef=sqrt(1.d0+Z_f**2)
  return
  end function kf_00
end subroutine Compt_Ave_mag00


subroutine Test_Compt_Ave()
implicit none
real*8::E_i,beta,theta_i,p_ave,E_ave
integer::s_i
  E_i=0.02d0
  beta=-0.d0
  theta_i=0.d0
  s_i=1
  do while(E_i.le.1.d0)
    call Compt_Ave_mag00(0.1d0,E_i,beta,theta_i,s_i,p_ave,E_ave)
     write(*,*)511*E_i,p_ave,511*E_ave
     E_i=E_i*1.03d0
  end do
return
end subroutine Test_Compt_Ave
!=================================================================================================




!=================================================================================================
! ???
!=================================================================================================
subroutine ResonansBreak()
implicit none
real*8::theta
real*8::int_simpson_new  !==function==!
real*8::res
  theta=0.d0
  do while(theta.lt.1.57d0)
    res=int_simpson_new(I_beamed,0.d0,theta,1.d-2)/int_simpson_new(I_beamed,0.d0,1.57d0,1.d-2)
    write(*,*)theta,res
    theta=theta+0.02d0
  end do
return
contains
  real*8 function I_beamed(t)
  implicit none
  real*8,intent(in)::t
  real*8::A,E2Ecyc
    A=0.05d0
    E2Ecyc=0.1d0
    I_beamed=cos(t)/(A*cos(t)+sin(t)*sin(t)+E2Ecyc**2) *cos(t)*sin(t)  !==the last 2 because of the geometrical factors==!
  return
  end function I_beamed
end subroutine ResonansBreak

