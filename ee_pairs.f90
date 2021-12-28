!====================================================================================
!====================================================================================
subroutine test_ee_pairs()
implicit none
real*8::b,n_e_ini30,T_keV
real*8::n_e30,n_p30,mu
real*8::lg_TkeV,lg_n_e_ini30,Q23,mu_keV
integer::n_sum_max,n_max_e,n_max_p,n_max
character(len=100):: file_name
real*8::Z_e_max(5000),Z_p_max(5000)

  201 format (7(es11.4,"   "),3(I6,"   "))
  202 format (1(A13),(es8.2))

  b=0.3d0
  write(file_name,202)"./res/res_nu_",b
  write(*,*)"# file_name=",file_name; write(*,*)
  open (unit = 20, file = file_name)
  write(20,*)"# b=",b
  write(20,*)"# Format: b,lg_TkeV,lg_n_e_ini30,mu_keV,Q23,Q23_app,res/Q23_app,n_sum_max"
  write(20,*)
  close(20)
  do while(b.le.1.d0)
    lg_TkeV=2.d0
    do while(lg_TkeV.le.3.d0)
      T_keV=10.d0**lg_TkeV
      lg_n_e_ini30=-5.d0

      do while(lg_n_e_ini30.le.2.d0)
        n_e_ini30=10.d0**lg_n_e_ini30
        call EE_pairs(b,n_e_ini30,T_keV,n_e30,n_p30,mu,n_max,n_max_e,n_max_p,Z_e_max,Z_p_max)
        open(unit = 20, file = file_name, status = 'old',form='formatted',position="append")
        write(*,201)b,lg_TkeV,lg_n_e_ini30,mu
        write(20,201)b,lg_TkeV,lg_n_e_ini30,mu!,n_max,n_max_e,n_max_p
        close(20)
        lg_n_e_ini30=lg_n_e_ini30+0.1d0
      end do
      lg_TkeV=lg_TkeV+0.02d0
    end do
    b=b+0.1d0
  end do
return
end subroutine test_ee_pairs


!====================================================================================================
! This module is used in the subroutne calculating number dencities of ee+pairs in non-magnetic case.
!====================================================================================================
module EE_pairs_help
  implicit none
  real*8::mu_mod,TkeV_mod
end module
!====================================================================================================



!=========================================================================================
! The module is used in the subroutne calculating number dencities of ee+pairs in B-field.
!=========================================================================================
module EE_pairs_mag_help
  implicit none
  real*8::mu_mod,TkeV_mod,b_mod
  integer::n_mod
end module
!=========================================================================================



!====================================================================================================================
! The subroutine calculates the maximal temperature in accretion column regulated by electron-positron pair creation.
! b - dimentionless magnetic field strength
! n_b30 - the number density of barions
! n the accurasy of temperature mesh.
!====================================================================================================================
subroutine T_max_PairInfluence(b,n_b30,n)
implicit none
real*8,intent(in)::b,n_b30
integer,intent(in)::n
real*8::TkeV_1,TkeV_2,TkeV_max,dT,TkeV,n_e30,n_p30,mu,beta_min,beta_max,dbeta,beta,gamma,n_30max
integer::i,n_max,n_max_e,n_max_p,n_beta
real*8::mas
dimension mas(n,2)
real*8::n_p30_max  !==function==!
real*8::Z_e_max(5000),Z_p_max(5000)
  TkeV_1=1.d2
  TkeV_2=2.d3
  dT=(TkeV_2-TkeV_1)/(n-1)
  i=1
  do while(i.le.n)
    TkeV=TkeV_1+(i-1)*dT
    call EE_pairs(b,n_b30,TkeV,n_e30,n_p30,mu,n_max,n_max_e,n_max_p,Z_e_max,Z_p_max)
    mas(i,1)=n_p30
    mas(i,2)=TkeV
    write(*,*)"# ",i,TkeV,n_p30
    i=i+1
  end do
  !==the data base is ready==!

  beta_min=0.1d0
  beta_max=0.7d0
  n_beta=100
  dbeta=(beta_max-beta_min)/(n_beta-1)
  i=1
  do while(i.le.n_beta)
    beta=beta_min+(i-1)*dbeta
    n_30max=n_p30_max(n_b30,beta)
    TkeV_max=find_H(n_30max,mas,n)
    write(*,*)b,n_b30,beta,TkeV_max
    i=i+1
  end do
return
contains

  real*8 function find_H(x,H,n)
  implicit none
  integer,intent(in)::n
  real*8,intent(in)::x,H
  dimension H(n,2)
  integer::i
  real*8::x_down,x_up,H_down,H_up
    x_down=-1313.d0; x_up=1313.d0
    H_down=-0.d0; H_up=0.d0
    ! going through the array of data
    i=1
    do while(i.le.n)
      if((H(i,1).gt.x_down).and.(H(i,1).le.x))then
        x_down=H(i,1); H_down=H(i,2)
      end if
      if((H(i,1).lt.x_up).and.(H(i,1).ge.x))then
        x_up=H(i,1); H_up=H(i,2)
      end if
      i=i+1
    end do

    if(x_down.eq.-1313.d0)then
      x_down=x_up; H_down=H_up
    end if
    if(x_up.eq.(1313.d0))then
      x_up=x_down; H_up=H_down
    end if
    if(x_up.ne.x_down)then
      find_H=H_down+(H_up-H_down)*(x-x_down)/(x_up-x_down)
    else
      find_H=H_down
    end if
  return
  end function find_H

end subroutine T_max_PairInfluence


!===================================================================================================
! The function calculates the maximal number density of positrons in accretion flow as a function 
! of local free-fall velosity and the number density of barions n_b30.
!===================================================================================================
real*8 function n_p30_max(n_b30,beta_ff)
implicit none
real*8,intent(in)::n_b30,beta_ff
real*8::gamma_ff
  gamma_ff=1.d0/sqrt(1.d0-beta_ff**2)
  n_p30_max=918.d0*n_b30*(gamma_ff-1.d0)
return
end function n_p30_max



!===================================================================================================
! The subroutine calculates the number dencities of electron and positron pairs 
! at a given temperature TkeV, B-field strength b, and initial number dencity of electrons n_e_ini30.
! Parameters: 
! n_e_ini30 - the initial nuber density of electrons
! TkeV - the temperature in keV.
! Finally we will get: n_e_ini30=n_e30-n_p30
! Chemical potential is in keV.
!===================================================================================================
subroutine EE_pairs(b,n_e_ini30,TkeV,n_e30,n_p30,mu,n_max,n_max_e,n_max_p,Z_e_max,Z_p_max)
implicit none
real*8,intent(in)::b,n_e_ini30,TkeV
real*8::n_e30,n_p30,mu1,mu2,mu
integer::det,n_max,n_max_e,n_max_p,i
real*8::Z_e_max(5000),Z_p_max(5000)

  n_max   = 0
  n_max_e = 0
  n_max_p = 0
  !==(1)getting the interval, where chemical potential is==!
  mu2=1.d-5
  mu1=mu2
  if(b.eq.0.d0)then
    call EE_pairs_mu_fix(TkeV,mu1,n_e30,n_p30)
  else
    call EE_pairs_mag_mu_fix(b,TkeV,mu1,n_e30,n_p30,n_max,n_max_e,n_max_p)
  end if

  if(n_e_ini30.ge.(n_e30-n_p30))then
    do while(n_e_ini30.ge.(n_e30-n_p30))
      mu1=mu1*2
      if(b.eq.0.d0)then
        call EE_pairs_mu_fix(TkeV,mu1,n_e30,n_p30)
      else
        call EE_pairs_mag_mu_fix(b,TkeV,mu1,n_e30,n_p30,n_max,n_max_e,n_max_p)
      end if
    end do
    mu2=mu1/2
  else
    do while(n_e_ini30.ge.(n_e30-n_p30))
      mu1=mu1/2
      if(b.eq.0.d0)then
        call EE_pairs_mu_fix(TkeV,mu1,n_e30,n_p30)
      else
        call EE_pairs_mag_mu_fix(b,TkeV,mu1,n_e30,n_p30,n_max,n_max_e,n_max_p)
      end if
    end do
    mu2=mu1*2
  end if
  !==(1)end==!

  if(mu1.ge.mu2)then
    mu=mu1
    mu1=mu2
    mu2=mu
  end if

  !=we know: the actual chemical potential is between mu2 and mu1==!
  !==(2)getting the actual chemical potential==!
  det=11
  do while(det.eq.11)
    mu=(mu1+mu2)/2
    !call EE_pairs_mu_fix(TkeV,mu,n_e30,n_p30)
    !call EE_pairs_mag_mu_fix(b,TkeV,mu,n_e30,n_p30)
    if(b.eq.0.d0)then
      call EE_pairs_mu_fix(TkeV,mu,n_e30,n_p30)
    else
      call EE_pairs_mag_mu_fix(b,TkeV,mu,n_e30,n_p30,n_max,n_max_e,n_max_p)
    end if
    if(n_e_ini30.ge.(n_e30-n_p30))then
      mu1=mu
    else
      mu2=mu
    end if
    if(((mu2-mu1)/mu1).le.1.d-3)then
      det=22
    end if
    !write(*,*)mu1,mu2
  end do
  !==(2)end: now we know the actual chemical potential==!

  !==(3)calculation of number densities==!
  if(b.eq.0.d0)then
    call EE_pairs_mu_fix(TkeV,mu1,n_e30,n_p30)
  else
    call EE_pairs_mag_mu_fix(b,TkeV,mu1,n_e30,n_p30,n_max,n_max_e,n_max_p)
    !==now we would like to get the maximal Z at the level==!
    i=0
    do while(i.le.n_max_e)
      Z_e_max(i+1)=get_Z_max_mag(b,TkeV,mu1,i,1,1.d-3)
      i=i+1
    end do
    i=0
    do while(i.le.n_max_p)
      Z_p_max(i+1)=get_Z_max_mag(b,TkeV,mu1,i,2,1.d-3)
      i=i+1
    end do
  end if
  !==(3)end==!

return
contains

  !==================================================================================================================
  ! Non-magnetic case.
  ! Subroutine calculates electron and positron number dencities for a fixed chemical potential and temperature
  !==================================================================================================================
  subroutine EE_pairs_mu_fix(TkeV,mu,n_e30,n_p30)
  use EE_pairs_help
  implicit none
  real*8,intent(in)::mu,TkeV
  real*8::n_e30,n_p30,eps,b_final
    mu_mod=mu
    TkeV_mod=TkeV
    eps=1.d-3
    call int_ImproperInt_simpson(Ie ,0.d0,1.d-13,eps,n_e30,b_final)
    n_e30=n_e30*3.56d23
    if(mu_mod.le.10.d0*TkeV_mod)then
      call int_ImproperInt_simpson(Ip ,0.d0,1.d-13,eps,n_p30,b_final)
      n_p30=n_p30*3.56d23
    else
      call int_ImproperInt_simpson(Ip2,0.d0,1.d-13,eps,n_p30,b_final)
      n_p30=n_p30*3.56d23/exp(mu_mod/TkeV_mod)
    end if
  return
  end subroutine EE_pairs_mu_fix


  !==================================================================================================================
  ! Magnetic case.
  ! Subroutine calculates electron and positron number dencities for a fixed chemical potential and temperature.
  ! n_max_e (n_max_p) - the maximal Landau level number, which is still useful for electron (positrons) distribution.
  !==================================================================================================================
  subroutine EE_pairs_mag_mu_fix(b,TkeV,mu,n_e30,n_p30,n_max,n_max_e,n_max_p)
  use EE_pairs_mag_help
  implicit none
  real*8,intent(in)::b,mu,TkeV
  real*8::n_e30,n_p30,eps,g_n,sum_e,sum_p,n_e30_add,n_p30_add,param,b_final
  integer::i,n_max,det
  integer::n_max_e,n_max_p         !==the maximal Landau level numbers for electrons and positrons==!

    eps=2.d-2
    mu_mod=mu
    TkeV_mod=TkeV
    b_mod=b
    !==rough estimation for the maximal possible Landau level==!
    if(TkeV/(b*511.d0).le.0.1d0)then
      n_max=1
    else
      n_max=7*TkeV/(b*511.d0)+15
    end if

    i=0; n_max_e=0; n_max_p=0
    param=1.d-7
    n_e30=0.d0
    n_p30=0.d0
    sum_e=n_e30
    sum_p=n_p30

    det=11
    do while(det.eq.11)
      do while(i.le.n_max)
        n_mod=i
        !==dealing with electrons==!
        if(i.eq.0)then; g_n=1.d0; else; g_n=2.d0; end if
        call int_ImproperInt_simpson(Ie_mag ,0.d0,1.d-3,eps,n_e30_add,b_final)
        n_e30_add=2*g_n*n_e30_add
        if(n_e30_add.lt.(param*n_e30))then
          n_max_e=min(n_max_e,i)
        else
          n_max_e=i
        end if
        n_e30     = n_e30+n_e30_add

        !==deaking with positrons==!
        if(mu_mod.le.10.d0*TkeV_mod)then
          call int_ImproperInt_simpson(Ip_mag ,0.d0,1.d-3,eps,n_p30_add,b_final)
          n_p30_add = 2*g_n*n_p30_add
          if(n_p30_add.lt.(param*n_p30))then
            n_max_p=min(n_max_p,i)
          else
            n_max_p=i
          end if
          n_p30     = n_p30+n_p30_add
        else
          call int_ImproperInt_simpson(Ip_mag2,0.d0,1.d-3,eps,n_p30_add,b_final)
          n_p30_add = 2*g_n*n_p30_add/exp(mu_mod/TkeV_mod)
          if(n_p30_add.lt.(param*n_p30))then
            n_max_p=min(n_max_p,i)
          else
            n_max_p=i
          end if
          n_p30     = n_p30+n_p30_add
        end if
        !write(*,*)"#a3",n_max,n_max_e,n_max_p
        i=i+1
      end do

      !write(*,*)"#a2 ",n_max,n_max_e,n_max_p

      !if( (((n_e30-sum_e)/n_e30).lt.eps) .and. (((n_p30-sum_p)/n_p30).lt.eps) )then
      if( (((n_e30-sum_e)/n_e30).lt.eps) )then
        !==we have obtained the necessary accuracy==!
!write(*,*)"!!!",mu,n_e30,sum_e,((n_e30-sum_e)/n_e30)
        det=1
      else
        !==we have to make the next step==!
        sum_e=n_e30
        sum_p=n_p30
        n_max=n_max*2
      end if
    end do
    n_e30=0.4415d0*n_e30*b
    n_p30=0.4415d0*n_p30*b
  122 return
  end subroutine EE_pairs_mag_mu_fix


  real*8 function Ie(p)
  use EE_pairs_help
  implicit none
  real*8,intent(in)::p
  real*8::E
    E=sqrt(511.d0**2+p**2*9.d20)
    Ie=p**2/(exp((E-mu_mod)/TkeV_mod)+1.d0)
  return
  end function Ie


  real*8 function Ie_mag(p_z)
  use EE_pairs_mag_help
  implicit none
  real*8,intent(in)::p_z
  real*8::E,t
    E=sqrt(1.d0+p_z**2+2*b_mod*n_mod)
    t=TkeV_mod*1.d3*11600.d0/5.93d9
    Ie_mag=1.d0/(exp((E-mu_mod/511)/t)+1.d0)
    !write(*,*)Ie_mag    
  return
  end function Ie_mag


  !==the accurate distribution function of pozitrons for a given chemical potential==!
  real*8 function Ip_mag(p_z)
  use EE_pairs_mag_help
  implicit none
  real*8,intent(in)::p_z
  real*8::E,t
    E=sqrt(1.d0+p_z**2+2*b_mod*n_mod)
    t=TkeV_mod*1.d3*11600.d0/5.93d9
    Ip_mag=1.d0/(exp((E+mu_mod/511)/t)+1.d0)   !!!!!!+1.d0 in Kaminker!!!
  return
  end function Ip_mag

  !==the approximate distribution function of pozitrons for a given chemical potential==!
  real*8 function Ip_mag2(p_z)
  use EE_pairs_mag_help
  implicit none
  real*8,intent(in)::p_z
  real*8::E,t
    E=sqrt(1.d0+p_z**2+2*b_mod*n_mod)
    t=TkeV_mod*1.d3*11600.d0/5.93d9
    Ip_mag2=1.d0/(exp((E)/t))   !!!!!!+1.d0 in Kaminker!!!
  return
  end function Ip_mag2


  real*8 function Ip(p)
  use EE_pairs_help
  implicit none
  real*8,intent(in)::p
  real*8::E
    E=sqrt(511.d0**2+p**2*9.d20)
    Ip=p**2/(exp((E+mu_mod)/TkeV_mod)+1.d0)
  return
  end function Ip

  real*8 function Ip2(p)
  use EE_pairs_help
  implicit none
  real*8,intent(in)::p
  real*8::E
    E=sqrt(511.d0**2+p**2*9.d20)
    Ip2=p**2/exp((E)/TkeV_mod)
  return
  end function Ip2

  !=======================================================================================
  ! The function returns Z_max, which covers (1-eps) particles at given Landau level.
  ! det_p=1 - electron
  ! det_p=2 - positron
  !=======================================================================================
  real*8 function get_Z_max_mag(b,T_keV,mu_keV,n,det_p,eps)
  use EE_pairs_mag_help
  implicit none
  real*8,intent(in)::b,T_keV,mu_keV,eps
  integer,intent(in)::n,det_p
  real*8::help,res
    TkeV_mod=T_keV
    b_mod=b
    n_mod=n
    mu_mod=mu_keV
    if(det_p.eq.1)then
      call int_ImproperInt_simpson(Ie_mag,0.d0,1.d-3,eps,help,res)
    else
      if(mu_mod.le.10.d0*TkeV_mod)then
        call int_ImproperInt_simpson(Ip_mag ,0.d0,1.d-3,eps,help,res)
      else
        call int_ImproperInt_simpson(Ip_mag2,0.d0,1.d-3,eps,help,res)
      end if
    end if
    get_Z_max_mag=res

  return
  end function get_Z_max_mag

end subroutine EE_pairs





!============================================================================================
! The function returnt the distribution function of electrons/positrons.
! p_z - particle momentum along B-field in units of [mc]
! n - Landau level number
! b - dimensionless magnetic field strength
! TkeV - temperature in keV
! mu - chemical potential in keV
! det_par=1 - electrons, det_par=2 - positrons
! The distribution normalised on the total number of paricles in units of 1.e30 cm^{-3}.
!============================================================================================
real*8 function EE_distribution(p_z,n,b,TkeV,mu,det_par)
implicit none
real*8,intent(in)::p_z,b,TkeV,mu
integer,intent(in)::n,det_par
integer::g_n
  if(n.eq.0)then; g_n=1; else; g_n=2; end if
  if(det_par.eq.1)then
    EE_distribution=0.4415d0*b*g_n*Ie_mag(p_z,b,n,mu,TkeV)
  else
    EE_distribution=0.4415d0*b*g_n*Ip_mag(p_z,b,n,mu,TkeV)
  end if
return
contains
  real*8 function Ie_mag(p_z,b,n,mu,TkeV)
  implicit none
  real*8,intent(in)::p_z,b,mu,TkeV
  integer,intent(in)::n
  real*8::E,t
    E=sqrt(1.d0+p_z**2+2*b*n)
    t=TkeV*1.d3*11600.d0/5.93d9
    Ie_mag=1.d0/(exp((E-mu/511.d0)/t)+1.d0)
  return
  end function Ie_mag

  !==the accurate distribution function of pozitrons for a given chemical potential==!
  real*8 function Ip_mag(p_z,b,n,mu,TkeV)
  implicit none
  real*8,intent(in)::p_z,b,mu,TkeV
  integer,intent(in)::n
  real*8::E,t
    E=sqrt(1.d0+p_z**2+2*b*n)
    t=TkeV*1.d3*11600.d0/5.93d9
    Ip_mag=1.d0/(exp((E+mu/511)/t)+1.d0)   !!!!!!+1.d0 in Kaminker!!!
  return
  end function Ip_mag

end function EE_distribution



subroutine Test_EE_distribution(n_e_ini30,b,TkeV)
implicit none
real*8,intent(in)::n_e_ini30,b,TkeV
real*8::EE_distribution,MaxwellRel1d_gen  !==functions==!
real*8::p_min,p_max,dp,p,res0,res1,res2,res3,n_e30,n_p30,mu,res_M
integer::nn,i,n_max,n_max_e,n_max_p
real*8::Z_e_max(5000),Z_p_max(5000)
  call EE_pairs(b,n_e_ini30,TkeV,n_e30,n_p30,mu,n_max,n_max_e,n_max_p,Z_e_max,Z_p_max)
  write(*,*)"#",mu
  !read(*,*)
  p_min=-6.d0
  p_max=+6.d0
  nn=1000
  dp=(p_max-p_min)/(nn-1)
  i=1
  do while(i.le.nn)
    p=p_min+dp*(i-1)
    res_M=MaxwellRel1d_gen(p,TkeV/511,0.d0)
    res0=EE_distribution(p,0,b,TkeV,mu,1)
    res1=EE_distribution(p,1,b,TkeV,mu,1)
    res2=EE_distribution(p,2,b,TkeV,mu,1)
    write(*,*)i,p,res0,res1,res2,res_M
    i=i+1
  end do
return
end subroutine Test_EE_distribution
!============================================================================================



!============================================================================================
! The ubroutine calculates the gas pressure along z-direction due to electron-positron pairs.
! b - dimensionless magnetic field
! n_e_ini30 - initial consentration of electrons
! TkeV - temperature
! mu - chemical potential in keV
! n_max - the maximal Landau level number, which is taken into account
!============================================================================================
subroutine EE_pairs_PressureZ(b,n_e_ini30,TkeV,mu,n_max,P23)
implicit none
real*8,intent(in)::b,n_e_ini30,TkeV,mu
integer,intent(in)::n_max
real*8::p_z_max,P23,res
integer::i,g_n
real*8::int_simpson_recursive_mas !==function==!
real*8::mas
dimension mas(4)
  p_z_max=sqrt((TkeV/511)**2+2*(TkeV/511))*5
  res=0.d0
  i=0
  do while(i.le.n_max)
    if(i.eq.0)then; g_n=1; else; g_n=2; end if
    mas(1)=b
    mas(2)=i*1.d0
    mas(3)=mu
    mas(4)=TkeV
    res=res+g_n*int_simpson_recursive_mas(fun,0.d0,p_z_max,mas,4,1.d-2)
    i=i+1
  end do
  P23=res*b*3.5d0
  !write(*,*)res,int_simpson_recursive_mas(fun,0.d0,p_z_max,mas,4,1.d-2)
return
contains
  !==the function which we will integrate==!
  real*8 function fun(p_z,mas,n_mas)
  implicit none
  real*8,intent(in)::p_z,mas
  integer,intent(in)::n_mas
  dimension mas(n_mas)
  real*8::b,n_real,mu,TkeV
  integer::n
    b   =mas(1)
    n_real=mas(2)
    mu  =mas(3)
    TkeV=mas(4)
    n=n_real
    fun=p_z**2*( Ie_mag(p_z,b,n,mu,TkeV) + Ip_mag(p_z,b,n,mu,TkeV) )/sqrt(1.d0+p_z**2)
    !write(*,*)fun,n,p_z**2,Ie_mag(p_z,b,n,mu,TkeV),sqrt(1.d0+p_z**2)
    !read(*,*)
  return
  end function fun

  real*8 function Ie_mag(p_z,b,n,mu,TkeV)
  implicit none
  real*8,intent(in)::p_z,b,mu,TkeV
  integer,intent(in)::n
  real*8::E,t
    E=sqrt(1.d0+p_z**2+2*b*n)
    t=TkeV*1.d3*11600.d0/5.93d9
    Ie_mag=1.d0/(exp((E-mu/511)/t)+1.d0)
!write(*,*)"aqd ",Ie_mag,p_z,b,n,mu,TkeV
!read(*,*)
  return
  end function Ie_mag


  !==the accurate distribution function of pozitrons for a given chemical potential==!
  real*8 function Ip_mag(p_z,b,n,mu,TkeV)
  implicit none
  real*8,intent(in)::p_z,b,mu,TkeV
  integer,intent(in)::n
  real*8::E,t
    E=sqrt(1.d0+p_z**2+2*b*n)
    t=TkeV*1.d3*11600.d0/5.93d9
    Ip_mag=1.d0/(exp((E+mu/511)/t)+1.d0)   !!!!!!+1.d0 in Kaminker!!!
  return
  end function Ip_mag

  !==the approximate distribution function of pozitrons for a given chemical potential==!
  real*8 function Ip_mag2(p_z,b,n,mu,TkeV)
  implicit none
  real*8,intent(in)::p_z,b,mu,TkeV
  integer,intent(in)::n
  real*8::E,t
    E=sqrt(1.d0+p_z**2+2*b*n)
    t=TkeV*1.d3*11600.d0/5.93d9
    Ip_mag2=1.d0/(exp((E)/t))   !!!!!!+1.d0 in Kaminker!!!
  return
  end function Ip_mag2

end subroutine EE_pairs_PressureZ



!==================================================================================
! Approximation for x=(kT/mc^2)<<1 and non-degenerate gas. See Zeldovich (8.3.11)
!==================================================================================
subroutine EE_pairs_app1(x,n30)
implicit none
real*8,intent(in)::x
real*8::n30
  n30=1.766d0*exp(-1.d0/x)*x**1.5d0
return
end subroutine EE_pairs_app1



!==================================================================================
! Approximation for x=(kT/mc^2)>>1 and non-degenerate gas. See Zeldovich (8.3.12)
!==================================================================================
subroutine EE_pairs_app2(x,n30)
implicit none
real*8,intent(in)::x
real*8::n30
  n30=1.766d0*2.d0*x**3*(1.d0-(0.5d0)**3+(0.3333d0)**3-(0.25)**3+(0.2d0)**3)
return
end subroutine EE_pairs_app2

