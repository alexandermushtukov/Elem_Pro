!====================================================================================================
! Subroutine constructs the Data Base of neutrino energy losses due to the annihilation of e-e-pairs.
!====================================================================================================
subroutine Test_NeutrinoEmAnnih()
implicit none
real*8::T_kev,b,T10,lg_TkeV,B12,Q23_app,n_e_ini30,lg_n_e_ini30,Q23,mu_keV
integer::n_sum_max,n_max_e,n_max_p,det_res
character(len=100):: file_name
real*8::Z_e_max(5000),Z_p_max(5000)
real*8::IntNeutrinoAnn,x,Fnn  !==function==!

  200 format (8(es11.4,"   "),1(I6,"   "))
  201 format (7(es11.4,"   "),4(I6,"   "))
  202 format (1(A21),(es8.2))
  b=0.9d0
  !write(file_name,202)"./res/res3_nu_MC2e6_b",b
  write(file_name,202)"./res/res4_nu_MCnew_d",b
  write(*,*)"# file_name=",file_name; write(*,*)
  open (unit = 20, file = file_name)
  write(20,*)"# b=",b
  write(20,*)"# Format: b,lg_TkeV,lg_n_e_ini30,mu_keV,Q23,Q23_app,res/Q23_app,n_sum_max,n_max_e,n_max_p,det_res"
  write(20,*)
  close(20)

  !lg_TkeV=2.d0
  lg_TkeV=2.6d0
  !do while(lg_TkeV.lt.2.5d0)
  do while(lg_TkeV.le.3.d0)
    T_keV=10.d0**lg_TkeV
    T10=T_keV*1.1602d-3 !*1.1602d-3
    B12=b*44.12d0

    !if(lg_TkeV.eq.2.58d0)then
    !  lg_n_e_ini30=1.8d0
    !else
      lg_n_e_ini30=-5.d0
    !end if

    do while(lg_n_e_ini30.le.2.d0)
      n_e_ini30=10.d0**lg_n_e_ini30
      !call EE_pairs(b,n_e_ini30,TkeV,n_e30,n_p30,mu_keV,n_max,Z_e_max,Z_p_max)
      call NeutrinoAnn_app(T10,B12,Q23_app)
      call NeutrinoAnn(b,T_keV,n_e_ini30,Q23,mu_keV,n_sum_max,n_max_e,n_max_p,det_res)
      open(unit = 20, file = file_name, status = 'old',form='formatted',position="append")
      !write(*,200)b,lg_TkeV,T_keV,lg_n_e_ini30,n_e_ini30,res,Q23_app,res/Q23_app,n_sum_max
      !write(*,201)b,T_keV,n_e_ini30,mu_keV,res,Q23_app,res/Q23_app,n_sum_max
      write(*,201)b,lg_TkeV,lg_n_e_ini30,mu_keV,Q23,Q23_app,Q23/Q23_app,n_sum_max,n_max_e,n_max_p,det_res
      write(20,201)b,lg_TkeV,lg_n_e_ini30,mu_keV,Q23,Q23_app,Q23/Q23_app,n_sum_max,n_max_e,n_max_p,det_res

      close(20)
      lg_n_e_ini30=lg_n_e_ini30+0.1d0
    end do
    lg_TkeV=lg_TkeV+2.d-2
    write(*,*)
  end do

return
end subroutine Test_NeutrinoEmAnnih
!==================================================================================================


!=======================================================================================
! The subroutine calculates approximate(!) neutrino emissivity due to pair annihilation.
! The calculations are based on the approximation given in Kaminker+ 1992,Ph.Rev.D,46,10
! See eq.(26) there. The temperature should be well above Fermi temperature.
!=======================================================================================
subroutine NeutrinoAnn_app(T10,B12,Q23)
implicit none
real*8,intent(in)::T10,B12
real*8::Q23,Qc23,C2V,C2A,t,b
real*8::C_Ve,C_Ae,C_Vo,C_Ao
real*8::pi=3.141592653589793d0

  t=T10/0.593d0  !==see Kaminker 1992==!
  b=B12/44.12d0
  Qc23=1.015d0
  C_Ve=2*0.223d0+0.5d0
  C_Ae=0.5d0
  C_Vo=2*0.223d0-0.5d0
  C_Ao=-0.5d0
  C2V= C_Ve**2 +2*C_Vo**2 !(2*0.223d0+0.5d0)**2 !0.93d0
  C2A= C_Ae**2 +2*C_Ao**2 !0.75d0

  Q23=Qc23/pi**4*( (t**3*C2V*(1.d0+3.75d0*t)+t**4*(C2V+C2A)*P(t))*F(t,b) &
      + t*b**2/6.d0/(1.d0+b)*(C2V+C2A*b/(1.d0+b))*S(t) )*exp(-2.d0/t)
return
contains
  real*8 function S(t)
  implicit none
  real*8,intent(in)::t
  real*8::b1,b2,b3
    b1=1.058d0; b2=0.6701d0; b3=0.9143d0
    S=1.d0+b1*t+b2*t**2+b3*t**3+0.472*t**4
  return
  end function S

  real*8 function P(t)
  implicit none
  real*8,intent(in)::t
  real*8::a1,a2,a3,a4
    a1=3.581d0; a2=39.64d0; a3=24.43d0; a4=36.49d0
    P=1.d0+a1*t+a2*t**2+a3*t**3+a4*t**4+18.75d0*t**5
  return
  end function P

  real*8 function F(t,b)
  implicit none
  real*8,intent(in)::t,b
  real*8::R,c
  dimension R(3),c(3)
    c(1)=3.106d-6
    c(2)=1.491d-3
    c(3)=4.839d-6
    R(1:3)=1.d0+c(1:3)*b/t**2*exp(sqrt(2*b)/3/t)
    F=1.d0/R(1)/R(2)/R(3)
  return
  end function F
end subroutine NeutrinoAnn_app



!====================================================================================================
! The subroutine calculates accurately Annihilation emission of neutrinos in a strong magnetic field.
! The calculations are based on the paper written by Kaminker+ 1992, PhRevD.
! b - dimensionless magnetic filed strength
! TkeV - temperature in [keV]
! res is represented in units of [1.e23 erg/s/cm^3]
! n_sum_max - the maximal n_e+n_p taken into accpount in calculations.
!====================================================================================================
subroutine NeutrinoAnn(b,TkeV,n_e_ini30,res,mu_keV,n_sum_max,n_max_e,n_max_p,det_n)
implicit none
real*8,intent(in)::b,TkeV,n_e_ini30
integer::n_e,n_p,n_max,n_max_e,n_max_p,det,n_sum,det_n,n_sum_max,i
real*8::res,res_add,Z1,Z2,T,mu,mu_keV,f1,f2,eps,f,n_e30,n_p30,mas(1000)
real*8::pi=3.141592653589793d0
integer::ifEven    !==function==!
real*8::sum_array  !==function==!
real*8::Z_e_max(5000),Z_p_max(5000),sum_array_5_last

  eps=2.d-2
  !==getting the actual chemical potential==!
  call EE_pairs(b,n_e_ini30,TkeV,n_e30,n_p30,mu_keV,n_max,n_max_e,n_max_p,Z_e_max,Z_p_max)
  T=TkeV/511
  mu=mu_keV/511
  !write(*,*)"# chemical potential: ",b,TkeV,mu,n_max,n_max_e,n_max_p
  !read(*,*)

  !n_max=1000   !==the actual n_max is found later authomatically==!

  res=0.d0

  !==limits of electron/pozitron momentum, where we integrate==!
  Z1= -sqrt(T**2+2*T)!*4 !-2.d0
  Z2=  sqrt(T**2+2*T)!*4 !+2.d0

  mas(1:1000)=0.d0; det_n=11

  n_sum=0
  !n_sum=77
  !write(*,*)"#",n_sum
  do while(det_n.eq.11)
    if(n_sum.le.8)then
      call sum_over_diagonal_par(n_sum,n_max_e,n_max_p,Z_e_max,Z_p_max,b,T,Z1,Z2,mu,res_add)
    else
      !call sum_over_diagonal_par(n_sum,n_max_e,n_max_p,Z_e_max,Z_p_max,b,T,Z1,Z2,mu,res_add)
      call sum_over_diagonal_approx_par(n_sum,n_max_e,n_max_p,Z_e_max,Z_p_max,b,T,Z1,Z2,mu,res_add)
    end if
    res=res+res_add
    mas(n_sum+1)=res_add

    if( ((n_sum.gt.2).and.(ifEven(n_sum+1).eq.0)) .or. (res_add.lt.(1.d-5*res)) )then
      if((sum_array(mas,1000,1,n_sum+1).eq.0.d0).or.&
         ((sum_array(mas,1000,(n_sum+1)/2+1,n_sum+1)/sum_array(mas,1000,1,n_sum+1)).lt.1.d-2))then
        det_n=13
        n_sum_max=n_sum
      end if
      !==additional criterio 1==!
      if(n_sum.gt.10)then
        i=0; sum_array_5_last=0.d0
        do while(i.le.4)
          sum_array_5_last=sum_array_5_last+mas(n_sum+1-i)
          i=i+1
        end do
        !write(*,*)"#criterio:",(sum_array_5_last/sum_array(mas,1000,1,n_sum+1)),(eps*5.d0/(n_sum+1))
        if( (sum_array_5_last/sum_array(mas,1000,1,n_sum+1)) .lt. (eps*5/(n_sum+1))  )then
          det_n=14
          n_sum_max=n_sum
        end if
        !==additional criterio 2==!
        if(n_sum.ge.n_max)then
          det_n=15
        end if
      end if
    end if
    !write(*,*)"###",n_sum,res,res_add
    n_sum=n_sum+1
  end do
  res=res*b/3/(2*pi)**5   !==10^(23)==!
return
contains

  subroutine sum_over_diagonal(n_sum,b,T,Z1,Z2,mu,res)
  implicit none
  integer,intent(in)::n_sum
  real*8,intent(in)::b,T,Z1,Z2,mu
  integer::n_e,n_p
  real*8::res
    res=0.d0
    n_e=0
    do while(n_e.le.n_sum)
      n_p=n_sum-n_e
      !==integration over electron momentum==!

      f=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu)
      res=res+f
      !write(*,*)" ",n_e,n_p,res,f,f/res
      n_e=n_e+1
    end do
  return
  end subroutine sum_over_diagonal


  !=========================================================================================================
  !=========================================================================================================
  subroutine sum_over_diagonal_par(n_sum,n_max_e,n_max_p,Z_e_max,Z_p_max,b,T,Z1,Z2,mu,res)
  use omp_lib 
  implicit none
  integer,intent(in)::n_sum,n_max_e,n_max_p
  real*8,intent(in)::b,T,Z1,Z2,mu,Z_e_max(5000),Z_p_max(5000)
  integer::n_e,n_p,n_e_task,n_p_task,n_MC
  real*8::res
  real*8::Ann_integration_MC    !==function==!

    !write(*,*)"##1"
    !n_e_task=38; n_p_task=39
    !n_e_task=50; n_p_task=50
    !f=Ann_integration_MC(b,T*511,mu*511,n_e_task,n_p_task,&
    !-Z_e_max(n_e_task+1),Z_e_max(n_e_task+1),-Z_p_max(n_p_task+1),Z_p_max(n_p_task+1),int(2.e5))
    !write(*,*)f
    !write(*,*)"##2"
    !n_e_task=38; n_p_task=41
    !n_e_task=50; n_p_task=50
    !f=Ann_integration_MC(b,T*511,mu*511,n_e_task,n_p_task,&
    !  -Z_e_max(n_e_task+1),Z_e_max(n_e_task+1),-Z_p_max(n_p_task+1),Z_p_max(n_p_task+1),int(2.e5))
    !write(*,*)f
    !read(*,*)

    n_MC=int(1.e4)
    res=0.d0
    n_e=0
    f=0.d0
    call omp_set_num_threads(n_sum+1)
    !$omp parallel private(n_e_task,n_p_task) firstprivate(n_sum) reduction(+: f)
       n_e_task=omp_get_thread_num()
       n_p_task=n_sum-n_e_task
       if((n_e_task.le.n_max_e).and.(n_p_task.le.n_max_p))then
         !f=NeutrinoSyn_intZe(b,T,n_e_task,n_p_task,Z1,Z2,mu)
         f=Ann_integration_MC(b,T*511,mu*511,n_e_task,n_p_task,&
           -Z_e_max(n_e_task+1),Z_e_max(n_e_task+1),-Z_p_max(n_p_task+1),Z_p_max(n_p_task+1),n_MC)
       else
         f=0.d0
       end if
       !write(*,*)n_e_task,n_p_task,f
    !$omp end parallel
    res=f
  return
  end subroutine sum_over_diagonal_par


  !=========================================================================================================
  !=========================================================================================================
  subroutine sum_over_diagonal_approx_par(n_sum,n_max_e,n_max_p,Z_e_max,Z_p_max,b,T,Z1,Z2,mu,res)
  use omp_lib
  implicit none
  integer,intent(in)::n_sum,n_max_e,n_max_p
  real*8,intent(in)::b,T,Z1,Z2,mu,Z_e_max(5000),Z_p_max(5000)
  integer::n_e,n_p,M1,n_e_task,n_p_task
  integer::num,num_new,add_num,added_num
  real*8::fun,fun_new
  dimension num(n_sum+1),fun(n_sum+1),num_new(n_sum+1),fun_new(n_sum+1),added_num(n_sum+1)
  integer::i,N1,step_num,n_MC
  real*8::res,sum_array,sum_array_new,eps,fun_task
  real*8::Ann_integration_MC    !==function==!
    n_MC=int(2.e5)
    eps=1.d-2

    !==get the first three points for the sum==!
    N1=3
    num(1)=1; n_e=num(1); n_p=n_sum-n_e 
    if((n_e.le.n_max_e).and.(n_p.le.n_max_p))then
      !fun(1)=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu)
      fun(1)=Ann_integration_MC(b,T*511,mu*511,n_e,n_p,&
            -Z_e_max(n_e+1),Z_e_max(n_e+1),-Z_p_max(n_p+1),Z_p_max(n_p+1),n_MC)
    else; fun(1)=0.d0; end if

    num(2)=(n_sum+1)/2; n_e=num(2); n_p=n_sum-n_e
    if((n_e.le.n_max_e).and.(n_p.le.n_max_p))then
      !fun(2)=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu)
      fun(2)=Ann_integration_MC(b,T*511,mu*511,n_e,n_p,&
            -Z_e_max(n_e+1),Z_e_max(n_e+1),-Z_p_max(n_p+1),Z_p_max(n_p+1),n_MC)
    else; fun(2)=0.d0; end if

    num(3)=n_sum-1; n_e=num(3); n_p=n_sum-n_e 
    if((n_e.le.n_max_e).and.(n_p.le.n_max_p))then
      !fun(3)=NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu)
      fun(3)=Ann_integration_MC(b,T*511,mu*511,n_e,n_p,&
            -Z_e_max(n_e+1),Z_e_max(n_e+1),-Z_p_max(n_p+1),Z_p_max(n_p+1),n_MC)
    else; fun(3)=0.d0; end if
    !=============================================!

    sum_array=13.d0
    sum_array_new=2*sum_array
    step_num=0
    M1=11
    do while((M1.gt.0).and.( (abs(sum_array_new-sum_array)/sum_array_new).gt.eps ) )
      !write(*,*)"#1 ",n_sum+1,N1,(abs(sum_array_new-sum_array)/sum_array_new)
      step_num=step_num+1
      sum_array=sum_array_new
      call get_new_points_for_sum(n_sum+1,N1,num,fun,num_new,fun_new,added_num,M1)
      !write(*,*)"#3 M1=",M1,N1

      if(M1.gt.0)then
        call omp_set_num_threads(M1)
        !$omp parallel private(i,n_e_task,n_p_task,fun_task) firstprivate(n_sum,b,T,Z1,Z2,mu)
        i=omp_get_thread_num()+1
        n_e_task=num_new(added_num(i)) 
        n_p_task=n_sum-n_e_task
        if((n_e_task.le.n_max_e).and.(n_p_task.le.n_max_p))then
          !fun_task=NeutrinoSyn_intZe(b,T,n_e_task,n_p_task,Z1,Z2,mu)
          fun_task=Ann_integration_MC(b,T*511,mu*511,n_e_task,n_p_task,&
                   -Z_e_max(n_e_task+1),Z_e_max(n_e_task+1),-Z_p_max(n_p_task+1),Z_p_max(n_p_task+1),n_MC)
        else 
          fun_task=0.d0
        end if
        fun_new(added_num(i))=fun_task
        !$omp end parallel
      end if
      num(1:n_sum+1)=num_new(1:n_sum+1)
      fun(1:n_sum+1)=fun_new(1:n_sum+1)
      N1=N1+M1
      call sum_over_incomplete_array2(n_sum+1,N1,num,fun,sum_array_new)
    end do
    res=sum_array_new

    !==add two "edge" points to the sum==!
    if(n_sum.le.n_max_e)then
      !res=res+NeutrinoSyn_intZe(b,T,n_sum,0,Z1,Z2,mu)
      res=res+Ann_integration_MC(b,T*511,mu*511,n_sum,0,&
             -Z_e_max(n_sum+1),Z_e_max(n_sum+1),-Z_p_max(1),Z_p_max(1),n_MC)
    end if
    if(n_sum.le.n_max_p)then
      !res=res+NeutrinoSyn_intZe(b,T,0,n_sum,Z1,Z2,mu)
      res=res+Ann_integration_MC(b,T*511,mu*511,0,n_sum,&
                  -Z_e_max(1),Z_e_max(1),-Z_p_max(n_sum+1),Z_p_max(n_sum+1),n_MC)
    end if
  return
  end subroutine sum_over_diagonal_approx_par


  !===================================================================================
  ! Updated adaptive integration over the electron momentum. The integration is improper
  ! and Z1,Z2 - the initial interval which we take to try the integration.
  !===================================================================================
  real*8 function NeutrinoSyn_intZe(b,T,n_e,n_p,Z1,Z2,mu)
  implicit none
  real*8,intent(in)::b,T,Z1,Z2,mu
  integer,intent(in)::n_e,n_p
  real*8::mas(7),res,eps,Zp1,Zp2
  integer::n_max_int
  real*8::Zp1_,Zp2_
    eps=1.d-2; n_max_int=64
    Zp1= -sqrt(T**2+2*T)
    Zp2=  sqrt(T**2+2*T)
    mas(1)=Zp1
    mas(2)=Zp2
    mas(3)=b
    mas(4)=n_e*1.d0
    mas(5)=n_p*1.d0
    mas(6)=mu
    mas(7)=T
    call mf_int_Simp_adap_improp_1(NeutrinoSyn_intZp,Z1,Z2,eps,mas,7,n_max_int,res,Zp1_,Zp2_)
    NeutrinoSyn_intZe=res
  return
  end function NeutrinoSyn_intZe


  real*8 function NeutrinoSyn_intZp(Z,mas,n)
  implicit none
  real*8,intent(in)::Z,mas
  integer,intent(in)::n
  dimension mas(n),mas_(6)
  real*8::Zp1,Zp2,b,mu,T,res,eps,mas_
  integer::n_e,n_p,n_max_int
    eps=1.d-2; n_max_int=32
    Zp1=mas(1)
    Zp2=mas(2)
    b=mas(3)
    n_e=int(mas(4))
    n_p=int(mas(5))
    mu=mas(6)
    T=mas(7)
    mas_(1)=Z; mas_(2)=b; mas_(3)=n_e*1.d0; mas_(4)=n_p*1.d0; mas_(5)=mu; mas_(6)=T
    call mf_int_Simpson_adaptive(fun_int_Zp,Zp1,Zp2,eps,mas_,6,n_max_int,res)
    NeutrinoSyn_intZp=res
  return
  end function NeutrinoSyn_intZp
  !====================================================================================



  !=================================================================================================
  ! This function is integrated over positron momentum.
  !=================================================================================================
  real*8 function fun_int_Zp(Zp,mas,n)
  implicit none
  real*8,intent(in)::Zp,mas
  integer,intent(in)::n
  dimension mas(n),mas_(8)
  real*8::eps,mas_,Z_i,b,mu,T,q_z,omega,q_p_min,q_p_max,res
  integer::n_e,n_p
    eps=1.d-2
    !==reading input array==!
    Z_i=mas(1)
    b=mas(2)
    n_e=int(mas(3))
    n_p=int(mas(4))
    mu=mas(5)
    T=mas(6)
    !==end reading the array==!

    !==get limits of integration==!
    q_z=Z_i+Zp
    omega=sqrt(1.d0+2*n_e*b+Z_i**2)+sqrt(1.d0+2*b*n_p+Zp**2)
    q_p_min=0.d0
    q_p_max=sqrt(omega**2-q_z**2)
    !=============================!
    mas_(1)=n_e*1.d0; mas_(2)=n_p*1.d0; mas_(3)=Z_i; mas_(4)=Zp; mas_(5)=q_z; mas_(6)=b; mas_(7)=mu; mas_(8)=T

    call mf_int_Simpson_adaptive(fun_int_q_p,q_p_min,q_p_max,eps,mas_,8,111,res)
    fun_int_Zp=res
  return
  end function fun_int_Zp


  !=======================================================================================
  ! This function is adopted for integration over q_p.
  !=======================================================================================
  real*8 function fun_int_q_p(q_p,mas,n)
  implicit none
  real*8,intent(in)::q_p,mas
  integer,intent(in)::n
  dimension mas(n)
  real*8::IntNeutrinoAnn  !==function==!
    !mas(1)=n_e*1.d0; mas(2)=n_p*1.d0; mas(3)=Z_e; mas(4)=Z_p; mas(5)=q_z; mas(6)=b; mas(7)=mu; mas(8)=T
    fun_int_q_p=IntNeutrinoAnn(int(mas(1)),int(mas(2)),mas(3),mas(4),mas(5),q_p,mas(6),mas(7),mas(8))
  return
  end function fun_int_q_p


end subroutine NeutrinoAnn


!=============================================================================================
! Monte Carlo 3d integration.
! n_max - number of runs in MC simulation.
!=============================================================================================
real*8 function Ann_integration_MC(b,T_keV,mu_keV,n_e,n_p,Ze_min,Ze_max,Zp_min,Zp_max,n_m)
implicit none
real*8,intent(in)::b,T_keV,mu_keV,Ze_min,Ze_max,Zp_min,Zp_max
integer,intent(in)::n_e,n_p,n_m
real*8::random,Ze,Zp,q,q_z,omega,q_p_max,res,res_,eps
integer::i,n_max,det,n_m_max
real*8::IntNeutrinoAnn  !==function==!
  eps=2.d-2
  n_m_max=int(1.e7)        !==the maximal possible number of point in MC integration==!
  call init_random_seed()

  n_max=n_m
  res=0.d0
  i=1
  do while(i.le.n_max)
    call RANDOM_NUMBER(random)
    Ze=Ze_min+random*(Ze_max-Ze_min)
    call RANDOM_NUMBER(random)
    Zp=Zp_min+random*(Zp_max-Zp_min)
    q_z=Ze+Zp
    omega=sqrt(1.d0+2*n_e*b+Ze**2)+sqrt(1.d0+2*b*n_p+Zp**2)
    q_p_max=sqrt(omega**2-q_z**2)
    call RANDOM_NUMBER(random)
    q=random*q_p_max
    res=res+IntNeutrinoAnn(n_e,n_p,Ze,Zp,q_z,q,b,mu_keV/511,T_keV/511)*q_p_max
    i=i+1
  end do
  res=res/(i-1)
  res_=res

  !==try to get the necessary accuracy==!
  det=11
  do while(det.eq.11)
    n_max=n_max*2
    res=res*(i-1)
    !i=1   !==we start from the last i==!
    do while(i.le.n_max)
      call RANDOM_NUMBER(random)
      Ze=Ze_min+random*(Ze_max-Ze_min)
      call RANDOM_NUMBER(random)
      Zp=Zp_min+random*(Zp_max-Zp_min)
      q_z=Ze+Zp
      omega=sqrt(1.d0+2*n_e*b+Ze**2)+sqrt(1.d0+2*b*n_p+Zp**2)
      q_p_max=sqrt(omega**2-q_z**2)
      call RANDOM_NUMBER(random)
      q=random*q_p_max
      res=res+IntNeutrinoAnn(n_e,n_p,Ze,Zp,q_z,q,b,mu_keV/511,T_keV/511)*q_p_max
      i=i+1
    end do
    res=res/(i-1)
    if( ( abs((res-res_)/res).lt.eps ).or.( (res.eq.0.d0).and.(res_.eq.0.d0) ).or.(n_max.gt.n_m_max) )then
      det=12
    else
      res_=res
    end if
  end do

  Ann_integration_MC=res*(Ze_max-Ze_min)*(Zp_max-Zp_min)!*omega     !==fix it==! : is it a mistake?
return
end function Ann_integration_MC


Subroutine test_IntAnn()
implicit none
real*8::q_p1,q_p2,q_z,Ze,Zp,omega,b,mu,T,dq_p,res,q_p,coeff,f_i
integer::nn,n_e,n_p,i
real*8::SimpsonCoeff,IntNeutrinoAnn
  q_p1=3.d0 !0.d0
  q_p2=6.d0 !12.727695838029952
  q_z= -3.4142980019670275d0
  Ze=-1.7071490009835137d0
  Zp=-1.7071490009835137d0
  omega=13.177696012264093d0
  b=0.5d0
  n_e=40
  n_p=39
  mu=2.0962818003913899d-4
  T=0.97847358121330719d0
  nn=2000
  dq_p=(q_p2-q_p1)/(nn-1)
  res=0.d0
  i=1
  do while(i.le.nn)
    q_p=q_p1+(i-1)*dq_p
    coeff=SimpsonCoeff(i,nn)
    f_i=IntNeutrinoAnn(n_e,n_p,Ze,Zp,q_z,q_p,b,mu,T)
    !if(nn.gt.4.e2)then
    !  write(*,*)i,q_p,f_i
    !end if
    res=res+coeff*f_i*dq_p/3
    i=i+1
  end do

return
end subroutine test_IntAnn


!=============================================================================================================
! The function which is used in integration to get the intencity if neutrino emission due to the annihilation.
! mu - chemical potential in units of [m_{e}c^2]
! T  - temperature in units of [m_{e}c^2]
!=============================================================================================================
real*8 function IntNeutrinoAnn(n_e,n_p,Z_e,Z_p,q_z,q_p,b,mu,T)
implicit none
integer,intent(in)::n_e,n_p
real*8,intent(in)::Z_e,Z_p,q_z,q_p,b,mu,T
!integer::n_e,n_p
!real*8::Z_e,Z_p,q_z,q_p,b,mu,T
real*8::E_e,E_p,omega
real*8::res

!n_e=12
!n_p=2
!Z_e=3.5456465591638029d0
!Z_p=-1.2384556508395681d0
!q_z=2.3071909083242348d0
!q_p=0.39062599493478073d0
!b=0.1d0
!mu=6.0023483365949132d-4
!T=0.19569471624266144d0
  if(q_p.ne.0.d0)then
    E_e=sqrt(1.d0+2*b*n_e+Z_e**2)
    E_p=sqrt(1.d0+2*b*n_p+Z_p**2)
    omega=E_e+E_p
    res=q_p*A(b,n_e,n_p,Z_e,Z_p,E_e,E_p,q_z,q_p,omega)*omega*fe(E_e,mu,T)*fp(E_p,mu,T)   !==the function under the integral==!
    !write(*,*)"res: ",res,A(b,n_e,n_p,Z_e,Z_p,E_e,E_p,q_z,q_p,omega)
  else
    res=0.d0
  end if

  if(isnan(res))then
    open (unit = 20, file = "./res/resNaN", status = 'old',form='formatted',position="append")
    write(20,*)"#IntNeutrinoAnn - NaN:",res,n_e,n_p,Z_e,Z_p,q_z,q_p,b,mu,T
    !write(*,*)"#IntNeutrinoAnn - NaN:",res,n_e,n_p,Z_e,Z_p,q_z,q_p,b,mu,T
    !write(*,*)"#### ",q_p,A(b,n_e,n_p,Z_e,Z_p,E_e,E_p,q_z,q_p,omega),omega,fe(E_e,mu,T),fp(E_p,mu,T)
    close(20)
    res=0.d0
    !read(*,*)
  !else
  !  write(*,*)"#IntNeutrinoAnn - nonNaN:",res,n_e,n_p,Z_e,Z_p,q_z,q_p,b,mu,T
  !  read(*,*)
  end if

  IntNeutrinoAnn=res
return
contains
  !==Electron distribution function for a given chemical potential mu, see eq (7) in Kaminker 1992==!
  real*8 function fe(E,mu,T)
  implicit none
  real*8,intent(in)::E,mu,T
    fe=1.d0/(exp((E-mu)/T)+1.d0)
  return
  end function fe

  !==Positron distribution function for a given chemical potential mu, see eq (7) in Kaminker 1992==!
  real*8 function fp(E,mu,T)
  implicit none
  real*8,intent(in)::E,mu,T
    fp=1.d0/(exp((E+mu)/T)+1.d0)
  return
  end function fp

  !================================================================================================
  ! Function A given by (24) in Kaminker 1992 ==!
  !================================================================================================
  real*8 function A(b,n_i,n_f,Z_i,Z_f,E_i,E_f,q_z,q_p,omega)
  implicit none
  real*8,intent(in)::b,Z_i,Z_f,E_i,E_f,q_z,q_p,omega
  integer,intent(in)::n_i,n_f
  real*8::funPhi,funPsi   !==functions==!
  real*8::C_V_e,C_A_e,C_V_o,C_A_o,u,Psi_min,Psi_plus,Phi_min,Phi_plus,p_p_i2,p_p_f2
  real*8::A1,A2,A3
    C_V_e=2*0.223d0+0.5d0
    C_A_e= 0.5d0
    C_V_o=2*0.223d0-0.5d0
    C_A_o=-0.5d0

    u=q_p**2/2/b
    p_p_i2=2*b*n_i
    p_p_f2=2*b*n_f
    Phi_plus=funPhi(n_f,n_i,u,1)
    Phi_min =funPhi(n_f,n_i,u,-1)
    Psi_plus=funPsi(n_f,n_i,u,1)
    Psi_min =funPsi(n_f,n_i,u,-1)

    A1=Psi_plus*( 2*(E_i*E_f-Z_i*Z_f)**2+(E_i*E_f-Z_i*Z_f)*(2.d0+p_p_i2+p_p_f2-2*q_p**2)+q_p**2/2*(q_p**2-p_p_i2-p_p_f2-3.d0) )
    A1=A1+Phi_plus*( (E_i*E_f-Z_i*Z_f)*(1.d0+p_p_i2+p_p_f2) +1.d0+ 0.5d0*(p_p_i2+p_p_f2)*(3.d0-q_p**2+p_p_i2+p_p_f2) )
    A1=A1*(C_V_e**2+C_A_e**2)/E_i/E_f +2*A1*(C_V_o**2+C_A_o**2)/E_i/E_f

    A2=( omega**2-q_z**2-q_p**2/2 )*Psi_plus + ( omega**2/2-q_z**2/2-q_p**2 )*Phi_plus
    A2=A2*(C_V_e**2-C_A_e**2)/E_i/E_f +2*A2*(C_V_o**2-C_A_o**2)/E_i/E_f

    A3=( 2*omega**2-2*q_z**2-3*q_p**2 )*Psi_min + (p_p_i2-p_p_f2)*Phi_min
    A3=A3*(E_i*Z_f-E_f*Z_i)*C_V_e*C_A_e/E_i/E_f +2*A3*(E_i*Z_f-E_f*Z_i)*C_V_o*C_A_o/E_i/E_f
    A=A1+A2-A3
  return
  end function A

end function IntNeutrinoAnn


!================================================================================================================
! N  - total number of points for summation
! N1 - number of points which is taken to calculate approximate summ
! num(N) - array containing numbers of taken pointa
! fun(N) - array containing function in a given points
! M1 - a new number of points to be considered to calculate approximate sum (output)
! num_new(N) - array containing numbers of new points, actually we use only first (M1-N1) positions of this array
!================================================================================================================
subroutine get_new_points_for_sum(N,N1,num,fun,num_new,fun_new,added_num,M1)
implicit none
integer,intent(in)::N,N1,num
real*8,intent(in)::fun
dimension num(N),fun(N),num_new(N),fun_new(N),added_num(N)
integer::M1,num_new,added_num,i,j,k,new_1,new_2,add_1,add_2,help1,help2
real*8::eps,sum_1,sum_2,sum_tot,fun_new
  
  eps=1.d-2
  M1=0
  i=1; j=0; k=1
  do while(i.lt.(N1-1))
    !write(*,*)"a1"
    call mf_sum_fn_approx_(num(i+1)-num(i)+1  ,fun(i)  ,fun(i+1),sum_1,new_1,add_1)
    new_1=(new_1-1)+num(i)
    call mf_sum_fn_approx_(num(i+2)-num(i+1)+1,fun(i+1),fun(i+2),sum_2,new_2,add_2)
    new_2=(new_2-1)+num(i+1)   
    call mf_sum_fn_approx_(num(i+2)-num(i)+1  ,fun(i)  ,fun(i+2),sum_tot,help1,help2) 
    if(( abs((sum_1+sum_2-sum_tot)/(sum_1+sum_2)) ).lt.eps)then
      !==we will not add points==!
      !write(*,*)"a1.1",abs((sum_1+sum_2-sum_tot)/(sum_1+sum_2))
      M1=M1+0
      num_new(k)  =num(i)  ; fun_new(k)  =fun(i)
      num_new(k+1)=num(i+1); fun_new(k+1)=fun(i+1)
      num_new(k+2)=num(i+2); fun_new(k+2)=fun(i+2)
      k=k+2
    else
      !==we have to add points==!
      M1=M1+add_1+add_2
      num_new(k)=num(i)  ; fun_new(k)=fun(i)
      if(add_1.eq.1)then; 
        j=j+1; added_num(j)=k+add_1
        num_new(k+add_1)=new_1
      end if
      num_new(k+add_1+1)=num(i+1)  ; fun_new(k+add_1+1)=fun(i+1)
      if(add_2.eq.1)then; 
        j=j+1; added_num(j)=k+add_1+1+add_2
        num_new(k+add_1+1+add_2)=new_2
      end if  
      num_new(k+add_1+add_2+2)=num(i+2); fun_new(k+add_1+add_2+2)=fun(i+2)
      k=k+2+add_1+add_2
      !==we have added a new point to the massinve==!
    end if  
    i=i+2
  end do
  !==if there is only one interval left==!
  if(i.eq.(N1-1))then
    !write(*,*)"a2"
    call mf_sum_fn_approx_(num(i)-num(i-1)+1,fun(i-1),fun(i)  ,sum_1,new_1,add_1) 
    call mf_sum_fn_approx_(num(i+1)-num(i)+1,fun(i)  ,fun(i+1),sum_2,new_2,add_2)   
    call mf_sum_fn_approx_(num(i+1)-num(i-1)+1,fun(i-1),fun(i+1),sum_tot,help1,help2) 
    if(( (sum_1+sum_2-sum_tot)/(sum_1+sum_2) ).lt.eps)then
      !we will not add points
      M1=M1+0
      num_new(k+1)=num(i+1); fun_new(k+1)=fun(i+1)
    else
      ! we have to add points
      M1=M1+add_2
      if(add_2.eq.1)then
        j=j+1; added_num(j)=k+add_2
        num_new(k+add_2)=new_2
      end if  
      num_new(k+add_2+1)=num(i+1); fun_new(k+add_2+1)=fun(i+1)
      !==we have added a new point to the massinve==!
    end if  
  end if
  !==now we have all numbers to add==!
  !write(*,*)"a_res: ",M1
return
end subroutine get_new_points_for_sum



!================================================================================================
! N - a number of elements in a set
! f_1=f(1)
! f_N=f(N)
! number_add - number of points to add
!================================================================================================
subroutine mf_sum_fn_approx_(N,f_1,f_N,sum_app,new_N,number_add)
implicit none
integer,intent(in)::N
real*8,intent(in)::f_1,f_N
real*8::sum_app,df,f_i
integer::new_N,i,number_add
  if(N.gt.2)then
    df=(f_N-f_1)/(N-1)
    i=1; sum_app=0.d0
    do while(i.le.N) 
      f_i=f_1+(i-1)*df
      sum_app=sum_app+f_i
      i=i+1
    end do
    new_N=1+N/2
    number_add=1
  else
    if(N.eq.2)then
      sum_app=f_1+f_N
    end if
    if(N.eq.1)then
      sum_app=f_1
    end if
    new_N=1
    number_add=0
  end if
return
end subroutine mf_sum_fn_approx_



!===================================================================================================
! The subroutine calculates the summ over the elements of array fun(N), which total length is N, but where
! only N1<N elements are known. The numbers of these elements and fun(i) there are reprennted in first N1
! elements of arrays num() and fun().
!===================================================================================================
subroutine sum_over_incomplete_array(N,N1,num,fun,sum_array)
implicit none
integer,intent(in)::N,N1,num
real*8,intent(in)::fun
dimension num(N),fun(N)
integer::i
real*8::sum_array,sum_add
integer::new_N,number_add  !==for help==!
  sum_array=0.d0
  i=1
  do while(i.le.(N1-1))
    call mf_sum_fn_approx_(num(i+1)-num(i)+1,fun(i),fun(i+1),sum_add,new_N,number_add)
    sum_add=sum_add-fun(i+1)
    sum_array=sum_array+sum_add
    i=i+1
  end do
  sum_array=sum_array+fun(N1)
return
end subroutine sum_over_incomplete_array



!========================================================================================================
! The subroutine calculates the summ over the elements of array fun(N), which total length is N, but where
! only N1<N elements are known. The numbers of these elements and fun(i) there are reprennted in first N1
! elements of arrays num() and fun().
! Advansed version of subroutine sum_over_incomplete_array.
! There we use parabolic approximation for point without the data.
!========================================================================================================
subroutine sum_over_incomplete_array2(N,N1,num,fun,sum_array)
implicit none
integer,intent(in)::N,N1,num(N)
real*8,intent(in)::fun(N)
integer::i,j
real*8::sum_array,sum_add,a,b,c
  sum_array=0.d0
  i=1
  do while(i.le.(N1-1))
    if((num(i+1)-num(i)).le.0)then
      sum_add=0.d0
    else
      if((num(i+1)-num(i)).eq.1)then
        sum_add=fun(i)
      else
        !==there are some points b/w num(i) and num(i+1)==!
        sum_add=0.d0
        if(i.eq.1)then
          !==we are looking at the first point in a row==!
          if(i.le.(N1-2))then
            !==there are two points after i==!
            call get_parabola_coeff(1.d0*num(i),fun(i),1.d0*num(i+1),fun(i+1),1.d0*num(i+2),fun(i+2),a,b,c)
            j=num(i)
            do while(j.lt.num(i+1))
              sum_add=sum_add+( a*j**2 + b*j + c )
              j=j+1
            end do
          else
            !==there are NO two points after i==!
            j=num(i)
            do while(j.lt.num(i+1))
              sum_add=sum_add+( fun(i) +(fun(i+1)-fun(i))/(num(i+1)-num(i))*(j-num(i)) )
              j=j+1
            end do
          end if
        else
          !==the point is not the first in a row==!
          call get_parabola_coeff(1.d0*num(i-1),fun(i-1),1.d0*num(i),fun(i),1.d0*num(i+1),fun(i+1),a,b,c)
          j=num(i)
          do while(j.lt.num(i+1))
            sum_add=sum_add+( a*j**2 + b*j + c )
            j=j+1
          end do
        end if
      end if
    end if
    sum_array=sum_array+sum_add
    i=i+1
  end do
  sum_array=sum_array+fun(N1)
return
contains
  !===========================================================================================
  ! The subroutine calculates coeffisients of parabola using coordinates of 3 points. 
  !===========================================================================================
  subroutine get_parabola_coeff(x1,y1,x2,y2,x3,y3,a,b,c)
  implicit none
  real*8,intent(in)::x1,y1,x2,y2,x3,y3
  real*8::a,b,c
    a=( y3 - (x3*(y2-y1)+x2*y1-x1*y2)/(x2-x1) )/( x3*(x3-x1-x2) + x1*x2 )
    b=(y2-y1)/(x2-x1) - a*(x1+x2)
    c=(x2*y1-x1*y2)/(x2-x1) + a*x1*x2
  return
  end subroutine get_parabola_coeff
end subroutine sum_over_incomplete_array2
!=============================================================================================




