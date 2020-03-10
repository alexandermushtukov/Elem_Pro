!==================================================================================
! Some approximations are taken from Kaminker 1992 (11)
!==================================================================================
real*8 function Fnn(n1,n2,u)
implicit none
integer,intent(in)::n1,n2
real*8,intent(in)::u
real*8::laguer_general,fact,fact2fact,ermit,laguer_gen_row_new  !==function==!
real*8::l,big_num,f,t,eps,ksi
real*8::pi=3.141592653589793d0
integer::nu
  big_num=20.d0
  eps=2.d-2
  if((n1.ge.0).and.(n2.ge.0))then
    l=n1-n2
    !write(*,*)"# ",u,( 1.d-1*(1+abs(l))/min(n1,n2) )
    if(u.lt.( eps*(1+abs(l))/min(n1,n2) ))then
      !==approximation (23) from Kaminker&Yakovlev 1981 ==!
      if((abs(l)).le.big_num)then 
        !write(*,*)"e"      
        Fnn=(-1)**((l+int(abs(l)))*0.5d0) / fact(int(abs(l))) *u**(0.5d0*int(abs(l))) *exp(-u/2) &
              * sqrt(fact2fact(max(n1,n2),min(n1,n2)))
      else
        f=max(n1,n2)
        t=min(n1,n2)
        if(min(n1,n2).le.big_num)then
          !write(*,*)"a"
          Fnn=(-1)**((l+int(abs(l)))*0.5d0) /(2*pi)**0.25d0 /sqrt(gamma(t+1.d0))&
              *exp( 0.25d0*log(f)*(1.d0+2*f) -0.5d0*log(abs(l))*(1.d0+2*abs(l)) +0.5d0*abs(l)*(log(u)+2.d0) -0.5d0*(f+u) )  
        else
          !write(*,*)"b"
          !==both max(n1,n2) and min(n1,n2) are large numbers==!
          Fnn=(-1)**((l+int(abs(l)))*0.5d0) /(2*pi)**0.50d0 &
              *exp( 0.25d0* (log(f)*(1.d0+2*f)-log(t)*(1.d0+2*t)) &
                   -0.5d0*log(abs(l))*(1.d0+2*abs(l)) +0.5d0*abs(l)*(log(u)+1.d0) -0.5d0*u )  
        end if
      end if
      !Fnn=1.d0
      !==NOTE: can be simplified further for large l ==!
    else
      if(u.le.( eps*(1+abs(l)) ))then
        !==approximation (24a) from Kaminker&Yakovlev 1981 ==!
        !write(*,*)"c_"
        nu=n1+n2+1
        f=max(n1,n2)
        t=min(n1,n2)
        if(f.ge.big_num)then
          if(t.lt.big_num)then
            !==f - big, t - not big; approximation (24a)==!
            Fnn=(2*pi*f)**0.25d0/sqrt(gamma(t+1)) *exp( f/2*(log(2*f/(f+t+1.d0))-1.d0) -t/2*log(2/(f+t+1.d0)) )&
                 *bessel_jn(n1-n2, -sqrt(2*nu*u))
          else
            !==f - big, t - big; approximation (24a)==!
            Fnn=(f/t)**0.25d0 *exp( f/2*(log(2*f/(f+t+1.d0))-1.d0) -t/2*(log(2*t/(f+t+1.d0)) -1.d0) )&
                 *bessel_jn(n1-n2, -sqrt(2*nu*u))
          end if
        else
          Fnn=sqrt(gamma(f+1.d0)/gamma(t+1.d0)) *(2.d0/nu)**(abs(l)/2) *bessel_jn(n1-n2, -sqrt(2*nu*u))
        end if
      else
        if(n2.gt.(10*n1))then
          !==use approximation (32) from Kaminker&Yakovlev==! 
          Fnn=1.d0/(2*pi*n2)**0.25d0 /sqrt(2.d0**n1*fact(n1)) *exp(-ksi**2/4) *ermit(n1,ksi/sqrt(2.d0))
        else
          !==direct formula==!problem is in calculations of Lagguerre polinoms
          !write(*,*)"c" 
          if((n1.gt.35).and.(n2.gt.35).and.(u.gt.50.d0))then
          !if((n1.gt.35).and.(n2.gt.35).and.(abs(laguer_general(n1,n2-n1,u)).gt.2.d4))then
            Fnn=Fnn_app_q_classic(n1,n2,u)
          else
            Fnn=(fact2fact(n1,n2))**0.5d0 *u**((n2-n1)/2.d0) *exp(-u/2) *laguer_general(n1,n2-n1,u) 
          end if
          !Fnn=1.d0/(2*pi*n2)**0.25d0 /sqrt(2.d0**n1*fact(n1)) *exp(-ksi**2/4) *ermit(n1,ksi/sqrt(2.d0))
        end if
      end if
      !!!
    end if
  else
    Fnn=0.d0
  end if
return
contains 
  real*8 function Fnn_app_q_classic(n1,n2,u)
  implicit none
  real*8::pi=3.141592653589793d0
  integer,intent(in)::n1,n2
  real*8,intent(in)::u
  real*8::u1,u2,Delta,f_u,p,Phi,res,fi,eps,ksi1,ksi2
  integer::l,nu
  real*8::funAi  !==functions==!
    eps=2.d-2
    l=n1-n2
    nu=n1+n2+1
    u1=(sqrt(1.d0*n2)-sqrt(1.d0*n1))**2
    u2=(sqrt(1.d0*n2)+sqrt(1.d0*n1))**2
    Delta=u2-u1
    f_u=(u-u1)*(u-u2)/4/u**2
    !write(*,*)u1,u,u2
    !write(*,*)"## ",eps*Delta,abs(u-u1),abs(u-u2)
    if(((abs(u-u1)).lt.(eps*Delta)).or.((abs(u-u2)).lt.(eps*Delta)))then
      !==&&==!
      if((abs(u-u1)).gt.(eps*Delta))then
        !==(29a)==!
        ksi1=-(u-u1)*sqrt(Delta/4/u1**2)
        res=(-1)**((l+int(abs(l)))*0.5d0) *sqrt(2.d0/pi) /(2*u1*Delta)**(1.d0/6.d0) *funAi(ksi1)
      else
        !==(29b)==!
        ksi2= (u-u2)*sqrt(Delta/4/u2**2)  
        res=(-1)**n1 *sqrt(2.d0/pi) /(2*u2*Delta)**(1.d0/6.d0) *funAi(ksi2)
      end if
    else
      !==&&==!
      if((u.lt.u1).or.(u.gt.u2))then
        p=sqrt(-f_u)
        Phi=abs(1.d0*l)/2*log( (sqrt(u1*abs(u2-u)) +sqrt(u2*(u1-u)))**2/u/Delta ) &
           + 0.5d0*nu*log( (sqrt(abs(u2-u)) - sqrt(abs(u1-u)))**2/Delta )
        if(u.lt.u1)then
          !==use (27a) from Kaminker&Yakovlev (1981) ==!
          res=(-1)**((l+int(abs(l)))*0.5d0) /sqrt(4*pi*abs(p)*u) *exp(abs(p)*u-Phi)
        else
          !==use (27b) from Kaminker&Yakovlev (1981) ==!
          res=(-1)**n1 /sqrt(4*pi*abs(p)*u) *exp(-abs(p)*u-Phi)
        end if
      else
        !==use (30) from Kaminker&Yakovlev (1981) ==!
        p=sqrt(-f_u)
        fi=0.5d0*nu *asin((u1+u2-2*u)/Delta) -abs(0.5d0*l)*asin((2*u1*u2-u*u1-u*u2)/u/Delta) -pi/4*(nu+abs(l)+2*l-1)
        res=1.d0/sqrt(pi*p*u)*cos(p*u-fi)
      end if
    end if

    Fnn_app_q_classic=res
  return
  end function Fnn_app_q_classic
end function Fnn


subroutine test_Fnn()
implicit none
real*8::Fnn,Fnn_approx01  !==functions==!
integer::i,n1,n2
real*8::x,laguer_general,laguer_gen_asymp_2
  x=1.d0
  n1=4000
  n2=3
  write(*,*)Fnn(n1,n2,x),laguer_general(n1,1,x),laguer_gen_asymp_2(n1,1,x)

return
end subroutine test_Fnn
!==============================================================================================
