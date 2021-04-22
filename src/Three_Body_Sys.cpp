#include "Three_Body_Sys.h"

Real* Three_Body_Sys::Get_Jacobi(const std::vector<Real>  &u) const
{
  Real C_1=u[1]*u[1]+u[2]*u[2]+(u[0]+miu-1)*(u[0]+miu-1);
  Real C_2=u[1]*u[1]+u[2]*u[2]+(u[0]+miu)*(u[0]+miu);
  Real* m=new Real [36];
  for (int i=0;i<36;++i) m[i]=0;
  m[18]=m[25]=m[32]=1;
  auto f1=[&](Real x,Real y)->Real {return (miu*sqrt(Pow(C_1,3))-3*sqrt(C_1)*x*x*miu)/Pow(C_1,3)+
				    ((1-miu)*sqrt(Pow(C_2,3))-3*sqrt(C_2)*y*y*(1-miu))/Pow(C_2,3);};
  auto f2=[&](Real x,Real y,Real z)->Real{return (3*sqrt(C_1)*x*y*miu)/Pow(C_1,3)+
					  (3*sqrt(C_2)*x*z*(1-miu))/Pow(C_2,3);};
  m[3]=1-f1(u[0]+miu-1,u[0]+miu);
  m[9]=f2(u[1],u[0]+miu-1,u[0]+miu);
  m[15]=f2(u[2],u[0]+miu-1,u[0]+miu);
  m[27]=2;

  m[4]=f2(u[1],u[0]+miu-1,u[0]+miu);
  m[10]=1-f1(u[1],u[1]);
  m[16]=f2(u[1],u[2],u[2]);
  m[22]=-2;

  m[5]=f2(u[2],u[0]+miu-1,u[0]+miu);
  m[11]=f2(u[1],u[2],u[2]);
  m[17]=-f1(u[2],u[2]);
  return m;
}
  

std::vector<Real> Three_Body_Sys::Get_diff(const std::vector<Real> & u) const
{
  Real C_1=sqrt(Pow(u[1]*u[1]+u[2]*u[2]+(u[0]+miu-1)*(u[0]+miu-1),3));
  Real C_2=sqrt(Pow(u[1]*u[1]+u[2]*u[2]+(u[0]+miu)*(u[0]+miu),3));
  std::vector<Real> temp(n);
  temp[0]=u[3];
  temp[1]=u[4];
  temp[2]=u[5];
  temp[3]=2*u[4]+u[0]-miu*(u[0]+miu-1)/C_1-(1-miu)*(u[0]+miu)/C_2;
  temp[4]=-2*u[3]+u[1]-miu*u[1]/C_1-(1-miu)*u[1]/C_2;
  temp[5]=-miu*u[2]/C_1-(1-miu)*u[2]/C_2;
  return temp;
}
