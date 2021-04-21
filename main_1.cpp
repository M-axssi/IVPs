#include <iostream>
#include "IVPs.h"
#include "Adams_Bashforth.h"
#include "TimeIntegratorFactory.h"
#include <cmath>
#include <fstream>
#include <functional>

Real Pow(Real x,int y)
{
  Real s=1;
  for (int i=1;i<=y;++i)
    {
      s*=x;
    }
  return s;
}

class Three_Body_Sys:public IVPs
{
public:
  Three_Body_Sys(int N,std::vector<Real> & u):IVPs(N,u){};
  std::vector<Real> Get_diff(const std::vector<Real> & u) const;
  Real * Get_Jacobi(const std::vector<Real> &u) const;
private:
  const Real miu=1/81.45;
};



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

int main()
{
  std::ifstream fin("main1_input.txt");

  std::string s;
  int p;
  int steps;
  Real T;
  
  fin>>s>>p>>steps;
  
  TimeIntegrator * TI=TimeIntegratorFactory::instance().CreateTimeIntegrator(s);

  std::vector<Real> u0;
  u0={0.994,0,0,0,-2.0015851063790825224,0};
  T=17.06521656015796;
  //u0={0.87978,0,0,0,-0.3797,0};
  
  std::vector<std::vector<Real>> A;
  Three_Body_Sys F(6,u0);
  TI->Solve_IVPs(steps,T,F,A,p);

  // for (auto x:A.back())
  //   std::cout<<x<<"  ";
    
  std::ofstream fout("main1_output.m");
  fout<<"X=[";
  for (int i=0;i<=steps;++i)
    {
      fout<<A[i][0]<<" ";
    }
  fout<<"]"<<std::endl;
  fout<<"Y=[";
  for (int i=0;i<=steps;++i)
    {
      fout<<A[i][1]<<" ";
    }
  fout<<"]"<<std::endl;
  fout<<"plot(X,Y);";
  fout.close();
}
 
