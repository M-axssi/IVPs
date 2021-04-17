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
  std::vector<std::vector<Real>> Get_Jacobi(const std::vector<Real> &u) const;
private:
  const Real miu=1/81.45;
};



std::vector<std::vector<Real>> Three_Body_Sys::Get_Jacobi(const std::vector<Real>  &u) const
{
  Real C_1=u[1]*u[1]+u[2]*u[2]+(u[0]+miu-1)*(u[0]+miu-1);
  Real C_2=u[1]*u[1]+u[2]*u[2]+(u[0]+miu)*(u[0]+miu);
  std::vector<std::vector<Real>> m;
  m.push_back({0,0,0,1,0,0});
  m.push_back({0,0,0,0,1,0});
  m.push_back({0,0,0,0,0,1});
  auto f1=[&](Real x,Real y)->Real {return (miu*sqrt(Pow(C_1,3))-3*sqrt(C_1)*x*x*miu)/Pow(C_1,3)+
				    ((1-miu)*sqrt(Pow(C_2,3))-3*sqrt(C_2)*y*y*(1-miu))/Pow(C_2,3);};
  auto f2=[&](Real x,Real y,Real z)->Real{return (3*sqrt(C_1)*x*y*miu)/Pow(C_1,3)+
					  (3*sqrt(C_2)*x*z*(1-miu))/Pow(C_2,3);};
  m.push_back({1-f1(u[0]+miu-1,u[0]+miu),f2(u[1],u[0]+miu-1,u[0]+miu),f2(u[2],u[0]+miu-1,u[0]+miu),0,2,0});
  m.push_back({f2(u[1],u[0]+miu-1,u[0]+miu),1-f1(u[1],u[1]),f2(u[1],u[2],u[2]),-2,0,0});
  m.push_back({f2(u[2],u[0]+miu-1,u[0]+miu),f2(u[1],u[2],u[2]),-f1(u[2],u[2]),0,0,0});
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
  std::ifstream fin("input.txt");

  std::string s;
  int p;
  int question_num;
  Real T;
  int steps;

  fin>>s>>p>>question_num>>T>>steps;
  
  TimeIntegrator * ABF=TimeIntegratorFactory::instance().CreateTimeIntegrator(s);

  std::vector<Real> u0;
  if (question_num==1)
     u0={0.994,0,0,0,-2.0015851063790825224,0};
  else
     u0={0.87978,0,0,0,-0.3797,0};
  
  std::vector<std::vector<Real>> A;
  Three_Body_Sys F(6,u0);
  ABF->Solve_IVPs(steps,T,F,A,p);

  for (auto x:A.back())
    std::cout<<x<<"  ";
    
  std::ofstream fout("output.m");
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
 
