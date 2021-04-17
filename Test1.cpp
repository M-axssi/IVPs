#include <iostream>
#include "IVPs.h"
#include "Adams_Bashforth.h"
#include "TimeIntegratorFactory.h"
#include <cmath>
#include <fstream>
#include <functional>
#include<ctime>

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

Real Get_max_norm(std::vector<Real> x,std::vector<Real> y)
{
  int n=x.size();
  Real temp=0;
  for (int i=0;i<n;++i)
    {
      if (fabs(x[i]-y[i])>temp)
      temp=fabs((x[i]-y[i]));
    }
  return temp;
}

int main()
{
  std::vector<Real> u0;
  u0={0.994,0,0,0,-2.0015851063790825224,0};
  Real T=17.06521656015796;
  Three_Body_Sys F(6,u0);
  TimeIntegrator * ABF=TimeIntegratorFactory::instance().CreateTimeIntegrator("Bdfs");
  //TimeIntegrator * RK=TimeIntegratorFactory::instance().CreateTimeIntegrator("Runge_Kutta");
  
  std::ofstream fout("Test1_Bdfs.out");

  clock_t startTime,endTime;
  
  for (int i=1;i<=4;++i)
    {
      fout<<"The solution errors,convergence rates,CPU time for Bdfs of "
	<<i<<" order of accuracy:"<<std::endl;		 
      int steps=3000000;
      std::vector<std::vector<Real>> temp;
      
      startTime=clock();
      ABF->Solve_IVPs(steps,T,F,temp,i);
      endTime=clock();
      std::vector<Real> S=temp.back();
      Real E1=Get_max_norm(S,u0);
      fout<<"Steps ="<<steps<<":   solution errors = "<<E1<<"  ,  CPU time = "<<" "<<(double)(endTime - startTime) / CLOCKS_PER_SEC<<" s "<<std::endl;
      
      for (int j=1;j<=5;++j)
	{
	  steps=steps+1000000;
	  std::vector<std::vector<Real>> A;
	  
	  startTime=clock();
	  ABF->Solve_IVPs(steps,T,F,A,i);
	  endTime=clock();
	  
	  std::vector<Real> S=A.back();
	  Real E2=Get_max_norm(S,u0);
	  fout<<"Steps ="<<steps<<":   solution errors = "<<E2<<"  ,  CPU time = "<<" "<<(double)(endTime - startTime) / CLOCKS_PER_SEC<<" s "<<"  ,  convergence rates = "<<log(E1/E2)/log(double(steps)/(steps-1000000))<<std::endl;
	  E1=E2;
	}
        fout<<std::endl<<std::endl<<std::endl;
    }

  // for (int i=4;i<=4;++i)
  //   {
  //     fout<<"The solution errors,convergence rates,CPU time for Runge_Kutta of "
  // 	<<i<<" order of accuracy:"<<std::endl;		 
  //     int steps=25600;
  //     std::vector<std::vector<Real>> temp;
      
  //     startTime=clock();
  //     RK->Solve_IVPs(steps,T,F,temp,i);
  //     endTime=clock();
  //     std::vector<Real> S=temp.back();
  //     Real E1=Get_max_norm(S,u0);
  //     fout<<"Steps ="<<steps<<":   solution errors = "<<E1<<"  ,  CPU time = "<<" "<<(double)(endTime - startTime) / CLOCKS_PER_SEC<<" s "<<std::endl;
      
  //     for (int j=1;j<=6;++j)
  // 	{
  // 	  steps=steps+10000;
  // 	  std::vector<std::vector<Real>> A;
	  
  // 	  startTime=clock();
  // 	  RK->Solve_IVPs(steps,T,F,A,i);
  // 	  endTime=clock();
	  
  // 	  std::vector<Real> S=A.back();
  // 	  Real E2=Get_max_norm(S,u0);
  // 	  fout<<"Steps ="<<steps<<":   solution errors = "<<E2<<"  ,  CPU time = "<<" "<<(double)(endTime - startTime) / CLOCKS_PER_SEC<<" s "<<"  ,  convergence rates = "<<log(E1/E2)/log(double(steps)/(steps-10000))<<std::endl;
  // 	  E1=E2;
  // 	}
  //       fout<<std::endl<<std::endl<<std::endl;
  //   }
  
}
 
