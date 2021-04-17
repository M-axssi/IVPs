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
  Real*  Get_Jacobi(const std::vector<Real> &u) const;
private:
  const Real miu=1/81.45;
};


Real*  Three_Body_Sys::Get_Jacobi(const std::vector<Real>  &u) const
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

Real Get_1_norm(std::vector<Real> x,std::vector<Real> y)
{
  int n=x.size();
  Real temp=0;
  for (int i=0;i<n;++i)
    {
      temp+=fabs((x[i]-y[i]));
    }
  return temp;
}

int main()
{
  std::string s;
  int p;
  long long  steps;
  int turns;

  std::ifstream fin("Test2_input.txt");
  fin>>s>>p>>steps>>turns;
  
  std::vector<Real> u0;
  u0={0.87978,0,0,0,-0.3797,0};
  Real T=19.14045706162071;
  Three_Body_Sys F(6,u0);
  TimeIntegrator * TI=TimeIntegratorFactory::instance().CreateTimeIntegrator(s);
  
  std::ofstream fout("Test2.out");

  clock_t startTime,endTime;
  fout<<"The solution errors,convergence rates,CPU time for "<<s<<" which of "<<p<<" order of accuracy:"<<std::endl;		 

  std::vector<std::vector<Real>>  temp;

  long long  st[turns+1];
  Real time[turns+1];
  Real error[turns+1];
  
  for (int j=0;j<=turns;++j)
    {
      st[j]=steps;
      startTime=clock();
      std::vector<std::vector<Real>> A;
      TI->Solve_IVPs(steps,T,F,A,p);
      endTime=clock();
      time[j]=(double)(endTime - startTime) / CLOCKS_PER_SEC;

      if (j!=0)
	{
	  Real h=2*T/steps;
	  for (int i=0;i<=steps;i+=2)
	    {
	      int t=i/2;
	      error[j-1]+=h*(Get_1_norm(temp[t],A[i]));
	    }
	}
      temp=A;
      steps=steps*2;
    }
  for (int j=0;j<=turns-1;++j)
    {
      fout<<"Steps = "<<st[j]<<": CPU time = "<<time[j];
      fout<<" , Solution errors = "<<error[j];
      if (j!=0)
	{
	  fout<<" , Convergence rates = "<<log(error[j-1]/error[j])/log(2.0);
	}
      fout<<std::endl;
    }
}
 
