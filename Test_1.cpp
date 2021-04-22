#include <iostream>
#include "Three_Body_Sys.h"
#include "Adams_Bashforth.h"
#include "TimeIntegratorFactory.h"
#include <cmath>
#include <fstream>
#include <functional>
#include<ctime>

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
  std::string s;
  int p;
  int steps;
  int turns;
  int increment;

  std::ifstream fin("Test1_input.txt");
  fin>>s>>p>>steps>>turns>>increment;
  
  std::vector<Real> u0;
  u0={0.994,0,0,0,-2.0015851063790825224,0};
  Real T=17.06521656015796;
  Three_Body_Sys F(6,u0);
  TimeIntegrator * TI=TimeIntegratorFactory::instance().CreateTimeIntegrator(s);
  
  std::ofstream fout("Test1.out");

  clock_t startTime,endTime;
  fout<<"The solution errors,convergence rates,CPU time for "<<s<<" which of "<<p<<" order of accuracy:"<<std::endl;		 
  std::vector<std::vector<Real>> temp;
      
  startTime=clock();
  TI->Solve_IVPs(steps,T,F,temp,p);
  endTime=clock();
  std::vector<Real> S=temp.back();
  Real E1=Get_max_norm(S,u0);
  fout<<"Steps ="<<steps<<":   Solution errors = "<<E1<<"  ,  CPU time = "<<" "<<(double)(endTime - startTime) / CLOCKS_PER_SEC<<" s "<<std::endl;
      
  for (int j=1;j<=turns;++j)
    {
      steps=steps+increment;
      std::vector<std::vector<Real>> A;
	  
      startTime=clock();
      TI->Solve_IVPs(steps,T,F,A,p);
      endTime=clock();
	  
      std::vector<Real> S=A.back();
      Real E2=Get_max_norm(S,u0);
      fout<<"Steps ="<<steps<<":   Solution errors = "<<E2<<"  ,  CPU time = "<<" "<<(double)(endTime - startTime) / CLOCKS_PER_SEC<<" s "<<"  ,  Convergence rates = "<<log(E1/E2)/log(double(steps)/(steps-increment))<<std::endl;
      E1=E2;
    }
  
}
 
