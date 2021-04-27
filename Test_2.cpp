#include <iostream>
#include "Three_Body_Sys.h"
#include "Adams_Bashforth.h"
#include "TimeIntegratorFactory.h"
#include <cmath>
#include <fstream>
#include <functional>
#include<ctime>

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
  std::ofstream fout("Test2.out");
  std::ifstream fin("Test2_input.txt");

  std::string s;
  int p;
  long long  steps;
  int turns;

  std::vector<Real> u0;
  u0={0.87978,0,0,0,-0.3797,0};
  Real T=19.14045706162071;
  Three_Body_Sys F(6,u0);
  
  while (fin>>s>>p>>steps>>turns)
    {
      TimeIntegrator * TI=TimeIntegratorFactory::instance().CreateTimeIntegrator(s);

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
      fout<<std::endl;
    }
}
 
