#include <iostream>
#include "Three_Body_Sys.h"
#include "Adams_Bashforth.h"
#include "TimeIntegratorFactory.h"
#include <cmath>
#include <fstream>
#include <functional>

int main()
{
  std::ifstream fin("main2_input.txt");

  std::string s;
  int p;
  int steps;
  Real T;
  
  fin>>s>>p>>steps;
  
  TimeIntegrator * TI=TimeIntegratorFactory::instance().CreateTimeIntegrator(s);

  std::vector<Real> u0;
  T=19.14045706162071;
  u0={0.87978,0,0,0,-0.3797,0};
  
  std::vector<std::vector<Real>> A;
  Three_Body_Sys F(6,u0);
  TI->Solve_IVPs(steps,T,F,A,p);

  // for (auto x:A.back())
  //   std::cout<<x<<"  ";
    
  std::ofstream fout("main2_output.m");
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
 
