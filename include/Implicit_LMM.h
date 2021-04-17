#ifndef IMPLICIT_LMM_H
#define IMPLICIT_LMM_H
#include "TimeIntegrator.h"

class Implicit_LMM:public TimeIntegrator
{
 public:
 Implicit_LMM():TimeIntegrator(){};
 virtual void  Solve_IVPs(const int steps,
			  const Real T,const IVPs & F,std::vector<std::vector<Real>>& S,int p=1)=0;
 void Nonlinear_Solve(std::vector<Real> &  initial,std::vector<Real> &C,const IVPs & F,Real K);
 virtual ~Implicit_LMM()=default;
};

#else
//Do nothing!

#endif
