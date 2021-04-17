#ifndef EXPLICIT_LMM_H
#define EXPLICIT_LMM_H
#include "TimeIntegrator.h"

class Explicit_LMM:public TimeIntegrator
{
 public:
 Explicit_LMM():TimeIntegrator(){};
 virtual void  Solve_IVPs(const int steps,
			  const Real T,const IVPs & F,std::vector<std::vector<Real>>& S,int p=1)=0;
 virtual ~Explicit_LMM()=default;
};

#else
//Do nothin!

#endif
