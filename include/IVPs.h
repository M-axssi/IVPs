#ifndef IVPS_H
#define IVPS_H
#include "Utility.h"

class IVPs
{
 protected:
  int n;
  std::vector<Real> u0; 
 public:
 IVPs(int N,std::vector<Real> & u):n(N),u0(u){};
  int Get_n()const {return n;};
  std::vector<Real>  Get_u0() const {return u0;};
  virtual std::vector<Real>  Get_diff(const std::vector<Real> &) const =0;
  virtual std::vector<std::vector<Real>> Get_Jacobi(const std::vector<Real> &) const=0;
};


#else
//Do nothing!
#endif
