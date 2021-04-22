#ifndef THREE_BODY_SYS_H
#define THREE_BODY_SYS_H
#include "IVPs.h"

class Three_Body_Sys:public IVPs
{
public:
 Three_Body_Sys():IVPs(){};
  Three_Body_Sys(int N,std::vector<Real> & u):IVPs(N,u){};
  virtual std::vector<Real> Get_diff(const std::vector<Real> & u) const;
  virtual Real * Get_Jacobi(const std::vector<Real> &u) const;
  ~Three_Body_Sys(){};
private:
  const Real miu=1/81.45;
};

#else
//Do nothing!
#endif
