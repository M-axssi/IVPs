#include "Bdfs.h"

const std::vector<std::vector<Real>> Bdfs::Alpha=
  {
    {-1.0},
    {1.0/3,-4.0/3},
    {-2.0/11,9.0/11,-18.0/11},
    {3.0/25,-16.0/25,36.0/25,-48.0/25}
  };

const std::vector<Real> Bdfs::Beta={1.0,2.0/3,6.0/11,12.0/25};

const int Bdfs::_s[4]={1,2,3,4};

void Bdfs::Solve_IVPs(const int steps,                                  
		      const Real T,const IVPs & F,std::vector<std::vector<Real>>& S, int p)
{
  int n=F.Get_n();
  S.push_back(F.Get_u0());
  int s=_s[p-1];
  if (steps==0) return;
  const Real k=T/steps;
  for (int i=1;i<=s-1;++i)
    {
      std::vector<Real> temp(n,0);
      std::vector<Real> initial(S[i-1]);
      //initial+=k*F.Get_diff(S[i-1]);
      for (int j=0;j<=i-1;++j)
        {
      	  temp+=Alpha[i-1][j]*S[j];
      	}
      Nonlinear_Solve(initial,temp,F,-k*Beta[i-1]);
      S.push_back(initial);
    }
   for (int i=s;i<=steps;++i)
    {
      std::vector<Real>temp(n,0);
      std::vector<Real> initial(S[i-1]);
      //initial+=k*F.Get_diff(S[i-1]);
      for (int j=i-s;j<=i-1;++j)
	{
	  temp+=Alpha[s-1][j-i+s]*S[j];
	}
      Nonlinear_Solve(initial,temp,F,-k*Beta[s-1]);
      S.push_back(initial);
    }
}
