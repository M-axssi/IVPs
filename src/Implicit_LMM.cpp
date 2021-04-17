#include "Implicit_LMM.h"

void Implicit_LMM::Nonlinear_Solve(std::vector<Real> &  initial,std::vector<Real> &C,const IVPs & F,Real K)
{
  int n=F.Get_n();
  Real * temp=Change(initial+C+K*F.Get_diff(initial));
  int count=0;
  Real error=cblas_dnrm2(n,temp,1);
  while (error>1e-14 && count<=10)
    {
      ++count;
      Real * mat=Change(F.Get_Jacobi(initial));
      for (int i=0;i<n;++i)
	for (int j=0;j<n;++j)
	  mat[i+j*n]=K*mat[i+j*n];
      for (int i=0;i<n;++i)
	mat[i+i*n]=1+mat[i+i*n];
      Real *Y=Change((-K*F.Get_diff(initial)-C-initial));
      int ipiv[n];
      LAPACKE_dgesv(LAPACK_COL_MAJOR, n, 1, mat, n, ipiv, Y, n);
      initial=initial+Change(Y,n);
      delete [] Y;
      delete [] mat;
      temp=Change(initial+C+K*F.Get_diff(initial));
      error=cblas_dnrm2(n,temp,1);
    }
  delete [] temp;
}

