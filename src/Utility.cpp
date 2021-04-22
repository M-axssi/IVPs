#include "Utility.h"

void operator+=(std::vector<Real> &a,const std::vector<Real> & b)
{
  int size=a.size();
  for (int i=0;i<size;++i) a[i]+=b[i];
}

std::vector<Real>  operator*(const Real a,const std::vector<Real> & c)
{
  std::vector<Real> temp;
  for (auto x:c)
    {
      temp.push_back(x*a);
    }
  return temp;
}

Real Pow(Real x,int y)
{
  Real s=1;
  for (int i=1;i<=y;++i)
    {
      s*=x;
    }
  return s;
}

std::vector<Real>  operator-(const std::vector<Real> & a,const std::vector<Real> & b)
{
  std::vector<Real> temp;
  int n=a.size();
  for (int i=0;i<n;++i)
    {
      temp.push_back(a[i]-b[i]);
    }
  return temp;
}

std::vector<Real>  operator+(const std::vector<Real> & a,const std::vector<Real> & b)
{
  std::vector<Real> temp;
  int n=a.size();
  for (int i=0;i<n;++i)
    {
      temp.push_back(a[i]+b[i]);
    }
  return temp;
}

Real * VV_ADD(const std::vector<Real> &a,const std::vector<Real> &b)
{
  int n=a.size();
  Real * temp=new Real[n];
  for (int i=0;i<n;++i)
    temp[i]=a[i]+b[i];
  return temp;
}

void  VP_ADD(std::vector<Real> &a,const Real * b)
{
  int n=a.size();
  for (int i=0;i<n;++i)
    a[i]+=b[i];
}
