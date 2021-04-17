#include "TimeIntegratorFactory.h"

bool TimeIntegratorFactory::RegisterTimeIntegrator(const std::string &s,TimeIntegratorCreator c)
{
  return Creator_Map.insert(std::make_pair(s,c)).second;
}

bool TimeIntegratorFactory::UnregisterTimeIntegrator(const std::string & s)
{
  return Creator_Map.erase(s)==1;
}

TimeIntegrator * TimeIntegratorFactory::CreateTimeIntegrator( const std::string & s)
{
  auto i=Creator_Map.find(s);
  if (i==Creator_Map.cend())
    {
      throw std::runtime_error("Unknow TimeIntegrator ID.");
    }else
    {
      return (i->second)();
    }
}
