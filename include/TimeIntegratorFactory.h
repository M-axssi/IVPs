#ifndef  TIMEINTEGRATOR_FACTORY_H
#define TIMEINTEGRATOR_FACTORY_H
#include "TimeIntegrator.h"
#include "Adams_Bashforth.h"
#include "Adams_Moulton.h"
#include "Bdfs.h"
#include "Runge_Kutta.h"

 class TimeIntegratorFactory
  {
  public:
    typedef TimeIntegrator * (*TimeIntegratorCreator) ();
    static TimeIntegratorFactory & instance()
    {
      static TimeIntegratorFactory instance;
      return instance;
    }
    bool RegisterTimeIntegrator(const std::string &,TimeIntegratorCreator);
    bool UnregisterTimeIntegrator(const std::string &);
    TimeIntegrator * CreateTimeIntegrator(const std::string &);
  private:
    TimeIntegratorFactory()=default;
    ~TimeIntegratorFactory()=default;
    std::map<std::string,TimeIntegratorCreator> Creator_Map;
  };

namespace 
{
  TimeIntegrator * Create_Adams_Bashforth()
  {
    return new Adams_Bashforth();
  };

  TimeIntegrator* Create_Adams_Moulton()
  {
    return new Adams_Moulton();
  };

  TimeIntegrator* Create_Bdfs()
  {
    return new Bdfs();
  };

  TimeIntegrator* Create_Runge_Kutta()
  {
    return new Runge_Kutta();
  };

  bool _Adams_Bashforth=TimeIntegratorFactory::instance().RegisterTimeIntegrator("Adams_Bashforth",Create_Adams_Bashforth);
  bool _Adams_Moulton =TimeIntegratorFactory::instance().RegisterTimeIntegrator("Adams_Moulton",Create_Adams_Moulton);
  bool _Bdfs=TimeIntegratorFactory::instance().RegisterTimeIntegrator("Bdfs",Create_Bdfs);
  bool _Runge_Kutta=TimeIntegratorFactory::instance().RegisterTimeIntegrator("Runge_Kutta",Create_Runge_Kutta);
};

#else
//Do nothing!
#endif
