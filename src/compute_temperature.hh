#ifndef __COMPUTE_TEMPERATURE__HH__
#define __COMPUTE_TEMPERATURE__HH__

/* -------------------------------------------------------------------------- */
#include "compute.hh"

//! Compute contact interaction between ping-pong balls
class ComputeTemperature : public Compute {

  // Virtual implementation
public:
  //! Penalty contact implementation
  void compute(System &system) override;
  void setDeltaT(Real dt);
  void setk(Real k);
private:
  Real dt;
  // conductivity,
  Real kappa = 1.;
};

/* -------------------------------------------------------------------------- */
#endif //__COMPUTE_TEMPERATURE__HH__
