#include "compute_temperature.hh"
#include "fft.hh"
#include "material_point.hh"
#include <cmath>

/* -------------------------------------------------------------------------- */

void ComputeTemperature::compute(System &system) {
    
  UInt size2 = system.getNbParticles();
  UInt size = sqrt(size2);
  
  Matrix<complex> T(size); // Temperature
  Matrix<std::complex<int>> K(size); // Frequency
  Matrix<complex> h(size); // Heat transfer
  Matrix<complex> dTdt(size); // Temperature
 
  for (auto &&entry : index(T)) {
    int i = std::get<0>(entry);
    int j = std::get<1>(entry);
    auto &val = std::get<2>(entry);
    auto &p = static_cast<MaterialPoint &>(system.getParticle(j * size + i));
    val = p.getTemperature();
  }
  for (auto &&entry : index(h)) {
    int i = std::get<0>(entry);
    int j = std::get<1>(entry);
    auto &val = std::get<2>(entry);
    auto &p = static_cast<MaterialPoint &>(system.getParticle(j * size + i));
    val = p.getHeatRate();
  }
  
  K = FFT::computeFrequencies(size);
  T = FFT::transform(T);
  h = FFT::transform(h);
  
  for (UInt j = 0; j<size; j++){
      for (UInt i = 0; i<size; i++){
        double k = std::abs(K(i,j));
        dTdt(i,j) = h(i,j) - kappa*M_PI*M_PI*T(i,j)*k*k;
      }
  }
  
  dTdt = FFT::itransform(dTdt);
  T = FFT::itransform(T);
  for (UInt j = 0; j<size; j++){
      for (UInt i = 0; i<size; i++){
        T(i,j) = T(i,j) + dt*dTdt(i,j);
        auto &p = static_cast<MaterialPoint &>(system.getParticle(j * size + i));
        p.getTemperature() = T(i,j).real();
      }
  }

}

void ComputeTemperature::setDeltaT(Real dt){ this->dt = dt;}
void ComputeTemperature::setk(Real k){ this->kappa = k;}

/* -------------------------------------------------------------------------- */
