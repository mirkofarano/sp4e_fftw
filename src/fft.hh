#ifndef FFT_HH
#define FFT_HH
/* ------------------------------------------------------ */
#include "matrix.hh"
#include "my_types.hh"
#include <fftw3.h>
/* ------------------------------------------------------ */

struct FFT {

  static Matrix<complex> transform(Matrix<complex>& m);
  static Matrix<complex> itransform(Matrix<complex>& m);

  static Matrix<std::complex<int>> computeFrequencies(int size);

};

/* ------------------------------------------------------ */

inline Matrix<complex> FFT::transform(Matrix<complex>& m_in) {
    
    fftw_plan p;
    fftw_complex *in, *out;
    UInt Ny = m_in.rows();
    UInt Nx = m_in.cols();
    Matrix<complex> m_out(m_in.size());
    
    in = reinterpret_cast<fftw_complex*>(m_in.data());
    out = reinterpret_cast<fftw_complex*>(m_out.data());
    p = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p);
    m_out/=Ny*Nx; // Normalization
    fftw_destroy_plan(p);
    return m_out;
    
}

/* ------------------------------------------------------ */

inline Matrix<complex> FFT::itransform(Matrix<complex>& m_in) {
        
    fftw_plan p;
    fftw_complex *in, *out;
    UInt Ny = m_in.rows();
    UInt Nx = m_in.cols();
    Matrix<complex> m_out(m_in.size());

    in = reinterpret_cast<fftw_complex*>(m_in.data());
    out = reinterpret_cast<fftw_complex*>(m_out.data());
    p = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(p);
    fftw_destroy_plan(p);
    return m_out;
}

/* ------------------------------------------------------ */


/* ------------------------------------------------------ */

inline Matrix<std::complex<int>> FFT::computeFrequencies(int size) {
    
    Matrix<std::complex<int>> m(size);
    
      for (auto&& entry : index(m)) {
        int i = std::get<0>(entry);
        int j = std::get<1>(entry);
        auto& val = std::get<2>(entry);
        int kx = (i <= size / 2) ? i : i - size;
        int ky = (j <= size / 2) ? j : j - size;
        val = std::complex<int>(kx,ky);
        std::cout << val << std::endl;
    }

    return m;
}

#endif  // FFT_HH
