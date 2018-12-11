#include "compute_temperature.hh"
#include "my_types.hh"
#include "csv_reader.hh"
#include "csv_writer.hh"
#include "material_point.hh"
#include "material_points_factory.hh"
#include "system.hh"
#include <gtest/gtest.h>

/*****************************************************************/
// Fixture class
class HomogeneousTemperature : public ::testing::Test {
protected:
  void SetUp() override {
    MaterialPointsFactory::getInstance();
    std::vector<MaterialPoint> points;
    n_points = 16;
    for (UInt i = 0; i < n_points*n_points; ++i) {
      MaterialPoint p;
      p.getTemperature() = 1.;
      p.getHeatRate() = 0.;
      points.push_back(p);
    }

    for (auto &p : points) {
      system.addParticle(std::make_shared<MaterialPoint>(p));
    }
    
     temerature = std::make_shared<ComputeTemperature>();
  }

  System system;
  UInt n_points;
  std::shared_ptr<ComputeTemperature> temerature;
};

/*****************************************************************/
TEST_F(HomogeneousTemperature, homogeneous) {
    
  temerature->setDeltaT(5e-4);
  
  for (UInt i = 0; i < 10; ++i) {
    temerature->compute(system);
  }
  
  for (auto &p : system){
    auto &mpoint = static_cast<MaterialPoint&>(p);
    ASSERT_NEAR(mpoint.getTemperature(), 1., 1e-15);
  }

}

// Fixture class
class HeatSource : public ::testing::Test {
protected:
  void SetUp() override {
    MaterialPointsFactory::getInstance();
    std::vector<MaterialPoint> points;
    n_points = 64;
    for (UInt i = 0; i < n_points*n_points; ++i) {
      MaterialPoint p;
      p.getTemperature() = 0.;
      if ((i%n_points)==n_points/4){
        p.getHeatRate() = -1.;
    }else if((i%n_points)==3*n_points/4){
        p.getHeatRate() = +1.;
    }else{
        p.getHeatRate() = 0.;
    }
      points.push_back(p);
    }

    for (auto &p : points) {
      // std::cout << p << std::endl;
      system.addParticle(std::make_shared<MaterialPoint>(p));
    }
    
     temerature = std::make_shared<ComputeTemperature>();
  }

  System system;
  UInt n_points;
  std::shared_ptr<ComputeTemperature> temerature;
};

/*****************************************************************/
TEST_F(HeatSource, homogeneous) {
    
  temerature->setDeltaT(5e-5);
  
  for (UInt i = 0; i < 100000; ++i) {
    temerature->compute(system);
  }
  
  CsvWriter writer("tmp_file");
  writer.compute(system);
  
  for (UInt i = 0; i < n_points*n_points; ++i) {
    auto& mpoint = static_cast<MaterialPoint&>(system.getParticle(i));
    if ((i%n_points)<=n_points/4){
        ASSERT_NEAR(mpoint.getTemperature(), -1 - (1.0/(n_points - 1))*mpoint.getPosition()[0], 1e-4);
    }else if ((i%n_points)>=n_points/4 && (i%n_points)<=3*n_points/4){
        ASSERT_NEAR(mpoint.getTemperature(), (1.0/(n_points - 1))*mpoint.getPosition()[0], 1e-4);
    }else{
        ASSERT_NEAR(mpoint.getTemperature(), 1 - (1.0/(n_points - 1))*mpoint.getPosition()[0], 1e-4);
    }
}

}
