#ifndef SK_SIMULATION_H
#define SK_SIMULATION_H

#include <vector>
#include <random>
#include <iostream>

class SKSimulation{
#ifdef DEVBUG
  friend class TestSKSimulation;  
#endif
 public:
  int a;
  SKSimulation(size_t n, double tau, double H_rate, std::ostream &out = std::cout,
                                                  std::ostream &avout = std::cout);
  SKSimulation(size_t n, double tau, double H_rate, std::vector<std::vector<double>> J,
               std::ostream &outstream, std::ostream &avout);
  void Initialize();
  void PrintJ();
  void HysteresisLoop(size_t loops);  
 private:
  unsigned int Seed();
  void RecalculateLocalFields();
  std::vector<std::vector<double>> generate_exchange_matrix(size_t n);
  std::vector<size_t> FindUnstable();
  void Advance(double H_end);
  void StableAdvance(double delta_t);
  void Avalanche();
  
  std::ostream &output_;
  std::ostream &avalanche_output_;
  std::exponential_distribution<double> exp_distr_;
  std::mt19937 engine_;
  const std::vector<std::vector<double>> J_;                          // exchange constants J_{ij} / 2
  std::vector<double> local_fields_;                            // on-site fields h_i
  std::vector<int> spins_;                                       // spin states
  const size_t n_;                                                    // number of spins
  double tau_;                                                  // time step
  double H_rate_;                                               // rate of magnetic field change
  double H_;                                                    // current magnetic field  
  double H_max_;                                                // "all up" magnetic field
  double H_min_;                                                // "all down" magnetic field
  int M_;                                                       // total magnetization
  double m_;                                                    // mean magnetization

  
};
#endif
