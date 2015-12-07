#include "simulation.h"
#include <fstream>

int main() {
  std::ofstream out("3.dat");
  SKSimulation sim(10, 0.001, 0.05, out);
  sim.Initialize();
  sim.HysteresisLoop(1);
} 
