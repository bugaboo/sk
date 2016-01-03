#include "simulation.h"
#include <fstream>


int main() {
  std::ofstream out("3.dat", std::ios::out | std::ios::app);
  std::ofstream avout("avalan", std::ios::out | std::ios::app);
  for (int i = 0; i < 1; ++i) {
    SKSimulation sim(1600, 0.001, 0.05, out, avout);
    sim.Initialize();
    sim.HysteresisLoop(1);
  }
} 
