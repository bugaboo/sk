#include "simulation.h"
#include <fstream>

int main() {
  std::ofstream out("2.dat");
  SKSimulation sim(30, 0.0001, 0.1, out);
//  SKSimulation sim(2, 0.01, 0.1, {{0.0, 0.1764}, {-0.4201, 0.0}}, out);
  sim.Initialize();
  sim.PrintJ();
  sim.HysteresisLoop(1);
}
