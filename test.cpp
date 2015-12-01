#ifndef DEVBUG
#define DEVBUG
#endif
#include "simulation.h"
#include "nullbuffer.h"

constexpr double epsilon = 1e-13;

bool compare(const std::vector<double> &a, const std::vector<double> &b) {
  if (a.size() != b.size()) {
    return false;
  }
  bool ret = true;
  for (size_t i = 0; i < a.size(); ++i) {
    if (fabs(a[i] - b[i]) > epsilon) {
      ret = false;
      break;
    }
  }
  return ret;
}

class TestSKSimulation {
 public:
  bool RunTests() {
    return TestAvalanche_();
  }
 private:
  bool TestAvalanche_(){
    bool ret = true;
    NullBuffer null_buf;
    std::ostream out(&null_buf);
    SKSimulation sim(2, 0.1, 0.001, {{0.0, 0.5}, {-0.5, 0.0}}, out);
    sim.spins_.resize(2, -1);
    sim.local_fields_.resize(2);
    sim.H_ = 0.0;
    sim.m_ = -1;
    sim.RecalculateLocalFields_();
    sim.Avalanche_();
    auto a = std::vector<double>{1.0002, 0.0002};
    if (!compare(sim.local_fields_, std::vector<double>{1.0002, 0.0002})) {
      std::cout << "Avalanche fail\n";
      ret = false;
    }
    return ret;
  }
};

int main() {
  TestSKSimulation test;
  return test.RunTests();
}
