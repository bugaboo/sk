#ifndef DEVBUG
#define DEVBUG
#endif
#include "simulation.h"
#include "nullbuffer.h"
#include <fstream>

constexpr double epsilon = 1e-13;
NullBuffer null_buf;
std::ostream out(&null_buf);

template <typename T>
bool compare(const std::vector<T> &a, const std::vector<T> &b) {
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
    bool ret = true;
    ret = ret && TestAvalanche_();
    ret = ret && TestAdvance_();
    return ret;
  }
 private:
  bool TestAdvance_() {
    constexpr size_t n = 10;
    std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));
    std::ifstream ifile("jmatrix");
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        ifile >> J[i][j];
      }
    }
    SKSimulation sim(n, 0.01, 0.1, J, out);
    sim.Initialize();
/*    sim.H_ = 0.634;
    sim.spins_ = {-1, 1, 1, -1, 1, 1, 1, 1, 1, 1};
    sim.m_ = 0.6;
    sim.RecalculateLocalFields_();
 */   std::cout << "H_=" << sim.H_ << "\tH_max=" << sim.H_max_ << '\n';
    sim.Advance_(sim.H_max_);
    if (!compare(sim.spins_, std::vector<int>(n, 1))) {
      std::cout << "Advance fail\n";
      for (auto elem: sim.spins_) {
        std::cout << elem << '\t';
      }
      std::cout << std::endl;
      return false;
    }
    return true;
  }

  bool TestAvalanche_(){
    bool ret = true;
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
