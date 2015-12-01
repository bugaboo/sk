#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>

constexpr double J_mean = 0.0;

class SKSimulation{
#ifdef DEVBUG
  friend class TestSKSimulation;  
#endif
 public:
  int a;
  SKSimulation(size_t n, double tau, double H_rate, std::ostream &out = std::cout);
  SKSimulation(size_t n, double tau, double H_rate, std::vector<std::vector<double>> J,
               std::ostream &outstream);
  void Initialize();
  void PrintJ();
  void HysteresisLoop(size_t loops);  
 private:
  unsigned int Seed_();
  void RecalculateLocalFields_();
  std::vector<std::vector<double>> generate_exchange_matrix_(size_t n);
  std::vector<size_t> FindUnstable_();
  void Advance_(double H_end);
  void StableAdvance_(int n_steps);
  void Avalanche_();
  
  std::ostream &output_;
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

SKSimulation::SKSimulation(size_t n, double tau, double H_rate, std::vector<std::vector<double>> J,
                           std::ostream &outstream) : output_(outstream),exp_distr_(tau), 
                            engine_(Seed_()),J_(J), n_(n), tau_(tau), H_rate_(H_rate) {}

SKSimulation::SKSimulation(size_t n, double tau, double H_rate, std::ostream &out) : 
                            SKSimulation(n, tau, H_rate, generate_exchange_matrix_(n), out) {}

unsigned int SKSimulation::Seed_() {
  std::random_device rd;
  return rd();
}

std::vector<std::vector<double>> SKSimulation::generate_exchange_matrix_(size_t n) {
  std::vector<std::vector<double>> ret(n, std::vector<double>(n, 0.0));
  std::random_device rd;
  std::mt19937 generator(rd());
  std::normal_distribution<double> distr(J_mean, sqrt(1.0 / static_cast<double>(n)));
  for (size_t i = 0; i < ret.size(); ++i) {
    for (size_t j = 0; j < ret.size(); ++j) {
      ret[i][j] = i == j ? 0.0 : distr(generator) / 2.0;
    }
  }
  return ret;
}

void SKSimulation::PrintJ() {
  for (size_t i = 0; i < J_.size(); ++i) {
    for (size_t j = 0; j < J_[i].size(); ++j) {
      std::cout << J_[i][j] << '\t';
    }
    std::cout << std::endl;
  }
  std::cout << H_min_ << '\t' << H_max_ << '\n';
}

void SKSimulation::Initialize() {
  spins_.resize(n_);
  local_fields_.resize(n_);
  std::fill(spins_.begin(), spins_.end(), 1);
  H_ = -1.0;
  do {
    H_ += 0.01;
    RecalculateLocalFields_();
  } while (FindUnstable_().size() > 0);
  H_max_ = H_;
  std::fill(spins_.begin(), spins_.end(), -1);
  H_ = 1.0;
  do {
    H_ -= 0.01;
    RecalculateLocalFields_();
  } while (FindUnstable_().size() > 0);
  m_ = -1;
  H_min_ = H_;
}

void SKSimulation::HysteresisLoop(const size_t loops_) {
  for (size_t loop = 0; loop < loops_; ++loop) {
    Advance_(H_max_);
//    Advance_(H_min_);
  }
}

void SKSimulation::Advance_(double H_end) {
  if ((H_end - H_) * H_rate_ < 0.0) {
    H_rate_ *= -1.0;
  }
  while ((H_end - H_) * H_rate_ > 0.0) {
    Avalanche_();
    int index = 0;
    while (index < n_ && spins_[index] * H_rate_ > 0.0) {
      ++index;
    }
    if (index < n_) {
      double min_field_to_flip = local_fields_[index] * spins_[index];
      while (index < n_) {
        if (spins_[index] * H_rate_ > 0.0 && 
            min_field_to_flip < local_fields_[index] * spins_[index]) {
          min_field_to_flip = local_fields_[index] * spins_[index];
        }
        ++index;
      }
#ifdef DEVBUG
      std::cout << min_field_to_flip << '\t' << H_ << '\t' << spins_[0] << '\t' << spins_[1] << '\n';
      if (min_field_to_flip < 0.0) {
        std::cout << "Negative Advance!!\n\n";
      }
#endif
      StableAdvance_((int)(min_field_to_flip / fabs(H_rate_) / tau_) + 1);
    }
  }
}

void SKSimulation::RecalculateLocalFields_() {
  if (local_fields_.size() != n_) {
    std::cerr << "Wrong size of local fields array!\n";
  }
  if (J_.size() != n_) {
    std::cerr << "Wrong size of exchange matrix\n";
  }
  for (size_t i = 0; i < n_; ++i) {
    local_fields_[i] = H_;
    for (size_t j = 0; j < n_; ++j) {
      local_fields_[i] += J_[i][j] * (double)(spins_[j]);
    }
  }
}

std::vector<size_t> SKSimulation::FindUnstable_() {
  std::vector<size_t> ret;
  if (spins_.size() != n_ || local_fields_.size() != n_) {
    std::cerr << "Array size problem in find unstable\n";
    return ret;
  }
  for (size_t i = 0; i < n_; ++i) {
    if ((double)(spins_[i]) * local_fields_[i] < 0) {
      ret.push_back(i);
    }
  }
  return ret;
}

void SKSimulation::StableAdvance_(int n_steps) {
  double delta_H = n_steps * tau_ * H_rate_;
  H_ += delta_H;
  for (size_t i = 0; i < n_; ++i) {
    local_fields_[i] += delta_H;
  }  
  output_ << H_ << '\t' << m_ << '\n';
}

void SKSimulation::Avalanche_() {
  auto unstable = FindUnstable_();
  if (unstable.size() == 0) return;
  do {
    double min = 1000 * tau_;
    size_t min_index = 0;
    for (size_t i = 0; i < unstable.size(); ++i) {
      double tmp = exp_distr_(engine_);
      if (tmp < min) {
        min = tmp;
        min_index = unstable[i];
      }
    }
    auto flipped = spins_[min_index];
    m_ -= 2.0 * (double)flipped / n_;
    spins_[min_index] *= -1;
    for (size_t i = 0; i < n_; ++i) {
      local_fields_[i] += tau_ * H_rate_ - 2 * flipped * J_[i][min_index];
    }
    H_ += tau_ * H_rate_;
    output_ << H_ << '\t' << m_ << '\n';
    unstable = FindUnstable_();
  } while (unstable.size() > 0);
}


int main() {
  std::ofstream out("1.dat");
  SKSimulation sim(10, 0.1, 0.001, out);
  sim.Initialize();
  sim.PrintJ();
  sim.HysteresisLoop(1);
}
