#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include "simulation.h"

constexpr double J_mean = 0.0;                                    // Mean value of J_ij
constexpr double eps = 1e-10;                                     // Finite float-point precision
constexpr double avalanche_statistics_bound = 0.3;                // Avalanche statistics is saved only when magnetization is lower than this bound


SKSimulation::SKSimulation(size_t n, double tau, double H_rate, std::vector<std::vector<double>> J,
                           std::ostream &outstream, std::ostream &avout) : output_(outstream), 
                           avalanche_output_(avout), exp_distr_(1.0 / tau), engine_(Seed()),J_(J), 
                           n_(n), tau_(tau), H_rate_(H_rate) {}

SKSimulation::SKSimulation(size_t n, double tau, double H_rate, std::ostream &out, std::ostream &avout) : 
                            SKSimulation(n, tau, H_rate, generate_exchange_matrix(n), out, avout) {}

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
  H_max_ = H_;
  std::fill(spins_.begin(), spins_.end(), -1);
  H_ = -2.99;
  do {
    H_ -= 0.01;
    RecalculateLocalFields();
  } while (FindUnstable().size() > 0);
  m_ = -1;
  H_min_ = H_;
  H_max_ = -H_;
}

void SKSimulation::HysteresisLoop(size_t loops_) {
  for (size_t loop = 0; loop < loops_; ++loop) {
    Advance(H_max_);
    Advance(H_min_);
  }
}

unsigned int SKSimulation::Seed() {
  std::random_device rd;
  return rd();
}

std::vector<std::vector<double>> SKSimulation::generate_exchange_matrix(size_t n) {
  std::vector<std::vector<double>> ret(n, std::vector<double>(n, 0.0));
  std::random_device rd;
  std::mt19937 generator(rd());
  std::normal_distribution<double> distr(J_mean, sqrt(1.0 / static_cast<double>(n)));
  for (size_t i = 0; i < ret.size(); ++i) {
    for (size_t j = i + 1; j < ret.size(); ++j) {
      ret[i][j] = distr(generator);
      ret[j][i] = ret[i][j];
    }
  }
  return ret;
}

void SKSimulation::RecalculateLocalFields() {
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

void SKSimulation::Advance(double H_end) {
  output_ << H_ << '\t' << m_ << '\n';
  if ((H_end - H_) * H_rate_ < 0.0) {
    H_rate_ *= -1.0;
  }
  while ((H_end - H_) * H_rate_ > 0.0) {
    Avalanche();
    if ((H_end - H_) * H_rate_ < 0.0) {
      break;
    }

    double min_field_to_flip = fabs((H_max_ - H_min_) / H_rate_);
    bool found = false;
    for (size_t i = 0; i < n_; ++i) {
      if (spins_[i] * H_rate_ < 0.0 && min_field_to_flip > -local_fields_[i] / H_rate_) {
        min_field_to_flip = -local_fields_[i] / H_rate_;
        found = true;
      }
    }
    if (!found) {
      min_field_to_flip = (H_end - H_) / H_rate_;
    }
#ifdef DEVBUG
    std::cout << min_field_to_flip << '\t' << H_ << '\t' << spins_[0] << '\t' << spins_[1] << '\n';
    if (min_field_to_flip < 0.0) {
      std::cout << "Negative Advance!!\nFields:\t" << H_end << '\t' << H_ << '\n';
    }
#endif
    StableAdvance(min_field_to_flip + eps);
  }
}

std::vector<size_t> SKSimulation::FindUnstable() {
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

void SKSimulation::StableAdvance(double delta_t) {
#ifdef DEVBUG
  output_ << "StableAdvance\t";
#endif
  double delta_H = delta_t * H_rate_;
  H_ += delta_H;
  for (size_t i = 0; i < n_; ++i) {
    local_fields_[i] += delta_H;
  }  
  output_ << H_ << '\t' << m_ << '\n';
}

void SKSimulation::Avalanche() {
  auto unstable = FindUnstable();
  if (unstable.size() == 0) return;
  int flips = 0;
  int signed_flips = 0;
  do {
    double min = exp_distr_(engine_);
    size_t min_index = unstable[0];
    for (size_t i = 1; i < unstable.size(); ++i) {
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
      local_fields_[i] += min * H_rate_ - 2 * flipped * J_[i][min_index];
    }
    H_ += min * H_rate_;
    output_ << H_ << '\t' << m_ << '\n';
    ++flips;
    signed_flips += 2 * flipped;
    unstable = FindUnstable();
  } while (unstable.size() > 0);
  if (fabs(m_) < avalanche_statistics_bound) {
    avalanche_output_ << flips << '\t' << signed_flips << '\n';
  }
}

/*
int main() {
  std::ofstream out("1.dat");
  SKSimulation sim(300, 0.01, 0.0005, out);
  sim.Initialize();
  sim.HysteresisLoop(1);
}
*/

