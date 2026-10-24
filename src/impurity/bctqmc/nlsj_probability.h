#ifndef BCTQMC_NLSJ_PROBABILITY_H
#define BCTQMC_NLSJ_PROBABILITY_H

// Optional angular-momentum labels for the already existing Prob_eig output.
// atomh.py writes <cix>.nlsj in exactly the direct basis/order of <cix>.
// This reader runs only on the master rank and never changes the solver's
// propagators, Hamiltonian eigenvectors, or Probability.dat.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace nlsj_probability {

struct Label {
  int n = 0;
  double energy = 0.0; // eigenvalue of the shifted solver Hamiltonian, before lambda0
  double mean_n = 0.0; // <N> in this Hamiltonian eigenstate
  double mean_pairs = 0.0; // <N(N-1)/2> in this eigenstate
  int twice_l = 0;
  int twice_s = 0;
  int twice_j = 0;
  double weight = 0.0; // overlap with the entire (N,L,S,J) subspace
  bool available = false;
};

struct Labels {
  std::string filename;
  std::vector<std::vector<Label>> block;
  bool has_sidecar = false;
};

inline std::string fnv1a_file(const std::string& filename)
{
  std::ifstream input(filename, std::ios::binary);
  if (!input) throw std::runtime_error("Cannot open CIX for NLSJ fingerprint: " + filename);
  std::uint64_t value = UINT64_C(14695981039346656037);
  char buffer[1 << 16];
  while (input) {
    input.read(buffer, sizeof(buffer));
    for (std::streamsize i = 0; i < input.gcount(); ++i) {
      value ^= static_cast<unsigned char>(buffer[i]);
      value *= UINT64_C(1099511628211);
    }
  }
  std::ostringstream result;
  result << std::hex << std::setfill('0') << std::setw(16) << value;
  return result.str();
}

inline void require(bool condition, const std::string& message)
{
  if (!condition) throw std::runtime_error("Invalid NLSJ sidecar: " + message);
}

inline std::string next_data_line(std::istream& input)
{
  // Readable sidecars follow CIX: numeric fields first, comments after '#'.
  std::string line;
  while (std::getline(input, line)) {
    if (const auto hash = line.find('#'); hash != std::string::npos) line.resize(hash);
    const auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos) continue;
    const auto last = line.find_last_not_of(" \t\r");
    return line.substr(first, last - first + 1);
  }
  return {};
}

inline int twice_quantum_number(double value)
{
  const int doubled = static_cast<int>(std::lround(2.0 * value));
  require(doubled >= 0 && std::abs(0.5 * doubled - value) < 1e-8,
          "L, S, or J is not an integer or half-integer");
  return doubled;
}

inline double energy_shift(const Cix& cix, int n)
{
  // boldc.cc applies this diagonal shift after Read_cix. In the atomic
  // fixed-N blocks used by V3 it is a scalar and cannot change eigenvectors.
  const double u = cix.Nsite == 1 ? cix.F0_input : 0.0;
  return u * n * (n - 1) / 2.0 - cix.mu_input * n;
}

inline void attach_block_labels(Labels& labels, int ib, int dimension,
                                const std::vector<double>& basis,
                                const std::vector<std::tuple<int,int,int,int>>& quantum_numbers,
                                const std::vector<double>& h_eigenvectors)
{
  // Each column of basis is one NLSJ vector in the CIX direct basis. Six-digit
  // decimal coefficients are intentionally approximate; the complete set of
  // vectors must nevertheless be nearly orthonormal.
  for (int a = 0; a < dimension; ++a) {
    for (int b = a; b < dimension; ++b) {
      double overlap = 0.0;
      for (int row = 0; row < dimension; ++row)
        overlap += basis[row * dimension + a] * basis[row * dimension + b];
      require(std::abs(overlap - (a == b ? 1.0 : 0.0)) < 1e-4,
              "NLSJ vectors are not orthonormal in block " + std::to_string(ib + 1));
    }
  }
  for (int eig = 0; eig < dimension; ++eig) {
    std::map<std::tuple<int,int,int,int>, double> weights;
    for (int col = 0; col < dimension; ++col) {
      double amplitude = 0.0;
      for (int row = 0; row < dimension; ++row)
        amplitude += h_eigenvectors[eig * dimension + row] * basis[row * dimension + col];
      weights[quantum_numbers[col]] += amplitude * amplitude;
    }
    double total = 0.0, largest = -1.0;
    std::tuple<int,int,int,int> best;
    for (const auto& [qn, weight] : weights) {
      total += weight;
      if (weight > largest) { largest = weight; best = qn; }
    }
    require(std::abs(total - 1.0) < 1e-4, "NLSJ weights do not sum to one");
    auto& target = labels.block[ib][eig];
    target.n = std::get<0>(best);
    target.twice_l = std::get<1>(best);
    target.twice_s = std::get<2>(best);
    target.twice_j = std::get<3>(best);
    target.weight = largest / total; // Pure states cannot print above one.
    target.available = true;
  }
}

inline void read_v2(std::istream& input, const Cix& cix, Labels& labels,
                    const std::vector<std::vector<double>>& h_eigenvectors)
{
  int natst = 0, nblock = 0, norbitals = 0;
  std::string fingerprint;
  std::istringstream header(next_data_line(input));
  header >> natst >> nblock >> norbitals >> fingerprint;
  require(bool(header) && natst == cix.Natst && nblock == cix.Nblock &&
          norbitals == cix.Npsi, "V2 counts differ from CIX");
  require(fingerprint == fnv1a_file(cix.cix_filename), "CIX fingerprint differs; copy the matching sidecar");

  std::istringstream dimensions(next_data_line(input));
  for (int ib = 0; ib < cix.Nblock; ++ib) {
    int dimension = 0;
    dimensions >> dimension;
    require(bool(dimensions) && dimension == cix.Block[ib].size(),
            "V2 block dimensions differ from CIX");
  }
  std::string extra;
  require(!(dimensions >> extra), "extra data on V2 block-dimensions line");
  for (int orbital = 0; orbital < norbitals; ++orbital) {
    int index = -1;
    double j = 0.0, mj = 0.0;
    std::istringstream row(next_data_line(input));
    row >> index >> j >> mj;
    require(bool(row) && index == orbital, "invalid V2 orbital table");
  }

  int next_state = 1;
  for (int ib = 0; ib < cix.Nblock; ++ib) {
    int file_block = 0, n = -1, dimension = 0;
    double jz = 0.0;
    std::istringstream block(next_data_line(input));
    block >> file_block >> n >> jz >> dimension;
    require(bool(block) && file_block == ib + 1 &&
            dimension == cix.Block[ib].size(), "invalid V2 block header");
    // Coefficient positions use the exact CIX block order. The human-readable
    // occupation annotations list only dominant components and are comments,
    // so they never replace the complete vectors needed for projection.
    for (int row = 0; row < dimension; ++row) {
      const int direct_state = cix.Block[ib](row);
      require(cix.iN(direct_state - 1) == n &&
              std::abs(cix.iSz(direct_state - 1) - jz) < 1e-8,
              "V2 N or Jz differs from CIX");
    }
    std::string line = next_data_line(input);
    require(!line.empty(), "truncated V2 block");
    if (line == "UNLABELED") continue;
    std::vector<double> basis(dimension * dimension);
    std::vector<std::tuple<int,int,int,int>> quantum_numbers(dimension);
    for (int col = 0; col < dimension; ++col) {
      if (col > 0) line = next_data_line(input);
      std::istringstream vector(line);
      int global_state = 0, vector_block = 0, local = -1, vector_n = -1, vector_dim = 0;
      double l = 0.0, s = 0.0, j = 0.0;
      vector >> global_state >> vector_block >> local >> vector_n >> l >> s >> j >> vector_dim;
      require(bool(vector) && global_state == next_state++ && vector_block == ib + 1 &&
              local == col && vector_n == n && vector_dim == dimension,
              "invalid V2 NLSJ state header");
      quantum_numbers[col] = {n, twice_quantum_number(l),
                              twice_quantum_number(s), twice_quantum_number(j)};
      for (int row = 0; row < dimension; ++row)
        vector >> basis[row * dimension + col];
      require(bool(vector), "truncated V2 NLSJ eigenvector");
      require(!(vector >> extra), "extra data in V2 NLSJ eigenvector");
    }
    attach_block_labels(labels, ib, dimension, basis, quantum_numbers, h_eigenvectors[ib]);
  }
  require(next_data_line(input).empty(), "unexpected data after final V2 block");
  labels.has_sidecar = true;
}

inline void read_v3(std::istream& input, const Cix& cix, Labels& labels,
                    const std::vector<std::vector<double>>& h_eigenvectors)
{
  int natst = 0, nblock = 0, norbitals = 0;
  std::string fingerprint;
  std::istringstream header(next_data_line(input));
  header >> natst >> nblock >> norbitals >> fingerprint;
  require(bool(header) && natst == cix.Natst && nblock == cix.Nblock &&
          norbitals == cix.Npsi, "V3 counts differ from CIX");
  require(fingerprint == fnv1a_file(cix.cix_filename),
          "CIX fingerprint differs; copy the matching sidecar");

  std::istringstream dimensions(next_data_line(input));
  for (int ib = 0; ib < cix.Nblock; ++ib) {
    int dimension = 0;
    dimensions >> dimension;
    require(bool(dimensions) && dimension == cix.Block[ib].size(),
            "V3 block dimensions differ from CIX");
  }
  std::string extra;
  require(!(dimensions >> extra), "extra data on V3 block-dimensions line");

  // The first section is intentionally a compact table of physical energy
  // eigenstates. A single dominant NLSJ label is printed for each row, with
  // its weight making an ambiguous or mixed assignment visible.
  std::vector<std::vector<Label>> saved(cix.Nblock);
  int global_state = 1;
  for (int ib = 0; ib < cix.Nblock; ++ib) {
    const int dimension = cix.Block[ib].size();
    saved[ib].resize(dimension);
    const int block_n = cix.iN(cix.Block[ib](0) - 1);
    for (int row = 1; row < dimension; ++row)
      require(cix.iN(cix.Block[ib](row) - 1) == block_n,
              "V3 sidecar requires a fixed-N Hamiltonian block");
    for (int eig = 0; eig < dimension; ++eig, ++global_state) {
      int state = 0, block = 0, n = -1;
      double l = 0.0, s = 0.0, j = 0.0, jz = 0.0, energy = 0.0, weight = 0.0;
      std::istringstream row(next_data_line(input));
      row >> state >> block >> n >> l >> s >> j >> jz >> energy >> weight;
      require(bool(row) && state == global_state && block == ib + 1 &&
              n == labels.block[ib][eig].n &&
              std::abs(jz - cix.iSz(cix.Block[ib](eig) - 1)) < 1e-8,
              "invalid V3 eigenstate row");
      require(!(row >> extra), "extra data on V3 eigenstate row");
      const double reference_energy = labels.block[ib][eig].energy -
                                      energy_shift(cix, n);
      require(std::abs(energy - reference_energy) <
              1e-6 * (1.0 + std::abs(reference_energy)),
              "V3 energy spectrum differs from CIX Hamiltonian");
      require(std::isfinite(weight) && weight >= 0.0 && weight <= 1.0001,
              "invalid V3 NLSJ label weight");
      auto& target = saved[ib][eig];
      target.n = n;
      target.energy = energy;
      target.weight = weight;
      if (l == -1.0 && s == -1.0 && j == -1.0) {
        require(weight == 0.0, "unlabeled V3 state has nonzero weight");
      } else {
        target.twice_l = twice_quantum_number(l);
        target.twice_s = twice_quantum_number(s);
        target.twice_j = twice_quantum_number(j);
        target.available = true;
      }
    }
  }

  global_state = 1;
  for (int ib = 0; ib < cix.Nblock; ++ib) {
    const int dimension = cix.Block[ib].size();
    std::vector<double> vectors(dimension * dimension);
    for (int eig = 0; eig < dimension; ++eig, ++global_state) {
      int state = 0;
      std::istringstream row(next_data_line(input));
      row >> state;
      require(bool(row) && state == global_state, "invalid V3 eigenvector state index");
      for (int component = 0; component < dimension; ++component)
        row >> vectors[eig * dimension + component];
      require(bool(row) && !(row >> extra), "invalid V3 eigenvector components");
    }
    for (int a = 0; a < dimension; ++a) {
      for (int b = a; b < dimension; ++b) {
        double overlap = 0.0;
        for (int row = 0; row < dimension; ++row)
          overlap += vectors[a * dimension + row] * vectors[b * dimension + row];
        require(std::abs(overlap - (a == b ? 1.0 : 0.0)) < 1e-4,
                "V3 eigenvectors are not orthonormal in block " + std::to_string(ib + 1));
      }
    }
    for (int eig = 0; eig < dimension; ++eig) {
      for (int row = 0; row < dimension; ++row) {
        double h_times_vector = 0.0;
        for (int col = 0; col < dimension; ++col)
          h_times_vector += cix.Ham[ib](row,col) * vectors[eig * dimension + col];
        const double shifted_energy = saved[ib][eig].energy +
                                      energy_shift(cix, saved[ib][eig].n);
        const double residual = h_times_vector -
          shifted_energy * vectors[eig * dimension + row];
        require(std::abs(residual) < 1e-5 * (1.0 + std::abs(shifted_energy)),
                "V3 eigenvector does not match its CIX energy in block " +
                std::to_string(ib + 1));
      }
    }
    // This projection also handles rotations within an exactly degenerate
    // eigenspace: vectors with the same label add before choosing the largest
    // weight. Information below the one dominant label per sidecar row is
    // intentionally unavailable, so mixed-state weights remain conservative.
    for (int eig = 0; eig < dimension; ++eig) {
      std::map<std::tuple<int,int,int,int>, double> weights;
      for (int col = 0; col < dimension; ++col) {
        const auto& source = saved[ib][col];
        if (!source.available) continue;
        double amplitude = 0.0;
        for (int row = 0; row < dimension; ++row)
          amplitude += h_eigenvectors[ib][eig * dimension + row] *
                       vectors[col * dimension + row];
        weights[{source.n, source.twice_l, source.twice_s, source.twice_j}] +=
          amplitude * amplitude * source.weight;
      }
      auto& target = labels.block[ib][eig];
      for (const auto& [qn, weight] : weights) {
        if (!target.available || weight > target.weight) {
          target.n = std::get<0>(qn);
          target.twice_l = std::get<1>(qn);
          target.twice_s = std::get<2>(qn);
          target.twice_j = std::get<3>(qn);
          target.weight = std::min(1.0, weight);
          target.available = true;
        }
      }
    }
  }
  require(next_data_line(input).empty(), "unexpected data after final V3 eigenvector");
  labels.has_sidecar = true;
}

inline Labels load(const Cix& cix)
{
  Labels labels;
  labels.filename = cix.cix_filename;
  labels.block.resize(cix.Nblock);
  std::vector<std::vector<double>> h_eigenvectors(cix.Nblock);
  for (int ib = 0; ib < cix.Nblock; ++ib) {
    const int dimension = cix.Block[ib].size();
    labels.block[ib].resize(dimension);
    for (int i = 0; i < dimension; ++i)
      labels.block[ib][i].n = cix.iN(cix.Block[ib](i) - 1);
    // Compute exactly the eigenvalue ordering used by the probability solver.
    // The cached vectors also permit NLSJ overlap when a sidecar is present.
    Array<double,2> h(cix.Ham[ib].copy());
    Array<double,1> energies(dimension), work(3 * dimension);
    require(xsyev(dimension, h.data(), dimension, energies.data(),
                  work.data(), work.size(), "V") == 0, "Hamiltonian diagonalization failed");
    h_eigenvectors[ib].resize(dimension * dimension);
    for (int eig = 0; eig < dimension; ++eig) {
      labels.block[ib][eig].energy = energies(eig);
      for (int row = 0; row < dimension; ++row) {
        h_eigenvectors[ib][eig * dimension + row] = h(eig,row);
        const int n = cix.iN(cix.Block[ib](row) - 1);
        const double weight = h(eig,row) * h(eig,row);
        labels.block[ib][eig].mean_n += weight * n;
        labels.block[ib][eig].mean_pairs += weight * n * (n - 1) / 2.0;
      }
    }
  }

  const std::string filename = cix.cix_filename + ".nlsj";
  if (!std::filesystem::exists(filename)) return labels; // N-only fallback
  std::ifstream input(filename);
  require(bool(input), "cannot open " + filename);
  std::string first_line;
  std::getline(input, first_line);
  if (first_line == "# NLSJ_EIGENSTATES_V3") {
    read_v3(input, cix, labels, h_eigenvectors);
    return labels;
  }
  if (first_line == "# NLSJ_BASIS_V2") {
    read_v2(input, cix, labels, h_eigenvectors);
    return labels;
  }
  input.clear();
  input.seekg(0);
  std::string token, fingerprint;
  int natst = 0, nblock = 0;
  input >> token >> natst >> nblock >> fingerprint;
  require(bool(input) && token == "NLSJ_BASIS_V1", "bad format/version in " + filename);
  require(natst == cix.Natst && nblock == cix.Nblock, "state or block count differs from CIX");
  require(fingerprint == fnv1a_file(cix.cix_filename), "CIX fingerprint differs; regenerate/copy the matching sidecar");

  for (int ib = 0; ib < cix.Nblock; ++ib) {
    int file_block = 0, dimension = 0, n = 0;
    input >> token >> file_block >> dimension >> n;
    require(bool(input) && token == "BLOCK" && file_block == ib + 1 &&
            dimension == cix.Block[ib].size(), "block identifier or dimension differs from CIX");
    for (int i = 0; i < dimension; ++i)
      require(labels.block[ib][i].n == n, "N differs from CIX in block " + std::to_string(ib + 1));

    input >> token;
    require(bool(input), "truncated block " + std::to_string(ib + 1));
    if (token == "UNLABELED") continue; // approximate-valence subspace

    // Matrix rows are CIX direct states and columns are orthonormal NLSJ
    // vectors. Components have six significant decimal digits. Do not match
    // separate diagonalizations by energy: overlap the *solver's* eigenvectors.
    std::vector<double> basis(dimension * dimension);
    std::vector<std::tuple<int,int,int,int>> quantum_numbers(dimension);
    for (int col = 0; col < dimension; ++col) {
      int file_column = -1, vector_n = -1, twice_l = -1, twice_s = -1, twice_j = -1;
      require(token == "VECTOR", "expected VECTOR in block " + std::to_string(ib + 1));
      input >> file_column >> vector_n >> twice_l >> twice_s >> twice_j;
      require(bool(input) && file_column == col && vector_n == n &&
              twice_l >= 0 && twice_s >= 0 && twice_j >= 0, "invalid vector label");
      quantum_numbers[col] = {vector_n, twice_l, twice_s, twice_j};
      for (int row = 0; row < dimension; ++row)
        input >> basis[row * dimension + col];
      require(bool(input), "truncated vector in block " + std::to_string(ib + 1));
      if (col + 1 < dimension) input >> token;
    }
    for (int a = 0; a < dimension; ++a) {
      for (int b = a; b < dimension; ++b) {
        double overlap = 0.0;
        for (int row = 0; row < dimension; ++row)
          overlap += basis[row * dimension + a] * basis[row * dimension + b];
        require(std::abs(overlap - (a == b ? 1.0 : 0.0)) < 1e-4,
                "NLSJ vectors are not orthonormal in block " + std::to_string(ib + 1));
      }
    }

    for (int eig = 0; eig < dimension; ++eig) {
      std::map<std::tuple<int,int,int,int>, double> weights;
      for (int col = 0; col < dimension; ++col) {
        double amplitude = 0.0;
        for (int row = 0; row < dimension; ++row)
          amplitude += h_eigenvectors[ib][eig * dimension + row] * basis[row * dimension + col];
        weights[quantum_numbers[col]] += amplitude * amplitude;
      }
      double total = 0.0, largest = -1.0;
      std::tuple<int,int,int,int> best;
      for (const auto& [qn, weight] : weights) {
        total += weight;
        if (weight > largest) { largest = weight; best = qn; }
      }
      require(std::abs(total - 1.0) < 1e-4, "NLSJ weights do not sum to one");
      auto& target = labels.block[ib][eig];
      target.n = std::get<0>(best);
      target.twice_l = std::get<1>(best);
      target.twice_s = std::get<2>(best);
      target.twice_j = std::get<3>(best);
      // Six-digit ASCII coefficients change the norm by a few parts in a
      // million. Normalize the weights so a pure state never prints >1.
      target.weight = largest / total;
      target.available = true;
    }
  }
  require(!(input >> token), "unexpected data after final block");
  labels.has_sidecar = true;
  return labels;
}

inline const Labels& cached(const Cix& cix)
{
  // The CIX and its sidecar do not change during a solver run. Reuse the
  // labels through all self-consistency iterations, after one validation.
  static std::optional<Labels> result;
  if (!result || result->filename != cix.cix_filename) result = load(cix);
  return *result;
}

inline void write(const Cix& cix, const Array<double,1>& probability, std::ostream& output)
{
  const Labels& labels = cached(cix);
  output << "# state block eigen P E_cix E_F0 E_gc N L S J label_weight status\n";
  output << "# eigen is zero-based, as in Probability.dat; state is one-based.\n";
  output << "# N is <N>; E_F0=E_cix+F0*<N(N-1)/2> for one site.\n";
  output << "# E_gc=E_F0-mu*<N> is the solver eigenvalue before lambda0.\n";
  output << "# L,S,J are dominant many-electron labels; weight measures their purity.\n";
  output << std::setprecision(12);
  int state = 0;
  for (int ib = 0; ib < cix.Nblock; ++ib) {
    for (int eig = 0; eig < cix.Block[ib].size(); ++eig) {
      ++state;
      const auto& label = labels.block[ib][eig];
      const double u = cix.Nsite == 1 ? cix.F0_input : 0.0;
      const double e_gc = label.energy;
      const double e_f0 = e_gc + cix.mu_input * label.mean_n;
      const double e_cix = e_f0 - u * label.mean_pairs;
      output << state << ' ' << ib + 1 << ' ' << eig << ' '
             << probability(cix.Block[ib](eig) - 1) << ' '
             << e_cix << ' ' << e_f0 << ' '
             << e_gc << ' ' << label.mean_n << ' ';
      if (label.available) {
        output << label.twice_l / 2.0 << ' ' << label.twice_s / 2.0 << ' '
               << label.twice_j / 2.0 << ' ' << label.weight << ' '
               << (label.weight < 0.7 ? "mixed" : "dominant");
      } else {
        output << "- - - 0 unavailable";
      }
      output << '\n';
    }
  }
}

} // namespace nlsj_probability
#endif
