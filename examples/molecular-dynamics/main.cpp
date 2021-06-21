// #include "vector.hpp"
#include "version.hpp"

#include <numeric>
#include <iostream>
#include <cmath>
#include <string>

//

#include <array>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include <thrust/copy.h>
#include <thrust/transform.h>
#include <thrust/sort.h>

#define MAX_NEIGHBOR_CAPACITY 100

struct Host {};
struct Device {};

struct HostToDevice {};
struct DeviceToHost {};

template<typename property_type>
class state_vector {

  using prop_host = thrust::host_vector<property_type>;
  using prop_device = thrust::device_vector<property_type>;

  prop_host prop_host_;
  prop_device prop_device_;

public:
  explicit state_vector(const std::size_t size) : prop_host_(size), prop_device_(size) {}

  explicit state_vector(const thrust::host_vector<property_type>& prop_host) :
      prop_host_(prop_host),
      prop_device_(prop_host) {}

  explicit state_vector(const std::vector<property_type>& prop_host) :
      prop_host_(prop_host.begin(), prop_host.end()),
      prop_device_(prop_host.begin(), prop_host.end()) {}

  void mirror(HostToDevice) {
    thrust::copy(prop_host_.begin(), prop_host_.end(), prop_device_.begin());
  }

  void mirror(DeviceToHost) {
    thrust::copy(prop_device_.begin(), prop_device_.end(), prop_host_.begin());
  }

  std::size_t size() const {
    return prop_host_.size();
  }

  prop_host& get(Host) {
    return prop_host_;
  }

  const prop_host& get(Host) const {
    return prop_host_;
  }

  prop_device& get(Device) {
    return prop_device_;
  }

  const prop_device& get(Device) const {
    return prop_device_;
  }

  void write(const std::string& file_name) const {

  }

  // void print() const {
  //   thrust::copy(vec.begin(), vec.end(), std::ostream_iterator<float>(std::cout, " "));
  // }
};

//

struct InterParticleForce {

  const double epsilon_;

  const double sigma_;
  const double sigma_6_;
  const double sigma_12_;

  InterParticleForce(const double epsilon, const double sigma) :
      epsilon_(epsilon),
      sigma_(sigma),
      sigma_6_(std::pow(sigma_, 6.0)),
      sigma_12_(std::pow(sigma_, 12.0)) {}

  double squared_distance(const std::array<double, 3>& a, const std::array<double, 3>& b) const {

    double distance_squared = 0.0;

    for (std::size_t d = 0; d < 3; ++d) {
      distance_squared += (a[d] - b[d]) * (a[d] - b[d]);
    }
    return distance_squared;
  }

  std::array<double, 3> operator()(const std::array<double, 3>& x_i,
                                   const std::array<std::array<double, 3>, MAX_NEIGHBOR_CAPACITY>& x_neighbor) const {

    std::array<double, 3> f_ij = {0.0, 0.0, 0.0};

    for (const auto x_j : x_neighbor) {

      // compute the squared distance between the ith and jth particle, where the latter is in the local neighborhood to
      // the former
      const double x_ij_squared = squared_distance(x_i, x_j);

      // compute the potential between the ith and jth particle
      const double phi_ij = 48.0 * epsilon_ *
                            (sigma_12_ / std::pow(x_ij_squared, 6.0) - sigma_6_ / std::pow(x_ij_squared, 3.0)) /
                            x_ij_squared;

      // compute the force from the potential
      for (std::size_t d = 0; d < 3; ++d) {
        f_ij.at(d) += x_i.at(d) * phi_ij;
      }
    }

    return f_ij;
  }
};

//

// neighborhood_matrix

struct MaskNeighborhood {

  const double r_cutoff_;

  explicit MaskNeighborhood(const double r_cutoff) :
    r_cutoff_(r_cutoff)
  {}

  double distance(const std::array<double, 3>& x, const std::array<double, 3>& y) const {
    double distance_squared = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
      distance_squared += (x[d] - y[d]) * (x[d] - y[d]);
    }
    return std::sqrt(distance_squared);
  }

  bool is_within_cutoff(const std::size_t distance) const {
    return (distance < r_cutoff_) ? true : false;
  }

  std::size_t operator() (const std::size_t index_i, const std::array<double, 3>& x_i,
                          const std::size_t index_j, const std::array<double, 3>& x_j) const {
    return is_within_cutoff(distance(x_i, x_j)) ? index_j : index_i;
  }
};

struct CellList {

  const double grid_delta_;

  const std::array<double, 3> grid_min_;
  const std::array<double, 3> grid_max_;

  const std::uint32_t n_cells;

  explicit CellList(const double grid_delta, const std::array<double, 3>& grid_min,
                    const std::array<double, 3>& grid_max) :
      grid_delta_(grid_delta),
      grid_min_(grid_min),
      grid_max_(grid_max),
      n_cells(1.0 / grid_delta) {}
};

class particle_to_cell_index : public thrust::unary_function<std::array<double, 3>, std::uint32_t> {

  const CellList& cell_list_;

public:
  explicit particle_to_cell_index(const CellList& cell_list) : cell_list_(cell_list) {}

  std::uint32_t operator() (const std::array<double, 3>& r) const {

    // index of the grid cell containing the particle with the coordinate `r`
    const std::uint32_t id_x = r.at(0) * cell_list_.n_cells;
    const std::uint32_t id_y = r.at(1) * cell_list_.n_cells;
    const std::uint32_t id_z = r.at(2) * cell_list_.n_cells;

    // https://stackoverflow.com/questions/10903149/how-do-i-compute-the-linear-index-of-a-3d-coordinate-and-vice-versa
    return id_x + (cell_list_.n_cells + 1) * id_y + (cell_list_.n_cells + 1) * (cell_list_.n_cells + 1) * id_z;
  }
};

//

using fprec = double;

using position_t = std::array<fprec, 3>;
using velocity_t = std::array<fprec, 3>;
using force_t = std::array<fprec, 3>;

int main(void) {

  const fprec epsilon = 1.0;
  const fprec sigma = 0.1;
  const fprec cutoff_radius = 0.1; //3.0 * sigma;
  const fprec time_step = 0.0005;
  const fprec particle_delta = 0.1;
  const fprec grid_delta = cutoff_radius;

  const std::array<fprec, 3> grid_min = {0.0, 0.0, 0.0};
  const std::array<fprec, 3> grid_max = {1.0, 1.0, 1.0};

  const auto cell_list = CellList(grid_delta, grid_min, grid_max);

  const auto _ = particle_to_cell_index(cell_list);

  // std::cout << cell_list.n_cells << std::endl;
  // std::cout << _({0.0,0.15,0.0}) << std::endl;

  // number of particles into each spatial dimension
  const std::array<std::size_t, 3> n_particles = {
      static_cast<std::size_t>((grid_max.at(0) - grid_min.at(0)) / particle_delta),
      static_cast<std::size_t>((grid_max.at(1) - grid_min.at(1)) / particle_delta),
      static_cast<std::size_t>((grid_max.at(2) - grid_min.at(2)) / particle_delta)};

  // number of particles 
  const auto n_total_particles =
      std::accumulate(n_particles.begin(), n_particles.end(), 1, std::multiplies<const std::size_t>());

  // std::cout << n_total_particles << std::endl;

  // initialize
  // std::vector<std::array<double, 3>> x_init;
  // std::vector<std::array<double, 3>> u_init;
  // std::vector<std::array<double, 3>> f_init;

  auto position = state_vector<position_t>(n_total_particles);
  auto velocity = state_vector<velocity_t>(n_total_particles);
  auto force = state_vector<force_t>(n_total_particles);

  auto& host_particles = position.get(Host());

  std::size_t l = 0;

  for (std::size_t k = 0; k < n_particles.at(0); ++k) {
    for (std::size_t j = 0; j < n_particles.at(1); ++j) {
      for (std::size_t i = 0; i < n_particles.at(2); ++i) {
        host_particles[l] = {(0.5 + static_cast<fprec>(i)) * particle_delta, (0.5 + static_cast<fprec>(j)) * particle_delta, (0.5 + static_cast<fprec>(k)) * particle_delta};
        l++;
      }
    }
  }
  position.mirror(HostToDevice());

  // the velocity is initialized to zero on the device memory and then mirrored to the host
  auto& velocity_d = velocity.get(Device());
  thrust::fill(velocity_d.begin(), velocity_d.end(), std::array<fprec, 3>({0.0, 0.0, 0.0}));

  velocity.mirror(DeviceToHost());

  // initialize the neighborhood matrix with dimension (number of particles, number of neighbors)
  // auto neighborhood_matrix = state_vector<std::uint32_t>(n_total_particles * 100);

  auto grid_index = state_vector<std::uint32_t>(cell_list.n_cells * cell_list.n_cells * cell_list.n_cells);

  thrust::transform(position.get(Device()).begin(), position.get(Device()).end(), grid_index.get(Device()).begin(), particle_to_cell_index(cell_list));

  // sort
  thrust::sort_by_key(grid_index.get(Device()).begin(), grid_index.get(Device()).end(), position.get(Device()).begin());

  // compute the size of each cell
  auto grid_size = state_vector<std::uint32_t>(cell_list.n_cells * cell_list.n_cells * cell_list.n_cells);

  //
  //
  //
  //
  //

  struct transform_neighbors : public thrust::unary_function<
    // the type of the neighboring keys
    std::array<std::uint32_t, MAX_NEIGHBOR_CAPACITY>,
    // the type of the neighboring positions
    std::array<std::array<double, 3>, MAX_NEIGHBOR_CAPACITY>> {

    std::array<std::array<double, 3>, MAX_NEIGHBOR_CAPACITY> operator() (const std::array<std::uint32_t, MAX_NEIGHBOR_CAPACITY>& keys) {

      std::array<std::array<double, 3>, MAX_NEIGHBOR_CAPACITY> x_neighbor;

      return x_neighbor;
    }

  };

  // thrust::make_zip_iterator(thrust::make_tuple(position.begin(),));

  // compute force
  // thrust::transform()

  // for (auto& x : particles.)
}