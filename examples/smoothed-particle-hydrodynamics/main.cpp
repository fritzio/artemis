#include "vector.hpp"

#include "kernel.hpp"
#include "particle_initializer.hpp"

#include <array>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>
#include <thrust/tuple.h>

#include <utils/execution_space_annotations.hpp>

using fprec = float;

using Host = artemis::Host;
using Device = artemis::Device;

using HostToDevice = artemis::HostToDevice;
using DeviceToHost = artemis::DeviceToHost;

// template<typename fprec, std::size_t dim>
// struct count_cells_in_each_dimension {

//     auto operator()(const fprec grid_delta, const std::array<fprec, dim>& grid_min,
//                     const std::array<fprec, dim>& grid_max) const -> const std::array<std::size_t, dim> {

//         std::array<std::size_t, dim> cell_count;

//         for (std::size_t d = 0; d < dim; ++d) {
//             cell_count.at(d) = static_cast<std::size_t>((grid_max.at(d) - grid_min.at(d)) / grid_delta) + 1;
//         }
//     }
// };

struct print {
    template<typename T>
    HOST_FUNCTION void operator()(T t) const {
        std::cout << t << std::endl;
    }
};

struct write {

    std::ofstream& file;

    explicit write(std::ofstream& file) : file(file) {}

    template<typename T>
    HOST_FUNCTION void operator()(T t) {
        file << thrust::get<0>(t) << "," << thrust::get<1>(t) << "," << thrust::get<2>(t) << "," << thrust::get<3>(t)
             << "," << thrust::get<4>(t) << std::endl;
    }
};

// struct write {

//     const std::string filename = "test.csv";
//     std::ofstream file;

//     explicit write() {
//         file.open(filename);
//     }

//     template<typename T>
//     HOST_FUNCTION void operator()(T t) {
//         file << t << std::endl;
//     }
// };

namespace {

template<typename fprec, std::size_t dim>
struct count_cells {

    const fprec dx;

    explicit count_cells(const fprec dx) : dx(dx){};

    // std::size_t operator()(const fprec x_min, const fprec x_max) const {
    //     return static_cast<std::size_t>((x_max - x_min) / dx) + 1;
    // }

    template<typename T>
    std::size_t operator()(const T t) const {

        const fprec x_min = thrust::get<0>(t);
        const fprec x_max = thrust::get<1>(t);

        return static_cast<std::size_t>((x_max - x_min) / dx);
    }
};

namespace initialize {

template<typename fprec, std::size_t dim>
std::array<std::size_t, dim> cell_count(const fprec grid_delta, const std::array<fprec, dim>& grid_min,
                                        const std::array<fprec, dim>& grid_max) {

    std::array<std::size_t, dim> cell_count;

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(grid_min.begin(), grid_max.begin())),
                      thrust::make_zip_iterator(thrust::make_tuple(grid_min.end(), grid_max.end())), cell_count.begin(),
                      count_cells<fprec, dim>(grid_delta));

    return cell_count;
}

template<std::size_t dim>
std::size_t total_number_of_cells(std::array<std::size_t, dim> cell_count) {
    return thrust::reduce(cell_count.begin(), cell_count.end(), 1, thrust::multiplies<std::size_t>());
}

} // namespace initialize

} // namespace

template<typename fprec, std::size_t dim>
struct CellList {
    const fprec grid_delta;

    const std::array<fprec, dim> grid_min;
    const std::array<fprec, dim> grid_max;

    const std::array<std::size_t, dim> cell_count;

    const std::size_t total_number_of_cells;

    explicit CellList(const fprec grid_delta, const std::array<fprec, dim>& grid_min,
                      const std::array<fprec, dim>& grid_max) :
        grid_delta(grid_delta),
        grid_min(grid_min),
        grid_max(grid_max),
        cell_count(initialize::cell_count(grid_delta, grid_min, grid_max)),
        total_number_of_cells(initialize::total_number_of_cells<dim>(cell_count)) {}
};

template<typename fprec, std::size_t dim>
class get_cell_id {

    const CellList<fprec, dim>& cell_list_;

  public:
    explicit get_cell_id(const CellList<fprec, dim>& cell_list) : cell_list_(cell_list) {}

    std::size_t operator()(const thrust::tuple<fprec, fprec, fprec> particle_position) const {

        const fprec x_particle = thrust::get<0>(particle_position);
        const fprec y_particle = thrust::get<1>(particle_position);
        const fprec z_particle = thrust::get<2>(particle_position);

        const std::size_t x_cell =
            static_cast<std::size_t>((x_particle - cell_list_.grid_min.at(0)) / cell_list_.grid_delta);
        const std::size_t y_cell =
            static_cast<std::size_t>((y_particle - cell_list_.grid_min.at(1)) / cell_list_.grid_delta);
        const std::size_t z_cell =
            static_cast<std::size_t>((z_particle - cell_list_.grid_min.at(2)) / cell_list_.grid_delta);

        return x_cell + y_cell * cell_list_.cell_count.at(0) +
               z_cell * cell_list_.cell_count.at(0) * cell_list_.cell_count.at(1);
    }
};

int main(void) {

    // particle spacing on cartesian lattice
    const fprec dx = static_cast<fprec>(0.1);

    // smoothing length of kernel function
    const fprec h = static_cast<fprec>(1.0) * dx;

    // initialize the kernel function
    const auto kernel = Kernel<fprec, 3>(h);

    // initialize grid
    const auto grid = CellList<fprec, 3>(kernel.cutoff_radius, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});

    // initial particle coordinates and particle identifier
    std::vector<fprec> x_init, y_init, z_init;
    std::vector<int> id_init;

    // initalize particle coordinates and identifier onto a cartesian lattice
    particle_initializer::cartesian_lattice(x_init, y_init, z_init, id_init, dx, grid);

    // initialize trajectory
    artemis::vector<fprec> x(x_init), y(y_init), z(z_init);
    artemis::vector<int> id(id_init);

    // cell_ids just needed for sorting

    artemis::vector<int> cell_ids(id.size());

    thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(x.begin(Device()), y.begin(Device()), z.begin(Device()))),
        thrust::make_zip_iterator(thrust::make_tuple(x.end(Device()), y.end(Device()), z.end(Device()))),
        cell_ids.begin(Device()), get_cell_id(grid));

    cell_ids.mirror(DeviceToHost());

    // count the number of particles in each cell

    // artemis::vector<int> neighbor_ids();

    //
    //
    //

    // thrust::for_each(cell_ids.begin(Host()), cell_ids.end(Host()), print());

    // std::ofstream file;

    // file.open("trajectory_0.csv");
    // thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(id.begin(Host()), x.begin(Host()), y.begin(Host()),
    //                                                               z.begin(Host()), cell_ids.begin(Host()))),
    //                  thrust::make_zip_iterator(thrust::make_tuple(id.end(Host()), x.end(Host()), y.end(Host()),
    //                                                               z.end(Host()), cell_ids.end(Host()))),
    //                  write(file));

    // file.close();
}