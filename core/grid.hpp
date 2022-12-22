#ifndef ARTEMIS_GRID_HPP
#define ARTEMIS_GRID_HPP

#include <cmath>
#include <tuple>

namespace artemis {

namespace grid {

using int3 = std::tuple<int, int, int>;

int periodic_mapping(const int x, const int nx);

int map_coordinates_to_identifier(const int x, const int y, const int z, const int nx, const int ny);

int3 map_identifier_to_coordinates(const int cell_id, const int nx, const int ny, const int nz);

int convert_to_cell_id(const int3 coordinate, const int3 cell_count);

int3 convert_to_coordinates(const int cell_id, const int3 cell_count);

std::array<int3, 27> get_neighboring_coordinates(const int3 coordinate, const int3 cell_count);

std::array<int, 27> get_neighboring_ids(const int cell_id, const int3 cell_count);

} // namespace grid

} // namespace artemis

#endif // ARTEMIS_GRID_HPP