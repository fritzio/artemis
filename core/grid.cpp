#include "grid.hpp"

#include <array>

int artemis::grid::periodic_mapping(const int x, const int nx) {
	return std::abs(x + nx) % nx;
}

int artemis::grid::map_coordinates_to_identifier(const int x, const int y, const int z, const int nx, const int ny) {
	return x + y * nx + z * nx * ny;
}

artemis::grid::int3 artemis::grid::map_identifier_to_coordinates(const int cell_id, const int nx, const int ny, const int nz) {
    
}

int artemis::grid::convert_to_cell_id(const artemis::grid::int3 coordinate, const artemis::grid::int3 cell_count) {

	const int x = periodic_mapping(std::get<0>(coordinate), std::get<0>(cell_count));
	const int y = periodic_mapping(std::get<1>(coordinate), std::get<1>(cell_count));
	const int z = periodic_mapping(std::get<2>(coordinate), std::get<2>(cell_count));

	const int nx = std::get<0>(cell_count);
	const int ny = std::get<1>(cell_count);

	return map_coordinates_to_identifier(x, y, z, nx, ny);
}

artemis::grid::int3 artemis::grid::convert_to_coordinates(const int cell_id, const artemis::grid::int3 cell_count) {

    const int x = cell_id % std::get<0>(cell_count);

    return artemis::grid::int3({x, x, x});
}

std::array<artemis::grid::int3, 28> artemis::grid::get_neighboring_coordinates(const artemis::grid::int3 coordinate, 
                                                                               const artemis::grid::int3 cell_count) {

    const int x = std::get<0>(coordinate);
    const int y = std::get<1>(coordinate);
    const int z = std::get<2>(coordinate);

    std::array<artemis::grid::int3, 28> neighboring_coordinates = {
        std::tuple<int, int, int>{x, y, z},     // 
        std::tuple<int, int, int>{x - 1, y, z}, // left
        std::tuple<int, int, int>{x + 1, y, z}, // right
        std::tuple<int, int, int>{x, y - 1, z}, // bottom
        std::tuple<int, int, int>{x, y + 1, z}, // top
        std::tuple<int, int, int>{x, y, z - 1}, // back
        std::tuple<int, int, int>{x, y, z + 1}, // front
    };

    return neighboring_coordinates;
}

std::array<int, 28> artemis::grid::get_neighboring_ids(const int cell_id, const artemis::grid::int3 cell_count) {

    const artemis::grid::int3 x = artemis::grid::convert_to_coordinates(cell_id, cell_count);

    const std::array<artemis::grid::int3, 28> neighboring_coordinates = get_neighboring_coordinates(x, cell_count);
    std::array<int, 28> neighboring_ids;

    for (std::size_t i = 0; i < 28; ++i) {
        neighboring_ids.at(i) = convert_to_cell_id(neighboring_coordinates.at(i), cell_count);
    }

    return neighboring_ids;
}