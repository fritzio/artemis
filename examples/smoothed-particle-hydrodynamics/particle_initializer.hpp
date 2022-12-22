#ifndef PARTICLE_INiTIALIZER_HPP
#define PARTICLE_INiTIALIZER_HPP

#include <vector>

namespace particle_initializer {

template<typename fprec, typename CellList>
void cartesian_lattice(std::vector<fprec>& x, std::vector<fprec>& y, std::vector<fprec>& z, std::vector<int>& id,
                       const fprec dx, const CellList& cell_list) {

    int counter = 0;

    for (int k = 0; k < (static_cast<int>((cell_list.grid_max.at(2) - cell_list.grid_min.at(2)) / dx)); ++k) {

        const fprec zp = static_cast<fprec>(0.5 + k) * dx;

        for (int j = 0; j < (static_cast<int>((cell_list.grid_max.at(1) - cell_list.grid_min.at(1)) / dx)); ++j) {

            const fprec yp = static_cast<fprec>(0.5 + j) * dx;

            for (int i = 0; i < (static_cast<int>((cell_list.grid_max.at(0) - cell_list.grid_min.at(0)) / dx)); ++i) {

                const fprec xp = static_cast<fprec>(0.5 + i) * dx;

                x.push_back(xp);
                y.push_back(yp);
                z.push_back(zp);

                id.push_back(counter++);
            }
        }
    }
}

} // namespace particle_initializer

#endif // PARTICLE_INiTIALIZER_HPP