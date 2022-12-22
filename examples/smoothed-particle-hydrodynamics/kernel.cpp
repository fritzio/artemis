#include <algorithm>

#include "kernel.hpp"

template<typename fprec, std::size_t dim>
fprec Kernel<fprec, dim>::w(const fprec r) const {
    const fprec q = r * h_inv;

    const fprec q1 = std::max(static_cast<fprec>(0.0), static_cast<fprec>(1.0) - q);
    const fprec q2 = std::max(static_cast<fprec>(0.0), static_cast<fprec>(2.0) - q);
    const fprec q3 = std::max(static_cast<fprec>(0.0), static_cast<fprec>(3.0) - q);

    return sigma * (q3 * q3 * q3 * q3 * q3 - c2 * q2 * q2 * q2 * q2 * q2 + c1 * q1 * q1 * q1 * q1 * q1);
}

template struct Kernel<double, 1>;
template struct Kernel<double, 2>;
template struct Kernel<double, 3>;

template struct Kernel<float, 1>;
template struct Kernel<float, 2>;
template struct Kernel<float, 3>;