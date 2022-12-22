#ifndef KERNEL_HPP
#define KERNEL_HPP

#include <math.h>

template<typename fprec, std::size_t dim>
struct Kernel {

    const fprec h;
    const fprec h_inv;
    const fprec cutoff_radius;

    static constexpr int ratio_cutoff_over_h = 3;

    const fprec c1 = static_cast<fprec>(15.0);
    const fprec c2 = static_cast<fprec>(6.0);

    const fprec sigma = (dim == 3) ? (static_cast<fprec>(1.0) / static_cast<fprec>(120.0) / M_PI)
                                   : (static_cast<fprec>(7.0) / static_cast<fprec>(478.0) / M_PI);

    Kernel(const fprec smoothing_length) :
        h(smoothing_length),
        h_inv(static_cast<fprec>(1.0) / smoothing_length),
        cutoff_radius(smoothing_length * static_cast<fprec>(ratio_cutoff_over_h)) {}

    fprec w(const fprec r) const;
};

#endif // KERNEL_HPP