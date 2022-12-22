#ifndef ARTEMIS_VERSION_HPP
#define ARTEMIS_VERSION_HPP

#include <string>
#include <thrust/version.h>

namespace artemis {

namespace version {

std::string thrust() {

    std::string major = std::to_string(THRUST_MAJOR_VERSION);
    std::string minor = std::to_string(THRUST_MINOR_VERSION);

    return major + "." + minor;
}

} // namespace version

} // namespace artemis

#endif // ARTEMIS_VERSION_HPP