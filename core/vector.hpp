#ifndef ARTEMIS_VECTOR_HPP
#define ARTEMIS_VECTOR_HPP

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

namespace artemis {

struct Host {};
struct Device {};

struct HostToDevice {};
struct DeviceToHost {};

template<typename property_type>
class vector {

    using prop_host = thrust::host_vector<property_type>;
    using prop_device = thrust::device_vector<property_type>;

    prop_host prop_host_;
    prop_device prop_device_;

  public:
    explicit vector(const std::size_t size) : prop_host_(size), prop_device_(size) {}

    explicit vector(const thrust::host_vector<property_type>& prop_host) :
        prop_host_(prop_host),
        prop_device_(prop_host) {}

    explicit vector(const std::vector<property_type>& prop_host) :
        prop_host_(prop_host.begin(), prop_host.end()),
        prop_device_(prop_host.begin(), prop_host.end()) {}

    void mirror(HostToDevice) {
        thrust::copy(prop_host_.begin(), prop_host_.end(), prop_device_.begin());
    }

    void mirror(DeviceToHost) {
        thrust::copy(prop_device_.begin(), prop_device_.end(), prop_host_.begin());
    }

    auto begin(Device) const -> decltype(prop_device_.begin()) {
        return prop_device_.begin();
    }

    auto begin(Device) -> decltype(prop_device_.begin()) {
        return prop_device_.begin();
    }

    auto end(Device) const -> decltype(prop_device_.end()) {
        return prop_device_.end();
    }

    auto end(Device) -> decltype(prop_device_.end()) {
        return prop_device_.end();
    }

    auto begin(Host) const -> decltype(prop_host_.begin()) {
        return prop_host_.begin();
    }

    auto begin(Host) -> decltype(prop_host_.begin()) {
        return prop_host_.begin();
    }

    auto end(Host) const -> decltype(prop_host_.end()) {
        return prop_host_.end();
    }

    auto end(Host) -> decltype(prop_host_.end()) {
        return prop_host_.end();
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

    // void write(const std::string& file_name) const {}

    // void print() const {
    //     thrust::copy(vec.begin(), vec.end(), std::ostream_iterator<float>(std::cout, " "));
    // }
};

} // namespace artemis

#endif // ARTEMIS_VECTOR_HPP