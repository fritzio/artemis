#ifndef ARTEMIS_UTILS_XECUTION_SPACE_ANNOTATIONS_HPP
#define ARTEMIS_UTILS_XECUTION_SPACE_ANNOTATIONS_HPP

// host and device execution-space annotations
#ifdef __CUDA_ARCH__
#define HOST_AND_DEVICE_FUNCTION __device__ __host__
#else
#define HOST_AND_DEVICE_FUNCTION /* empty */
#endif

// host execution-space annotations
#ifdef __CUDA_ARCH__
#define HOST_FUNCTION __host__
#else
#define HOST_FUNCTION /* empty */
#endif

// device execution-space annotations
#ifdef __CUDA_ARCH__
#define DEVICE_FUNCTION __device__
#else
#define DEVICE_FUNCTION /* empty */
#endif

#endif // ARTEMIS_UTILS_XECUTION_SPACE_ANNOTATIONS_HPP