#ifndef EMA_TERMS_AVERAGING_H
#define EMA_TERMS_AVERAGING_H

#include <complex>
#include <type_traits>


namespace ema::averaging {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> isotropic3D(const std::complex<T> &axialSum) {
    return axialSum / static_cast<T>(3.0);
}

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> isotropic2D(const std::complex<T> &axialSum) {
    return axialSum / static_cast<T>(2.0);
}

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> aligned(const std::complex<T> &axialSum) {
    return axialSum;
}

}   // namespace ema::averaging

#endif
