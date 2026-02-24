#ifndef EMA_3D_TERMS_HXX
#define EMA_3D_TERMS_HXX

#include <complex>
#include <type_traits>
#include "ema-data.hxx"


//  Naming: Par - particles, Dist - random
namespace ema::bruggeman::v3D {
namespace term {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParRandDist(const ema::data::MaterialNode<T> &node,
                                           const std::complex<T> &effectiveParam) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    numerator{ node.getParam() - effectiveParam };
    for (const auto L : node.getDepolarizationFV())
        termsSum += ( numerator / (effectiveParam + L * numerator) );
    return (termsSum * node.getVolumeFraction() / 3.0);
}

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParRandDist(const ema::data::MaterialNode<T> &node,
                                           const std::complex<T> &effectiveParam,
                                           const T volumeFraction) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    numerator{ node.getParam() - effectiveParam };
    for (const auto L : node.getDepolarizationFV())
        termsSum += ( numerator / (effectiveParam + L * numerator) );
    return (termsSum * volumeFraction / 3.0);
}


template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParRandDist(const ema::data::MaterialNode<T> &node,
                                           const std::complex<T> &effectiveParam,
                                           const std::complex<T> &param) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    numerator{ param - effectiveParam };
    for (const auto L : node.getDepolarizationFV())
        termsSum += ( numerator / (effectiveParam + L * numerator) );
    return (termsSum * node.getVolumeFraction() / 3.0);
}

}   // namespace term

namespace derivative_term {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParRandDist(const ema::data::MaterialNode<T> &node,
                                           const std::complex<T> &effectiveParam) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    param{ node.getParam() },
                    denominator_root{ 0.0, 0.0 };
    for (const auto L : node.getDepolarizationFV())
    {
        denominator_root = effectiveParam + L * (param - effectiveParam);
        termsSum += -param / (denominator_root * denominator_root);
    }
    return (termsSum * node.getVolumeFraction() / 3.0);
}

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParRandDist(const ema::data::MaterialNode<T> &node,
                                           const std::complex<T> &effectiveParam,
                                           const T volumeFraction) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    param{ node.getParam() },
                    denominator_root{ 0.0, 0.0 };
    for (const auto L : node.getDepolarizationFV())
    {
        denominator_root = effectiveParam + L * (param - effectiveParam);
        termsSum += -param / (denominator_root * denominator_root);
    }
    return (termsSum * volumeFraction / 3.0);
}

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParRandDist(const ema::data::MaterialNode<T> &node,
                                           const std::complex<T> &effectiveParam,
                                           const std::complex<T> &param) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    denominator_root{ 0.0, 0.0 };
    for (const auto L : node.getDepolarizationFV())
    {
        denominator_root = effectiveParam + L * (param - effectiveParam);
        termsSum += -param / (denominator_root * denominator_root);
    }
    return (termsSum * node.getVolumeFraction() / 3.0);
}

}   // namespace derivative_term
}   // namespace ema::bruggeman::v3D


#endif  // EMA_3D_TERMS_HXX
