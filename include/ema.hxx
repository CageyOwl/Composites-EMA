#ifndef EMA_HXX
#define EMA_HXX

#include <complex>
#include <type_traits>
#include "data.hxx"


//  Naming: Par - particles, Dist - distribution
namespace bruggeman {
namespace base {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDistTerm(const data::MaterialNode<T> &node, const std::complex<T> &effectiveParam) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    numerator{ node.getParam() - effectiveParam };
    for (const auto L : node.getDepolarizationF())
        termsSum += ( numerator / (effectiveParam + L * numerator) );
    return (termsSum * node.getVolumeFraction() / 3.0);
}

}   // namespace base

namespace base_variable {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDistTerm(const data::MaterialNode<T> &node, const std::complex<T> &effectiveParam, const T volumeFraction) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    numerator{ node.getParam() - effectiveParam };
    for (const auto L : node.getDepolarizationF())
        termsSum += ( numerator / (effectiveParam + L * numerator) );
    return (termsSum * volumeFraction / 3.0);
}

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDistTerm(const data::MaterialNode<T> &node, const std::complex<T> &effectiveParam, const std::complex<T> &param) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    numerator{ param - effectiveParam };
    for (const auto L : node.getDepolarizationF())
        termsSum += ( numerator / (effectiveParam + L * numerator) );
    return (termsSum * node.getVolumeFraction() / 3.0);
}

}   // namespace base_controlled

namespace derivative {

//  TODO

}   // namespace derivative
}   // namespace bruggema


#endif  // EMA_HXX
