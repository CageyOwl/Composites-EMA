#ifndef EMA_BRUGGEMAN_TERMS_H
#define EMA_BRUGGEMAN_TERMS_H

#include <complex>
#include <type_traits>
#include "ema-data.hxx"


namespace ema::bruggeman::func::term {
namespace cmp {

    template <typename T> requires std::is_floating_point_v<T>
    inline std::complex<T> axialTerm(const std::complex<T> &param,
                                     const std::complex<T> &effectiveParam,
                                     const T L) {
        std::complex<T> numerator{param - effectiveParam};
        return (numerator / (effectiveParam + L * numerator));
    }

}   // namespace cmp

namespace derivative_cmp {

    template <typename T> requires std::is_floating_point_v<T>
    inline std::complex<T> axialTerm(const std::complex<T> &param,
                                     const std::complex<T> &effectiveParam,
                                     const T L) {
        std::complex<T> denominator{effectiveParam + L * (param - effectiveParam)};
        return -param / (denominator * denominator);
    }

}   //  derivative_cmp

template <typename T, typename AxialTermFuncT> requires std::is_floating_point_v<T>
                                                     && std::is_invocable_r_v<std::complex<T>, AxialTermFuncT, const std::complex<T>&, const std::complex<T>&, const T>
inline std::complex<T> componentAxialTermsSum(AxialTermFuncT func,
                                              const ema::data::MaterialNode<T> &node,
                                              const std::complex<T> &effectiveParam) {
    std::complex<T> axialTermsSum{0.0, 0.0};
    for (const auto &L : node.getDepolarizationFV())
        axialTermsSum += func(node.getParam(), effectiveParam, L);
    return axialTermsSum;
}

template <typename T, typename AxialTermFuncT> requires std::is_floating_point_v<T>
                                                     && std::is_invocable_r_v<std::complex<T>, AxialTermFuncT, const std::complex<T>&, const std::complex<T>&, const T>
std::complex<T> random3D(AxialTermFuncT func,
                         const ema::data::MaterialNode<T> &node,
                         const std::complex<T> &effectiveParam) {
    return componentAxialTermsSum(func, node, effectiveParam) * node.getVolumeFraction() / static_cast<T>(3.0);
}

// template <typename T> requires std::is_floating_point_v<T>
// std::complex<T> random2D(const ema::data::MaterialNode<T> &node,
//                          const std::complex<T> &effectiveParam) {
//     return (cmp::axialTerm(node.getParam(), effectiveParam, node.getDepolarizationFV()[ema::data::Axis::x])
//           + cmp::axialTerm(node.getParam(), effectiveParam, node.getDepolarizationFV()[ema::data::Axis::y]))
//           * node.getVolumeFraction() / static_cast<T>(2.0);
// }

// template <typename T> requires std::is_floating_point_v<T>
// std::complex<T> alignedZ(const ema::data::MaterialNode<T> &node,
//                          const std::complex<T> &effectiveParam) {
//     return cmp::axialTerm(node.getParam(), effectiveParam, node.getDepolarizationFV()[ema::data::Axis::z]) * node.getVolumeFraction();
// }

}   // namespace ema::bruggeman::func::term


#endif
