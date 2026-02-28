#ifndef EMA_BRUGGEMAN_TERMS_H
#define EMA_BRUGGEMAN_TERMS_H

#include <complex>
#include <type_traits>
#include "ema-data.hxx"


namespace ema::bruggeman::func::term {
namespace axial {

    template <typename T> requires std::is_floating_point_v<T>
    inline std::complex<T> contribution(const std::complex<T> &param,
                                        const std::complex<T> &effectiveParam,
                                        const T L) {
        std::complex<T> numerator{param - effectiveParam};
        return (numerator / (effectiveParam + L * numerator));
    }

    template <typename T> requires std::is_floating_point_v<T>
    inline std::complex<T> derivativeContribution(const std::complex<T> &param,
                                                  const std::complex<T> &effectiveParam,
                                                  const T L) {
        std::complex<T> denominator{effectiveParam + L * (param - effectiveParam)};
        return -param / (denominator * denominator);
    }

}   //  axial

template <typename T, typename AxialContributionFuncT> requires std::is_floating_point_v<T>
                                                    && std::is_invocable_r_v<std::complex<T>, AxialContributionFuncT, const std::complex<T>&, const std::complex<T>&, const T>
std::complex<T> axialContributionsSum(const ema::data::MaterialNode<T> &node,
                                 const std::complex<T> &effectiveParam,
                                 AxialContributionFuncT func) {
    std::complex<T> sum{0.0, 0.0};
    const std::complex<T> &param{node.getParam()};
    for (const auto &L : node.getDepolarizationFV())
        sum += func(param, effectiveParam, L);
    return sum;
}

template <typename T, typename AveragingFuncT> requires std::is_floating_point_v<T>
                                                     && std::is_invocable_r_v<std::complex<T>, AveragingFuncT, const std::complex<T>&>
std::complex<T> calcTerm(const T fraction,
                         const std::complex<T> &axialContributionsSumVal,
                         AveragingFuncT averFunc) {
    return fraction * averFunc(axialContributionsSumVal);
}

}   // namespace ema::bruggeman::func::term


#endif
