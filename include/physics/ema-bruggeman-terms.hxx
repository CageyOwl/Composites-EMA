#ifndef EMA_BRUGGEMAN_TERMS_H
#define EMA_BRUGGEMAN_TERMS_H

#include <complex>
#include <concepts>
#include <type_traits>
#include "ema-data.hxx"


namespace ema::bruggeman::func::term {

template <typename Func, typename T>
concept AxialContributionCalcFunc = std::floating_point<T> && requires(Func func, const std::complex<T> &param, const std::complex<T> &effectiveParam, const T L) {
    { func(param, effectiveParam, L) } -> std::convertible_to<std::complex<T>>;
};

template <typename Func, typename T>
concept AveragingFunc = std::floating_point<T> && requires(Func func, const std::complex<T> &axialSum) {
    { func(axialSum) } -> std::convertible_to<std::complex<T>>;
};

namespace axial {

    template <std::floating_point T>
    inline std::complex<T> contribution(const std::complex<T> &param,
                                        const std::complex<T> &effectiveParam,
                                        const T L) {
        std::complex<T> numerator{param - effectiveParam};
        return (numerator / (effectiveParam + L * numerator));
    }

    template <std::floating_point T>
    inline std::complex<T> derivativeContribution(const std::complex<T> &param,
                                                  const std::complex<T> &effectiveParam,
                                                  const T L) {
        std::complex<T> denominator{effectiveParam + L * (param - effectiveParam)};
        return -param / (denominator * denominator);
    }

}   //  namespace axial

template <std::floating_point T>
std::complex<T> axialContributionsSum(const ema::data::MaterialNode<T> &node,
                                      const std::complex<T> &effectiveParam,
                                      AxialContributionCalcFunc<T> auto axConFunc) {
    std::complex<T> sum{0.0, 0.0};
    const std::complex<T> &param{node.getParam()};
    for (const auto &L : node.getDepolarizationFV())
        sum += axConFunc(param, effectiveParam, L);
    return sum;
}

template <std::floating_point T>
std::complex<T> calcTerm(const T fraction,
                         const std::complex<T> &axialContributionsSumVal,
                         AveragingFunc<T> auto avgFunc) {
    return fraction * avgFunc(axialContributionsSumVal);
}

}   // namespace ema::bruggeman::func::term


#endif
