#ifndef EMA_FUNC_H
#define EMA_FUNC_H

#include <complex>
#include <type_traits>
#include "ema-averaging.hxx"
#include "ema-bruggeman-terms.hxx"
#include "ema-data.hxx"


namespace ema::bruggeman::func {

template <typename T, typename AxialContributionFuncT> requires std::is_floating_point_v<T>
                                                             && std::is_invocable_r_v<std::complex<T>, AxialContributionFuncT, const std::complex<T>&, const std::complex<T>&, const T>
std::complex<T> isotropic3D(const std::vector<ema::data::MaterialNode<T>> &composite,
                            const std::complex<T> &effectiveParam,
                            AxialContributionFuncT axConFunc) {
    std::complex<T> F{0.0, 0.0};
    for (const auto &c : composite) {
        F += term::calcTerm(c.getVolumeFraction(),
                            term::axialContributionsSum(c, effectiveParam, axConFunc),
                            ema::averaging::isotropic3D<T>);
    }
    return F;
}

}   // namespace ema::bruggeman::funcEMA


#endif
