#ifndef EMA_FUNC_H
#define EMA_FUNC_H

#include <complex>
#include <concepts>
#include <type_traits>
#include "ema-averaging.hxx"
#include "ema-bruggeman-terms.hxx"
#include "ema-data.hxx"


namespace ema::bruggeman::func {

template <std::floating_point T>
std::complex<T> isotropic3D(const std::vector<ema::data::MaterialNode<T>> &composite,
                            const std::complex<T> &effectiveParam,
                            ema::bruggeman::func::term::AxialContributionCalcFunc<T> auto axConFunc) {
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
