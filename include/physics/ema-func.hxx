#ifndef EMA_FUNC_HXX
#define EMA_FUNC_HXX

#include <complex>
#include <type_traits>
#include "ema-data.hxx"
#include "ema-3d-terms.hxx"


namespace ema::bruggeman::funcEMA {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParRandDist(const std::vector<ema::data::MaterialNode<T>> &composite,
                                       const std::complex<T> &effectiveParam) {
    std::complex<T> F{ 0.0, 0.0 };
    for (const auto &c : composite)
        F += ema::bruggeman::v3D::term::anisotropicParRandDist(c, effectiveParam);
    return F;
}

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParRandDist(const std::vector<ema::data::MaterialNode<T>> &composite,
                                       const std::complex<T> &effectiveParam,
                                       const data::VolumeFractions<T> &volumeFractions) {
    std::complex<T> F{ 0.0, 0.0 };
    size_t i{ 0 };
    for (const auto &c : composite)
        F += ema::bruggeman::v3D::term::anisotropicParRandDist(c, effectiveParam, volumeFractions[i++]);
    return F;
}

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParRandDist(const std::vector<ema::data::MaterialNode<T>> &composite,
                                       const std::complex<T> &effectiveParam,
                                       const data::FillersParams<T> &fillersParams) {
    std::complex<T> F{ 0.0, 0.0 };
    // TODO
    for (const auto &c : composite)
        F += ;
    return F;
}

}   // namespace ema::bruggeman::funcEMA


#endif
