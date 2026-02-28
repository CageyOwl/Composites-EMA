#ifndef EMA_FUNC_H
#define EMA_FUNC_H

#include <complex>
#include <type_traits>
#include "ema-data.hxx"
#include "ema-bruggeman-terms.hxx"


namespace ema::bruggeman::func {

// template <typename T> requires std::is_floating_point_v<T>
// std::complex<T> anisotropicParRandDist(const std::vector<ema::data::MaterialNode<T>> &composite,
//                                        const std::complex<T> &effectiveParam) {
//     std::complex<T> F{ 0.0, 0.0 };
//     for (const auto &c : composite)
//         F += ema::bruggeman::term::v3D::anisotropicParRandDist(c, effectiveParam);
//     return F;
// }

// template <typename T> requires std::is_floating_point_v<T>
// std::complex<T> anisotropicParRandDist(const std::vector<ema::data::MaterialNode<T>> &composite,
//                                        const std::complex<T> &effectiveParam,
//                                        const data::VolumeFractions<T> &volumeFractions) {
//     std::complex<T> F{ 0.0, 0.0 };
//     size_t i{ 0 };
//     for (const auto &c : composite)
//         F += ema::bruggeman::term::v3D::anisotropicParRandDist(c, effectiveParam, volumeFractions[i++]);
//     return F;
// }

// template <typename T> requires std::is_floating_point_v<T>
// std::complex<T> anisotropicParRandDist(const std::vector<ema::data::MaterialNode<T>> &composite,
//                                        const std::complex<T> &effectiveParam,
//                                        const data::FillersParams<T> &fillersParams) {
//     std::complex<T> F{ 0.0, 0.0 };
//     // TODO
//     for (size_t pos{0}; pos < composite.size(); ++pos)
//         F += (fillersParams.contains(pos)
//                 ? ema::bruggeman::term::v3D::anisotropicParRandDist(composite[pos], effectiveParam, fillersParams[pos].value())
//                 : ema::bruggeman::term::v3D::anisotropicParRandDist(composite[pos], effectiveParam));
//     return F;
// }

}   // namespace ema::bruggeman::funcEMA


#endif
