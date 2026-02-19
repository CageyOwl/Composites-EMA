#ifndef EMA_HXX
#define EMA_HXX

#include <complex>
#include <type_traits>
#include "data.hxx"


//  Naming: Par - particles, Dist - distribution
namespace bruggeman {
namespace term {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDistTerm(const data::MaterialNode<T> &node,
                                              const std::complex<T> &effectiveParam) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    numerator{ node.getParam() - effectiveParam };
    for (const auto L : node.getDepolarizationF())
        termsSum += ( numerator / (effectiveParam + L * numerator) );
    return (termsSum * node.getVolumeFraction() / 3.0);
}

}   // namespace term

namespace term_variable {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDistTerm(const data::MaterialNode<T> &node,
                                              const std::complex<T> &effectiveParam,
                                              const T volumeFraction) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    numerator{ node.getParam() - effectiveParam };
    for (const auto L : node.getDepolarizationF())
        termsSum += ( numerator / (effectiveParam + L * numerator) );
    return (termsSum * volumeFraction / 3.0);
}

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDistTerm(const data::MaterialNode<T> &node,
                                              const std::complex<T> &effectiveParam,
                                              const std::complex<T> &param) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    numerator{ param - effectiveParam };
    for (const auto L : node.getDepolarizationF())
        termsSum += ( numerator / (effectiveParam + L * numerator) );
    return (termsSum * node.getVolumeFraction() / 3.0);
}

}   // namespace term_variable

namespace derivative_term {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDistTerm(const data::MaterialNode<T> &node,
                                              const std::complex<T> &effectiveParam) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    param{ node.getParam() },
                    denominator_root{ 0.0, 0.0 };
    for (const auto L : node.getDepolarizationF())
    {
        denominator_root = effectiveParam + L * (param - effectiveParam);
        termsSum += -param / (denominator_root * denominator_root);
    }
    return (termsSum * node.getVolumeFraction() / 3.0);
}

}   // namespace derivative_term

namespace derivative_term_variable {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDistTerm(const data::MaterialNode<T> &node,
                                              const std::complex<T> &effectiveParam,
                                              const T volumeFraction) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    param{ node.getParam() },
                    denominator_root{ 0.0, 0.0 };
    for (const auto L : node.getDepolarizationF())
    {
        denominator_root = effectiveParam + L * (param - effectiveParam);
        termsSum += -param / (denominator_root * denominator_root);
    }
    return (termsSum * volumeFraction / 3.0);
}

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDistTerm(const data::MaterialNode<T> &node,
                                              const std::complex<T> &effectiveParam,
                                              const std::complex<T> &param) {
    std::complex<T> termsSum{ 0.0, 0.0 },
                    denominator_root{ 0.0, 0.0 };
    for (const auto L : node.getDepolarizationF())
    {
        denominator_root = effectiveParam + L * (param - effectiveParam);
        termsSum += -param / (denominator_root * denominator_root);
    }
    return (termsSum * node.getVolumeFraction() / 3.0);
}

}   // namespace derivative_term_variable

namespace funcEMA {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDist(const std::vector<data::MaterialNode<T>> composite,
                                          const std::complex<T> &effectiveParam) {
    std::complex<T> F{ 0.0, 0.0 };
    for (const auto &c : composite)
        F += term::anisotropicParChaoticDistTerm(c, effectiveParam);
    return F;
}

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDist(const std::vector<data::MaterialNode<T>> composite,
                                          const std::complex<T> &effectiveParam,
                                          const std::vector<T> &volumeFractions) {
    std::complex<T> F{ 0.0, 0.0 };
    size_t i{ 0 };
    for (const auto &c : composite)
        F += term_variable::anisotropicParChaoticDistTerm(c, effectiveParam, volumeFractions[i++]);
    return F;
}

// TODO: check the correctness of the parameters injection!!!!!!!
template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParChaoticDist(const std::vector<data::MaterialNode<T>> composite,
                                          const std::complex<T> &effectiveParam,
                                          const std::complex<T> &fillerParam) {
    std::complex<T> F{ 0.0, 0.0 };

    size_t fillerNodePosition{ 0 };
    std::complex<T> fillerNodeParam{ composite[fillerNodePosition].getParam() };

    for (size_t i{ 1 }; i < composite.size(); ++i) {
        if (composite[i].getParam().real() > fillerNodeParam.real()) {
            fillerNodePosition = i;
            fillerNodeParam = composite[i].getParam();
        }
    }

    for (size_t j{ 0 }; j < composite.size(); ++j) {
        F += (
            j == fillerNodePosition
                ? term_variable::anisotropicParChaoticDistTerm(composite[j], effectiveParam, fillerParam)
                : term::anisotropicParChaoticDistTerm(composite[j], effectiveParam)
        );
    }

    return F;
}

}   // namespace funcEMA

namespace derivative_funcEMA {

// TODO:

}   // namespace derivative_funcEMA
}   // namespace bruggeman


#endif  // EMA_HXX
