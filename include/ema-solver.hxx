#ifndef EMA_SOLVER_HXX
#define EMA_SOLVER_HXX

#include <type_traits>


namespace bruggeman {

template <typename T> requires std::is_floating_point_v<T>
class SolverEMA {
public:

private:

};




namespace funcEMA {

template <typename T> requires std::is_floating_point_v<T>
std::complex<T> anisotropicParRandDist(const std::vector<data::MaterialNode<T>> &composite,
                                          const std::complex<T> &effectiveParam) {
    std::complex<T> F{ 0.0, 0.0 };
    for (const auto &c : composite)
        F += term::anisotropicParRandDistTerm(c, effectiveParam);
    return F;
}

namespace var_vol_fraction {

    template <typename T> requires std::is_floating_point_v<T>
    std::complex<T> anisotropicParRandDist(const std::vector<data::MaterialNode<T>> &composite,
                                              const std::complex<T> &effectiveParam,
                                              const std::vector<T> &volumeFractions) {
        std::complex<T> F{ 0.0, 0.0 };
        size_t i{ 0 };
        for (const auto &c : composite)
            F += term::var_vol_fraction::anisotropicParRandDistTerm(c, effectiveParam, volumeFractions[i++]);
        return F;
    }

}   // namespace var_vol_fraction


}


#endif
