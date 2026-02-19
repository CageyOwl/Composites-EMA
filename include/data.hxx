#ifndef DATA_H
#define DATA_H

#include <complex>
#include <type_traits>
#include <unordered_map>
#include <vector>


namespace data {

template <typename T> requires std::is_floating_point_v<T>
class MaterialNode {
public:
    MaterialNode(const T volumeFraction,
                 const std::complex<T> &param,
                 const std::vector<T> &depolarizationF = { 1.0/3.0, 1.0/3.0, 1.0/3.0 })
        : volumeFraction(volumeFraction), param(param), depolarizationF(depolarizationF) {}

    T getVolumeFraction() const { return this->volumeFraction; }
    std::complex<T> getParam() const { return this->param; }
    std::vector<T> getDepolarizationF() const { return this->depolarizationF; }

    // TODO: Implement setters
    // TODO: Implement the volume fractions sum check

private:
    T volumeFraction;
    std::complex<T> param;
    std::vector<T> depolarizationF;
};

enum Mode {
    E_VS_FILLER_FRACTION,
    E_VS_E_FILLER
};

extern const std::unordered_map<unsigned int, std::string> modeCaption;

}   // namespace data


#endif
