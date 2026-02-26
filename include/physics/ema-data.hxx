#ifndef EMA_DATA_H
#define EMA_DATA_H

#include <complex>
#include <type_traits>
#include <unordered_map>
#include <vector>


namespace ema::data {

template <typename T> requires std::is_floating_point_v<T>
class MaterialNode {
public:
    MaterialNode(const T volumeFraction,
                 const std::complex<T> &param,
                 const std::vector<T> &depolarizationFV = { 1.0/3.0, 1.0/3.0, 1.0/3.0 })
        : volumeFraction(volumeFraction), param(param), depolarizationFV(depolarizationFV) {
    }

    T getVolumeFraction() const { return this->volumeFraction; }
    std::complex<T> getParam() const { return this->param; }
    std::vector<T> getDepolarizationFV() const { return this->depolarizationFV; }

    void setVolumeFraction(const T volumeFraction) { this->volumeFraction = volumeFraction; }
    void setParam(const std::complex<T> &param) { this->param = param; }
    void setDepolarizationFV(const std::vector<T> &depolarizationFV) { this->depolarizationFV = depolarizationFV; }
    void setSingleDepolarizationF(const T sDepolarizationF, const size_t index) {
        if (index < this->depolarizationFV.size()) {
            this->depolarizationFV[index] = sDepolarizationF;
        }
    }

private:
    T volumeFraction;
    std::complex<T> param;
    std::vector<T> depolarizationFV;
};

template <typename T> requires std::is_floating_point_v<T>
struct VolumeFractions {
    std::vector<T> volumeFractions;
};

template <typename T> requires std::is_floating_point_v<T>
struct FillersParams {
    std::unordered_map<size_t, std::complex<T>> fillersParams;
};

enum Axis {
    x,y,z
};

}   // namespace ema::data


#endif
