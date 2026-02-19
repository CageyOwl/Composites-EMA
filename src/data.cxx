#include <unordered_map>
#include "../include/data.hxx"


namespace data {

const std::unordered_map<unsigned int, std::string> modeCaption {
    {E_VS_FILLER_FRACTION, "Effective parameter vs volume fraction of the filler"},
    {E_VS_E_FILLER, "Effective parameter vs dielectric constant of the filler"}
};

}   // namespace data
