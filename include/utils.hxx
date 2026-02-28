#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <unordered_map>


namespace utils {

enum Mode {
    E_VS_FILLER_FRACTION,
    E_VS_E_FILLER
};

extern const std::unordered_map<unsigned int, std::string> modeCaption;

}


#endif
