#include <iostream>
#include "physics/ema-data.hxx"
#include "physics/ema-bruggeman-terms.hxx"


int main() {
    // TODO
    std::vector<ema::data::MaterialNode<double>> composite{
        {0.33, {60000.0, 0.0}, {1.0/3.0, 1.0/3.0, 1.0/3.0}},
        {0.67, {2.6, 0.0}, {1.0/3.0, 1.0/3.0, 1.0/3.0}}
    };

    std::cout << ema::bruggeman::func::term::random3D<double>(ema::bruggeman::func::term::cmp::axialTerm<double>, composite[0], 4.0) << '\n';

    return 0;
}
