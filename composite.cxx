#include <iostream>
#include "physics/ema-bruggeman-func.hxx"
#include "physics/ema-bruggeman-terms.hxx"
#include "physics/ema-data.hxx"


int main() {
    // TODO
    namespace ebf = ema::bruggeman::func;

    std::vector<ema::data::MaterialNode<double>> composite{
        {0.33, {60000.0, 0.0}, {1.0/3.0, 1.0/3.0, 1.0/3.0}},
        {0.67, {2.6, 0.0}, {1.0/3.0, 1.0/3.0, 1.0/3.0}}
    };

    std::cout << ebf::isotropic3D(composite, std::complex<double>{6.0, 0.0}, ebf::term::axial::contribution<double>) << '\n';

    return 0;
}
