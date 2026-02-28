#include <iostream>
#include "physics/ema-data.hxx"
#include "physics/ema-averaging.hxx"
#include "physics/ema-bruggeman-terms.hxx"


int main() {
    // TODO
    namespace ebf = ema::bruggeman::func;

    std::vector<ema::data::MaterialNode<double>> composite{
        {0.33, {60000.0, 0.0}, {1.0/3.0, 1.0/3.0, 1.0/3.0}},
        {0.67, {2.6, 0.0}, {1.0/3.0, 1.0/3.0, 1.0/3.0}}
    };

    std::cout << ebf::term::calcTerm(composite[0].getVolumeFraction(),
                                     ebf::term::axialContributionsSum(composite[0], {7.0,0.0}, ebf::term::axial::contribution<double>),
                                     ema::averaging::isotropic3D<double>) << '\n';

    std::cout << ebf::term::calcTerm(composite[0].getVolumeFraction(),
                                     ebf::term::axialContributionsSum(composite[0], {7.0,0.0}, ebf::term::axial::contribution<double>),
                                     ema::averaging::isotropic2D<double>) << '\n';

    return 0;
}
