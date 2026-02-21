#include "ema-data.hxx"


int main() {
    // TODO
    std::vector<ema_data::MaterialNode<double>> composite{
        {0.2, {60000.0, 0.0}, {1.0/3.0, 1.0/3.0, 1.0/3.0}},
        {0.2, {2.6, 0.0}, {1.0/3.0, 1.0/3.0, 1.0/3.0}}
    };

    for (double i{ 0.0 }; i < 1.01; i += 0.01)
        ;

    return 0;
}
