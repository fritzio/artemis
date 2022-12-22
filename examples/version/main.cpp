#include "version.hpp"
#include <iostream>

int main(void) {
    std::cout << "THRUST VERSION : " + artemis::version::thrust() << std::endl;
}