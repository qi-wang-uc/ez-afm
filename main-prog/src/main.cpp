#include <iostream>
#include <string>
#include "../include/config.hpp"
#include "../include/util.hpp"

int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cout << "\nEZAFM> Easy Atomic Force Miscroscopy Simulations (v2.0)\n";
    if(argc < 2) return 1;
	
	/* Read configure file. */
    std::string config_name = argv[1];
    Config config;
    if(!config.read_config(config_name)) return 1;

    /* Execute configure file. */
    auto tick = TimeStamp::now();
    if(!config.exec_config()) return 1;
    std::cout << "\nEZAFM> Normal termination of program.\n";
    auto tock = TimeStamp::now();
    time_elapsed(tick, tock);

    return 0;
}
