#include <iostream>
#include <string>
#include "../include/config.hpp"
#include "../include/util.hpp"

int main(int argc, char* argv[]) {
	std::cout << std::endl
			  << "EZAFM> Easy Atomic Force Miscroscopy Simulations (v1.0)" 
	          << std::endl;
	if(argc < 2) return 1;
	
	/* Read configure file. */
	std::string config_name = argv[1];
	Config config;
	if(!config.read_config(config_name)) return 1;

	/* Execute configure file. */
	auto tick = TimeStamp::now();
	if(!config.exec_config()) 
		return 1;
	else
		std::cout << std::endl 
                  << "EZAFM> Normal termination of program." 
				  << std::endl;
	auto tock = TimeStamp::now();
	time_elapsed(tick, tock);

	return 0;
}
