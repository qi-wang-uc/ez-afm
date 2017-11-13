#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include "../include/afm.h"
#include "../include/config.h"
#include "../include/read.h"
#include "../include/uservar.h"
#include "../include/mini.h"
#include "../include/dyna.h"
#include "../include/util.h"

int main(int argc, char* argv[]) {
	/* Initialize timer and check if configure file name is passed*/
	print_banner();
	/* Simple sanity check*/
	if(!has_input(argc)) return 1;

	/* Read configure file. */
	std::string config_name = argv[1];
	auto cmds_set = read_config(config_name);
	if(cmds_set.empty()) return 1;

	/* Now parse commands set*/
	for(auto it_cmds=cmds_set.begin(); it_cmds!=cmds_set.end(); ++it_cmds) {
		/* Processing system variables */
		if((*it_cmds).front()=="set")
			set_uservar(*it_cmds);

		/* Read RTF PRM PSF and COR */
		if((*it_cmds).front()=="read")
			//read_data(*it_cmds);
			if(!read_data(*it_cmds))
				return 1;

		/* Molecular mechanics part */
		if((*it_cmds).front()=="mini") 
			prep_mini(*it_cmds);

		if((*it_cmds).front()=="afm")
			if(!prep_afm(*it_cmds))
				return 1;

		if((*it_cmds).front()=="dyna")
			if(!prep_dyna(*it_cmds)) 
				return 1;

		/* Normal termination */
		if((*it_cmds).front()=="stop") 
			break;
	}
	//print_uservar();
    /* Finalize the program and print time statistics.*/
	return 0;
}
