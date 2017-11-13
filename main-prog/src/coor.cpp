#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <iomanip>
#include "../include/coor.h"
#include "../include/util.h"
#include "../include/psf.h"

/*******************************************************************
                         COORDINATE DATA  
  For historical reasons, coordinates in dcd file are populated like 
  [natom*X, natom*Y, natom*Z]. Thus perhaps it's not a good idea to 
  decalre coordinates storage as a 3*Natom array.
*******************************************************************/
static std::vector<real> xcoor;
static std::vector<real> ycoor;
static std::vector<real> zcoor;

static unsigned int natom_from_cor;
/* This function reads the CHARMM coordinate file as string streams.
   Pros: 1. Can read both NORMAL or EXPAND format wihout modifications.
         2. Make the coding easier.
   Cons: It assumes there are no comments in the coordinates.*/
bool read_cor(const std::string &inp_name) {
	std::ifstream inp_file(inp_name);
	unsigned int natom_cor = 0;
	if(!inp_file.is_open()) {
		std::cout << "ERROR> Could not open file [" << inp_name 
			<< "]" << std::endl;
		return false;
	} else {
		std::cout << "READCOR> Reading coordinate file from [" 
			<< inp_name << "]" << std::endl; 
	}
	std::string each_line;
	std::stringstream each_stream;
	/* For temporary variable holders */
	std::string tmp_etc;  // direct all unwanted string streams here.
	real tmp_xcoor = 0.0;
	real tmp_ycoor = 0.0;
	real tmp_zcoor = 0.0;
	while(std::getline(inp_file, each_line)) {
		if(each_line.empty()) continue;
		/* The title may be printed out but not the comments.*/
		if(each_line.substr(0,1).compare("!")==0) continue;
		if(each_line.substr(0,1).compare("*")==0) {
			std::cout << each_line << std::endl;
			continue;
		}
		each_stream.str(each_line);
		/* The NATOM line: Here we assume the single word line is the NATOM line, 
		   but this does not work with any space separated comments. */
		if(1==n_of_words(each_line)) {
			each_stream >> natom_cor;
		} else {
			each_stream >> tmp_etc >> tmp_etc >> tmp_etc >> tmp_etc
					>> tmp_xcoor >> tmp_ycoor >> tmp_zcoor;
			xcoor.push_back(tmp_xcoor);
			ycoor.push_back(tmp_ycoor);
			zcoor.push_back(tmp_zcoor);
		}
		each_stream.clear();
	}
	inp_file.close();
	std::cout << "READCOR> After reading, (" << natom_cor 
		<< ") atoms found in coordinate file." << std::endl;
	natom_from_cor = natom_cor;
	return true;
}

unsigned int get_natom_cor() {
	return natom_from_cor;
}

bool is_good_cor() {
	return natom_from_cor >= 1;
}

real* p_xcoor() {
	return xcoor.data();
}

real* p_ycoor() {
	return ycoor.data();
}

real* p_zcoor() {
	return zcoor.data();
}

void print_coor() {
	if(0==natom_from_cor) return;
	for(unsigned int i=0; i<natom_from_cor; i++) {
		std::cout << std::setw(8)  << i+1
				  << std::setw(12) << xcoor.at(i)
				  << std::setw(12) << ycoor.at(i)
				  << std::setw(12) << zcoor.at(i)
				  << std::endl;
	}
}


const real& get_xcoor(const unsigned int atomid) {
	return xcoor[atomid];
}

const real& get_ycoor(const unsigned int atomid) {
	return ycoor[atomid];
}

const real& get_zcoor(const unsigned int atomid) {
	return zcoor[atomid];
}

void add_to_coor(const unsigned int atomid, const Vector& coor) {
	if(is_movable(atomid+1)) {
		xcoor[atomid] += coor.x;
		ycoor[atomid] += coor.y;
		zcoor[atomid] += coor.z;
	}
}