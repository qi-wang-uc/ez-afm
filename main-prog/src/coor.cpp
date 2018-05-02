#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <iomanip>
#include "../include/coor.hpp"
#include "../include/util.hpp"
#include "../include/psf.hpp"

/* This function reads the CHARMM coordinate file as string streams.
   Pros: 1. Can read both NORMAL or EXPAND format wihout modifications.
         2. Make the coding easier.
   Cons: It assumes there are no comments in the coordinates.*/
bool CorData::read_cor(const std::string &inp_name) {
	std::ifstream inp_file(inp_name);
	this->_sizeof_coor = 0;
	if(!inp_file.is_open()) {
		std::cout << "ERROR> Could not open file [" << inp_name << "]" 
                  << std::endl;
		return false;
	} else {
		std::cout << "ReadCOR> Reading coordinate file from [" << inp_name << "]"
                  << std::endl; 
	}
	std::string each_line;
	std::stringstream each_stream;
	/* For temporary variable holders */
	std::string tmp_etc;  // direct all unwanted string streams here.
	double tmp_xcoor = 0.0;
	double tmp_ycoor = 0.0;
	double tmp_zcoor = 0.0;
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
			each_stream >> this->_sizeof_coor;
		} else {
			each_stream >> tmp_etc >> tmp_etc >> tmp_etc >> tmp_etc
					>> tmp_xcoor >> tmp_ycoor >> tmp_zcoor;
			this->_xcoor.push_back(tmp_xcoor);
			this->_ycoor.push_back(tmp_ycoor);
			this->_zcoor.push_back(tmp_zcoor);
		}
		each_stream.clear();
	}
	inp_file.close();
	std::cout << "ReadCOR> After reading, (" << this->_sizeof_coor
              << ") atoms found in coordinate file." << std::endl;
	return true;
}

void CorData::move_cor(Vec3d& dr, const size_t& atomid) {
	this->_xcoor.at(atomid) += dr.x;
	this->_ycoor.at(atomid) += dr.y;
	this->_zcoor.at(atomid) += dr.z;
}

const float* CorData::p_xcoor(void) {
	return this->_xcoor.data();
}

const float* CorData::p_ycoor(void) {
	return this->_ycoor.data();
}

const float* CorData::p_zcoor(void) {
	return this->_zcoor.data();
}

const Vec3d CorData::get_atom_coor(const size_t& atomid) const {
	return Vec3d(static_cast<double>(this->_xcoor.at(atomid)),
                 static_cast<double>(this->_ycoor.at(atomid)),
                 static_cast<double>(this->_zcoor.at(atomid)));
}

void CorData::print(void) const {
	for(size_t i=0; i<this->_sizeof_coor; ++i) {
		std::cout << std::setw(8)  << i+1
				  << std::setw(12) << this->_xcoor.at(i)
				  << std::setw(12) << this->_ycoor.at(i)
				  << std::setw(12) << this->_zcoor.at(i)
				  << std::endl;
	}
}