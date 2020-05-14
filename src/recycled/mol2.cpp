#include "mol2.h"
#include "../utils/trim.h"
#include "../utils/string2vector.h"

// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <string>
// #include <map>
// #include <cstring>

// #include <sstream>
// #include <algorithm>

namespace mol{
	// const std::regex MOL2Parser::r_MOLECULE = std::regex("^(HEADER)[ ]{4}([^[:space:]]{0,})");
	
	// const std::regex MOL2Parser::r_ATOM = std::regex("^(ATOM)");
	// const std::regex MOL2Parser::r_ENDMDL = std::regex("^(ENDMDL)");
	// const std::regex MOL2Parser::r_MODEL = std::regex("^(MODEL)\\s{1,}(\\d{1,})");


	const MOL2Parser::BLOCK_TYPE MOL2Parser::parseTYPE(const std::string &line){
        if(util::Trim(line)=="@<TRIPOS>MOLECULE") return MOLECULE_;
		else if(util::Trim(line)=="@<TRIPOS>ATOM") return ATOM_;
        else if(util::Trim(line)=="@<TRIPOS>BOND") return BOND_;
        else if(util::Trim(line)=="@<TRIPOS>SUBSTRUCTURE") return SUBSTRUCTURE_;
        else if(util::Trim(line).substr(0,9)=="@<TRIPOS>") return UNKNOWN_;
        else return NONE_;
    }

	// const int MOL2Parser::parseMODEL(const std::string &sline){
	// 	std::smatch results;
	// 	std::regex_search(sline,results,PDBParser::r_MODEL);
	// 	std::cout << results.str(0) << std::endl;
	// 	std::cout << "record: " << results.str(1) << std::endl;
	// 	std::cout << "serial: " << results.str(2) << std::endl;
	// 	// exit(2);
	// 	return std::stoi(results.str(2));
	// }
	// void MOL2Parser::parseENDMDL(const std::string &sline){
	// 	// return 1;
	// }
    // void MOL2Parser::parseMOLECULE(const std::string &sline){
    //     std::vector<std::string> atom_items = string2vector(sline);
	// 	Atom atom;
		
	// 	return atom;
	// }
	const Atom MOL2Parser::parseATOM(const std::string &sline){
        std::vector<std::string> atom_items = util::String2Vector(sline);
		Atom atom;
		// atom.recordName = util::Trim( sline.substr(0,6) );
		atom.serial = std::stoi(atom_items[0]);
		atom.name = util::Trim(atom_items[1]);
		// atom.altLoc = util::Trim(sline.substr(16,1));
		atom.res_name = util::Trim(atom_items[7]);
		// atom.chainID = sline.substr(21,1);
		atom.res_serial = std::stoi(atom_items[6]);
		// atom.iCode = sline.substr(26,1);
		atom.x = std::stod(atom_items[2]);
		atom.y = std::stod(atom_items[3]);
		atom.z = std::stod(atom_items[4]);
		// atom.occupancy = std::stod( util::Trim( sline.substr(54,6) ) );
		// atom.tempFactor = std::stod( util::Trim( sline.substr(60,6) ) );
		atom.sybyl_atom_type = util::Trim(atom_items[5]);
		atom.element = util::Trim(atom_items[5].substr(0,atom_items[5].find(".")));
		atom.charge = util::Trim(atom_items[8]);
        atom.partial_charge = std::stod(atom_items[8]);

        
		return atom;
	}
}//namespace bio
