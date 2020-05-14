#include "pdb.h"
#include "../utils/trim.h"

// #include <iostream>
// #include <fstream>
// #include <vector>
#include <string>
// #include <map>
// #include <cstring>

// #include <sstream>
// #include <algorithm>

namespace mol{
	// int PDB::PDB_number = 0;
	
	// const std::regex PDBParser::r_CRYST1_ = std::regex("");
	// const std::regex PDBParser::r_END_ = std::regex("");
	const std::regex PDBParser::r_HEADER_ = std::regex("^(HEADER)[ ]{4}([^[:space:]]{0,})");
	// const std::regex PDBParser::r_NUMMDL_ = std::regex("");
	// const std::regex PDBParser::r_MASTER_ = std::regex("");
	// const std::regex PDBParser::r_ORIGXn_ = std::regex("");
	// const std::regex PDBParser::r_SCALEn_ = std::regex("");
	// const std::regex PDBParser::r_AUTHOR_ = std::regex("");
	// const std::regex PDBParser::r_CAVEAT_ = std::regex("");
	// const std::regex PDBParser::r_COMPND_ = std::regex("");
	// const std::regex PDBParser::r_EXPDTA_ = std::regex("");
	// const std::regex PDBParser::r_MDLTYP_ = std::regex("");
	// const std::regex PDBParser::r_KEYWDS_ = std::regex("");
	// const std::regex PDBParser::r_OBSLTE_ = std::regex("");
	// const std::regex PDBParser::r_SOURCE_ = std::regex("");
	// const std::regex PDBParser::r_SPLIT_ = std::regex("");
	// const std::regex PDBParser::r_SPRSDE_ = std::regex("");
	// const std::regex PDBParser::r_TITLE_ = std::regex("");
	// const std::regex PDBParser::r_ANISOU_ = std::regex("(^ANISOU)");
	const std::regex PDBParser::r_ATOM_ = std::regex("^(ATOM)");
	// const std::regex PDBParser::r_ATOM_ = std::regex("(ATOM)\\s{2,}(\\d{1,5})(\\s{2,}(\\S{1,3})|\\s{1,}(\\S{1,4}))");
	const std::regex PDBParser::r_HETATM_ = std::regex("^(HETATM)");
	// const std::regex PDBParser::r_CISPEP_ = std::regex("");
	// const std::regex PDBParser::r_CONECT_ = std::regex("");
	// const std::regex PDBParser::r_DBREF_ = std::regex("");
	// const std::regex PDBParser::r_HELIX_ = std::regex("");
	// const std::regex PDBParser::r_HET_ = std::regex("");
	// const std::regex PDBParser::r_LINK_ = std::regex("");
	// const std::regex PDBParser::r_MODRES_ = std::regex("");
	// const std::regex PDBParser::r_MTRIXn_ = std::regex("");
	// const std::regex PDBParser::r_REVDAT_ = std::regex("");
	// const std::regex PDBParser::r_SEQADV_ = std::regex("");
	// const std::regex PDBParser::r_SHEET_ = std::regex("");
	// const std::regex PDBParser::r_SSBOND_ = std::regex("");
	// const std::regex PDBParser::r_FORMUL_ = std::regex("");
	// const std::regex PDBParser::r_HETNAM_ = std::regex("");
	// const std::regex PDBParser::r_HETSYN_ = std::regex("");
	// const std::regex PDBParser::r_SEQRES_ = std::regex("");
	// const std::regex PDBParser::r_SITE_ = std::regex("");
	const std::regex PDBParser::r_ENDMDL_ = std::regex("^(ENDMDL)");
	const std::regex PDBParser::r_MODEL_ = std::regex("^(MODEL)\\s{1,}(\\d{1,})");
	// const std::regex PDBParser::r_TER_ = std::regex("");
	// const std::regex PDBParser::r_JRNL_ = std::regex("");
	// const std::regex PDBParser::r_REMARK_ = std::regex("");
	// try{
	// }catch(std::regex_error e){
	// 	std::cout << e.what() << "\ncode: " << e.code() << std::endl;
	// }


	const PDBParser::RECORD_TYPE PDBParser::parseRECORD(const std::string &line){
        if(util::Trim(line.substr(0,6))=="ATOM") return ATOM_;
		if(util::Trim(line.substr(0,6))=="HETATM") return HETATM_;
		if(util::Trim(line.substr(0,6))=="MODEL") return MODEL_;
		if(util::Trim(line.substr(0,6))=="ENDMDL") return ENDMDL_;
    }

	const int PDBParser::parseMODEL(const std::string &sline){
		std::smatch results;
		std::regex_search(sline,results,PDBParser::r_MODEL_);
		// std::cout << results.str(0) << std::endl;
		// std::cout << "record: " << results.str(1) << std::endl;
		// std::cout << "serial: " << results.str(2) << std::endl;
		// exit(2);
		return std::stoi(results.str(2));
	}
	void PDBParser::parseENDMDL(const std::string &sline){
		// return 1;
	}
	const Atom PDBParser::parseATOM(const std::string &sline){
		Atom atom;
		atom.record_name = util::Trim(sline.substr(0,6));
		atom.serial = std::stoi(util::Trim(sline.substr(6,5)));
		atom.name = util::Trim(sline.substr(12,4));
		atom.alt_loc = util::Trim(sline.substr(16,1));
		atom.res_name = util::Trim(sline.substr(17,3));
		atom.chain_name = sline.substr(21,1);
		atom.res_serial = std::stoi(util::Trim( sline.substr(22,4)));
		atom.i_code = sline.substr(26,1);
		atom.x = std::stod(util::Trim( sline.substr(30,8)));
		atom.y = std::stod(util::Trim( sline.substr(38,8)));
		atom.z = std::stod(util::Trim( sline.substr(46,8)));
		atom.occupancy = std::stod(util::Trim( sline.substr(54,6)));
		atom.temp_factor = std::stod(util::Trim( sline.substr(60,6)));
		atom.element = util::Trim(sline.substr(76,2));
		atom.charge = util::Trim(sline.substr(78,2));
		// atom.partial_charge = std::stod(atom.charge);
		return atom;
	}
	const Atom PDBParser::parseHETATM(const std::string &sline){
		Atom atom;
		atom.record_name = util::Trim(sline.substr(0,6));
		atom.serial = std::stoi(util::Trim( sline.substr(6,5)));
		atom.name = util::Trim(sline.substr(12,4) );
		atom.alt_loc = util::Trim(sline.substr(16,1));
		atom.res_name = util::Trim(sline.substr(17,3));
		atom.chain_name = sline.substr(21,1);
		atom.res_serial = std::stoi(util::Trim( sline.substr(22,4)));
		atom.i_code = sline.substr(26,1);
		atom.x = std::stod(util::Trim(sline.substr(30,8)));
		atom.y = std::stod(util::Trim(sline.substr(38,8)));
		atom.z = std::stod(util::Trim(sline.substr(46,8)));
		atom.occupancy = std::stod(util::Trim( sline.substr(54,6)));
		atom.temp_factor = std::stod(util::Trim( sline.substr(60,6)));
		atom.element = util::Trim(sline.substr(76,2));
		atom.charge = util::Trim(sline.substr(78,2));
		// atom.partial_charge = std::stod(atom.charge);
		return atom;
	}
}//namespace bio
