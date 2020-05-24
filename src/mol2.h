#pragma once

// #include <iostream>
#include <fstream>
#include <string>
// #include <vector>
#include <regex>
#include "atom.h"
#include "molecule.h"
#include "str.h"


namespace tmd {

	enum MOL2_BLOCK_TYPE {MOL2_MOLECULE,MOL2_ATOM,MOL2_BOND,MOL2_SUBSTRUCTURE,MOL2_UNKNOWN,MOL2_NONE};

	inline const MOL2_BLOCK_TYPE to_mol2_block_type(const std::string s) {
		if(trim(s)=="@<TRIPOS>MOLECULE") return MOL2_MOLECULE;
		else if(trim(s)=="@<TRIPOS>ATOM") return MOL2_ATOM;
        else if(trim(s)=="@<TRIPOS>BOND") return MOL2_BOND;
        else if(trim(s)=="@<TRIPOS>SUBSTRUCTURE") return MOL2_SUBSTRUCTURE;
        else if(trim(s).substr(0,9)=="@<TRIPOS>") return MOL2_UNKNOWN;
        else return MOL2_NONE;
	}

	inline const Atom parse_mol2_atom(const int& ci, const std::string &sline) {
		std::vector<std::string> atom_items = string2vector(sline);
		const std::string sybyl_atom_type = trim(atom_items[5]);
		Atom atom(sybyl_atom_type,SYBYL);
		// std::cout << sline << std::endl;
		atom.set_charge(std::stod(atom_items[8]));
		const Float x = std::stod(atom_items[2]);
		const Float y = std::stod(atom_items[3]);
		const Float z = std::stod(atom_items[4]);
		atom.set_coord(Vec3d(x,y,z));
		atom.set_name(trim(atom_items[1]));
		atom.set_res_name(trim(atom_items[7]));
		atom.set_chain_name("");
		atom.set_serial(std::stoi(atom_items[0]));

		atom.set_res_serial(std::stoi(atom_items[6]));
		atom.set_context_index(ci);
		return atom;
	}




	inline void read_rna_mol2(const std::string rna_path, RNA& r, std::ostream& tee) {
		std::ifstream in_rna(rna_path);
		assert(in_rna);
		std::string sline;

		int model_num = 0;
		MOL2_BLOCK_TYPE Mol2_Block_Type_Flag = MOL2_NONE;

		while(std::getline(in_rna,sline)) {
			int context_index = r.contexts.size();
			r.contexts.push_back(Context(context_index,sline));
			// tee << sline << std::endl;
			if(sline=="") continue;

			MOL2_BLOCK_TYPE Temp_Mol2_Block_Flag = to_mol2_block_type(sline);
			if(Temp_Mol2_Block_Flag != MOL2_NONE) {
				Mol2_Block_Type_Flag = Temp_Mol2_Block_Flag;
				std::getline(in_rna,sline);
			}
			switch(Mol2_Block_Type_Flag) {
				case MOL2_ATOM: {
					const Atom a = parse_mol2_atom(context_index,sline);
					int a_i = r.add_atom(a);
					int r_a_i = r.add_ref_atom(a);
					assert(a_i == r_a_i);
					break;
				}
				case MOL2_MOLECULE:
					++model_num;
					break;
				case MOL2_BOND:
				case MOL2_SUBSTRUCTURE:
				case MOL2_UNKNOWN:
					;
					break;
				case MOL2_NONE:
				default:
					assert(false);
					break;
			}
		}
		in_rna.close();

	}

}