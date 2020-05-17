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

	enum PDBQT_PARSE_TYPE {PDBQT_ROOT,PDBQT_ENDROOT,PDBQT_BRANCH,PDBQT_ENDBRANCH,PDBQT_TORSDOF,PDBQT_REMARK,PDBQT_ATOM,PDBQT_HETATM,PDBQT_UNKNOWN};

	//store each atom's information while parsing the input file
	struct Atom_Tag {
		Atom_Index atom_index;
		Node_Id node_id;
		Atom_Tag(const Atom_Index ai, const Node_Id& ni) : atom_index(ai), node_id(ni) {}
	};
	//store each node's connection information
	struct Node_Info {
		Atom_Index parent_bond_atom_index;
		Atom_Index self_bond_atom_index;
		Node_Id node_id;
		Node_Info(const Atom_Index pbai, const Atom_Index sbai, const Node_Id& ni) : parent_bond_atom_index(pbai), self_bond_atom_index(sbai), node_id(ni) {}
	};

	inline const PDBQT_PARSE_TYPE to_pdbqt_parse_type(const std::string s) {
		const std::vector<std::string> vs = string2vector(s);
		const std::string& r = vs[0];
		if(r=="ATOM") return PDBQT_ATOM;
		else if(r=="HETATM") return PDBQT_HETATM;
		else if(r=="REMARK") return PDBQT_REMARK;
		else if(r=="BRANCH") return PDBQT_BRANCH;
		else if(r=="ENDBRANCH") return PDBQT_ENDBRANCH;
		else if(r=="ROOT") return PDBQT_ROOT;
		else if(r=="ENDROOT") return PDBQT_ENDROOT;
		else if(r=="TORSDOF") return PDBQT_TORSDOF;
		else return PDBQT_UNKNOWN;
	}

	// const std::regex r_REMARK_ = std::regex("^(REMARK)");
	// const std::regex r_ATOM_ = std::regex("^(ATOM)");
	// const std::regex r_HETATM_ = std::regex("^(HETATM)");

	// const std::regex PDBParser::r_ATOM_ = std::regex("(ATOM)\\s{2,}(\\d{1,5})(\\s{2,}(\\S{1,3})|\\s{1,}(\\S{1,4}))");

	// const std::regex r_ROOT_ = std::regex("^(ROOT)");
	// const std::regex r_ENDROOT_ = std::regex("^(ENDROOT)");
	// const std::regex r_BRANCH_ = std::regex("^(BRANCH)");
	// const std::regex r_ENDBRANCH_ = std::regex("^(ENDBRANCH)");
	// const std::regex r_TORSDOF_ = std::regex("^(TORSDOF)");

	inline const Node_Info parse_pdbqt_branch(const std::string &sline, const Node_Id& ni) {
		// std::smatch results;
		// std::regex_search(sline,results,r_BRANCH_);
		// std::cout << results.str(0) << std::endl;
		// std::cout << "record: " << results.str(1) << std::endl;
		// std::cout << "serial: " << results.str(2) << std::endl;
		// exit(2);
		// return std::stoi(results.str(2));
		const std::vector<std::string> vs = string2vector(sline);
		return Node_Info(std::stoi(vs[1])-1, std::stoi(vs[2])-1, ni);
	}

	inline const Atom parse_pdbqt_atom(const Context_Index& ci, const std::string &sline) {
		const std::string autodock_atom_type = trim(sline.substr(77,2));
		const std::string sybyl_atom_type = trim(sline.substr(82,5));
		Atom atom(sybyl_atom_type,autodock_atom_type,SYBYL_AD4);
		atom.set_charge(std::stod(sline.substr(70,6)));
		const Float x = std::stod(trim(sline.substr(30,8)));
		const Float y = std::stod(trim(sline.substr(38,8)));
		const Float z = std::stod(trim(sline.substr(46,8)));
		atom.set_coord(Vec3d(x,y,z));
		atom.set_name(trim(sline.substr(12,4)));
		atom.set_res_name(trim(sline.substr(17,3)));
		atom.set_chain_name(sline.substr(21,1));
		// log << trim(sline.substr(6,5)) << " " << trim( sline.substr(22,4)) << std::endl;
		atom.set_serial(std::stoi(trim(sline.substr(6,5))));

		const std::string res_serial_str = trim( sline.substr(22,4));
		if(res_serial_str!="") {
			atom.set_res_serial(std::stoi(res_serial_str));
		} else {
			atom.set_res_serial(0);
		}
		atom.set_context_index(ci);
		return atom;
	}




	inline void read_ligand_pdbqt(const std::string ligand_path, Ligand& ligand, std::ofstream& log) {
		std::ifstream in_ligand(ligand_path);
		assert(in_ligand);
		std::string sline;

		//put atom into ligand and gather tree information
		Node_Index depth_level = 0;
		std::vector<Node_Index> width_in_each_depth;
		std::vector<Atom_Tag> atom_tags;
		std::vector<Node_Info> node_infos;
		std::map<std::string,Node_Index> node_id_to_node_index;
		while(std::getline(in_ligand,sline)) {
			Context_Index line_index = ligand.contexts.size();
			ligand.contexts.push_back(Context(line_index,sline));
			// log << sline << std::endl;
			if(sline=="") continue;
			PDBQT_PARSE_TYPE ppt = to_pdbqt_parse_type(sline);
			switch(ppt) {
				case PDBQT_ATOM:
				case PDBQT_HETATM: {
					const Atom a = parse_pdbqt_atom(line_index,sline);
					Atom_Index a_i = ligand.add_atom(a);
					if(a_i==0) {
						//make sure ligand atom start from 1
						assert(a.get_serial()==1);
					}
					Atom_Index r_a_i = ligand.add_ref_atom(a);
					assert(a_i == r_a_i);
					Node_Id ni = Node_Id(depth_level,width_in_each_depth[depth_level]);
					atom_tags.push_back(Atom_Tag(a_i,ni));
					assert(a_i==atom_tags.size()-1);
					break;
				}
				case PDBQT_BRANCH: {
					++depth_level;
					if(width_in_each_depth.size() <= depth_level) {
						width_in_each_depth.push_back(0);
					} else {
						++width_in_each_depth[depth_level];
					}
					Node_Id ni = Node_Id(depth_level,width_in_each_depth[depth_level]);
					Node_Info bni = parse_pdbqt_branch(sline, ni);
					node_infos.push_back(bni);
					node_id_to_node_index.insert({std::to_string(depth_level)+"-"+std::to_string(width_in_each_depth[depth_level]), node_infos.size()-1});
					break;
				}
				case PDBQT_ENDBRANCH:
					--depth_level;
					assert(depth_level>=0);
					break;
				case PDBQT_ROOT:
					assert(depth_level==0);
					assert(width_in_each_depth.size()==0);
					width_in_each_depth.push_back(0);
					//root atom node info push back, root should start from 1
					node_infos.push_back({0, 0, {depth_level, width_in_each_depth[depth_level]}});
					node_id_to_node_index.insert({"0-0", node_infos.size()-1});
				case PDBQT_ENDROOT:
					;
					break;
				case PDBQT_TORSDOF:
					;
					break;
				case PDBQT_REMARK:
					;
					break;
				case PDBQT_UNKNOWN:
				default:
					assert(false);
					break;
			}
		}
		in_ligand.close();

		assert(ligand.atom_num()==atom_tags.size());

		log << "Initialize ligand's nodes..." << std::endl;
		//Initialize ligand's tree i.e all nodes, this must be done first since the vector will allocate new memory when push_back new nodes, thus old pointer will fail
		const Node_Index max_depth = width_in_each_depth.size();
		log << "total max depth--->" << max_depth << std::endl;
		for(Node_Index depth_index = 0; depth_index < width_in_each_depth.size(); ++depth_index) {
			const Node_Index num_of_nodes = width_in_each_depth[depth_index]+1;
			log << "depth " << depth_index << " has " << num_of_nodes << " nodes" << std::endl;
		}

		for(const Node_Info& ni : node_infos) {
			Node tmp_node = Node();
			tmp_node.set_node_id(ni.node_id.depth_index, ni.node_id.width_index);
			tmp_node.set_parent_bond_atom_index(ni.parent_bond_atom_index);
			tmp_node.set_self_bond_atom_index(ni.self_bond_atom_index);
			ligand.nodes.push_back(tmp_node);
		}

		//assign atom range to each node
		log << "assign atom range to each node..." << std::endl;
		for(Atom_Index i = 0; i < atom_tags.size(); ++i) {
			const Node_Index depth_index = atom_tags.at(i).node_id.depth_index;
			const Node_Index width_index = atom_tags.at(i).node_id.width_index;
			Node_Index ni = node_id_to_node_index.at(std::to_string(depth_index)+"-"+std::to_string(width_index));
			ligand.nodes[ni].add_atom_index(atom_tags.at(i).atom_index);
		}

		log << "Establish parent and children's connection..." << std::endl;
		//Establish parent and children's connection, and parent bond atom index, self bond atom index
		//for root
		ligand.nodes[0].set_parent_pointer(nullptr);
		//for non root
		for(auto ni = 1; ni < node_infos.size(); ++ni) {
			Node_Id parent_node_id = atom_tags.at(node_infos[ni].parent_bond_atom_index).node_id;
			Node_Index parent_node_index = node_id_to_node_index.at(std::to_string(parent_node_id.depth_index)+"-"+std::to_string(parent_node_id.width_index));
			ligand.nodes[ni].set_parent_pointer(ligand.nodes[parent_node_index].get_self_pointer());
			ligand.nodes[parent_node_index].add_child_pointer(ligand.nodes[ni].get_self_pointer());
		}

		//for root, set origin and rot_axis, rot_angle to (0,0,0) (0,0,1) 0, because they will reset afterwards
		ligand.nodes[0].set_origin(Vec3d(0,0,0));
		ligand.nodes[0].set_rot_axis(Vec3d(0,0,1));
		ligand.nodes[0].set_rot_angle(0.0);
		// set non root node origin(parent bond atom coord), rot axis(unit vector from parent bond atom to self bond atom) and set all rot angle to 0
		for(Node_Index ni = 1; ni < ligand.nodes.size(); ++ni) {
			Node& n = ligand.nodes[ni];
			const Atom_Index parent_bond_atom_index = n.get_parent_bond_atom_index();
			const Atom_Index self_bond_atom_index = n.get_self_bond_atom_index();
			const Vec3d origin = ligand.atoms[parent_bond_atom_index].get_coord();
			Vec3d rot_axis = ligand.atoms[self_bond_atom_index].get_coord() - origin;
			rot_axis.normalize();
			n.set_origin(origin);
			n.set_rot_axis(rot_axis);
			n.set_rot_angle(0.0);
		}
	}

}