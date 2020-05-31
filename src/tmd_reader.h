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

	enum PDBQT_PARSE_TYPE {PDBQT_ROOT,PDBQT_ENDROOT,PDBQT_BRANCH,PDBQT_ENDBRANCH,PDBQT_TORSDOF,PDBQT_REMARK,PDBQT_ATOM,PDBQT_HETATM,PDBQT_BOND,PDBQT_UNKNOWN};

	//store each atom's information while parsing the input file
	struct Atom_Tag {
		int atom_index;
		Node_Id node_id;
		Atom_Tag(const int ai, const Node_Id& ni) : atom_index(ai), node_id(ni) {}
	};
	//store each node's connection information
	struct Node_Info {
		int parent_bond_atom_index;
		int self_bond_atom_index;
		Node_Id node_id;
		Node_Info(const int pbai, const int sbai, const Node_Id& ni) : parent_bond_atom_index(pbai), self_bond_atom_index(sbai), node_id(ni) {}
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
		else if(r=="<TRIPOS>BOND") return PDBQT_BOND;
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

	struct Bond_Info {
		int a1;
		int a2;
		BOND_TYPE bt;
		Bond_Info(const int aa1, const int aa2, const BOND_TYPE bbt) : a1(aa1), a2(aa2), bt(bbt) {}
	};

	inline const Bond_Info parse_pdbqt_bond_and_add_bond(const std::string& sline) {
		const std::vector<std::string> vs = string2vector(sline);
		assert(vs.size()==5 || vs.size()==6);
		int a1 = std::stoi(vs[2])-1;
		int a2 = std::stoi(vs[3])-1;
		BOND_TYPE bt;
		if(vs[4]=="1") {
			bt = SINGLE_BOND;
		} else if(vs[4]=="2") {
			bt = DOUBLE_BOND;
		} else if(vs[4]=="3") {
			bt = TRIPLE_BOND;
		} else if(vs[4]=="am") {
			bt = AMIDE_BOND;
		} else if(vs[4]=="ar") {
			bt = AROMATIC_BOND;
		} else if(vs[4]=="du") {
			bt = DUMMY_BOND;
		} else if(vs[4]=="un") {
			bt = UNKNOWN_BOND;
		} else if(vs[4]=="nc") {
			bt = NOT_CONNECTED_BOND;
		} else {
			assert(false);
		}
		return Bond_Info(a1,a2,bt);
	}

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

	inline const Atom parse_pdbqt_atom(const std::string &sline) {
		const std::string autodock_atom_type_name = trim(sline.substr(77,2));
		const std::string sybyl_atom_type_name = trim(sline.substr(82,5));
		Atom atom;
		atom.set_ad4_type(string_to_ad4_type(autodock_atom_type_name));
		atom.set_sybyl_type(string_to_sybyl_type(sybyl_atom_type_name));
		atom.set_element_type(sybyl_type_to_element_type(atom.get_sybyl_type()));
		assert(atom.get_sybyl_type() != SYBYL_Unk && atom.get_ad4_type() != AD4_Unk);
		assert(atom.get_element_type() == ad4_type_to_element_type(atom.get_ad4_type()));
		atom.set_charge(std::stod(sline.substr(70,6)));
		const Float x = std::stod(trim(sline.substr(30,8)));
		const Float y = std::stod(trim(sline.substr(38,8)));
		const Float z = std::stod(trim(sline.substr(46,8)));
		atom.set_coord(Vec3d(x,y,z));
		atom.set_name(trim(sline.substr(12,4)));
		atom.set_res_name(trim(sline.substr(17,3)));
		atom.set_chain_name(sline.substr(21,1));
		// tee << trim(sline.substr(6,5)) << " " << trim( sline.substr(22,4)) << std::endl;
		atom.set_serial(std::stoi(trim(sline.substr(6,5))));

		const std::string res_serial_str = trim(sline.substr(22,4));
		if(res_serial_str!="") {
			atom.set_res_serial(std::stoi(res_serial_str));
		} else {
			atom.set_res_serial(0);
		}
		// atom.set_context_index(ci);
		return atom;
	}




	inline void read_ligand_pdbqt(const std::string ligand_path, Ligand& ligand, std::ofstream& tee) {
		std::ifstream in_ligand(ligand_path);
		assert(in_ligand);
		std::string sline;

		//put atom into ligand and gather tree information
		int depth_level = 0;
		std::vector<int> width_in_each_depth;
		std::vector<Atom_Tag> atom_tags;
		std::vector<Node_Info> node_infos;
		std::map<std::string,int> node_id_to_node_index;
		while(std::getline(in_ligand,sline)) {
			if(sline=="") continue;
			//to_pdbqt_parse_type do not accept empty string
			PDBQT_PARSE_TYPE ppt = to_pdbqt_parse_type(sline);
			int line_index = ligand.contexts.size();
			if(ppt!=PDBQT_BOND) {
				ligand.contexts.push_back(Context(-1,sline));
			}
			// tee << sline << std::endl;
			switch(ppt) {
				case PDBQT_ATOM:
				case PDBQT_HETATM: {
						const Atom a = parse_pdbqt_atom(sline);
						int a_i = ligand.add_atom(a);
						if(a_i==0) {
							//make sure ligand atom start from 1
							assert(a.get_serial()==1);
						}
						int r_a_i = ligand.add_ref_atom(a);
						assert(a_i == r_a_i);
						ligand.contexts[ligand.contexts.size()-1].atom_index = a_i;
						Node_Id ni = Node_Id(depth_level,width_in_each_depth[depth_level]);
						atom_tags.push_back(Atom_Tag(a_i,ni));
						assert(a_i==atom_tags.size()-1);
					}
					break;
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
					}
					break;
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
				case PDBQT_BOND: {
						tee << "atom size: " << ligand.atoms.size() << " bond: " << sline << std::endl;
						Bond_Info bi = parse_pdbqt_bond_and_add_bond(sline);
						if(bi.bt != NOT_CONNECTED_BOND && bi.bt != NONE_BOND) {
							Bond b1(bi.a2,bi.bt);
							Bond b2(bi.a1,bi.bt);
							assert(bi.a1 < ligand.atoms.size());
							assert(bi.a2 < ligand.atoms.size());
							ligand.atoms[bi.a1].add_bond(b1);
							ligand.atoms[bi.a2].add_bond(b2);
							ligand.ref_atoms[bi.a1].add_bond(b1);
							ligand.ref_atoms[bi.a2].add_bond(b2);
						}
					}
					break;
				case PDBQT_UNKNOWN:
				default:
					assert(false);
					break;
			}
		}
		in_ligand.close();

		assert(ligand.atom_num()==atom_tags.size());

		tee << "Initialize ligand's nodes..." << std::endl;
		//Initialize ligand's tree i.e all nodes, this must be done first since the vector will allocate new memory when push_back new nodes, thus old pointer will fail
		const int max_depth = width_in_each_depth.size();
		tee << "total max depth--->" << max_depth << std::endl;
		for(int depth_index = 0; depth_index < width_in_each_depth.size(); ++depth_index) {
			const int num_of_nodes = width_in_each_depth[depth_index]+1;
			tee << "depth " << depth_index << " has " << num_of_nodes << " nodes" << std::endl;
		}

		for(const Node_Info& ni : node_infos) {
			Node tmp_node = Node();
			tmp_node.set_node_id(ni.node_id.depth_index, ni.node_id.width_index);
			tmp_node.set_parent_bond_atom_index(ni.parent_bond_atom_index);
			tmp_node.set_self_bond_atom_index(ni.self_bond_atom_index);
			ligand.nodes.push_back(tmp_node);
		}

		//assign atom range to each node
		tee << "assign atom range to each node..." << std::endl;
		for(int i = 0; i < atom_tags.size(); ++i) {
			const int depth_index = atom_tags.at(i).node_id.depth_index;
			const int width_index = atom_tags.at(i).node_id.width_index;
			int ni = node_id_to_node_index.at(std::to_string(depth_index)+"-"+std::to_string(width_index));
			ligand.nodes[ni].add_atom_index(atom_tags.at(i).atom_index);
		}

		tee << "Establish parent and children's connection..." << std::endl;
		//Establish parent and children's connection, and parent bond atom index, self bond atom index
		//for root
		ligand.nodes[0].set_parent_pointer(nullptr);
		//for non root
		for(auto ni = 1; ni < node_infos.size(); ++ni) {
			Node_Id parent_node_id = atom_tags.at(node_infos[ni].parent_bond_atom_index).node_id;
			int parent_node_index = node_id_to_node_index.at(std::to_string(parent_node_id.depth_index)+"-"+std::to_string(parent_node_id.width_index));
			ligand.nodes[ni].set_parent_pointer(ligand.nodes[parent_node_index].get_self_pointer());
			ligand.nodes[parent_node_index].add_child_pointer(ligand.nodes[ni].get_self_pointer());
		}

		//for root, set origin and rot_axis, rot_angle to (0,0,0) (0,0,1) 0, because they will reset afterwards
		ligand.nodes[0].set_origin(Vec3d(0,0,0));
		ligand.nodes[0].set_rot_axis(Vec3d(0,0,1));
		ligand.nodes[0].set_rot_angle(0.0);
		// set non root node origin(parent bond atom coord), rot axis(unit vector from parent bond atom to self bond atom) and set all rot angle to 0
		for(int ni = 1; ni < ligand.nodes.size(); ++ni) {
			Node& n = ligand.nodes[ni];
			const int parent_bond_atom_index = n.get_parent_bond_atom_index();
			const int self_bond_atom_index = n.get_self_bond_atom_index();
			const Vec3d origin = ligand.atoms[parent_bond_atom_index].get_coord();
			Vec3d rot_axis = ligand.atoms[self_bond_atom_index].get_coord() - origin;
			rot_axis.normalize();
			n.set_origin(origin);
			n.set_rot_axis(rot_axis);
			n.set_rot_angle(0.0);
		}

		//set ref heavy atom center
		Atoms ref_heavy_atoms;
		for(const auto& ha : ligand.ref_atoms) {
			if(ha.get_element_type() != EL_H) {
				ref_heavy_atoms.push_back(ha);
			}
		}
		ligand.ref_heavy_atoms_center = (1.0/Float(ref_heavy_atoms.size())) * std::accumulate( ref_heavy_atoms.begin(),ref_heavy_atoms.end(),Vec3d(0,0,0),[](const Vec3d& a, const Atom& b) {return a + b.get_coord();} );
		// ligand.ref_heavy_atoms_center = (1.0/Float(ligand.ref_atoms.size())) * std::accumulate( ligand.ref_atoms.begin(),ligand.ref_atoms.end(),Vec3d(0,0,0),[](const Vec3d& a, const Atom& b) {return a + b.get_coord();} );
	}

}