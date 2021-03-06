#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "atom.h"
#include "vec3d.h"
#include "mat3x3.h"
#include "quaternion.h"
#include "context.h"
// #include "conformer.h"

namespace tmd {

using Atom_Ranges = std::vector<int>;

//for non root nodes, Conf defines its rot_axis as the axis along the rotatable bond from parent node to itself, and rot_angle is the angle that will be used to rotate the bond (smapled), the origin is the coordinate of the parent atom of the rotatable bond
//for root nodes, rot_axis and rot_angle are sampled and origin is defined as geometric center
struct Conf {
    Vec3d origin;
    Vec3d rot_axis;
    Float rot_angle = 0;
};

//Each node's Node_Id tells us its position in the ligand torsional tree
struct Node_Id {
    int depth_index = -1;
    int width_index = -1;
    // Node_Id() {}
    Node_Id(const int di, const int wi) : depth_index(di), width_index(wi) {
        assert(depth_index >= 0 && width_index >= 0);
    }
};

using Node_Ids = std::vector<Node_Id>;

//each Node is linked to its parent and children through pointers and the transformation operations performed on one Node will be transmitted to all its children and children's children util reach the end
class Node {
private:
    Node * parent = nullptr;//if it is root Node, parent will be nullptr
    std::vector<Node *> children;

    Node_Id node_id;

    int parent_bond_atom_index = -1;//store parent bond atom, if parent does not exist store the root atom
    int self_bond_atom_index = -1;
    Atom_Ranges atom_ranges;

    Conf conf;

    // transitive_transform will propagate the transformation from current node to all of its children via a branch-wise manner
    void transitive_transform(const Vec3d o, const Vec3d rax, const Float ra, Atoms& atoms) {
        // std::cout << "NodeId[" << this->node_id.depth_index << "," << this->node_id.width_index << "]" << std::endl;
        //rotate node's atoms
        for(const int i : this->atom_ranges) {
            Atom& sa = atoms[i];
            sa.set_coord(o + (sa.get_coord() - o)*Quaternion(rax,ra).to_mat3x3());
        }
        //rotate children's atoms in the same way
        for(int i = 0; i < children.size(); ++i) {
            children[i]->transitive_transform(o,rax,ra,atoms);
        }
    }
public:
    friend class Ligand;
    friend void print(const Node& n, std::ostream& out);
    Node() : parent(nullptr), node_id(0,0) {}
    //move assignment operator//std::sort need this
    // Node& operator=(Node&& n) {
    //     //check if it is itself
    //     if(this!=(&n)) {
    //         //when move a node's content to another,
    //         this->parent = n.parent;
    //         this->children = n.children;



    //     }
    //     return *this;
    // }
    //copy constructor, be carefull about the parent and children pointers
    //this copy constructor can only be used in Tree's copy constructor, because it assgin all pointer to nullptr and does not retain the connnections
    // Node(const Node& n) {
    //     this->parent = nullptr;
    //     this->children = n.children;
    //     for(int c = 0; c < this->children.size(); ++c) {
    //         this->children[c] = nullptr;
    //     }
    //     this->node_id = n.node_id;
    //     this->parent_bond_atom_index = n.parent_bond_atom_index;
    //     this->self_bond_atom_index = n.self_bond_atom_index;
    //     this->atom_ranges = n.atom_ranges;
    //     this->conf = n.conf;
    // }
    // Node(Node *p, std::vector<Node*> c, Atom_Index pbai, Atom_Index sbai) : parent(p), children(c), parent_bond_atom_index(pbai), self_bond_atom_index(sbai) {
    //     this->origin = this->parent_bond_atom->get_coord();
    //     this->rot_axis = (this->bond_atom->get_coord() - this->parent_bond_atom->get_coord());
    //     this->rot_axis.normalize();
    // }
    // Node(const std::vector<Node*> c, const Vec3d o) : parent(nullptr), children(c), parent_bond_atom(nullptr), bond_atom(nullptr) {
    //     this->rot_axis = zero_vec3d;
    //     this->origin = o;
    // }//root node initialization

    void set_origin(const Vec3d& o) {
        this->conf.origin = o;
    }
    void set_rot_axis(const Vec3d& rax) {
        this->conf.rot_axis = rax;
    }
    void set_rot_angle(const Float rag) {
        this->conf.rot_angle = rag;
    }
    void set_parent_bond_atom_index(const int pbai) {
        this->parent_bond_atom_index = pbai;
    }
    void set_self_bond_atom_index(const int sbai) {
        this->self_bond_atom_index = sbai;
    }
    void set_node_id(const int di, const int wi) {
        this->node_id.depth_index = di;
        this->node_id.width_index = wi;
    }
    void set_parent_pointer(Node *const p) {
        this->parent = p;
    }
    void add_child_pointer(Node *const p) {
        this->children.push_back(p);
    }
    void add_atom_index(const int ai) {
        this->atom_ranges.push_back(ai);
    }
    void transform(Atoms& atoms) {
        // root rot_axis is sampled, so it should be check every time
        // but still check for all nodes
        if(!eq(this->conf.rot_axis.norm(),1.0)) {
            std::cout << "rot_axis not normalized! node_id: [" << this->node_id.depth_index << "," << this->node_id.width_index << "]" << std::endl;
            std::cout << "rot axis: " << this->conf.rot_axis[0] << " " << this->conf.rot_axis[1] << " " << this->conf.rot_axis[2] << std::endl;
            assert(false);
        }

        // transform all the nodes that stem from this node in the same way as this node
        if(this->parent != nullptr) {
            //update node conf: origin, rot_axis
            this->conf.origin = atoms[this->parent_bond_atom_index].get_coord();
            this->conf.rot_axis = (atoms[this->self_bond_atom_index].get_coord() - atoms[this->parent_bond_atom_index].get_coord());
            this->conf.rot_axis.normalize();
            //transform
            //rot_angle is sampled from other place
            this->transitive_transform(this->conf.origin, this->conf.rot_axis, this->conf.rot_angle, atoms);
        }

        //transform children
        for(int i = 0; i < this->children.size(); ++i) {
            this->children[i]->transform(atoms);
        }

    }

    Node *const get_self_pointer() {
        return this;
    }
    const int get_parent_bond_atom_index() const {
        return this->parent_bond_atom_index;
    }
    const int get_self_bond_atom_index() const {
        return this->self_bond_atom_index;
    }
    const Node_Id get_node_id() const {
        return this->node_id;
    }
    const Conf get_conf() const {
        return this->conf;
    }
};
using Nodes = std::vector<Node>;

inline const Float rmsd_between_atoms(const Atoms& as1, const Atoms& as2) {
    Float rmsd = 0.0;
    int heavy_atom_size = 0;
    for(int i = 0; i < as1.size(); ++i) {
        const Atom& a2 = as2[i];
        if(a2.get_element_type() == EL_H) continue;
        const Atom& a1 = as1[i];
        ++heavy_atom_size;
        rmsd += (a1.get_coord() - a2.get_coord()).norm_square();
    }
    return std::sqrt(rmsd/Float(heavy_atom_size));
}

enum LIGAND_SAMPLE_MODE {LIGAND_SAMPLE_RIGID,LIGAND_SAMPLE_FLEXIBLE};
class Ligand {
    Atoms ref_atoms;//store the input atoms from file
    // Coords ref_coords;
    Atoms atoms;//actual atoms reflect different sampled conformations
    // Coords coords;
    // Conformers conformers;//store the important conformations
    Nodes nodes;
    // DOFs dofs;//degrees of freedom: translational dof x,y,z (3) + rotational dof rot_axis,rot_angle (3) + ligand rotational dof len(nodes)-1 (3)// root rot_axis has a constraint which is normal
    Contexts contexts;

    Vec3d ref_heavy_atoms_center;//default constructor should be not_a_num

    const int add_atom(const Atom& a) {
        this->atoms.push_back(a);
        // this->coords.push_back(a.get_coord());
        return (this->atoms.size()-1);
    }
    const int add_ref_atom(const Atom& a) {
        this->ref_atoms.push_back(a);
        // this->ref_coords.push_back(a.get_coord());
        return (this->ref_atoms.size()-1);
    }

public:
    // friend class Sampling;
    friend void print(const Ligand& l, std::ostream& out);
    friend void read_ligand_pdbqt(const std::string ligand_path, Ligand& lig, std::ofstream& tee);
    //default constructor because we defined copy constructor
    Ligand() {}
    //move assignment operator//std::sort need this
    // Ligand& operator=(Ligand&& l) : {}
    //copy constructor, need to update Nodes' parent, children pointers and Ligand's dof_pointers
    // Ligand(const Ligand& l) {
    //     this->ref_atoms = l.ref_atoms;
    //     this->atoms = l.atoms;
    //     this->contexts = l.contexts;
    //     // std::cout << "finish copy normal things in Ligand..." << std::endl;
    //     this->nodes = l.nodes;
    //     std::map<std::string,Node&> node_map;

    //     //update this->tree connnection between all nodes through updating parent and children
    //     for(const Node_Id& ti : copied_ligand_node_ids) {
    //         // std::cout << ti.depth_index << " " << ti.width_index << std::endl;
    //         const Node& copied_node = l.tree(ti);
    //         Node& this_node = this->tree(ti);
    //         // std::cout << "access current int this tree and copied tree" << std::endl;
    //         //first set all parent pointers, root is set to nullptr
    //         if(copied_node.parent!=nullptr) {
    //             const Node_Id copied_node_parent_id = (*(copied_node.parent)).get_node_id();
    //             // std::cout << "get copied_node_parent_id" << std::endl;
    //             this_node.parent = this->tree(copied_node_parent_id).get_self_pointer();
    //         }
    //         //then set all children's pointers
    //         for(Node_Index i = 0; i < copied_node.children.size(); ++i) {
    //             const Node_Id copied_node_child_id = (*(copied_node.children[i])).get_node_id();
    //             this_node.children[i] = this->tree(copied_node_child_id).get_self_pointer();
    //         }
    //     }
    //     //for root
    //     this->nodes[0].parent = nullptr;
    //     // std::cout << "finish update Nodes' connections in this tree..." << std::endl;
    // }
    const int atom_num() const {
        return this->atoms.size();
    }
    const int heavy_atom_num() const {
        int heavy_num = 0;
        for(const Atom& a : atoms) {
            if(a.get_element_type()!=EL_H) {
                ++heavy_num;
            }
        }
        return heavy_num;
    }
    const Atoms get_atoms() const {
        return this->atoms;
    }
    const Atoms& get_atoms_reference() const {
        return this->atoms;
    }
    const Vec3d get_center_coord() const {
        const Vec3d sum = std::accumulate(this->atoms.begin(),this->atoms.end(),Vec3d(0,0,0),[](const Vec3d& a, const Atom& b){return a + b.get_coord();});
        assert(this->atoms.size()!=0);
        return sum*(1.0/Float(this->atoms.size()));
    }
    const Vec3d get_ref_center_coord() const {
        const Vec3d sum = std::accumulate(this->ref_atoms.begin(),this->ref_atoms.end(),Vec3d(0,0,0),[](const Vec3d& a, const Atom& b){return a + b.get_coord();});
        assert(this->ref_atoms.size()!=0);
        return sum*(1.0/Float(this->ref_atoms.size()));
    }
    const Contexts& get_contexts() const {
        return this->contexts;
    }

    const Vec3d& get_ref_heavy_center_coord() const {
        return this->ref_heavy_atoms_center;
    }

    void translate(const Vec3d& shift) {
        //first, ligand translation
        for(auto& a : this->atoms) {
            a.set_coord(a.get_coord()+shift);
        }
    }

    void rotate(const Vec3d& rot_origin, const Vec3d& rot_axis, const Float& rot_angle) {

        this->nodes[0].conf.rot_angle = rot_angle;
        this->nodes[0].conf.rot_axis = rot_axis;

        // this->nodes[0].conf.origin = std::accumulate(this->atoms.begin(),this->atoms.end(),Vec3d(0,0,0),[](const Vec3d& a, const Atom& b){return a + b.get_coord();})*(1.0/Float(this->atoms.size()));
        this->nodes[0].conf.origin = rot_origin;
        this->nodes[0].transitive_transform(this->nodes[0].conf.origin, this->nodes[0].conf.rot_axis, this->nodes[0].conf.rot_angle, this->atoms);
    }

    void transform(const Floats& torsional_dofs) {
        #ifdef DEBUG
        //make sure dofs and nodes has the same number
        assert((this->nodes.size()-1) == torsional_dofs.size());
        #endif

        //set node conf
        for(int i = 0; i < torsional_dofs.size(); ++i) {
            this->nodes[i+1].conf.rot_angle = torsional_dofs[i];
        }

        //second, ligand rotation and branch rotation begins, nodes[0] is the root node
        nodes[0].transform(this->atoms);

    }

    void reset_atoms_coord_to_ref_atoms_coord() {
        //reassign the inital coord to atoms
        for(int i = 0; i < this->atoms.size(); ++i) {
            this->atoms[i].set_coord(this->ref_atoms[i].get_coord());
        }
    }

    void sample(const Floats& dofs, const LIGAND_SAMPLE_MODE& l_mode) {
        this->reset_atoms_coord_to_ref_atoms_coord();
        #ifdef DEBUG
        assert(dofs.size() >= 6);
        #endif
        const Vec3d& shift = Vec3d(dofs[0],dofs[1],dofs[2]);
        this->translate(shift);

        const Vec3d& rot_origin = this->ref_heavy_atoms_center + shift;
        const Vec3d rot_vector(dofs[3],dofs[4],dofs[5]);
        const Float& rot_angle = rot_vector.norm();
        if(eq(std::abs(rot_angle),0.0)) {
            const Vec3d rot_axis(0,0,1);
            this->rotate(rot_origin, rot_axis, 0.0);
        } else {
            const Vec3d rot_axis = rot_vector * (1.0/rot_angle);
            // #ifdef DEBUG
            // if(!eq(rot_axis.norm(),1.0)) {
            //     std::cout << rot_angle << " " << rot_axis[0] << " " << rot_axis[1] << " " << rot_axis[2] << " " << dofs[3] << " " << dofs[4] << " " << dofs[5] << std::endl;
            //     for(int i = 0; i < dofs.size(); ++i) {
            //         std::cout << dofs[i] << " ";
            //     }
            //     std::cout << std::endl;
            //     assert(false);
            // }
            // assert(eq(rot_axis.norm(),1.0));
            // #endif
            this->rotate(rot_origin, rot_axis, rot_angle);
        }

        if(l_mode == LIGAND_SAMPLE_FLEXIBLE) {
            #ifdef DEBUG
            assert(dofs.size() == this->dof_num());
            #endif
            Floats torsional_dofs(dofs.size()-6,0.0);
            for(int i = 6; i < dofs.size(); ++i) {
                torsional_dofs[i-6] = dofs[i];
            }
            this->transform(torsional_dofs);
        }
    }

    const int dof_num() const {
        return (6+this->nodes.size()-1);//3 translation 3 overall rotation and nodes().size()-1 torsion
    }

    const Float rmsd_with_respect_to_ref_atoms() const {
        return rmsd_between_atoms(this->ref_atoms, this->atoms);
    }

    void write(std::ostream& out = std::cout) const {
        write_contexts(this->contexts,this->atoms,out);
    }

    void write_ref(std::ostream& out = std::cout) const {
        write_contexts(this->contexts,this->ref_atoms,out);
    }
};


class RNA {
    Atoms ref_atoms;
    Atoms atoms;
    // Coords ref_coords;
    // Coords coords;
    Contexts contexts;

    const int add_atom(const Atom& a) {
        this->atoms.push_back(a);
        // this->coords.push_back(a.get_coord());
        return (this->atoms.size()-1);
    }
    const int add_ref_atom(const Atom& a) {
        this->ref_atoms.push_back(a);
        // this->ref_coords.push_back(a.get_coord());
        return (this->ref_atoms.size()-1);
    }

public:
    friend void print(const RNA& r, std::ostream& out);
    friend void read_rna_mol2(const std::string rna_path, RNA& r, std::ostream& tee);
    const int atom_num() const {
        return this->atoms.size();
    }
    const Atoms& get_atoms_reference() const {
        return this->atoms;
    }
};


}
