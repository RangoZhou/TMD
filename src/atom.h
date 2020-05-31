#pragma once

#include <vector>
#include <cassert>
#include <string>
#include "atom_constant.h"
#include "vec3d.h"

namespace tmd {

class Atom_Type {
    ELEMENT_TYPE el_type = EL_Unk;
    AD4_TYPE ad4_type = AD4_Unk;
    SYBYL_TYPE sybyl_type = SYBYL_Unk;
public:
    // explicit Atom_Type(const AD4_TYPE& at, const SYBYL_TYPE& st) : ad4_type(at), sybyl_type(st) {
    //     this->el_type = ad4_infos.at(this->ad4_type).element_type;//ad4_type_to_element_type
    //     assert(this->sybyl_type != SYBYL_Unk && this->ad4_type != AD4_Unk);
    //     assert(this->el_type == sybyl_infos.at(this->sybyl_type).element_type);//sybyl_type_to_element_type
    // }

    // explicit Atom_Type(const SYBYL_TYPE& st) : sybyl_type(st) {
    //     this->el_type = sybyl_infos.at(this->sybyl_type).element_type;
    // }

    void set_element_type(const ELEMENT_TYPE& et) {
        this->el_type = et;
    }
    void set_sybyl_type(const SYBYL_TYPE& st) {
        this->sybyl_type = st;
    }
    void set_ad4_type(const AD4_TYPE& at) {
        this->ad4_type = at;
    }
    //Atom constants
    // element parameter
    const Float& get_covalent_radius() const {
        return element_infos.at(this->el_type).covalent_radius;
    }
    const Float& get_atomic_radius() const {
        return element_infos.at(this->el_type).atomic_radius;
    }
    const Float& get_vdw_radius() const {
        return element_infos.at(this->el_type).vdw_radius;
    }
    const int& get_atomic_number() const {
        return element_infos.at(this->el_type).atomic_number;
    }
    const std::string& get_element_type_name() const {
        return element_infos.at(this->el_type).name;
    }
    const std::string& get_sybyl_type_name() const {
        return sybyl_infos.at(this->sybyl_type).sybyl_name;
    }
    const std::string& get_ad4_type_name() const {
        return ad4_infos.at(this->ad4_type).ad4_name;
    }
    const ELEMENT_TYPE& get_element_type() const {
        return this->el_type;
    }
    const SYBYL_TYPE& get_sybyl_type() const {
        return this->sybyl_type;
    }
    const AD4_TYPE& get_ad4_type() const {
        return this->ad4_type;
    }

    //Atom chemical physical properties
    const bool is_metal() const {
        switch(this->el_type) {
            case EL_Ac:
            case EL_Ag:
            case EL_Al:
            case EL_Am:
            case EL_Au:
            case EL_Ba:
            case EL_Be:
            case EL_Bi:
            case EL_Bk:
            case EL_Ca:
            case EL_Cd:
            case EL_Ce:
            case EL_Cf:
            case EL_Cm:
            case EL_Co:
            case EL_Cr:
            case EL_Cs:
            case EL_Cu:
            case EL_Db:
            case EL_Dy:
            case EL_Er:
            case EL_Es:
            case EL_Eu:
            case EL_Fe:
            case EL_Fm:
            case EL_Fr:
            case EL_Ga:
            case EL_Gd:
            case EL_Ge:
            case EL_Hf:
            case EL_Hg:
            case EL_Ho:
            case EL_In:
            case EL_Ir:
            case  EL_K:
            case EL_La:
            case EL_Li:
            case EL_Lr:
            case EL_Lu:
            case EL_Md:
            case EL_Mg:
            case EL_Mn:
            case EL_Mo:
            case EL_Na:
            case EL_Nb:
            case EL_Nd:
            case EL_Ni:
            case EL_No:
            case EL_Np:
            case EL_Os:
            case EL_Pa:
            case EL_Pb:
            case EL_Pd:
            case EL_Pm:
            case EL_Po:
            case EL_Pr:
            case EL_Pt:
            case EL_Pu:
            case EL_Ra:
            case EL_Rb:
            case EL_Re:
            case EL_Rf:
            case EL_Rh:
            case EL_Ru:
            case EL_Sb:
            case EL_Sc:
            case EL_Sg:
            case EL_Sm:
            case EL_Sn:
            case EL_Sr:
            case EL_Ta:
            case EL_Tb:
            case EL_Tc:
            case EL_Th:
            case EL_Ti:
            case EL_Tl:
            case EL_Tm:
            case  EL_U:
            case  EL_V:
            case  EL_W:
            case  EL_Y:
            case EL_Yb:
            case EL_Zn:
            case EL_Zr:
                return true;
                break;
            default:
                return false;
                break;
        }
    }

    const bool is_hydrogen() const {
        return (this->el_type==EL_H);
    }
    // const bool is_acceptor() const {}
    // const bool is_donor() const {}
    // const bool is_hydrophobic() const {}
};

enum BOND_TYPE {SINGLE_BOND,DOUBLE_BOND,TRIPLE_BOND,AMIDE_BOND,AROMATIC_BOND,DUMMY_BOND,UNKNOWN_BOND,NOT_CONNECTED_BOND,NONE_BOND};
class Bond {
    int bond_atom_index;
    BOND_TYPE Bond_Type;
public:
    // Bond() : to_index(-1) {}
    explicit Bond(const int ai, const BOND_TYPE bt) : bond_atom_index(ai), Bond_Type(bt) {}
    const int& get_bonded_atom_index() const {
        return this->bond_atom_index;
    }
    const BOND_TYPE& get_bond_type() const {
        return this->Bond_Type;
    }
};
using Bonds = std::vector<Bond>;

class Atom_Base : public Atom_Type {
    Vec3d coord = k_nan_vec3d;
    Bonds bonds;
    Float charge = 0.0;
public:
    //Atom class properties
    void set_charge(const Float& c) {
        this->charge = c;
    }
    void set_coord(const Vec3d& cd) {
        this->coord = cd;
    }
    void add_bond(const Bond& b) {
        this->bonds.push_back(b);
    }
    const Float& get_charge() const {
        return this->charge;
    }
    const Vec3d& get_coord() const {
        return this->coord;
    }
    const Bonds& get_bonds() const {
        return this->bonds;
    }
};

class Atom : public Atom_Base {
    int serial;
    std::string name = "";
    int res_serial;
    std::string res_name = "";
    std::string chain_name = "";
public:
    void set_name(const std::string& n) {
        this->name = n;
    }
    void set_res_name(const std::string& rn) {
        this->res_name = rn;
    }
    void set_chain_name(const std::string& cn) {
        this->chain_name = cn;
    }
    void set_serial(const int& s) {
        this->serial = s;
    }
    void set_res_serial(const int& rs) {
        this->res_serial = rs;
    }
    const std::string& get_name() const {
        return this->name;
    }
    const std::string& get_res_name() const {
        return this->res_name;
    }
    const std::string& get_chain_name() const {
        return this->chain_name;
    }
    const int& get_serial() const {
        return this->serial;
    }
    const int& get_res_serial() const {
        return this->res_serial;
    }
};

using Atoms = std::vector<Atom>;

}
