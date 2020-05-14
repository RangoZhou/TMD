#pragma once

#include <vector>
#include <cassert>
#include <string>
#include "atom_constant.h"
#include "vec3d.h"

namespace tmd {

class Atom;//forward declaration

// using Ref_Atom = std::reference_wrapper<Atom>;
using Atom_Index = std::vector<Atom>::size_type;
using Atoms = std::vector<Atom>;
// using Ref_Atoms = std::vector<Atom>;
// using Input_Atoms std::vector<Atom>;
//Context stores the atom line and its position in the file content while reading
struct Context;//forward declaration
using Contexts = std::vector<Context>;
using Context_Index = Contexts::size_type;
struct Context {
    Context_Index index = -1;
    std::string c = "";
    Context() {}
    Context(const Context_Index ci, const std::string ct) : index(ci), c(ct) {}
};


enum ATOM_TYPE {SYBYL,AD4,SYBYL_AD4,AD4_SYBYL};

class Bond {
    Atom_Index to_index;
public:
    Bond() : to_index(-1) {}
    Bond(const Atom_Index ai) : to_index(ai) {}
    const Atom_Index get_bonded_atom_index() const {
        return this->to_index;
    }
};
using Bonds = std::vector<Bond>;

class Atom {
    ELEMENT_TYPE el_type = EL_UNK;
    AD4_TYPE ad4_type = AD4_UNK;
    SYBYL_TYPE sybyl_type = SYBYL_UNK;
    Vec3d coord;
    Bonds bonds;
    Atom_Index serial = -1;
    std::string name = "";
    std::string res_name = "";
    std::string chain_name = "";
    Atom_Index res_serial = -1;
    Context context;
    Float charge = 0.0;
    // Float occupancy = 0.0;
    // Float temp_factor = 0.0;
    // unsigned int model_serial = 1;
public:
    friend void print(const Atom& a, std::ostream& out);
    //constructors
    // Atom() {}
    explicit Atom(const std::string s, const ATOM_TYPE at) {
        switch(at) {
            case SYBYL:
                this->sybyl_type = (sybyl_type_lookup.find(s)==sybyl_type_lookup.end()) ? SYBYL_UNK : sybyl_type_lookup.at(s);//string_to_sybyl_type
                this->el_type = sybyl_infos.at(this->sybyl_type).EL_TYPE;//sybyl_type_to_element_type
                break;
            case AD4:
                this->ad4_type = (ad4_type_lookup.find(s)==ad4_type_lookup.end()) ? AD4_UNK : ad4_type_lookup.at(s);//string_to_ad4_type
                this->el_type = ad4_infos.at(this->ad4_type).EL_TYPE;//ad4_type_to_element_type
                break;
            default:
                assert(false);
                break;
        }
    }
    explicit Atom(const std::string s1, const std::string s2, const ATOM_TYPE at) {
        switch(at) {
            case AD4_SYBYL:
                this->ad4_type = (ad4_type_lookup.find(s1)==ad4_type_lookup.end()) ? AD4_UNK : ad4_type_lookup.at(s1);//string_to_ad4_type
                this->el_type = ad4_infos.at(this->ad4_type).EL_TYPE;//ad4_type_to_element_type
                this->sybyl_type = (sybyl_type_lookup.find(s2)==sybyl_type_lookup.end()) ? SYBYL_UNK : sybyl_type_lookup.at(s2);//string_to_sybyl_type
                assert(this->el_type == sybyl_infos.at(this->sybyl_type).EL_TYPE);//sybyl_type_to_element_type
                break;
            case SYBYL_AD4:
                this->sybyl_type = (sybyl_type_lookup.find(s1)==sybyl_type_lookup.end()) ? SYBYL_UNK : sybyl_type_lookup.at(s1);//string_to_sybyl_type
                this->el_type = sybyl_infos.at(this->sybyl_type).EL_TYPE;//sybyl_type_to_element_type
                this->ad4_type = (ad4_type_lookup.find(s2)==ad4_type_lookup.end()) ? AD4_UNK : ad4_type_lookup.at(s2);//string_to_ad4_type
                assert(this->el_type == ad4_infos.at(this->ad4_type).EL_TYPE);//ad4_type_to_element_type
                break;
            default:
                assert(false);
                break;
        }
    }
    // Atom(const SYBYL_TYPE st) : sybyl_type(st) {
    //     this->el_type = sybyl_infos.at(this->sybyl_type).EL_TYPE;//sybyl_type_to_element_type
    // }

    //Atom constants
    // element parameter
    const Float get_element_covalent_radius() const {
        return element_infos.at(this->el_type).covalent_radius;
    }
    // const Float get_atomic_radius() const;
    const Float get_element_vdw_radius() const {
        return element_infos.at(this->el_type).vdw_radius;
    }
    const std::string get_element_type_name() const {
        return element_infos.at(this->el_type).name;
    }

    // sybyl atom parameter
    const Float get_sybyl_covalent_radius() const {
        return sybyl_infos.at(this->sybyl_type).covalent_radius;
    }
    // const Float get_atomic_radius() const;
    const Float get_sybyl_vdw_radius() const {
        return sybyl_infos.at(this->sybyl_type).vdw_radius;
    }
    const Float get_sybyl_well_depth() const {
        return sybyl_infos.at(this->sybyl_type).well_depth;
    }
    const Float get_sybyl_vdw_volume() const {
        return sybyl_infos.at(this->sybyl_type).vdw_volume;
    }
    const Float get_sybyl_solvation() const {
        return sybyl_infos.at(this->sybyl_type).solvation;
    }
    const std::string get_sybyl_type_name() const {
        return sybyl_infos.at(this->sybyl_type).name;
    }

    //Atom class properties
    void set_charge(const Float c) {
        this->charge = c;
    }
    void set_element_type(const ELEMENT_TYPE et) {
        this->el_type = et;
    }
    void set_sybyl_type(const SYBYL_TYPE st) {
        this->sybyl_type = st;
    }
    void set_ad4_type(const AD4_TYPE at) {
        this->ad4_type = at;
    }
    void set_coord(const Vec3d cd) {
        this->coord = cd;
    }
    void set_name(const std::string n) {
        this->name = n;
    }
    void set_res_name(const std::string rn) {
        this->res_name = rn;
    }
    void set_chain_name(const std::string cn) {
        this->chain_name = cn;
    }
    void set_serial(const Atom_Index s) {
        this->serial = s;
    }
    void set_res_serial(const Atom_Index rs) {
        this->res_serial = rs;
    }
    void set_context(const Context_Index ci, const std::string c) {
        this->context.index = ci;
        this->context.c = c;
    }

    void add_bond(const Bond b) {
        this->bonds.push_back(b);
    }

    const Float get_charge() const {
        return this->charge;
    }
    const ELEMENT_TYPE get_element_type() const {
        return this->el_type;
    }
    const SYBYL_TYPE get_sybyl_type() const {
        return this->sybyl_type;
    }
    const AD4_TYPE get_ad4_type() const {
        return this->ad4_type;
    }
    const std::string get_name() const {
        return this->name;
    }
    const std::string get_res_name() const {
        return this->res_name;
    }
    const std::string get_chain_name() const {
        return this->chain_name;
    }
    const Atom_Index get_serial() const {
        return this->serial;
    }
    const Atom_Index get_res_serial() const {
        return this->res_serial;
    }
    const Context get_context() const {
        return this->context;
    }
    const Vec3d get_coord() const {
        return this->coord;
    }
    const Bonds get_bonds() const {
        return this->bonds;
    }

    //Atom chemical physical properties
    const bool is_metal() const {
        switch(this->el_type) {
            case EL_AC:
            case EL_AG:
            case EL_AL:
            case EL_AM:
            case EL_AU:
            case EL_BA:
            case EL_BE:
            case EL_BI:
            case EL_BK:
            case EL_CA:
            case EL_CD:
            case EL_CE:
            case EL_CF:
            case EL_CM:
            case EL_CO:
            case EL_CR:
            case EL_CS:
            case EL_CU:
            case EL_DB:
            case EL_DY:
            case EL_ER:
            case EL_ES:
            case EL_EU:
            case EL_FE:
            case EL_FM:
            case EL_FR:
            case EL_GA:
            case EL_GD:
            case EL_GE:
            case EL_HF:
            case EL_HG:
            case EL_HO:
            case EL_IN:
            case EL_IR:
            case  EL_K:
            case EL_LA:
            case EL_LI:
            case EL_LR:
            case EL_LU:
            case EL_MD:
            case EL_MG:
            case EL_MN:
            case EL_MO:
            case EL_NA:
            case EL_NB:
            case EL_ND:
            case EL_NI:
            case EL_NO:
            case EL_NP:
            case EL_OS:
            case EL_PA:
            case EL_PB:
            case EL_PD:
            case EL_PM:
            case EL_PO:
            case EL_PR:
            case EL_PT:
            case EL_PU:
            case EL_RA:
            case EL_RB:
            case EL_RE:
            case EL_RF:
            case EL_RH:
            case EL_RU:
            case EL_SB:
            case EL_SC:
            case EL_SG:
            case EL_SM:
            case EL_SN:
            case EL_SR:
            case EL_TA:
            case EL_TB:
            case EL_TC:
            case EL_TH:
            case EL_TI:
            case EL_TL:
            case EL_TM:
            case  EL_U:
            case  EL_V:
            case  EL_W:
            case  EL_Y:
            case EL_YB:
            case EL_ZN:
            case EL_ZR:
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

}
