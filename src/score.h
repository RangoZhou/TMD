#pragma once

#include <vector>
#include <set>
#include <thread>
#include <numeric>
#include <mutex>
#include <cmath>
#include <unordered_map>

#include "atom.h"
#include "matrix.h"
#include "molecule.h"
// #include "progress_bar.h"

namespace tmd {

// class RNA;
// class Ligand;
// class Matrix;

struct Grid {
    int x;
    int y;
    int z;
    std::vector<Atom_Index> rigid_atom_mapping;
    Grid(const int xx, const int yy, const int zz, const Atom_Index ai) : x(xx), y(yy), z(zz) {
        this->rigid_atom_mapping = std::vector<Atom_Index>(1,ai);
    }
};
using Grid_Map = std::unordered_map<std::string,Grid>;
using Grid_Index = Grid_Map::size_type;
struct Grids {
    Float width = 2.0;
    Grid_Map grid_map;
};

struct Interacting_Pair {
    Atom_Index i;
    Atom_Index j;
    Float vdw_well_depth;
    Float vdw_characteristic_dis;
    Float vdw_ACOEFF;
    Float vdw_BCOEFF;
    Interacting_Pair() : i(-1), j(-1), vdw_well_depth(k_not_a_num), vdw_characteristic_dis(k_not_a_num), vdw_ACOEFF(k_not_a_num), vdw_BCOEFF(k_not_a_num) {}
    Interacting_Pair(const Atom_Index ii, const Atom_Index jj, const Float vwd, const Float vcd) : i(ii), j(jj), vdw_well_depth(vwd), vdw_characteristic_dis(vcd) {
        this->vdw_ACOEFF = vwd*int_pow<12>(vcd);
        this->vdw_BCOEFF = 2.0*vwd*int_pow<6>(vcd);
    }
};

using Interacting_Pairs = std::vector<Interacting_Pair>;

struct thread_range {
    unsigned int begin;
    unsigned int end;
};

class VDW_Score {
    const Float lower_cutoff = 0.01;
    Float ligand_internal_vdw = k_not_a_num;
public:
    void init(const RNA& rna, const Ligand& lig, const Interacting_Pairs& ips, std::ostream& log) {
        //subtract the starting vdw interaction energy
        Float score = 0.0;
        const Atoms& lig_atoms = lig.get_atoms_reference();
        for(const Interacting_Pair& ip : ips) {
            const Float dis = (lig_atoms[ip.i].get_coord()-lig_atoms[ip.j].get_coord()).norm();
            if(dis < this->lower_cutoff) {
                std::cout << "dis " << dis << " atom_index: " << ip.i << " " << ip.j << std::endl;
                // std::cout << dis << " " << lig_atoms[ip.i].get_name() << " " << lig_atoms[ip.j].get_name() << std::endl;
                assert(false);
                score = k_max_float;
            }
            const Float dis_6 = int_pow<6>(dis);
            score += ip.vdw_ACOEFF/dis_6*dis_6 - ip.vdw_BCOEFF/dis_6;
        }
        this->ligand_internal_vdw = score;
    }
    const Float evaluate(const std::vector<thread_range>& tr, const RNA& rna, const Ligand& lig, const Grids& grids) const;
    const Float evaluate(const std::vector<thread_range>& tr, const Ligand& lig, const Interacting_Pairs& ips) const;
};

class YW_Score {
public:
    int contact_number = 0;
    std::set<std::string> contact_type;
    const Float rmax = 10.0;
    const Float dr = 0.2;
    const int usize = int(rmax/dr);
    // Float umax = 5.06356;//in the unit of kBT = 3 kcal/mol ;1kcal/mol = 1.68787kBT
    // Float gmin = exp(-umax);
    // Float kcal2kbT = 1.68787;
    std::map<std::string,std::vector<Float>> u;

    void init(const RNA& rna, const Ligand& lig, std::ostream& log) {
        this->contact_type = {
        "C.3-C.ar",
        "C.3-O.3",
        "C.2-C.3",
        "C.3-C.3",
        "C.3-N.ar",
        "C.ar-O.3",
        "C.2-C.ar",
        "C.ar-N.ar",
        "C.ar-C.ar",
        "C.3-N.pl3",
        "C.2-O.3",
        "C.3-O.co2",
        "N.ar-O.3",
        "C.3-O.2",
        "C.2-N.ar",
        "C.3-N.am",
        "C.ar-N.pl3",
        "C.ar-O.2",
        "C.ar-N.am",
        "O.3-O.3",
        "C.ar-O.co2",
        "N.pl3-O.3",
        "O.3-O.co2",
        "C.2-N.pl3",
        "C.3-N.2",
        "N.ar-N.pl3",
        "O.2-O.3",
        "C.3-P.3",
        "N.ar-N.ar",
        "N.ar-O.2",
        "C.3-N.4",
        "C.ar-N.2",
        "N.am-O.3",
        "C.2-C.2",
        "N.am-N.ar",
        "C.2-O.co2",
        "C.2-N.4",
        "N.ar-O.co2",
        "C.2-O.2",
        "N.4-O.3",
        "C.ar-N.4",
        "C.ar-P.3",
        "C.2-N.am",
        "N.2-O.3",
        "N.pl3-O.2",
        "N.pl3-O.co2",
        "O.3-P.3",
        "N.am-N.pl3",
        "N.4-N.ar",
        "N.2-N.ar",
        "C.2-N.2",
        "N.4-O.co2",
        "O.2-O.co2",
        "N.pl3-N.pl3",
        "N.4-O.2",
        "N.am-O.co2",
        "N.4-N.am",
        "N.ar-P.3",
        "N.2-N.pl3",
        "C.2-P.3",
        "N.4-N.pl3",
        "O.co2-O.co2",
        "N.am-O.2",
        "N.pl3-P.3",
        "O.co2-P.3",
        "N.2-O.co2",
        "N.2-N.4",
        "N.4-P.3",
        "N.2-O.2",
        "O.2-O.2",
        "C.3-S.3",
        "N.2-N.am",
        "O.2-P.3",
        "O.3-S.3",
        "N.am-N.am",
        "C.2-S.3",
        "C.2-C.cat",
        "C.3-C.cat",
        "N.am-P.3",
        "C.ar-S.3",
        "C.3-N.3",
        "C.2-N.3",
        "C.ar-N.3",
        "C.ar-C.cat",
        "C.3-Ha",
        "C.cat-O.3",
        "C.ar-Ha",
        "N.ar-S.3",
        "N.3-O.3",
        "N.3-N.ar",
        "C.cat-N.ar",
        "C.3-S.2",
        "Ha-O.3",
        "C.2-Ha",
        "N.2-P.3"
        };
        this->contact_number = this->contact_type.size();
        log << "number of contact " << this->contact_number << std::endl;

        ////////////////////////////////////////////////////
        ////////read u_now.list //////////////////////////// u_now.list record the statistical potential
        ////////////////////////////////////////////////////
        std::string sl;
        std::ifstream u_nowfile("/home/yuanzhe/TMD/src/u_now_VI.list",std::ios::in);
        assert(u_nowfile);
        while(getline(u_nowfile,sl)) {
            std::istringstream ss(sl);
            std::string buf;
            std::vector<std::string> token;
            while(ss >> buf) token.push_back(buf);
            std::string name=token[0];
            std::vector<Float> u_now_tmp;
            for(int j=1;j!=token.size();j++) {
                u_now_tmp.push_back(std::stod(token[j]));
            }
            this->u.insert(std::pair<std::string, std::vector<Float> >(name,u_now_tmp));
        }
        log << " u size "<< this->u.size() << std::endl;
        // for(const auto& u_pair : this->u) {
        //     log << u_pair.first << " " << u_pair.second.size() << std::endl;
        // }
    }

    const Float evaluate(const std::vector<thread_range>& tr, const RNA& rna, const Ligand& lig, const Grids& grids) const;
};

enum SCORE_TYPE {YW_SCORE,VDW_LIGAND,VDW_RNA_LIGAND,ALL};
class Scoring_Function {
    friend YW_Score;
    friend VDW_Score;
    const RNA& rna;
    const Ligand& lig;
    Grids grids;//stores the RNA atom mapping with the grids
    Interacting_Pairs ligand_pairs;
    Interacting_Pairs rna_pairs;
    // Matrix<Float> ligand_rna_dis_matrix;
    SCORE_TYPE Score_Type = YW_SCORE;
    // Matrix<Float> ligand_rna_dis6_matrix;
    // Matrix<Float> ligand_rna_dis12_matrix;
    YW_Score yw_score;
    VDW_Score vdw_score;
    unsigned int num_thread = 0;
    std::vector<thread_range> thread_ranges;
public:
    Scoring_Function(const RNA& r, const Ligand& l, std::ostream& log);//initial thread information, initialize distance matrix, and init score function parameters, and initialize mapping grids
    void set_score_type(const SCORE_TYPE st ) {
        this->Score_Type = st;
    }

    const Float evaluate() {

        switch(this->Score_Type) {
            case YW_SCORE:
                return this->yw_score.evaluate(this->thread_ranges,this->rna,this->lig,this->grids);
                break;
            case VDW_LIGAND:
                return this->vdw_score.evaluate(this->thread_ranges,this->lig,this->ligand_pairs);
                break;
            case VDW_RNA_LIGAND:
                return this->vdw_score.evaluate(this->thread_ranges,this->rna,this->lig,this->grids);
                break;
            case ALL: {
                Float score = 0.0;
                score += this->yw_score.evaluate(this->thread_ranges,this->rna,this->lig,this->grids);
                score += this->vdw_score.evaluate(this->thread_ranges,this->lig,this->ligand_pairs);
                score += this->vdw_score.evaluate(this->thread_ranges,this->rna,this->lig,this->grids);
                return score;
                break;
            }
            default:
                assert(false);
                break;
        }
    }
};

}
