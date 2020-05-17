#pragma once

#include <iostream>
#include <sstream>
#include <istream>
#include <ostream>
#include <fstream>
#include <map>
#include <string>
#include <iomanip>
#include <algorithm>
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


// struct Interacting_Pair {
//     const Atom_Index i;
//     const Atom_Index j;
//     Interacting_Pair(const Atom_Index ii, const Atom_Index jj) : i(ii), j(jj) {}
// };
// using Interacting_Pairs = std::vector<Interacting_Pair>;
// //precaluate interacting pairs
// Atom_Index num_of_ligand_interacing_pairs = 0;
// const Atom_Index num_lig_atom = lig.atom_num();
// const Atoms& lig_atoms = lig.get_atoms_reference();
// this->ligand_pairs.resize((num_lig_atom*(num_lig_atom-1))/2);
// for(Atom_Index i = 0; i < num_lig_atom; ++i) {
//     for(Atom_Index j = i+1; j < num_lig_atom; ++j) {
//         ++num_of_ligand_interacing_pairs;
//         const Float rij = lig_atoms[i].get_sybyl_vdw_radius() + lig_atoms[j].get_sybyl_vdw_radius();
//         const Float eij = std::sqrt(lig_atoms[i].get_sybyl_well_depth()*lig_atoms[j].get_sybyl_well_depth());
//         this->ligand_pairs[num_of_ligand_interacing_pairs-1] = std::move(Interacting_Pair(i,j,eij,rij));
//     }
// }


// class Coefficient_Matrix {

//     Triangular_Matrix<Float> vdw_ACOEFF_tri_matrix;
//     Triangular_Matrix<Float> vdw_BCOEFF_tri_matrix;
//     Triangular_Matrix<Float> vdw_characteristic_dis_tri_matrix;
//     Triangular_Matrix<Float> vdw_well_depth_tri_matrix;

// public:
//     Coefficient_Matrix(const RNA& r, const Ligand& l) : vdw_ACOEFF_tri_matrix(r.atom_num()+l.atom_num(),0.0), vdw_BCOEFF_tri_matrix(r.atom_num()+l.atom_num(),0.0), vdw_characteristic_dis_tri_matrix(r.atom_num()+l.atom_num(),0.0), vdw_well_depth_tri_matrix(r.atom_num()+l.atom_num(),0.0) {
//         Atom_Index rna_atom_num = r.atom_num();
//         Atom_Index lig_atom_num = l.atom_num();
//         for(Atom_Index i = 0; i < rna_atom_num; ++i) {
//             for(Atom_Index j = 0; j < rna_atom_num; ++j) {
//                 if(i>j) continue;
//                 const Atom& ai = r.get_atoms_reference()[i];
//                 const Atom& aj = r.get_atoms_reference()[j];
//                 this->vdw_well_depth_tri_matrix(i,j) = std::sqrt(ai.get_sybyl_well_depth() * aj.get_sybyl_well_depth());
//                 this->vdw_characteristic_dis_tri_matrix(i,j) = ai.get_sybyl_vdw_radius() + aj.get_sybyl_vdw_radius();
//                 this->vdw_ACOEFF_tri_matrix(i,j) = this->vdw_well_depth_tri_matrix(i,j)*int_pow<12>(this->vdw_characteristic_dis_tri_matrix(i,j));
//                 this->vdw_BCOEFF_tri_matrix(i,j) = 2.0*this->vdw_well_depth_tri_matrix(i,j)*int_pow<6>(this->vdw_characteristic_dis_tri_matrix(i,j));
//             }
//         }

//         for(Atom_Index i = rna_atom_num; i < lig_atom_num+rna_atom_num; ++i) {
//             for(Atom_Index j = rna_atom_num; j < lig_atom_num+rna_atom_num; ++j) {
//                 if(i>j) continue;
//                 const Atom& ai = l.get_atoms_reference()[i-rna_atom_num];
//                 const Atom& aj = l.get_atoms_reference()[j-rna_atom_num];
//                 this->vdw_well_depth_tri_matrix(i,j) = std::sqrt(ai.get_sybyl_well_depth() * aj.get_sybyl_well_depth());
//                 this->vdw_characteristic_dis_tri_matrix(i,j) = ai.get_sybyl_vdw_radius() + aj.get_sybyl_vdw_radius();
//                 this->vdw_ACOEFF_tri_matrix(i,j) = this->vdw_well_depth_tri_matrix(i,j)*int_pow<12>(this->vdw_characteristic_dis_tri_matrix(i,j));
//                 this->vdw_BCOEFF_tri_matrix(i,j) = 2.0*this->vdw_well_depth_tri_matrix(i,j)*int_pow<6>(this->vdw_characteristic_dis_tri_matrix(i,j));
//             }
//         }

//         for(Atom_Index i = 0; i < rna_atom_num; ++i) {
//             for(Atom_Index j = rna_atom_num; j < lig_atom_num+rna_atom_num; ++j) {
//                 if(i>j) continue;
//                 const Atom& ai = r.get_atoms_reference()[i];
//                 const Atom& aj = l.get_atoms_reference()[j-rna_atom_num];
//                 this->vdw_well_depth_tri_matrix(i,j) = std::sqrt(ai.get_sybyl_well_depth() * aj.get_sybyl_well_depth());
//                 this->vdw_characteristic_dis_tri_matrix(i,j) = ai.get_sybyl_vdw_radius() + aj.get_sybyl_vdw_radius();
//                 this->vdw_ACOEFF_tri_matrix(i,j) = this->vdw_well_depth_tri_matrix(i,j)*int_pow<12>(this->vdw_characteristic_dis_tri_matrix(i,j));
//                 this->vdw_BCOEFF_tri_matrix(i,j) = 2.0*this->vdw_well_depth_tri_matrix(i,j)*int_pow<6>(this->vdw_characteristic_dis_tri_matrix(i,j));
//             }
//         }


//     }

// };


struct thread_range {
    unsigned int begin;
    unsigned int end;
};

struct Grid {
    bool flag = false;
    std::vector<Atom_Index> rigid_atoms;
};
using Grid_Index = std::vector<Grid>::size_type;
struct Grids {
    Float width;
    std::vector<Grid> g_data;
    Grid_Index g_i, g_j, g_k;
    Grid_Index shift_i, shift_j, shift_k;
    Grid false_grid;
public:
    //column-major
    Grids() : g_i(0), g_j(0), g_k(0), shift_i(0), shift_j(0), shift_k(0), width(0) {}
    Grids(const Grid_Index i, const Grid_Index j, const Grid_Index k, const Grid_Index si, const Grid_Index sj, const Grid_Index sk, const Float w, const Grid& filler_val) : g_data(i*j*k, filler_val), g_i(i), g_j(j), g_k(k), shift_i(si), shift_j(sj), shift_k(sk), width(w) {}
    const Grid& operator()(Grid_Index i, Grid_Index j, Grid_Index k) const {
        Grid_Index gi = this->index(i,j,k);
        if(gi == -1) {
            return this->false_grid;
        }
        else {
            return g_data[this->index(i,j,k)];
        }
    }
    Grid& operator()(Grid_Index i, Grid_Index j, Grid_Index k) {
        Grid_Index gi = this->index(i,j,k);
        if(gi == -1) {
            return this->false_grid;
        }
        else {
            return g_data[this->index(i,j,k)];
        }
    }

    Grid_Index index(Grid_Index i, Grid_Index j, Grid_Index k) const {
        const Grid_Index& ii = i - this->shift_i;
        const Grid_Index& jj = j - this->shift_j;
        const Grid_Index& kk = k - this->shift_k;
        if(kk >= g_k || kk < 0 || jj >= g_j || jj < 0 || ii >= g_i || ii < 0) {
            return -1;
        }
		return ii + g_i*jj + g_i*g_j*kk; //column-major
    }
    Grid_Index dim_1() const { return g_i; }
	Grid_Index dim_2() const { return g_j; }
    Grid_Index dim_3() const { return g_k; }
};


class YW_Score {
public:
    Float rmax;
    Float dr;
    int u_bin_num;
    // Float umax = 5.06356;//in the unit of kBT = 3 kcal/mol ;1kcal/mol = 1.68787kBT
    // Float gmin = exp(-umax);
    // Float kcal2kbT = 1.68787;
    // store atom i,j 's u index
    Matrix<Size_Type> u_index_matrix;
    Matrix<Float> u_matrix;
    std::ostream& log;

    YW_Score(std::ostream& lg) : log(lg) {}

    YW_Score operator=(const YW_Score& rhs) {
        this->rmax = rhs.rmax;
        this->dr = rhs.dr;
        this->u_bin_num = rhs.u_bin_num;
        this->u_index_matrix = rhs.u_index_matrix;
        this->u_matrix = rhs.u_matrix;
        return *this;
    }

    YW_Score(const RNA& rna, const Ligand& lig, const Float& radius, const Float& delr, std::ostream& lg) : rmax(radius), dr(delr), u_bin_num(int(radius/delr)), log(lg) {
        this->log << "rmax: " << this->rmax << " dr: " << this->dr << " u_bin_num: " << this->u_bin_num << std::endl;
        ////////////////////////////////////////////////////
        ////////read u_now.list ////////////////////////////
        // u_now.list record the statistical potential
        ////////////////////////////////////////////////////
        struct U_Info {
            Size_Type index;
            Floats values;
            U_Info(const Size_Type& i, const Floats& v) : index(i), values(v) {}
        };
        std::map<std::string,U_Info> u_name_index_map;
        std::string sl;
        std::ifstream u_nowfile("/home/yuanzhe/TMD/src/u_now.list",std::ios::in);
        assert(u_nowfile);
        unsigned int u_file_line_count = 0;
        while(getline(u_nowfile,sl)) {
            std::istringstream ss(sl);
            std::string buf;
            std::vector<std::string> token;
            while(ss >> buf) token.push_back(buf);
            std::string name = token[0];
            Floats tmp_values;
            for(int j=1;j!=token.size();j++) {
                tmp_values.push_back(std::stod(token[j]));
            }
            u_name_index_map.insert({name,U_Info(u_file_line_count,tmp_values)});
            ++u_file_line_count;
        }
        const Size_Type num_u_name = u_name_index_map.size();
        const Size_Type num_dis_bin = u_name_index_map.begin()->second.values.size();
        for(const auto& uni : u_name_index_map) {
            assert(uni.second.values.size() == num_dis_bin);
        }
        assert(u_file_line_count == num_u_name);
        this->log << " u size ["<< num_u_name << "," << num_dis_bin << "]" << std::endl;

        u_matrix.resize(num_dis_bin,num_u_name,0.0);
        for(const auto& uni : u_name_index_map) {
            for(Size_Type i = 0; i < uni.second.values.size(); ++i) {
                u_matrix(i,uni.second.index) = uni.second.values[i];
            }
        }

        Atom_Index rna_atom_num = rna.atom_num();
        Atom_Index lig_atom_num = lig.atom_num();
        const Atoms& rna_atoms_ref = rna.get_atoms_reference();
        const Atoms& lig_atoms_ref = lig.get_atoms_reference();
        std::vector<std::string> rna_sybyl_types(rna_atom_num,"");
        for(Atom_Index i = 0; i < rna_atom_num; ++i) {
            const Atom& a = rna_atoms_ref[i];
            // const std::string& elemetn_type = a.get_element_type_name();
            // if(elemetn_type == "H") continue;
            const std::string& tmp_sybyl_type = a.get_sybyl_type_name();
            rna_sybyl_types[i] = (tmp_sybyl_type=="F" || tmp_sybyl_type=="Cl" || tmp_sybyl_type=="Br" || tmp_sybyl_type=="I") ? "Ha" : tmp_sybyl_type;
        }
        std::vector<std::string> lig_sybyl_types(lig_atom_num,"");
        for(Atom_Index i = 0; i < lig_atom_num; ++i) {
            const Atom& a = lig_atoms_ref[i];
            // const std::string& elemetn_type = a.get_element_type_name();
            // if(elemetn_type == "H") continue;
            const std::string& tmp_sybyl_type = a.get_sybyl_type_name();
            lig_sybyl_types[i] = (tmp_sybyl_type=="F" || tmp_sybyl_type=="Cl" || tmp_sybyl_type=="Br" || tmp_sybyl_type=="I") ? "Ha" : tmp_sybyl_type;
        }
        this->u_index_matrix.resize(rna_atom_num,lig_atom_num,-1);
        for(Atom_Index i = 0; i < rna_atom_num; ++i) {
            for(Atom_Index j = 0; j < lig_atom_num; ++j) {
                const std::string& rna_type = rna_sybyl_types[i];
                const std::string& lig_type = lig_sybyl_types[j];
                const std::string& atom_pair_type = (lig_type<rna_type) ? lig_type+"-"+rna_type : rna_type+"-"+lig_type;

                const auto& iter = u_name_index_map.find(atom_pair_type);
                if(iter != u_name_index_map.end()) {
                    this->u_index_matrix(i,j) = iter->second.index;
                    // std::cout << iter->second.index << " here " << this->u_index_matrix(i,j) << std::endl;
                }
            }
        }
    }

    const Float evaluate(const std::vector<thread_range>& tr, const RNA& rna, const Ligand& lig, const Grids& grids) const {
        // std::cout << " YW thread num: " << tr.size() << " begin: " << tr[0].begin << " end: " << tr[0].end << std::endl;
        Float score = 0.0;
        const unsigned int begin = tr[0].begin;
        const unsigned int end = tr[0].end;

        // std::mutex thread_lock;
        // thread_lock.lock();
        // std::cout << "evaluate YW thread " << std::this_thread::get_id() << " : " << begin << " " << end << std::endl;
        // thread_lock.unlock();
        const Atoms& rna_atoms_ref = rna.get_atoms_reference();
        const Atoms& lig_atoms_ref = lig.get_atoms_reference();
        for(Atom_Index l_i = begin; l_i < end; ++l_i) {
            const Atom& lig_atom = lig_atoms_ref[l_i];
            const Vec3d& lig_atom_coord = lig_atom.get_coord();
            const int atom_x_grid = static_cast<int>(lig_atom_coord[0]/grids.width);
            const int atom_y_grid = static_cast<int>(lig_atom_coord[1]/grids.width);
            const int atom_z_grid = static_cast<int>(lig_atom_coord[2]/grids.width);

            const Grid& grid = grids(atom_x_grid,atom_y_grid,atom_z_grid);
            if(grid.flag) {
                for(const Atom_Index& r_j : grid.rigid_atoms) {
                    const Atom& rna_atom = rna_atoms_ref[r_j];
                    const Size_Type u_index = this->u_index_matrix(r_j, l_i);
                    if(u_index == -1) continue;
                    const Float& dis = (lig_atom.get_coord() - rna_atom.get_coord()).norm();
                    if(dis > this->rmax) continue;
                    const unsigned int bin_index = int(dis/this->dr);
                    assert(bin_index <= this->u_bin_num);
                    score += this->u_matrix(bin_index,u_index);
                    // std::cout << atom_pair_type << std::endl;
                    // assert(false);
                }
            }
        }
        return score/static_cast<Float>(lig.heavy_atom_num());
    }
};

enum SCORE_MODE {YW_SCORE,VDW_LIGAND,VDW_RNA_LIGAND,ALL};
class Scoring_Function {
    const RNA& rna;
    const Ligand& lig;
    SCORE_MODE Score_Mode;
    // VDW_Score vdw_score;
    YW_Score yw_score;
    unsigned int num_thread = 0;
    std::vector<thread_range> thread_ranges;
    std::ostream& log;

    Grids grids;//stores the RNA atom mapping with the grids
public:
    friend YW_Score;
    // friend VDW_Score;
    //initiate
    Scoring_Function(const RNA& r, const Ligand& l, const SCORE_MODE& sm, std::ostream& lg) : rna(r), lig(l), Score_Mode(sm), log(lg), yw_score(lg) {
    //thread information
        this->log << "rna atom num: " << r.atom_num() << " ligand atom num: " << lig.atom_num() << std::endl;
        this->num_thread = std::thread::hardware_concurrency();
        if(this->num_thread>0) {
            this->log << "std::thread::hardware_concurrency() return <=0, so using 1 thread!" << std::endl;
            this->num_thread = 1;
            thread_range tmp_tr;
            tmp_tr.begin = 0;
            tmp_tr.end = lig.atom_num();
            this->thread_ranges.push_back(tmp_tr);
        } else {
            assert(false);
            int num_per_thread = -1;
            if(this->lig.atom_num()<this->num_thread) {
                this->log << "this->lig.atom_num() is smaller than num thread, using lig.atom_num() as thread num!" << std::endl;
                this->num_thread = this->lig.atom_num();
                num_per_thread = 1;
                for(auto i = 0; i < this->num_thread; ++i) {
                    thread_range tmp_tr;
                    tmp_tr.begin = i;
                    tmp_tr.end = i+1;
                    this->thread_ranges.push_back(tmp_tr);
                }
            } else {
                num_per_thread = this->lig.atom_num()/this->num_thread + 1;
                int total_atom_num_count = 0;
                while(total_atom_num_count < this->lig.atom_num()) {
                    thread_range tmp_tr;
                    tmp_tr.begin = total_atom_num_count;
                    total_atom_num_count += num_per_thread;
                    if(total_atom_num_count > this->lig.atom_num()) total_atom_num_count = this->lig.atom_num();
                    tmp_tr.end = total_atom_num_count;
                    this->thread_ranges.push_back(tmp_tr);
                }
                this->num_thread = this->thread_ranges.size();
            }
        }
        assert(this->thread_ranges.size()==this->num_thread);
        this->log << "thread num: " << this->thread_ranges.size() << " ligand atom num: " << this->lig.atom_num() << std::endl;
        /////////////////////////////////////////////////////
        switch(this->Score_Mode) {
            case YW_SCORE:
                this->yw_score = YW_Score(r,l,10.0,0.2,lg);
                this->initialize_grids(r,this->yw_score.rmax,2.0);
                break;
            case VDW_LIGAND:
                // this->
                break;
            case VDW_RNA_LIGAND:
                // this->
                break;
            case ALL:
                // this->
                break;
            default:
                assert(false);
                break;
        }
    }

    void initialize_grids(const RNA& rna, const Float interation_radius, const Float grid_width) {
        this->log << "start initializing grids..." << std::endl;
        struct Grid_Shift {
            const int x;
            const int y;
            const int z;
            Grid_Shift(const int xx, const int yy, const int zz) : x(xx), y(yy), z(zz) {}
        };
        //get all the grids shifts within the interacting radius(interation_radius)
        const int upper_cubic_shift = static_cast<int>(std::ceil((0.0+interation_radius)/grid_width))+1;//plus extra one to ensure it covers the area
        const int lower_cubic_shift = static_cast<int>(std::floor((0.0-interation_radius)/grid_width))-1;
        std::vector<Grid_Shift> grid_shifts;
        for(int x = lower_cubic_shift; x <= upper_cubic_shift; ++x) {
            for(int y = lower_cubic_shift; y <= upper_cubic_shift; ++y) {
                for(int z = lower_cubic_shift; z <= upper_cubic_shift; ++z) {
                    const Float& x_center = (static_cast<Float>(x)+0.5)*grid_width;
                    const Float& y_center = (static_cast<Float>(y)+0.5)*grid_width;
                    const Float& z_center = (static_cast<Float>(z)+0.5)*grid_width;
                    const Float& dis = std::sqrt((x_center*x_center+y_center*y_center+z_center*z_center));
                    // consider some grids partially covered by sphere
                    const Float& effective_radius = interation_radius + grid_width * std::sqrt(3);
                    if(dis<=effective_radius) {
                        grid_shifts.push_back(Grid_Shift(x,y,z));
                    }
                }
            }
        }
        this->log << "grid_shifts size -> " << grid_shifts.size() << std::endl;
        //initialize grid map
        this->log << "process rna atom grid_map..." << std::endl;
        std::map<std::string,Grid> grid_map;
        std::set<int> grid_x_set;
        std::set<int> grid_y_set;
        std::set<int> grid_z_set;
        Atom_Index rna_atom_index = 0;
        assert(grid_width>k_epsilon);
        const Atoms& rna_atoms_ref = rna.get_atoms_reference();
        for(const Atom& ra : rna_atoms_ref) {
            // progress_bar(rna_atom_index+1,rna.atom_num());
            const Vec3d& atom_coord = ra.get_coord();
            const int atom_x_grid = static_cast<int>(atom_coord[0]/grid_width);
            const int atom_y_grid = static_cast<int>(atom_coord[1]/grid_width);
            const int atom_z_grid = static_cast<int>(atom_coord[2]/grid_width);
            for(const Grid_Shift& gs : grid_shifts) {
                //get absolute grid index of this RNA atom
                const int actual_x_grid = atom_x_grid + gs.x;
                const int actual_y_grid = atom_y_grid + gs.y;
                const int actual_z_grid = atom_z_grid + gs.z;

                grid_x_set.insert(actual_x_grid);
                grid_y_set.insert(actual_y_grid);
                grid_z_set.insert(actual_z_grid);

                const std::string& grid_name = std::to_string(actual_x_grid)+"-"+std::to_string(actual_y_grid)+"-"+std::to_string(actual_z_grid);
                if(grid_map.find(grid_name) != grid_map.end()) {
                    grid_map.at(grid_name).rigid_atoms.push_back(rna_atom_index);
                } else {
                    Grid grid;
                    grid.flag = true;
                    grid.rigid_atoms.push_back(rna_atom_index);
                    grid_map.insert({grid_name,grid});
                }
            }
            rna_atom_index++;
        }
        const int min_grid_x = (*grid_x_set.begin());
        const int min_grid_y = (*grid_y_set.begin());
        const int min_grid_z = (*grid_z_set.begin());
        const int max_grid_x = (*grid_x_set.end());
        const int max_grid_y = (*grid_y_set.end());
        const int max_grid_z = (*grid_z_set.end());
        const unsigned int size_x = max_grid_x - min_grid_x;
        const unsigned int size_y = max_grid_y - min_grid_y;
        const unsigned int size_z = max_grid_z - min_grid_z;
        this->grids = Grids(size_x,size_y,size_z,min_grid_x,min_grid_y,min_grid_z,grid_width,Grid());

        for(int i = min_grid_x; i <= max_grid_x; ++i) {
            for(int j = min_grid_y; j <= max_grid_y; ++j) {
                for(int k = min_grid_z; k <= max_grid_z; ++k) {
                    const std::string& grid_name = std::to_string(i)+"-"+std::to_string(j)+"-"+std::to_string(k);
                    if(grid_map.find(grid_name) != grid_map.end()) {
                        this->grids(i,j,k) = grid_map.at(grid_name);
                    }
                }
            }
        }
        this->log << "finish initializing grids... size: " << this->grids.dim_1() << " " << this->grids.dim_2() << " " << this->grids.dim_3() << std::endl;
    }
    // void set_score_type(const SCORE_MODE st ) {
    //     this->Score_Mode = st;
    // }

    const Float evaluate() {

        switch(this->Score_Mode) {
            case YW_SCORE:
                return this->yw_score.evaluate(this->thread_ranges,this->rna,this->lig,this->grids);
                break;
            case VDW_LIGAND:
                // return this->vdw_score.evaluate(this->thread_ranges,this->lig,this->ligand_pairs);
                return 0.0;
                break;
            case VDW_RNA_LIGAND:
                // return this->vdw_score.evaluate(this->thread_ranges,this->rna,this->lig,this->grids);
                return 0.0;
                break;
            case ALL: {
                // Float score = 0.0;
                // score += this->yw_score.evaluate(this->thread_ranges,this->rna,this->lig,this->grids);
                // score += this->vdw_score.evaluate(this->thread_ranges,this->lig,this->ligand_pairs);
                // score += this->vdw_score.evaluate(this->thread_ranges,this->rna,this->lig,this->grids);
                // return score;
                return 0.0;
                break;
            }
            default:
                assert(false);
                break;
        }

        // const unsigned int num_thread = tr.size();
        // std::thread t[num_thread];
        // std::vector<Float> score_thread(num_thread,0);
        // for(unsigned int i = 0; i != num_thread; ++i) {
        //     t[i] = std::thread(yw_score_calculation, tr[i].begin, tr[i].end, std::ref(rna), std::ref(lig), std::ref(grids), std::ref(*this), std::ref(score_thread[i]));
        // }
        // for(unsigned int i = 0; i != num_thread; ++i) {
        //     t[i].join();
        // }
        // return std::accumulate(score_thread.begin(), score_thread.end(), 0.0)/static_cast<Float>(lig.heavy_atom_num());
    }
};

}
