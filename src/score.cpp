#include <iostream>
#include <sstream>
#include <istream>
#include <ostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include "score.h"
#include "molecule.h"
#include "matrix.h"

namespace tmd {

const Float VDW_Score::evaluate(const std::vector<thread_range>& tr, const RNA& rna, const Ligand& lig, const Grids& grids) const {
    ;
    return 0;
}
const Float VDW_Score::evaluate(const std::vector<thread_range>& tr, const Ligand& lig, const Interacting_Pairs& ips) const {
    Float score = 0.0;
    const Atoms& lig_atoms = lig.get_atoms_reference();
    for(const Interacting_Pair& ip : ips) {
        const Float dis = (lig_atoms[ip.i].get_coord()-lig_atoms[ip.j].get_coord()).norm();
        if(dis < this->lower_cutoff) {
            std::cout << "dis " << dis << " atom_index: " << ip.i << " " << ip.j << std::endl;
            // std::cout << dis << " " << lig_atoms[ip.i].get_name() << " " << lig_atoms[ip.j].get_name() << std::endl;
            assert(false);
            return k_max_float;
        }
        const Float dis_6 = int_pow<6>(dis);
        score += ip.vdw_ACOEFF/dis_6*dis_6 - ip.vdw_BCOEFF/dis_6;
    }
    return score;
}

inline void yw_score_calculation(unsigned int begin, unsigned int end, const RNA& rna, const Ligand& lig, const Grids& grids, const YW_Score& yw_score, Float& score) {
    // std::mutex thread_lock;
    // thread_lock.lock();
    // std::cout << "evaluate YW thread " << std::this_thread::get_id() << " : " << begin << " " << end << std::endl;
    // thread_lock.unlock();
    const Atoms& rna_atoms = rna.get_atoms_reference();
    const Atoms& lig_atoms = lig.get_atoms_reference();
    for(unsigned int i = begin; i < end; ++i) {
        const Atom& lig_atom = lig_atoms[i];
        const std::string lig_atom_elemetn_type = lig_atom.get_element_type_name();
        if(lig_atom_elemetn_type == "H") continue;
        const std::string tmp_lig_atom_sybyl_type = lig_atom.get_sybyl_type_name();
        const std::string lig_atom_sybyl_type = (tmp_lig_atom_sybyl_type=="F" || tmp_lig_atom_sybyl_type=="Cl" || tmp_lig_atom_sybyl_type=="Br" || tmp_lig_atom_sybyl_type=="I") ? "Ha" : tmp_lig_atom_sybyl_type;

        const Vec3d lig_atom_coord = lig_atom.get_coord();
        const int atom_x_grid = static_cast<int>(lig_atom_coord[0]/grids.width);
        const int atom_y_grid = static_cast<int>(lig_atom_coord[1]/grids.width);
        const int atom_z_grid = static_cast<int>(lig_atom_coord[2]/grids.width);
        const std::string grid_name = std::to_string(atom_x_grid)+"-"+std::to_string(atom_y_grid)+"-"+std::to_string(atom_z_grid);

        if(grids.grid_map.find(grid_name)!=grids.grid_map.end()) {
            for(const Atom_Index& j : grids.grid_map.at(grid_name).rigid_atom_mapping) {
                const Atom& rna_atom = rna_atoms[j];
                const std::string rna_atom_elemetn_type = rna_atom.get_element_type_name();
                if(rna_atom_elemetn_type == "H") continue;
                const std::string tmp_rna_atom_sybyl_type = rna_atom.get_sybyl_type_name();
                const std::string rna_atom_sybyl_type = (tmp_rna_atom_sybyl_type=="F" || tmp_rna_atom_sybyl_type=="Cl" || tmp_rna_atom_sybyl_type=="Br" || tmp_rna_atom_sybyl_type=="I") ? "Ha" : tmp_rna_atom_sybyl_type;
                const std::string atom_pair_type = (lig_atom_sybyl_type<rna_atom_sybyl_type) ? lig_atom_sybyl_type+"-"+rna_atom_sybyl_type : rna_atom_sybyl_type+"-"+lig_atom_sybyl_type;
                if(yw_score.contact_type.count(atom_pair_type) == 0) continue;
                const Float& dis = (lig_atom_coord - rna_atom.get_coord()).norm();
                if(dis > yw_score.rmax) continue;
                const unsigned int bin_index = int(dis/yw_score.dr);
                assert(bin_index <= yw_score.usize);
                if(yw_score.u.count(atom_pair_type)!=0) {
                    score += yw_score.u.at(atom_pair_type).at(bin_index);
                } else {
                    std::cout << atom_pair_type << std::endl;
                    assert(false);
                }
                // score += yw_score.u.at()[bin_index];
            }
        }
    }
}


// void yw_score_calculation(unsigned int begin, unsigned int end, const RNA& rna, const Ligand& lig, const YW_Score& yw_score, Float& score) {
//     // std::mutex thread_lock;
//     // thread_lock.lock();
//     // std::cout << "evaluate YW thread " << std::this_thread::get_id() << " : " << begin << " " << end << std::endl;
//     // thread_lock.unlock();
//     // const Atoms& rna_atoms = rna.get_atoms_reference();
//     const Atoms& lig_atoms = lig.get_atoms_reference();
//     for(unsigned int i = begin; i < end; ++i) {
//         const Atom& lig_atom = lig_atoms[i];
//         const std::string lig_atom_elemetn_type = lig_atom.get_element_type_name();
//         if(lig_atom_elemetn_type == "H") continue;
//         const std::string tmp_lig_atom_sybyl_type = lig_atom.get_sybyl_type_name();
//         const std::string lig_atom_sybyl_type = (tmp_lig_atom_sybyl_type=="F" || tmp_lig_atom_sybyl_type=="Cl" || tmp_lig_atom_sybyl_type=="Br" || tmp_lig_atom_sybyl_type=="I") ? "Ha" : tmp_lig_atom_sybyl_type;
//         const Vec3d atom_coord = lig_atom.get_coord();
//         const int atom_x_grid = static_cast<int>(atom_coord[0]);
//         const int atom_y_grid = static_cast<int>(atom_coord[1]);
//         const int atom_z_grid = static_cast<int>(atom_coord[2]);
//         const std::string grid_name = tmp_lig_atom_sybyl_type+"-"+std::to_string(atom_x_grid)+"-"+std::to_string(atom_y_grid)+"-"+std::to_string(atom_z_grid);
//         if(yw_score.grids.grid_map.find(grid_name)!=yw_score.grids.grid_map.end()) {
//             score += yw_score.grids.grid_map.at(grid_name).score;
//         }
//     }
// }

const Float YW_Score::evaluate(const std::vector<thread_range>& tr, const RNA& rna, const Ligand& lig, const Grids& grids) const {
    // std::cout << " YW thread num: " << tr.size() << " begin: " << tr[0].begin << " end: " << tr[0].end << std::endl;
    Float score = 0.0;
    yw_score_calculation(tr[0].begin, tr[0].end, rna, lig, grids, *this, score);
    return score/static_cast<Float>(lig.heavy_atom_num());
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

Scoring_Function::Scoring_Function(const RNA& r, const Ligand& l, std::ostream& log) : rna(r), lig(l) {
    //thread information
    log << "rna atom num: " << r.atom_num() << " ligand atom num: " << lig.atom_num() << std::endl;
    this->num_thread = std::thread::hardware_concurrency();
    if(this->num_thread>0) {
        log << "std::thread::hardware_concurrency() return <=0, so using 1 thread!" << std::endl;
        this->num_thread = 1;
        thread_range tmp_tr;
        tmp_tr.begin = 0;
        tmp_tr.end = lig.atom_num();
        this->thread_ranges.push_back(tmp_tr);
    } else {
        assert(false);
        int num_per_thread = -1;
        if(this->lig.atom_num()<this->num_thread) {
            log << "this->lig.atom_num() is smaller than num thread, using lig.atom_num() as thread num!" << std::endl;
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
    log << "thread num: " << this->thread_ranges.size() << " ligand atom num: " << this->lig.atom_num() << std::endl;
    /////////////////////////////////////////////////////
    //get the largest interaction radius
    const Float interation_radius = yw_score.rmax; //vdw_score_ligand_self; vdw_score_rna_ligand;

    log << "start initializing grids..." << std::endl;
    struct Grid_Shift {
        const int x;
        const int y;
        const int z;
        Grid_Shift(const int xx, const int yy, const int zz) : x(xx), y(yy), z(zz) {}
    };
    //get all the grids shifts within the interacting radius(interation_radius)
    const int upper_cubic_shift = static_cast<int>(std::ceil((0.0+interation_radius)/this->grids.width))+1;//plus extra one to ensure it covers the area
    const int lower_cubic_shift = static_cast<int>(std::floor((0.0-interation_radius)/this->grids.width))-1;
    std::vector<Grid_Shift> grid_shifts;
    for(int x = lower_cubic_shift; x <= upper_cubic_shift; ++x) {
        for(int y = lower_cubic_shift; y <= upper_cubic_shift; ++y) {
            for(int z = lower_cubic_shift; z <= upper_cubic_shift; ++z) {
                // const Float x_center = (static_cast<Float>(x)+0.5)*this->grids.width;
                // const Float y_center = (static_cast<Float>(y)+0.5)*this->grids.width;
                // const Float z_center = (static_cast<Float>(z)+0.5)*this->grids.width;
                // const Float dis = std::sqrt((x_center*x_center+y_center*y_center+z_center*z_center));
                // if(dis<=interation_radius) {
                    grid_shifts.push_back(Grid_Shift(x,y,z));
                // }
            }
        }
    }
    log << "grid_shifts size -> " << grid_shifts.size() << std::endl;
    //initialize grid map
    log << "process rna atom grid_map..." << std::endl;
    Atom_Index rna_atom_index = 0;
    assert(this->grids.width>k_epsilon);
    const Atoms& rna_atoms = rna.get_atoms_reference();
    for(const Atom& ra : rna_atoms) {
        // progress_bar(rna_atom_index+1,rna.atom_num());
        const Vec3d atom_coord = ra.get_coord();
        const int atom_x_grid = static_cast<int>(atom_coord[0]/this->grids.width);
        const int atom_y_grid = static_cast<int>(atom_coord[1]/this->grids.width);
        const int atom_z_grid = static_cast<int>(atom_coord[2]/this->grids.width);
        for(const Grid_Shift& gs : grid_shifts) {
            //get absolute grid index of this RNA atom
            const int actual_x_grid = atom_x_grid + gs.x;
            const int actual_y_grid = atom_y_grid + gs.y;
            const int actual_z_grid = atom_z_grid + gs.z;
            const std::string grid_name = std::to_string(actual_x_grid)+"-"+std::to_string(actual_y_grid)+"-"+std::to_string(actual_z_grid);
            if(this->grids.grid_map.find(grid_name)!=this->grids.grid_map.end()) {
                this->grids.grid_map.at(grid_name).rigid_atom_mapping.push_back(rna_atom_index);
            } else {
                this->grids.grid_map.insert({grid_name,Grid(actual_x_grid,actual_y_grid,actual_z_grid,rna_atom_index)});
            }
        }
        rna_atom_index++;
    }
    log << "finish initializing grids..." << this->grids.grid_map.size() << std::endl;


    //precaluate interacting pairs
    Atom_Index num_of_ligand_interacing_pairs = 0;
    const Atom_Index num_lig_atom = lig.atom_num();
    const Atoms& lig_atoms = lig.get_atoms_reference();
    this->ligand_pairs.resize((num_lig_atom*(num_lig_atom-1))/2);
    for(Atom_Index i = 0; i < num_lig_atom; ++i) {
        for(Atom_Index j = i+1; j < num_lig_atom; ++j) {
            ++num_of_ligand_interacing_pairs;
            const Float rij = lig_atoms[i].get_sybyl_vdw_radius() + lig_atoms[j].get_sybyl_vdw_radius();
            const Float eij = std::sqrt(lig_atoms[i].get_sybyl_well_depth()*lig_atoms[j].get_sybyl_well_depth());
            this->ligand_pairs[num_of_ligand_interacing_pairs-1] = std::move(Interacting_Pair(i,j,eij,rij));
        }
    }
    //init score function parameters
    this->yw_score.init(r,l,log);
    this->vdw_score.init(r,l,this->ligand_pairs,log);
}

}
