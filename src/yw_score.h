#pragma once

#include <ostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <fstream>

#include "common.h"
#include "matrix.h"
#include "molecule.h"
#include "atom.h"
#include "grid.h"


namespace tmd {

class YW_Score {
public:
    const Float rmax = 10.0;
    const Float dr = 0.2;
    const int u_bin_num = static_cast<int>(rmax/dr);
    // Float umax = 5.06356;//in the unit of kBT = 3 kcal/mol ;1kcal/mol = 1.68787kBT
    // Float gmin = exp(-umax);
    // Float kcal2kbT = 1.68787;
    // store atom i,j 's u index
    Matrix<int> u_index_matrix;
    Matrix<Float> u_matrix;
    std::ostream& tee;

    Grids grids;//stores the RNA atom mapping with the grids

    YW_Score(std::ostream& lg) : tee(lg) {
        this->tee << "YW_Score init" << std::endl;
    }

    // YW_Score operator=(const YW_Score& rhs) {
    //     this->rmax = rhs.rmax;
    //     this->dr = rhs.dr;
    //     this->u_bin_num = rhs.u_bin_num;
    //     this->u_index_matrix = rhs.u_index_matrix;
    //     this->u_matrix = rhs.u_matrix;
    //     return *this;
    // }

    void init(const RNA& rna, const Ligand& lig) {
        // processing input parameter and files
        this->tee << "rmax: " << this->rmax << " dr: " << this->dr << " u_bin_num: " << this->u_bin_num << std::endl;
        ////////////////////////////////////////////////////
        ////////read u_now.list ////////////////////////////
        // u_now.list record the statistical potential
        ////////////////////////////////////////////////////
        struct U_Info {
            int index;
            Floats values;
            U_Info(const int& i, const Floats& v) : index(i), values(v) {}
        };
        std::map<std::string,U_Info> u_name_index_map;
        std::string sl;
        std::ifstream u_file("/home/yuanzhe/TMD/src/res/u_now.list",std::ios::in);
        assert(u_file);
        int u_file_line_count = 0;
        while(std::getline(u_file,sl)) {
            if(sl[0]=='#') continue;
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
        const int num_u_name = u_name_index_map.size();
        const int num_dis_bin = u_name_index_map.begin()->second.values.size();
        for(const auto& uni : u_name_index_map) {
            assert(uni.second.values.size() == num_dis_bin);
        }
        assert(u_file_line_count == num_u_name);
        this->tee << " u size ["<< num_u_name << "," << num_dis_bin << "]" << std::endl;

        u_matrix.resize(num_dis_bin,num_u_name,0.0);
        for(const auto& uni : u_name_index_map) {
            for(int i = 0; i < uni.second.values.size(); ++i) {
                u_matrix(i,uni.second.index) = uni.second.values[i];
            }
        }

        int rna_atom_num = rna.atom_num();
        int lig_atom_num = lig.atom_num();
        const Atoms& rna_atoms_ref = rna.get_atoms_reference();
        const Atoms& lig_atoms_ref = lig.get_atoms_reference();
        std::vector<std::string> rna_sybyl_types(rna_atom_num,"");
        for(int i = 0; i < rna_atom_num; ++i) {
            const Atom& a = rna_atoms_ref[i];
            // const std::string& elemetn_type = a.get_element_type_name();
            // if(elemetn_type == "H") continue;
            const std::string& tmp_sybyl_type = a.get_sybyl_type_name();
            rna_sybyl_types[i] = (tmp_sybyl_type=="F" || tmp_sybyl_type=="Cl" || tmp_sybyl_type=="Br" || tmp_sybyl_type=="I") ? "Ha" : tmp_sybyl_type;
        }
        std::vector<std::string> lig_sybyl_types(lig_atom_num,"");
        for(int i = 0; i < lig_atom_num; ++i) {
            const Atom& a = lig_atoms_ref[i];
            // const std::string& elemetn_type = a.get_element_type_name();
            // if(elemetn_type == "H") continue;
            const std::string& tmp_sybyl_type = a.get_sybyl_type_name();
            lig_sybyl_types[i] = (tmp_sybyl_type=="F" || tmp_sybyl_type=="Cl" || tmp_sybyl_type=="Br" || tmp_sybyl_type=="I") ? "Ha" : tmp_sybyl_type;
        }
        this->u_index_matrix.resize(rna_atom_num,lig_atom_num,-1);
        for(int i = 0; i < rna_atom_num; ++i) {
            for(int j = 0; j < lig_atom_num; ++j) {
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
        // find the grids
        this->tee << "start initializing grids..." << std::endl;
        struct Grid_Shift {
            const int x;
            const int y;
            const int z;
            Grid_Shift(const int xx, const int yy, const int zz) : x(xx), y(yy), z(zz) {}
        };
        const double interation_radius = this->rmax;
        const double grid_width = 2.0;
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
        this->tee << "grid_shifts size -> " << grid_shifts.size() << std::endl;
        //initialize grid map
        this->tee << "process rna atom grid_map..." << std::endl;
        std::map<std::string,Grid> grid_map;
        std::set<int> grid_x_set;
        std::set<int> grid_y_set;
        std::set<int> grid_z_set;
        int rna_atom_index = 0;
        assert(grid_width>k_epsilon);
        // const Atoms& rna_atoms_ref = rna.get_atoms_reference();
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
        const int size_x = max_grid_x - min_grid_x + 1;
        const int size_y = max_grid_y - min_grid_y + 1;
        const int size_z = max_grid_z - min_grid_z + 1;
        this->grids = Grids(size_x,size_y,size_z,min_grid_x,min_grid_y,min_grid_z,grid_width,Grid(false));

        for(int i = min_grid_x; i <= max_grid_x; ++i) {
            for(int j = min_grid_y; j <= max_grid_y; ++j) {
                for(int k = min_grid_z; k <= max_grid_z; ++k) {
                    const std::string& grid_name = std::to_string(i)+"-"+std::to_string(j)+"-"+std::to_string(k);
                    if(grid_map.find(grid_name) != grid_map.end()) {
                        this->grids.assign(i,j,k,grid_map.at(grid_name));
                    }
                }
            }
        }
        this->tee << "finish initializing grids... size: " << this->grids.dim_1() << " " << this->grids.dim_2() << " " << this->grids.dim_3() << std::endl;
    }

    const Float evaluate(const RNA& rna, const Ligand& lig) const {
        // std::cout << " YW thread num: " << tr.size() << " begin: " << tr[0].begin << " end: " << tr[0].end << std::endl;
        Float score = 0.0;
        // const int begin = tr[0].begin;
        // const int end = tr[0].end;

        // std::mutex thread_lock;
        // thread_lock.lock();
        // std::cout << "evaluate YW thread " << std::this_thread::get_id() << " : " << begin << " " << end << std::endl;
        // thread_lock.unlock();
        const Atoms& rna_atoms_ref = rna.get_atoms_reference();
        const Atoms& lig_atoms_ref = lig.get_atoms_reference();
        for(int l_i = 0; l_i < lig_atoms_ref.size(); ++l_i) {
            const Atom& lig_atom = lig_atoms_ref[l_i];
            const Vec3d& lig_atom_coord = lig_atom.get_coord();
            const int atom_x_grid = static_cast<int>(lig_atom_coord[0]/this->grids.width);
            const int atom_y_grid = static_cast<int>(lig_atom_coord[1]/this->grids.width);
            const int atom_z_grid = static_cast<int>(lig_atom_coord[2]/this->grids.width);

            const Grid& grid = this->grids.at(atom_x_grid,atom_y_grid,atom_z_grid);
            if(grid.flag) {
                for(const int& r_j : grid.rigid_atoms) {
                    const Atom& rna_atom = rna_atoms_ref[r_j];
                    const int u_index = this->u_index_matrix(r_j, l_i);
                    if(u_index == -1) continue;
                    const Float& dis = (lig_atom.get_coord() - rna_atom.get_coord()).norm();
                    if(dis > this->rmax) continue;
                    const int bin_index = static_cast<int>(dis/this->dr);
                    assert(bin_index <= this->u_bin_num);
                    score += this->u_matrix(bin_index,u_index);
                    // std::cout << atom_pair_type << std::endl;
                    // assert(false);
                }
            }
        }
        // return score/static_cast<Float>(lig.heavy_atom_num());
        return score/static_cast<Float>(lig.heavy_atom_num());
    }
};

}