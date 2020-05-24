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
#include <thread>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <time.h>
#include <functional>

#include "atom.h"
#include "molecule.h"
#include "grid.h"
#include "vec3d.h"
#include "print.h"

namespace tmd {

struct Pocket_Site {
    double x;
    double y;
    double z;
};

class Pocket {
    std::ostream& tee;
public:
    std::vector<Pocket_Site> sites;
    Grids grids;
    Pocket(std::ostream& lg) : tee(lg) {}
    void find_pocket(const tmd::RNA& yz_rna) {

        //some constants
        const double r_water = 1.4;
        const double radius_small_sphere = 1.5;
        const double radius_large_sphere = 6.0;
        // const double Max_pocket_site = 100000;

        const Atoms& rna_atoms_ref = yz_rna.get_atoms_reference();

        tee << "pocket->rna_atoms_ref.size(): " << rna_atoms_ref.size() << std::endl;

        Atoms rna_heavy_atoms;

        Vec3d rec_min(rna_atoms_ref[0].get_coord());
        Vec3d rec_max(rec_min);
        Vec3d rec_center(0,0,0);

        int index = 0;
        for(const Atom& a : rna_atoms_ref) {
            if(a.get_element_type() != EL_H) {
                //reset radius
                // Atom b = a;
                // b.set
                //should I reset index/serial?
                rna_heavy_atoms.push_back(a);
                index++;

                const Vec3d& a_coord = a.get_coord();

                if(rec_max[0] < a_coord[0]) { rec_max[0] = a_coord[0]; }
                if(rec_max[1] < a_coord[1]) { rec_max[1] = a_coord[1]; }
                if(rec_max[2] < a_coord[2]) { rec_max[2] = a_coord[2]; }
                if(rec_min[0] > a_coord[0]) { rec_min[0] = a_coord[0]; }
                if(rec_min[1] > a_coord[1]) { rec_min[1] = a_coord[1]; }
                if(rec_min[2] > a_coord[2]) { rec_min[2] = a_coord[2]; }

                // print(rec_center,tee); tee << std::endl;

                rec_center += a_coord;
            }
        }
        rec_center = rec_center * (1.0/rna_heavy_atoms.size());

        tee << "pocket->rec_center: ";
        print(rec_center,tee);
        tee << " pocket->rec_min: ";
        print(rec_min,tee);
        tee << " pocket->rec_max: ";
        print(rec_max,tee);
        tee << std::endl;

        //////////Grid for binding site
        const double grid_width = 0.5;
        //grids size
        const Vec3d size_vec = (rec_max - rec_min + Vec3d(1,1,1)*4*r_water )*(1.0/grid_width) + Vec3d(1,1,1)*0.1;
        const int size_x = int(size_vec[0]);
        const int size_y = int(size_vec[1]);
        const int size_z = int(size_vec[2]);
        const Vec3d shift_vec = (rec_min - Vec3d(1,1,1)*2*r_water + Vec3d(1,1,1)*grid_width)*(1.0/grid_width);
        const int shift_x = int(shift_vec[0]);
        const int shift_y = int(shift_vec[1]);
        const int shift_z = int(shift_vec[2]);

        tee << "pocket->size_vec: ";
        print(size_vec,tee);
        tee << " pocket->shift_vec: ";
        print(shift_vec,tee);
        tee << " pocket->size_xyz: ";
        print(std::vector<int>{size_x,size_y,size_z},tee);
        tee << " pocket->shift_xyz: ";
        print(std::vector<int>{shift_x,shift_y,shift_z},tee);
        tee << std::endl;

        this->grids = Grids(size_x,size_y,size_z,shift_x,shift_y,shift_z,grid_width,Grid(true));

        int grids_true_count = 0;
        for (int i = shift_x; i < size_x+shift_x; i++) {
            for (int j = shift_y; j < size_y+shift_y; j++) {
                for (int k = shift_z; k < size_z+shift_z; k++) {
                    if(this->grids(i,j,k).flag) {
                        grids_true_count++;
                    }
                }
            }
        }
        tee << "pocket->true_grids_count: " << grids_true_count << std::endl;

        //////kick out girds on RNA
        for(const Atom& rr : rna_heavy_atoms) {
            const Vec3d& grid_min_vec = (rr.get_coord()-Vec3d(1,1,1)*(rr.get_element_vdw_radius()+radius_small_sphere)) * (1.0/grid_width);
            const Vec3d& grid_max_vec = (rr.get_coord()+Vec3d(1,1,1)*(rr.get_element_vdw_radius()+radius_small_sphere)) * (1.0/grid_width) + Vec3d(1,1,1);

            int xmin = grid_min_vec[0];
            int xmax = grid_max_vec[0];
            int ymin = grid_min_vec[1];
            int ymax = grid_max_vec[1];
            int zmin = grid_min_vec[2];
            int zmax = grid_max_vec[2];

            if(xmin<shift_x) xmin=shift_x;
            if(xmax>=size_x+shift_x) xmax=size_x+shift_x-1;

            if(ymin<shift_y) ymin=shift_y;
            if(ymax>=size_y+shift_y) ymax=size_y+shift_y-1;

            if(zmin<shift_z) zmin=shift_z;
            if(zmax>=size_z+shift_z) zmax=size_z+shift_z-1;

            for(int i=xmin; i<=xmax; i++)
            for(int j=ymin; j<=ymax; j++)
            for(int k=zmin; k<=zmax; k++)
            {
                if(this->grids(i,j,k).flag == true)
                {
                    const double dis = (Vec3d(grid_width*i,grid_width*j,grid_width*k) - rr.get_coord()).norm();
                    if( dis < rr.get_element_covalent_radius() + radius_small_sphere) {
                        this->grids(i,j,k).flag = false;
                        grids_true_count--;
                    }

                }
            }
        }

        tee << "pocket->after kick out rna, true_grids_count remain: " << grids_true_count << std::endl;


        for (int i = shift_x; i < size_x+shift_x; i++)
        for (int j = shift_y; j < size_y+shift_y; j++)
        for (int k = shift_z; k < size_z+shift_z; k++)
        {
            double xx = i*grid_width;
            double yy = j*grid_width;
            double zz = k*grid_width;

            int totalblock,blockrx=0,blocklx=0,blockry=0,blockly=0,blockrz=0,blocklz=0;

            if(this->grids(i,j,k).flag == true) {

                for(const Atom& a : rna_heavy_atoms){
                    const Vec3d& a_coord = a.get_coord();
                    double disxx = (xx - a_coord[0])*(xx - a_coord[0]);
                    double disyy = (yy - a_coord[1])*(yy - a_coord[1]);
                    double diszz = (zz - a_coord[2])*(zz - a_coord[2]);

                    if(blockrx==0){
                        if(a_coord[0] < xx && xx < a_coord[0] + radius_large_sphere) {
                            double dd = sqrt(disyy+diszz);
                            if(dd < a.get_element_covalent_radius() + radius_small_sphere) blockrx=1;
                        }
                    }

                    if(blocklx==0){
                        if(a_coord[0] > xx && xx > a_coord[0] - radius_large_sphere) {
                            double dd = sqrt(disyy+diszz);
                            if(dd < a.get_element_covalent_radius() + radius_small_sphere) blocklx=1;
                        }
                    }

                    if(blockry==0){
                        if(a_coord[1] < yy && yy < a_coord[1] + radius_large_sphere) {
                            double dd = sqrt(disxx+diszz);
                            if(dd < a.get_element_covalent_radius() + radius_small_sphere) blockry=1;
                        }
                    }

                    if(blockly==0){
                        if(a_coord[1] > yy && yy > a_coord[1] - radius_large_sphere) {
                            double dd = sqrt(disxx+diszz);
                            if(dd < a.get_element_covalent_radius() + radius_small_sphere) blockly=1;
                        }
                    }

                    if(blockrz==0){
                        if(a_coord[2] < zz && zz < a_coord[2] + radius_large_sphere) {
                            double dd = sqrt(disyy+disxx);
                            if(dd < a.get_element_covalent_radius() + radius_small_sphere) blockrz=1;
                        }
                    }

                    if(blocklz==0){
                        if(a_coord[2] > zz && zz > a_coord[2] - radius_large_sphere) {
                            double dd = sqrt(disyy+disxx);
                            if(dd < a.get_element_covalent_radius() + radius_small_sphere) blocklz=1;
                        }

                    }

                    totalblock = blocklx+blockrx+blockly+blockry+blocklz+blockrz;
                    if(totalblock>=4) break;
                }

                if(totalblock>=4) {
                    Pocket_Site p_tmp;
                    p_tmp.x = xx;
                    p_tmp.y = yy;
                    p_tmp.z = zz;
                    this->sites.push_back(p_tmp);
                    // if(this->sites.size() > Max_pocket_site) {
                    //     std::cout << "Too much pocket sites" << std::endl;
                    //     exit(2);
                    // }
                } else {
                    this->grids(i,j,k).flag == false;
                    grids_true_count--;
                    assert(grids_true_count>=0);
                }
            }
        }
        tee << "pocket->after pocket finding, true_grids_count remain: " << grids_true_count << std::endl;

        // for(int i = this->grids) {}
    }
};


}