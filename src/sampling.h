#pragma once

#include <vector>
#include "molecule.h"
#include "score.h"
#include "pocket.h"
#include "conformer.h"
//headers from numerical recipes
#include "nr3.h"
#include "ran.h"
#include "amebsa.h"
#include "ludcmp.h"
#include "qrdcmp.h"
#include "roots_multidim.h"
#include "quasinewton.h"

namespace tmd {

//this box defines the docking spaces for the center of the current ligand
//all coordinates stored in box use the same axis frame as the RNA input file
struct Box {
    Vec3d center;
    Vec3d corner1;
    Vec3d corner2;
    Box(const Vec3d c, const Vec3d c1, const Vec3d c2) : center(c), corner1(c1), corner2(c2) {}
};

enum SAMPLE_MODE {SAMPLE_RIGID,SAMPLE_FLEXIBLE};
enum OPTIMIZATION_MODE {SIMULATED_ANNEALING,QUASI_NEWTON};

// class Optimizer {
//     OPTIMIZATION_MODE Optimization_Mode = SIMULATED_ANNEALING;
// };

class Sampling {
    const RNA& rna;
    Ligand& lig;
    Conformers conformers;
    Scoring_Function* score_func_ptr;
    SAMPLE_MODE Sample_Mode = SAMPLE_RIGID;
    OPTIMIZATION_MODE Optimization_Mode = SIMULATED_ANNEALING;
    Box docking_box;
    const Pocket& pocket;
    //random generator
    RNGType& generator;
    std::ostream& tee;

    Float min_score = k_max_float;

public:
    Sampling(const RNA& r, Ligand& l, Scoring_Function* sf_ptr, const Box& b, const Pocket& p, SAMPLE_MODE s_mode, OPTIMIZATION_MODE o_mode, RNGType& g, std::ostream& lg) : rna(r), lig(l), conformers(l.get_contexts()), score_func_ptr(sf_ptr), docking_box(b), pocket(p), Sample_Mode(s_mode), Optimization_Mode(o_mode), generator(g), tee(lg) {}

    // void clear_sampled_ligands() {
    //     this->sampled_ligands.clear();
    // }
    // const Box get_box() const {
    //     return this->docking_box;
    // }
    // //when set box, the current ligand center will be placed at the center of the box
    // void set_box(const Vec3d c, const Vec3d c1, const Vec3d c2) {
    //     this->docking_box = Box(c,c1,c2);
    // }
    // void set_mode(const std::string m = "rigid") {
    //     if(m=="rigid") {
    //         this->MODE = SAMPLE_RIGID;
    //     } else if(m=="flexible") {
    //         this->MODE = SAMPLE_FLEXIBLE;
    //     } else {
    //         assert(false);
    //     }
    // }

    //this operator is overloaded for simulated annealing in nr3
    const Float operator()(const VecDoub& cf) {
        Floats dofs(cf.size(),0.0);
        for(int i = 0; i < dofs.size(); ++i) {
            dofs[i] = cf[i];
        }
        return (*this)(dofs);
    }
    const Float operator()(const Floats& cf) {
        // for(int i = 0; i < cf.size(); ++i) {
        //     if(std::isnan(cf[i])) {
        //         std::cout << " cf[" << i << "] is nan " << std::endl;
        //     }
        // }
        Float constraint_energy = 0.0;
        if(this->Optimization_Mode == SIMULATED_ANNEALING) {
            constraint_energy = this->simulated_check_dof_range(cf);
        } else if(this->Optimization_Mode == QUASI_NEWTON) {
            constraint_energy = this->quasi_check_dof_range(cf);
        } else {
            assert(false);
        }

        if(this->Sample_Mode == SAMPLE_FLEXIBLE) {
            this->lig.sample(cf,LIGAND_SAMPLE_FLEXIBLE);
        } else if(this->Sample_Mode == SAMPLE_RIGID) {
            this->lig.sample(cf,LIGAND_SAMPLE_RIGID);
        } else {
            assert(false);
        }
        const Float& s = score_func_ptr->evaluate();
        if(!eq(constraint_energy,0)) {
            return s+constraint_energy;
        }

        this->conformers.add_conformer(s,this->lig.rmsd_with_respect_to_ref_atoms(),this->lig.get_atoms(),cf);
        int ci = this->conformers.conformer_num() - 1;
        if(s < this->min_score) {
            this->min_score = s;
        }

        // if(std::isnan(s)) {
        //     std::cout << "sampled score is also nan" << std::endl;
        //     for(int i = 0; i < cf.size(); ++i) {
        //         std::cout << cf[i] << " ";
        //     }
        //     std::cout << std::endl;
        //     this->lig.write(std::cout);
        // }
        return s;
    }

    const Float simulated_check_dof_range(const Floats& cf) {
        const Float& x = cf[0];
        const Float& y = cf[1];
        const Float& z = cf[2];
        const Box& b = this->docking_box;
        const Vec3d& ref_center = this->lig.get_ref_heavy_center_coord();
        //corner1 is the lower corner
        const Vec3d& pos = Vec3d(x+ref_center[0],y+ref_center[1],z+ref_center[2]);
        //if outside box return a harmonic potential to get it back
        Float constraint = 0.0;
        if(pos[0] < b.corner1[0] || pos[0] > b.corner2[0]) {
            constraint += (pos[0]-b.corner1[0])*(pos[0]-b.corner1[0]);
        } else if(pos[0] > b.corner2[0]) {
            constraint += (pos[0]-b.corner2[0])*(pos[0]-b.corner2[0]);
        }
        if(pos[1] < b.corner1[1]) {
            constraint += (pos[1]-b.corner1[1])*(pos[1]-b.corner1[1]);
        } else if(pos[1] > b.corner2[1]) {
            constraint += (pos[1]-b.corner2[1])*(pos[1]-b.corner2[1]);
        }
        if(pos[2] < b.corner1[2]) {
            constraint += (pos[2]-b.corner1[2])*(pos[2]-b.corner1[2]);
        } else if(pos[2] > b.corner2[2]) {
            constraint += (pos[2]-b.corner2[2])*(pos[2]-b.corner2[2]);
        }
        return 100*constraint;

        //conver sample conf to ligand conf, complete the rot axis
        const Vec3d rot_vector(cf[3], cf[4], cf[5]);
        const Float rot_angle = rot_vector.norm();
        if(std::abs(rot_angle) > 2*k_pi) return 100*(std::abs(rot_angle) - 2*k_pi)*(std::abs(rot_angle) - 2*k_pi);

        if(this->Sample_Mode == SAMPLE_FLEXIBLE) {
            for(Int i = 6; i < cf.size(); ++i) {
                if(std::abs(cf[i]) > 2*k_pi) return 100*(std::abs(cf[i]) - 2*k_pi)*(std::abs(cf[i]) - 2*k_pi);
            }
        }

        // //not in pocket
        // const int& pocket_grid_x = (x+ref_center[0])/this->pocket.grids.width;
        // const int& pocket_grid_y = (y+ref_center[1])/this->pocket.grids.width;
        // const int& pocket_grid_z = (z+ref_center[2])/this->pocket.grids.width;
        // if(!this->pocket.grids(pocket_grid_x,pocket_grid_y,pocket_grid_z).flag) {
        //     return k_max_float;
        // }
        return 0;
    }

    const Float quasi_check_dof_range(const Floats& cf) {
        const Float& x = cf[0];
        const Float& y = cf[1];
        const Float& z = cf[2];
        const Box& b = this->docking_box;
        const Vec3d& ref_center = this->lig.get_ref_heavy_center_coord();
        //corner1 is the lower corner
        const Vec3d& pos = Vec3d(x+ref_center[0],y+ref_center[1],z+ref_center[2]);
        //if outside box return a harmonic potential to get it back
        Float constraint = 0.0;
        if(pos[0] < b.corner1[0] || pos[0] > b.corner2[0]) {
            constraint += (pos[0]-b.corner1[0])*(pos[0]-b.corner1[0]);
        } else if(pos[0] > b.corner2[0]) {
            constraint += (pos[0]-b.corner2[0])*(pos[0]-b.corner2[0]);
        }
        if(pos[1] < b.corner1[1]) {
            constraint += (pos[1]-b.corner1[1])*(pos[1]-b.corner1[1]);
        } else if(pos[1] > b.corner2[1]) {
            constraint += (pos[1]-b.corner2[1])*(pos[1]-b.corner2[1]);
        }
        if(pos[2] < b.corner1[2]) {
            constraint += (pos[2]-b.corner1[2])*(pos[2]-b.corner1[2]);
        } else if(pos[2] > b.corner2[2]) {
            constraint += (pos[2]-b.corner2[2])*(pos[2]-b.corner2[2]);
        }
        return 100*constraint;

        //conver sample conf to ligand conf, complete the rot axis
        const Vec3d rot_vector(cf[3], cf[4], cf[5]);
        const Float rot_angle = rot_vector.norm();
        if(std::abs(rot_angle) > 2*k_pi) return 100*(std::abs(rot_angle) - 2*k_pi)*(std::abs(rot_angle) - 2*k_pi);

        if(this->Sample_Mode == SAMPLE_FLEXIBLE) {
            for(Int i = 6; i < cf.size(); ++i) {
                if(std::abs(cf[i]) > 2*k_pi) return 100*(std::abs(cf[i]) - 2*k_pi)*(std::abs(cf[i]) - 2*k_pi);
            }
        }

        // //not in pocket
        // const int& pocket_grid_x = (x+ref_center[0])/this->pocket.grids.width;
        // const int& pocket_grid_y = (y+ref_center[1])/this->pocket.grids.width;
        // const int& pocket_grid_z = (z+ref_center[2])/this->pocket.grids.width;
        // if(!this->pocket.grids(pocket_grid_x,pocket_grid_y,pocket_grid_z).flag) {
        //     return k_max_float;
        // }
    }


    // const int conformer_num() const {
    //     return this->conformers.conformer_num();
    // }

    // const Conformers& cluster() const {
    //     return this->conformers;
    // }

    const Conformers& get_conformers() const {
        return this->conformers;
    }

    const Floats get_randomized_dofs() {
        // this->score_func.rl_score.set_mode(lz::RL_LJ_ONLY);
        Floats start_dofs;
        if(this->Sample_Mode == SAMPLE_RIGID) {
            start_dofs.resize(6,0.0);
        } else if(this->Sample_Mode == SAMPLE_FLEXIBLE) {
            start_dofs.resize(this->lig.dof_num(),0.0);
        } else {
            assert(false);
        }
        Float score_of_randomized = k_max_float;
        while(score_of_randomized>500) {
            static int randomized_times = 0;
            randomized_times++;
            const tmd::Vec3d random_xyz = tmd::random_in_box(this->docking_box.corner1, this->docking_box.corner2, generator) - this->lig.get_ref_heavy_center_coord();
            const tmd::Vec3d random_root_rot_vector = tmd::random_unit_vec3d(generator) * tmd::random_double(-tmd::k_pi, tmd::k_pi, generator);
            start_dofs[0] = random_xyz[0];
            start_dofs[1] = random_xyz[1];
            start_dofs[2] = random_xyz[2];
            start_dofs[3] = random_root_rot_vector[0];
            start_dofs[4] = random_root_rot_vector[1];
            start_dofs[5] = random_root_rot_vector[2];

            if(this->Sample_Mode == SAMPLE_FLEXIBLE) {
                for(int i = 6; i < start_dofs.size(); ++i) {
                    start_dofs[i] = tmd::random_double(-tmd::k_pi, tmd::k_pi, generator);
                }
            }

            score_of_randomized = (*this)(start_dofs);

            #ifdef DEBUG
            if(std::isnan(score_of_randomized)) {
                this->tee << "get nan randomized socre" << std::endl;
                this->lig.write(this->tee);
                assert(false);
            }
            #endif
            // this->lig.reset_atoms_coord_to_ref_atoms_coord();
            this->tee << "randomized_times: " << randomized_times << " score+constraints: " << score_of_randomized << std::endl;
            if(randomized_times>=10000) {
                this->tee << "can not be successfully randomized!" << std::endl;
                // assert(false);
            }
        }
        this->tee << "randomized initialized dofs --> ";
        for(const tmd::Float dof : start_dofs) {
            this->tee << dof << " ";
        }
        this->tee << std::endl;
        // this->score_func.rl_score.set_mode(lz::RL_ALL);
        return start_dofs;
    }

    void docking() {
        if(this->Optimization_Mode == SIMULATED_ANNEALING) {
            this->simulated_annealing();
        } else if(this->Optimization_Mode == QUASI_NEWTON) {
            this->quasi_newton();
        } else {
            assert(false);
        }
    }

    void quasi_newton() {
        for(int docking_time = 1; docking_time <= 10; ++docking_time) {
            const Floats start_dofs = this->get_randomized_dofs();
            VecDoub dofs(start_dofs.size(),0.0);
            for(int i = 0; i < start_dofs.size(); ++i) {
                dofs[i] = start_dofs[i];
            }

            Funcd<tmd::Sampling> quasi_func((*this));
            const Doub gtol = 0.001;
            Int iter = 0;
            Doub min_score = 10000000;
            dfpmin(dofs, gtol, iter, min_score, quasi_func);
            this->tee << "finish " << docking_time << "th qusi-newton..." << std::endl;
            this->tee << "use " << iter << " iterations" << std::endl;
            // this->tee << "min dofs: ";
            // for(Int i = 0; i < dofs.size(); ++i) {
            //     this->tee << dofs[i] << " ";
            // }
            // this->tee << std::endl;
            // const tmd::Float final_score = min_score;
            // this->tee << "min score: " << final_score << std::endl;
            // this->tee << "min score's rmsd: " << this->lig.rmsd_with_respect_to_ref_atoms() << std::endl;
            this->conformers.sort_by_score();
            this->conformers.prune_by_num(1000);
        }

        if(this->lig.rmsd_with_respect_to_ref_atoms()<=2.0) {
            std::cout << "successful" << std::endl;
        }
    }

    void simulated_annealing() {
        int docking_times = 10;
        const Floats start_dofs = this->get_randomized_dofs();

        // tmd::Float characteristic_score = (*this)(start_dofs);

        VecDoub dels(start_dofs.size(),0.001);
        // dels[0] = 0.1;
        // dels[1] = 0.1;
        // dels[2] = 0.1;
        Doub ftol = 0.1;
        VecDoub dofs(start_dofs.size(),0.0);
        for(int i = 0; i < start_dofs.size(); ++i) {
            dofs[i] = start_dofs[i];
        }
        Amebsa<tmd::Sampling> amebsa(dofs, dels, (*this), ftol);

        for(int docking_time = 1; docking_time <= docking_times; ++docking_time) {
            Doub temperature = 10000;
            Bool converged = false;
            while(!converged) {
                Int Iter = 1000;
                // if(docking_times >= 3) {
                //     Iter = 1000;
                // }
                // if(temperature>=std::abs(characteristic_score/2.0) && temperature<=std::abs(characteristic_score*1.5)) {
                //     Iter = 10000;
                // }
                converged = amebsa.anneal(Iter, temperature);
                temperature = 0.8*temperature;

                this->tee << "converged: " << converged <<  " score: " << amebsa.yb << " temperature: " << temperature << " already_saved_ligands -> " << this->conformers.conformer_num() << std::endl;
            }

            this->tee << "finish " << docking_time << "th simulated annealing..." << std::endl;
            this->tee << "min dofs: ";
            for(Int i = 0; i < amebsa.pb.size(); ++i) {
                this->tee << amebsa.pb[i] << " ";
            }
            this->tee << std::endl;
            const tmd::Float final_score = (*this)(amebsa.pb);
            assert(tmd::eq(amebsa.yb,final_score));
            this->tee << "min score: " << amebsa.yb << std::endl;
            this->tee << "min score's rmsd: " << this->lig.rmsd_with_respect_to_ref_atoms() << std::endl;
            // characteristic_score = final_score;
            this->conformers.sort_by_score();
            this->conformers.prune_by_num(1000);
        }
        if(this->lig.rmsd_with_respect_to_ref_atoms()<=2.0) {
            std::cout << "successful" << std::endl;
        }
    }
};


}
