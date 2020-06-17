#pragma once

#include <algorithm>
#include "common.h"
#include "atom.h"
#include "context.h"

namespace tmd {

// using Coords = std::vector<Vec3d>;

struct Conformer {
    Float score;
    Float rmsd;
    Atoms atoms;
    // Coords coords;
    Floats dofs;
    // Conformer() {}
    Conformer(const Float s, const Float r, const Atoms& a, const Floats& d) : score(s), rmsd(r), atoms(a), dofs(d) {}
    void write(const Contexts& cs, std::ostream& out = std::cout) {
        write_contexts(cs,this->atoms,out);
    }
};
// using Conformers = std::vector<Conformer>;

class Conformers {
    std::vector<Conformer> confs;
    std::vector<int> heavy_atom_indices;
    std::vector<int> cluster_indices;
    Contexts contexts;
public:
    Conformers(const Contexts& cs) : contexts(cs) {}
    //external refers to shift and overall orientation
    bool external_too_close(const Floats& new_dofs, const Floats& dofs) {
        for(int i = 0; i < 3; ++i) {
            const Float& diff = new_dofs[i] - dofs[i];
            if(diff > 0.25) {
                return false;
            }
        }
        const Vec3d new_orientaion(new_dofs[3],new_dofs[4],new_dofs[5]);
        const Vec3d old_orientation(dofs[3],dofs[4],dofs[5]);
        const Float& dot_norm = new_orientaion * old_orientation;
        const Float& new_norm = new_orientaion.norm();
        const Float& old_norm = old_orientation.norm();
        const Float& cos_theta = dot_norm / (new_norm*old_norm);
        if(std::abs(cos_theta-1) > 0.01 ) {
            return false;
        } else {
            const Float& diff = (new_norm-old_norm);
            const int& diff_n = diff/(2*k_pi);
            if(std::abs(diff-diff_n*2*k_pi) > 0.1) {
                return false;
            }
        }
        return true;
    }
    //internal refers to torsions
    bool internal_too_close(const Floats& new_dofs, const Floats& dofs) {
        for(int i = 6; i < dofs.size(); ++i) {
            const Float& diff = new_dofs[i] - dofs[i];
            const int& diff_n = diff / (2*k_pi);
            if(std::abs(diff-diff_n*2*k_pi) > 0.1) {
                return false;
            }
        }
        return true;
    }
    bool too_close(const Floats& new_dofs) {
        for(int i = 0; i < this->confs.size(); ++i) {
            const Floats& dofs = this->confs[i].dofs;
            #ifdef DEBUG
            assert(dofs.size() >= 6);
            #endif
            if(dofs.size() == 6) {
                if(this->external_too_close(new_dofs,dofs)) {
                    return true;
                }
            } else if(dofs.size() > 6) {
                if(this->external_too_close(new_dofs,dofs) && this->internal_too_close(new_dofs,dofs)) {
                    return true;
                }
            } else {
                assert(false);
            }
        }
        return false;
    }
    void add_conformer(const Float s, const Float r, const Atoms& a, const Floats& d) {
        if(!this->too_close(d)) {
            this->confs.push_back(Conformer(s,r,a,d));
        }
    }

    void prune_by_cutoff(const Float& cutoff) {
        for(std::vector<Conformer>::iterator c = this->confs.begin(); c != this->confs.end(); ++c) {
            if((c->score - this->confs[0].score) > cutoff) {
                this->confs.erase(c,this->confs.end());
                break;
            }
        }
    }

    void prune_by_num(const int& topN) {
        if(this->confs.size() > topN) {
            this->confs.erase(this->confs.begin()+topN,this->confs.end());
        }
    }

    void sort_by_score() {
        std::sort(this->confs.begin(), this->confs.end(), [](const Conformer& c1, const Conformer& c2){ return c1.score < c2.score; });
    }
    void get_heavy_atom_indices() {
        this->heavy_atom_indices.clear();
        #ifdef DEBUG
        assert(this->confs.size() > 0);
        #endif
        const Conformer& c = this->confs[0];
        int heavy_atom_index = 0;
        for(int i = 0; i < c.atoms.size(); ++i) {
            const Atom& a = c.atoms[i];
            if(a.get_element_type() != EL_H) {
                this->heavy_atom_indices.push_back(heavy_atom_index);
            }
            ++heavy_atom_index;
        }
    }

    const Float rmsd_between_conformer_fast(int ci, int cj) const {
        const Conformer& c1 = this->confs[ci];
        const Conformer& c2 = this->confs[cj];
        Float rmsd = 0.0;
        #ifdef DEBUG
        assert(c1.atoms.size() == c2.atoms.size());
        #endif
        int heavy_atom_size = this->heavy_atom_indices.size();
        for(int i = 0; i < this->heavy_atom_indices.size(); ++i) {
            const Atom& a2 = c2.atoms[i];
            const Atom& a1 = c1.atoms[i];
            rmsd += (a1.get_coord() - a2.get_coord()).norm_square();
        }
        return std::sqrt(rmsd/Float(heavy_atom_size));
    }

    const Float rmsd_between_conformer(int ci, int cj) const {
        const Conformer& c1 = this->confs[ci];
        const Conformer& c2 = this->confs[cj];
        Float rmsd = 0.0;
        #ifdef DEBUG
        assert(c1.atoms.size() == c2.atoms.size());
        #endif
        int heavy_atom_size = 0;
        for(int i = 0; i < c1.atoms.size(); ++i) {
            const Atom& a2 = c2.atoms[i];
            if(a2.get_element_type() == EL_H) continue;
            const Atom& a1 = c1.atoms[i];
            ++heavy_atom_size;
            rmsd += (a1.get_coord() - a2.get_coord()).norm_square();
        }
        return std::sqrt(rmsd/Float(heavy_atom_size));
    }

    const int conformer_num() const {
        return this->confs.size();
    }
    const std::vector<Conformer>& get_conformers() const {
        return this->confs;
    }

    const std::vector<int>& get_clusters() const {
        return this->cluster_indices;
    }
    void clustering() {
        //first update heavy atom size
        this->get_heavy_atom_indices();
        //suppose to use some clustering method here to cluster the conformers
        //update clusters
    }
    const int cluster_num() const {
        return this->cluster_indices.size();
    }

    void output_conformer(const int ci, const int model_num, std::ofstream& out) {
        #ifdef DEBUG
        assert(ci < this->confs.size());
        #endif
        out << "MODEL " << model_num << std::endl;
        out << "REMARK VINA RESULT:      " << this->confs[ci].score << "      " << this->confs[ci].rmsd << "      " << this->confs[ci].rmsd << std::endl;
        this->confs[ci].write(this->contexts, out);
        out << "ENDMDL" << std::endl;
    }

};

}