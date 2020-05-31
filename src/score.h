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
// #include <unordered_map>

#include "atom.h"
#include "matrix.h"
#include "molecule.h"
// #include "rl_find_pocket.h"
// #include "progress_bar.h"
#include "grid.h"
#include "rl_score.h"
#include "yw_score.h"

namespace tmd {

struct thread_range {
    int begin;
    int end;
};

enum SCORE_MODE {YW_SCORE,RL_SCORE,VDW_LIGAND,VDW_RNA_LIGAND,ALL};


class Scoring_Function {
public:
    const RNA& rna;
    const Ligand& lig;
    SCORE_MODE Score_Mode;
    // VDW_Score vdw_score;
    YW_Score yw_score;
    lz::RL_Score rl_score;
    int num_thread = 0;
    std::vector<thread_range> thread_ranges;
    std::ostream& tee;
public:
    friend YW_Score;
    // friend VDW_Score;
    //initiate
    Scoring_Function(const RNA& r, const Ligand& l, const SCORE_MODE& sm, std::ostream& lg) : rna(r), lig(l), Score_Mode(sm), tee(lg), yw_score(lg), rl_score(lg) {
    //thread information
        this->tee << "rna atom num: " << r.atom_num() << " ligand atom num: " << lig.atom_num() << std::endl;
        this->num_thread = std::thread::hardware_concurrency();
        if(this->num_thread>0) {
            this->tee << "std::thread::hardware_concurrency() return <=0, so using 1 thread!" << std::endl;
            this->num_thread = 1;
            thread_range tmp_tr;
            tmp_tr.begin = 0;
            tmp_tr.end = lig.atom_num();
            this->thread_ranges.push_back(tmp_tr);
        }
        // else {
        //     assert(false);
        //     int num_per_thread = -1;
        //     if(this->lig.atom_num()<this->num_thread) {
        //         this->tee << "this->lig.atom_num() is smaller than num thread, using lig.atom_num() as thread num!" << std::endl;
        //         this->num_thread = this->lig.atom_num();
        //         num_per_thread = 1;
        //         for(auto i = 0; i < this->num_thread; ++i) {
        //             thread_range tmp_tr;
        //             tmp_tr.begin = i;
        //             tmp_tr.end = i+1;
        //             this->thread_ranges.push_back(tmp_tr);
        //         }
        //     } else {
        //         num_per_thread = this->lig.atom_num()/this->num_thread + 1;
        //         int total_atom_num_count = 0;
        //         while(total_atom_num_count < this->lig.atom_num()) {
        //             thread_range tmp_tr;
        //             tmp_tr.begin = total_atom_num_count;
        //             total_atom_num_count += num_per_thread;
        //             if(total_atom_num_count > this->lig.atom_num()) total_atom_num_count = this->lig.atom_num();
        //             tmp_tr.end = total_atom_num_count;
        //             this->thread_ranges.push_back(tmp_tr);
        //         }
        //         this->num_thread = this->thread_ranges.size();
        //     }
        // }
        assert(this->thread_ranges.size()==this->num_thread);
        this->tee << "thread num: " << this->thread_ranges.size() << " ligand atom num: " << this->lig.atom_num() << std::endl;
        /////////////////////////////////////////////////////
        switch(this->Score_Mode) {
            case YW_SCORE:
                this->yw_score.init(r,l);
                break;
            case RL_SCORE:
                this->rl_score.init(r,l);
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

    const Float evaluate() {

        switch(this->Score_Mode) {
            case YW_SCORE:
                return this->yw_score.evaluate(this->rna,this->lig);
                break;
            case RL_SCORE:
                return this->rl_score.evaluate(this->rna,this->lig);
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

        // const int num_thread = tr.size();
        // std::thread t[num_thread];
        // std::vector<Float> score_thread(num_thread,0);
        // for(int i = 0; i != num_thread; ++i) {
        //     t[i] = std::thread(yw_score_calculation, tr[i].begin, tr[i].end, std::ref(rna), std::ref(lig), std::ref(grids), std::ref(*this), std::ref(score_thread[i]));
        // }
        // for(int i = 0; i != num_thread; ++i) {
        //     t[i].join();
        // }
        // return std::accumulate(score_thread.begin(), score_thread.end(), 0.0)/static_cast<Float>(lig.heavy_atom_num());
    }
};

}
