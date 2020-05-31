// 我能量还得重写，我之前是把能量赋予格点就想下面这个程序算LJ那样的，我发现和正经算出入还是
// 有点大。不过下面这个快是真的快。
// 我准备把rna分割到一个相对大的cell，每次遍历cell与ligand上的atom

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <thread>
#include <algorithm>
#include <time.h>
#include <functional>
#include "atom.h"
#include "molecule.h"
#include "grid.h"
#include "matrix.h"

namespace lz {

using namespace std;

inline void print_vector(const std::string s, const std::vector<double>& v) {
    std::cout << s << ":" << std::endl;
    std::cout << "****";
    for(const auto& a : v) {
        std::cout << a << " ";
    }
    std::cout << "****" << std::endl;
}

// struct thread_parameter
// {
//     int ibegin = 0;
//     int iend;
// };

// struct struct_MC_thread_info
// {
//     int seed1 = 25245822;
//     int seed2 = 25311542;
//     int sample_num = 0;
// };

struct atom_info
{
    int index;
    string name;
    double x;
    double y;
    double z;
    string sybyl_type;
    string resname;
    // int reseq;
    double charge;
    string element;
    double radius;
    double born_scale;
    int type_idx;
    int numH = 0;
    int yz_atom_index;
};

struct constant {
    const double radius_O = 1.5;
    const double radius_P = 1.9;
    const double radius_H = 1.0;
    const double radius_C = 1.7;
    const double radius_N = 1.65;
    const double radius_S = 1.80;
    const double born_scale_O = 0.85;
    const double born_scale_P = 0.86;
    const double born_scale_H = 0.85;
    const double born_scale_C = 0.72;
    const double born_scale_N = 0.79;
    const double born_scale_S = 0.80;
} CC;

class case_info
{
public:
    std::vector< atom_info > noH;
    std::vector< atom_info > wH;
    std::vector< atom_info > HH;
    void SetValue (const constant& CC, const tmd::Atoms& atoms)
    {
        int index=0,index_H=0;
        for(const tmd::Atom& atom : atoms) {
            atom_info a;
            a.index = index_H;
            a.yz_atom_index = index_H;
            index_H++;
            a.name = atom.get_name();
            a.x = atom.get_coord()[0];
            a.y = atom.get_coord()[1];
            a.z = atom.get_coord()[2];
            a.sybyl_type = atom.get_sybyl_type_name();
            a.resname = atom.get_res_name();
            a.charge = atom.get_charge();

            if(a.sybyl_type.find(".")!=std::string::npos)
            {
                int position;
                position = a.sybyl_type.find(".");
                a.element = a.sybyl_type.substr(0,position);
            }
            else
            {
                a.element = a.sybyl_type;
            }

            if(a.element == "H" || a.element == "h")
            {
                a.radius = CC.radius_H;
                a.born_scale = CC.born_scale_H;
                a.type_idx = 0;
            }
            else if(a.element == "C" || a.element == "c")
            {
                a.radius = CC.radius_C;
                a.born_scale = CC.born_scale_C;
                a.type_idx = 1;
            }
            else if(a.element == "N" || a.element == "n")
            {
                a.radius = CC.radius_N;
                a.born_scale = CC.born_scale_N;
                a.type_idx = 2;
            }
            else if(a.element == "O" || a.element == "o")
            {
                a.radius = CC.radius_O;
                a.born_scale = CC.born_scale_O;
                a.type_idx = 3;
            }
            else if(a.element == "P" || a.element == "p")
            {
                a.radius = CC.radius_P;
                a.born_scale = CC.born_scale_P;
                a.type_idx = 4;
            }
            else if(a.element == "S" || a.element == "s")
            {
                a.radius = CC.radius_S;
                a.born_scale = CC.born_scale_S;
                a.type_idx = 5;
            }
            else
            {
                a.radius = CC.radius_C;
                a.born_scale = CC.born_scale_C;
                a.type_idx = 6;
            }
            wH.push_back(a);

            if(a.element != "H" && a.element != "h")
            {
                a.index = index;
                index++;
                noH.push_back(a);
                // xcen+=a.x;
                // ycen+=a.y;
                // zcen+=a.z;
            }
            else
            {
                HH.push_back(a);
            }
        }
        // xcen/=(1.0*noH.size());
        // ycen/=(1.0*noH.size());
        // zcen/=(1.0*noH.size());
    }

    void H_charge_transfer()
    {
        for(auto& a : HH)
        {
            int idx;
            double dmin = 10000000;
            for(auto& b : noH)
            {
                double dis = (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z);
                if(dis<dmin) {
                    dmin=dis;
                    idx=b.index;
                }
            }
            noH[idx].charge += a.charge;
            noH[idx].numH++;
        }
    };

    void add_case(const case_info& ci) {
        const int old_noH_size = this->noH.size();
        this->noH.resize(old_noH_size+ci.noH.size());
        for(int i = old_noH_size; i < old_noH_size+ci.noH.size(); ++i) {
            this->noH[i] = ci.noH[i-old_noH_size];
        }
        const int old_wH_size = this->wH.size();
        this->wH.resize(old_wH_size+ci.wH.size());
        for(int i = old_wH_size; i < old_wH_size+ci.wH.size(); ++i) {
            this->wH[i] = ci.wH[i-old_wH_size];
        }
        const int old_HH_size = this->HH.size();
        this->HH.resize(old_HH_size+ci.HH.size());
        for(int i = old_HH_size; i < old_HH_size+ci.HH.size(); ++i) {
            this->HH[i] = ci.HH[i-old_HH_size];
        }
    }
    // std::vector< atom_info > center0;
    // std::vector< atom_info > center;
    // void Move_to_Center_Easy()
    // {
    //     center = noH;
    //     xcen=0,ycen=0,zcen=0;
    //     for(auto&a:noH){xcen+=a.x;ycen+=a.y;zcen+=a.z;}
    //     xcen/=1.0*noH.size();
    //     ycen/=1.0*noH.size();
    //     zcen/=1.0*noH.size();
    //     for(auto&a:center){a.x-=xcen;a.y-=ycen;a.z-=zcen;}
    // }

    void Find_Max_and_Min()
    {
        xmax=-1000000.0;
        for(auto&a:noH){if(xmax<a.x)xmax=a.x;}
        ymax=-1000000.0;
        for(auto&a:noH){if(ymax<a.y)ymax=a.y;}
        zmax=-1000000.0;
        for(auto&a:noH){if(zmax<a.z)zmax=a.z;}
        xmin=1000000.0;
        for(auto&a:noH){if(xmin>a.x)xmin=a.x;}
        ymin=1000000.0;
        for(auto&a:noH){if(ymin>a.y)ymin=a.y;}
        zmin=1000000.0;
        for(auto&a:noH){if(zmin>a.z)zmin=a.z;}
    }

    const double x_max() const {return xmax;}
    const double y_max() const {return ymax;}
    const double z_max() const {return zmax;}
    const double x_min() const {return xmin;}
    const double y_min() const {return ymin;}
    const double z_min() const {return zmin;}

private:
    // double xcen=0;
    // double ycen=0;
    // double zcen=0;

    double xmax;
    double ymax;
    double zmax;

    double xmin;
    double ymin;
    double zmin;
};

struct sasa_polar {
    // std::vector<double> phis;
    // std::vector<double> sitas;
    double dphi;
    double dsita;
    double phi;
    double sita;
    double sin_sita;
    double cos_sita;
    double sin_phi;
    double cos_phi;
    double Rx;
    double Ry;
    double Rz;
    // std::vector<double> sin_sitas;
    // std::vector<double> sin_phis;
    // std::vector<double> cos_phis;
    // std::vector<double> cos_sitas;
    // std::vector<double> Rxs;
    // std::vector<double> Rys;
    // std::vector<double> Rzs;

    // sasa_polar(int i) : phis(i), sitas(i), sin_sitas(i), sin_phis(i), cos_phis(i), cos_sitas(i), Rxs(i), Rys(i), Rzs(i) {
    //     assert(i>=0);
    // }
};


class parameter
{
public:
    // parameter() {
    //     std::cout << "parameter init" << std::endl;
    // }
    // int num_threads = 1;
    const double pi = 3.1415926;
    // int num_atom_type = 9;
    const double radius_small_sphere = 1.5;
    const double radius_large_sphere = 6.0;
    const int LJstep = 1000;
    const double LJ_cut = 2.50;
    const double equ_LJ = 0.8;
    // double excludLJ = 0.6;
    // double Num_in_sphere = 1000;
    // int Nats_RNA = 5500;
    // int Nats_ligand = 200;
    // int Nats_RNA_addH = 11000;
    // int Nats_ligand_addH = 400;
    const double gamma_SASA = 0.005;
    const double step_SASA = 0.25;
    const double radius_water = 1.4;
    const double Hbond_max = 1.3;
    const double Hbond_min = 0.8;
    // double step_pocket = 0.5;
    // int MaxN = 50;
    // int Maxpose = 20000;
    // double gamma_step1 = 10.0;
    // double Max_pocket_site = 10000;
    // double rmsd_cut_step1 = 0.5;
    // double ligand_cut_step1 = 2.0;
    // int TOP_step1 = 2;

    // double LJ_cut_new =8.0;
    // double lj_dr = LJ_cut_new/LJstep;
    // double iseed1 = 25478553;
    // double iseed2 = 25212112;
    const double Tc = 25.0;
    const double e1 = 20.0;
    const double radius_O = 1.5;
    const double radius_P = 1.9;
    const double radius_H = 1.0;
    const double radius_C = 1.7;
    const double radius_N = 1.65;
    const double radius_S = 1.80;
    const std::vector<double> radius_type = {radius_H, radius_C, radius_N, radius_O, radius_P, radius_S, radius_C, radius_small_sphere, radius_large_sphere};
    const double born_scale_O = 0.85;
    const double born_scale_P = 0.86;
    const double born_scale_H = 0.85;
    const double born_scale_C = 0.72;
    const double born_scale_N = 0.79;
    const double born_scale_S = 0.80;
    const std::vector<double> born_type = {born_scale_H, born_scale_C, born_scale_N, born_scale_O, born_scale_P, born_scale_S, born_scale_C};
    const std::vector<double> radius_with_water = {radius_H+radius_water, radius_C+radius_water, radius_N+radius_water, radius_O+radius_water, radius_P+radius_water, radius_S+radius_water, radius_C+radius_water, radius_small_sphere+radius_water, radius_large_sphere+radius_water};

    const double e2 = (87.740-0.4008*Tc+9.398*1e-4*Tc*Tc-1.41*1e-6*Tc*Tc*Tc);   // the dielctric consant of water in temperature T
    const double e25 = 87.740-0.4008*25+9.398*1e-4*25*25-1.41*1e-6*25*25*25;    // the dielctric consant of water in temperature 25
    const double lB = 7.15*(273+25)*e25/((273+Tc)*e2);   // e^2/(ebs*kB*Tc)
    const double lB0=e2*lB;        // in A lb0=e^2/(4*pi*e0*kB*T)
    const double rate_kcal_to_kt=0.593*(273+Tc)/(273+25);


    case_info lig;
    case_info rna;
    case_info complex;
    int rna_heavy_size;
    int lig_heavy_size;
    int complex_heavy_size;
    int rna_wH_size;
    int rna_HH_size;
    int lig_wH_size;
    int lig_HH_size;
    int complex_wH_size;
    int complex_HH_size;

    std::vector<double> re_weight = {3.300,1.320,0.780,2.520,0.920,0.120,4.960};


    double LJenergy = 0.;
    double ligand_LJenergy = 0.0;
    // double ligand_ELEenergy = 0.0;
    double ELEenergy = 0.;
    double HBenergy = 0.;
    double POLenergy = 0.0;
    double RNA_POLenergy = 0.;
    double ligand_POLenergy = 0.;
    double complex_POLenergy = 0.;
    double RNA_SELFenergy = 0.;
    double ligand_SELFenergy = 0.;
    double SASAenergy = 0.0;
    double sasa = 0.0;

    // tmd::Strictly_Triangular_Matrix<tmd::DISTANCE_TYPE> dis_type_mat;
    tmd::Strictly_Triangular_Matrix<tmd::BOND_TYPE> bond_type_mat;
    tmd::Grids sasa_grids;
    tmd::Grids lj_grids;

    tmd::Strictly_Triangular_Matrix<double> ele_prefix_mat;
    tmd::Triangular_Matrix<double> pol_prefix_mat;

    std::vector<std::vector<sasa_polar>> sasa_polar_table;

    vector<vector<vector<double>>> LJ_table;
    vector<vector<vector<double>>> SASA_table;
    vector<vector<vector<double>>> Hbond_table;
    vector<double> rna_born_radius_table;
    vector<double> lig_born_radius_table;
    vector<double> complex_born_radius_table;
    // vector<double> complex_fixed_born_part_table;//fixed part comes from rna rna atom pairs
    // vector<double> complex_after_born_radius_table;
    vector<double> rna_fixed_born_part_table;
    // vector<double> rna_after_born_radius_table;
    // vector<double> lig_after_born_radius_table;
    // tmd::Strictly_Triangular_Matrix<double> rna_dis_mat;
    // tmd::Strictly_Triangular_Matrix<double> lig_dis_mat;
    tmd::Strictly_Triangular_Matrix<double> complex_dis_mat;
    tmd::Strictly_Triangular_Matrix<double> complex_dis_square_mat;
    // double SumSASA_ligand;
    // double SumSASA_RNA;
    // double SumSASA_complex;

    // void print_lig_after_born_radius_table() {
    //     print_vector("lig_after_born_radius_table",this->lig_after_born_radius_table);
    // }
    // void print_rna_after_born_radius_table() {
    //     print_vector("rna_after_born_radius_table",this->rna_after_born_radius_table);
    // }
    void print_lig_born_radius_table() {
        print_vector("lig_born_radius_table",this->lig_born_radius_table);
    }void print_rna_born_radius_table() {
        print_vector("rna_born_radius_table",this->rna_born_radius_table);
    }
    // void print_complex_born_radius_table() {
    //     print_vector("complex_born_radius_table",this->complex_born_radius_table);
    // }
};

inline const double max(const double a, const double b)
{
    if (a > b)
        return (a);
    if (a < b)
        return (b);
    if (a == b)
        return (a);
}
inline const double min(const double a, const double b)
{
    if (a > b)
        return (b);
    if (a < b)
        return (a);
    if (a == b)
        return (b);
}

// enum RL_SCORE_MODE {RL_LJ_ONLY,RL_ALL};
class RL_Score {
    // std::vector<lz::pocket_info> pockets;
    parameter P;
    std::ostream& tee;
    // RL_SCORE_MODE RL_Score_Mode = RL_ALL;
public:
    RL_Score(std::ostream& lg) : tee(lg) {
        tee << "RL_Score init" << std::endl;
    }
    void init(const tmd::RNA& rna, const tmd::Ligand& lig);
    const double evaluate(const tmd::RNA& rna, const tmd::Ligand& lig);

    // void set_mode(const RL_SCORE_MODE& rl_sm) {
    //     this->RL_Score_Mode = rl_sm;
    // }
};


inline void RL_Score::init(const tmd::RNA& yz_rna, const tmd::Ligand& yz_lig) {
    // string filename = "input/"+argv_str[1]+"/"+argv_str[1]+"_ligand.mol2";
    P.lig.SetValue(CC,yz_lig.get_atoms_reference());
    if(P.lig.HH.size()>0) P.lig.H_charge_transfer();
    // lig.Move_to_Center_Easy();
    // cout<<"native Ligand input finish "<<lig.wH.size()<<" "<<P.lig_heavy_size<<endl;

    // filename = "input/"+argv_str[1]+"/"+argv_str[1]+"_RNA.mol2";
    P.rna.SetValue(CC,yz_rna.get_atoms_reference());
    P.rna.Find_Max_and_Min();
    if(P.rna.HH.size()>0) P.rna.H_charge_transfer();
    // cout<<"Receptor input finish "<<rna.wH.size()<<" "<<P.rna_heavy_size<<endl;

    P.complex = P.rna;
    P.complex.add_case(P.lig);

    P.rna_heavy_size = P.rna.noH.size();
    P.lig_heavy_size = P.lig.noH.size();
    P.complex_heavy_size = P.complex.noH.size();
    P.rna_wH_size = P.rna.wH.size();
    P.rna_HH_size = P.rna.HH.size();
    P.lig_wH_size = P.lig.wH.size();
    P.lig_HH_size = P.lig.HH.size();
    P.complex_wH_size = P.complex.wH.size();
    P.complex_HH_size = P.complex.HH.size();

    //cal LJ table
    // std::cout << "start init LJ" << std::endl;
    vector<double> tmp_lj(P.LJstep+1,0.0);
    vector<vector<double>> itmp_lj(P.radius_type.size(),tmp_lj);
    P.LJ_table = vector<vector<vector<double>>>(P.radius_type.size(),itmp_lj);
    ////////LJ index start from 1 in calculation 0 can not be divined.
    for(int i=0;i!=P.radius_type.size();i++)
    {
        for(int j=0;j!=P.radius_type.size();j++)
        {
            const double rdis = P.equ_LJ*(P.radius_type[i] + P.radius_type[j]);
            const double dr = P.LJ_cut*(P.radius_type[i] + P.radius_type[j])/static_cast<tmd::Float>(P.LJstep);
            for(int istep=1; istep<=P.LJstep; istep++)
            {
                double dis = istep*dr;
                double dx = rdis/dis;
                double x6 = dx*dx*dx*dx*dx*dx;
                P.LJ_table[i][j][istep] = x6*x6 - x6;
            }
        }
    }

    P.sasa_polar_table = vector<vector<sasa_polar>>(P.radius_type.size());
    for (int i=0;i!=P.radius_type.size();i++) {
        const double& r = P.radius_type[i] + P.radius_water;
        const double& dsita = P.step_SASA / r;
        const int& num_sita = static_cast<int>(P.pi / dsita);
        for (int si = 1; si <= num_sita; si++) {
            const double& sita = si * dsita;
            const double& sin_sita = sin(sita);
            const double& dphi = P.step_SASA / r * sin_sita;
            const int& num_phi = static_cast<int>(2. * P.pi / dphi + 1);
            for (int pk = 1; pk <= num_phi; pk++) {
                const double& phi = pk * dphi;
                const double& cos_sita = cos(sita);
                const double& sin_phi = sin(phi);
                const double& cos_phi = cos(phi);
                const double& Rx = r * sin(sita) * cos(phi);
                const double& Ry = r * sin(sita) * sin(phi);
                const double& Rz = r * cos(sita);
                sasa_polar sa;
                sa.cos_phi = cos_phi;
                sa.cos_sita = cos_sita;
                sa.dphi = dphi;
                sa.dsita = dsita;
                sa.phi = phi;
                sa.Rx = Rx;
                sa.Ry = Ry;
                sa.Rz = Rz;
                sa.sin_phi = sin_phi;
                sa.sin_sita = sin_sita;
                sa.sita = sita;
                P.sasa_polar_table[i].push_back(sa);
            }
        }
    }

    //cal Hbond table
    // std::cout << "start init Hbond" << std::endl;
    vector<double> tmp_hbond(P.LJstep+1,0.0);
    vector<vector<double>> itmp_bond(P.radius_type.size(),tmp_hbond);
    P.Hbond_table = vector<vector<vector<double>>>(P.radius_type.size(),itmp_bond);

    for (int i = 0; i < P.radius_type.size(); i++) {
        for (int j = 0; j < P.radius_type.size(); j++) {
            const double dismax = P.Hbond_max * (P.radius_type[i] + P.radius_type[j]);
            const double dismin = P.Hbond_min * (P.radius_type[i] + P.radius_type[j]);
            const double dr = dismax / static_cast<double>(P.LJstep);
            for(int istep = 1; istep <= P.LJstep; istep++)
            {
                const double dis = dr * static_cast<double>(istep);
                if (dis > dismin)
                {
                    P.Hbond_table[i][j][istep] = 1. - (dis - dismin) / (dismax - dismin);
                }
                else
                {
                    P.Hbond_table[i][j][istep] = 1.;
                }
            }
        }
    }

    //initialize bond type for ligand atoms
    // P.dis_type_mat.resize(P.complex_heavy_size,tmd::DISTANCE_FLEXIBLE);
    P.bond_type_mat.resize(P.complex_heavy_size,tmd::NONE_BOND);

    for (int i = 0; i < P.lig_heavy_size; i++)
    {
        for (int j = 0; j < P.lig_heavy_size; j++)
        {
            int real_i = i + P.rna_heavy_size;
            int real_j = j + P.rna_heavy_size;
            if(real_i < real_j) {
                const tmd::Atom& ai = yz_lig.get_atoms_reference()[P.complex.noH[real_i].yz_atom_index];
                const int& aj = P.complex.noH[real_j].yz_atom_index;
                for(const tmd::Bond& b : ai.get_bonds()) {
                    if(b.get_bonded_atom_index() == aj && b.get_bond_type() != tmd::NONE_BOND && b.get_bond_type() != tmd::NOT_CONNECTED_BOND) {
                        P.bond_type_mat(real_i,real_j) = tmd::DUMMY_BOND;
                        break;
                    }
                }
            }
        }
    }

    //init complex_dis_mat to 0.0
    P.complex_dis_mat.resize(P.complex_heavy_size,0.0);
    P.complex_dis_square_mat.resize(P.complex_heavy_size,0.0);
    //init born radius for rna
    ///////////////////////////////////////////////
    /// RNA //////////////////////////////////////
    P.rna_born_radius_table = vector<double>(P.rna_heavy_size,0.0);
    P.rna_fixed_born_part_table = vector<double>(P.rna_heavy_size,0.0);
    for (int i = 0; i < P.rna_heavy_size; i++)
    { // start to calculate rb0
        double Lij, Uij;
        double intb = 0;
        double radi = P.rna.noH[i].radius;
        for (int j = 0; j < P.rna_heavy_size; j++)
        {
            if (j!= i)
            {
                double dis = (P.rna.noH[j].x - P.rna.noH[i].x) * (P.rna.noH[j].x - P.rna.noH[i].x) +
                      (P.rna.noH[j].y - P.rna.noH[i].y) * (P.rna.noH[j].y - P.rna.noH[i].y) +
                      (P.rna.noH[j].z - P.rna.noH[i].z) * (P.rna.noH[j].z - P.rna.noH[i].z);
                dis = sqrt(dis);
                double radj = P.rna.noH[j].radius;
                double radc = dis + P.rna.noH[j].born_scale * radj;
                if (radi >= radc)
                {
                    Lij = 1.;
                    Uij = 1.;
                }
                else
                {
                    if (radi > (dis - P.rna.noH[j].born_scale * radj))
                    {
                        Lij = radi;
                    }
                    else
                    {
                        Lij = dis - P.rna.noH[j].born_scale * radj;
                    }
                    Uij = dis + P.rna.noH[j].born_scale * radj;
                }
                intb = intb + 0.5 * ((1. / Lij - 1. / Uij) + (P.rna.noH[j].born_scale * P.rna.noH[j].born_scale * radj * radj / (4. * dis) - dis / 4.) * (1. / (Lij * Lij) - 1. / (Uij * Uij)) + (1. / (2. * dis)) * log(Lij / Uij));
            }
        }
        P.rna_fixed_born_part_table[i] = intb;
        P.rna_born_radius_table[i] = 1. / (1. / radi - intb);
        //  printf("%d %lf %lf\n",i,rb0[i],rrr_RNA[i]);
    }

    ////////////////////////////////////////////////////////
    //// RNA alone ////////////////////////////////////////
    ///////////////////////////////////////////////////////
    // POL
    // std::cout << "cal RNA_POLenergy" << std::endl;
    for (int i = 0; i < P.rna_heavy_size; i++)
    {
        for (int j = 0; j < P.rna_heavy_size; j++)
        {
            if (i < j)
            {
                double dis =
                (P.rna.noH[i].x - P.rna.noH[j].x) * (P.rna.noH[i].x - P.rna.noH[j].x) +
                (P.rna.noH[i].y - P.rna.noH[j].y) * (P.rna.noH[i].y - P.rna.noH[j].y) +
                (P.rna.noH[i].z - P.rna.noH[j].z) * (P.rna.noH[i].z - P.rna.noH[j].z);
                dis = sqrt(dis);
                P.RNA_POLenergy = P.RNA_POLenergy + P.lB0 * (1. / P.e2 - 1. / P.e1) * P.rna.noH[i].charge * P.rna.noH[j].charge / sqrt(dis * dis + P.rna_born_radius_table[i] * P.rna_born_radius_table[j] * exp(-dis * dis / (4. * P.rna_born_radius_table[i] * P.rna_born_radius_table[j])));
                // P.rna_dis_mat[i][j] = dis;
                // P.rna_dis_mat[j][i] = dis;
                P.complex_dis_mat(i,j) = dis;
                P.complex_dis_square_mat(i,j) = dis*dis;
            }
        }
    }


    // std::cout << "cal sasa grids" << std::endl;
    struct Grid_Shift {
        const int x;
        const int y;
        const int z;
        Grid_Shift(const int xx, const int yy, const int zz) : x(xx), y(yy), z(zz) {}
    };
    //initalize sasa grids
    const double sasa_interaction_radius = (3.0 + P.radius_water)*2.0;
    const double sasa_grid_width = 1.0;
    //get all the grids shifts within the interacting radius(interation_radius)
    const int sasa_upper_cubic_shift = static_cast<int>(std::ceil((0.0+sasa_interaction_radius)/sasa_grid_width))+1;//plus extra one to ensure it covers the area
    const int sasa_lower_cubic_shift = static_cast<int>(std::floor((0.0-sasa_interaction_radius)/sasa_grid_width))-1;
    std::vector<Grid_Shift> sasa_grid_shifts;
    for(int x = sasa_lower_cubic_shift; x <= sasa_upper_cubic_shift; ++x) {
        for(int y = sasa_lower_cubic_shift; y <= sasa_upper_cubic_shift; ++y) {
            for(int z = sasa_lower_cubic_shift; z <= sasa_upper_cubic_shift; ++z) {
                const tmd::Float& x_center = (static_cast<tmd::Float>(x)+0.5)*sasa_grid_width;
                const tmd::Float& y_center = (static_cast<tmd::Float>(y)+0.5)*sasa_grid_width;
                const tmd::Float& z_center = (static_cast<tmd::Float>(z)+0.5)*sasa_grid_width;
                const tmd::Float& dis = std::sqrt((x_center*x_center+y_center*y_center+z_center*z_center));
                // consider some grids partially covered by sphere
                const tmd::Float& effective_radius = sasa_interaction_radius + sasa_grid_width * std::sqrt(3);
                if(dis<=effective_radius) {
                    sasa_grid_shifts.push_back(Grid_Shift(x,y,z));
                }
            }
        }
    }
    // std::cout << "grid_shifts size -> " << sasa_grid_shifts.size() << std::endl;
    //initialize grid map
    // std::cout << "process rna atom grid_map..." << std::endl;
    std::map<std::string,tmd::Grid> sasa_grid_map;
    std::set<int> sasa_grid_x_set;
    std::set<int> sasa_grid_y_set;
    std::set<int> sasa_grid_z_set;
    int sasa_rna_atom_index = 0;
    assert(sasa_grid_width>tmd::k_epsilon);
    for(const atom_info& ra : P.rna.noH) {
        const int atom_x_grid = static_cast<int>(ra.x/sasa_grid_width);
        const int atom_y_grid = static_cast<int>(ra.y/sasa_grid_width);
        const int atom_z_grid = static_cast<int>(ra.z/sasa_grid_width);
        for(const Grid_Shift& gs : sasa_grid_shifts) {
            //get absolute grid index of this RNA atom
            const int actual_x_grid = atom_x_grid + gs.x;
            const int actual_y_grid = atom_y_grid + gs.y;
            const int actual_z_grid = atom_z_grid + gs.z;

            sasa_grid_x_set.insert(actual_x_grid);
            sasa_grid_y_set.insert(actual_y_grid);
            sasa_grid_z_set.insert(actual_z_grid);

            const std::string& grid_name = std::to_string(actual_x_grid)+"-"+std::to_string(actual_y_grid)+"-"+std::to_string(actual_z_grid);
            if(sasa_grid_map.find(grid_name) != sasa_grid_map.end()) {
                sasa_grid_map.at(grid_name).rigid_atoms.push_back(sasa_rna_atom_index);
            } else {
                tmd::Grid grid;
                grid.flag = true;
                grid.rigid_atoms.push_back(sasa_rna_atom_index);
                sasa_grid_map.insert({grid_name,grid});
            }
        }
        sasa_rna_atom_index++;
    }
    // std::cout << "start insert grids" << std::endl;
    const int sasa_min_grid_x = (*sasa_grid_x_set.begin());
    const int sasa_min_grid_y = (*sasa_grid_y_set.begin());
    const int sasa_min_grid_z = (*sasa_grid_z_set.begin());
    const int sasa_max_grid_x = (*sasa_grid_x_set.end());
    const int sasa_max_grid_y = (*sasa_grid_y_set.end());
    const int sasa_max_grid_z = (*sasa_grid_z_set.end());
    const int sasa_size_x = sasa_max_grid_x - sasa_min_grid_x + 1;
    const int sasa_size_y = sasa_max_grid_y - sasa_min_grid_y + 1;
    const int sasa_size_z = sasa_max_grid_z - sasa_min_grid_z + 1;
    P.sasa_grids = tmd::Grids(sasa_size_x,sasa_size_y,sasa_size_z,sasa_min_grid_x,sasa_min_grid_y,sasa_min_grid_z,sasa_grid_width,tmd::Grid(false));

    for(int i = sasa_min_grid_x; i <= sasa_max_grid_x; ++i) {
        for(int j = sasa_min_grid_y; j <= sasa_max_grid_y; ++j) {
            for(int k = sasa_min_grid_z; k <= sasa_max_grid_z; ++k) {
                const std::string& grid_name = std::to_string(i)+"-"+std::to_string(j)+"-"+std::to_string(k);
                if(sasa_grid_map.find(grid_name) != sasa_grid_map.end()) {
                    P.sasa_grids.assign(i,j,k,sasa_grid_map.at(grid_name));
                }
            }
        }
    }
    // this->tee << "finish initializing grids... size: " << this->grids.dim_1() << " " << this->grids.dim_2() << " " << this->grids.dim_3() << std::endl;


    // std::cout << "cal lj grids" << std::endl;
    //initalize lj grids
    //2.5 should be larger than any atom vdw radius
    const double lj_interaction_radius = P.LJ_cut * 3.0;
    const double lj_grid_width = 1.0;
    //get all the grids shifts within the interacting radius(interation_radius)
    const int lj_upper_cubic_shift = static_cast<int>(std::ceil((0.0+lj_interaction_radius)/lj_grid_width))+1;//plus extra one to ensure it covers the area
    const int lj_lower_cubic_shift = static_cast<int>(std::floor((0.0-lj_interaction_radius)/lj_grid_width))-1;
    std::vector<Grid_Shift> lj_grid_shifts;
    for(int x = lj_lower_cubic_shift; x <= lj_upper_cubic_shift; ++x) {
        for(int y = lj_lower_cubic_shift; y <= lj_upper_cubic_shift; ++y) {
            for(int z = lj_lower_cubic_shift; z <= lj_upper_cubic_shift; ++z) {
                const tmd::Float& x_center = (static_cast<tmd::Float>(x)+0.5)*lj_grid_width;
                const tmd::Float& y_center = (static_cast<tmd::Float>(y)+0.5)*lj_grid_width;
                const tmd::Float& z_center = (static_cast<tmd::Float>(z)+0.5)*lj_grid_width;
                const tmd::Float& dis = std::sqrt((x_center*x_center+y_center*y_center+z_center*z_center));
                // consider some grids partially covered by sphere
                const tmd::Float& effective_radius = lj_interaction_radius + lj_grid_width * std::sqrt(3);
                if(dis<=effective_radius) {
                    lj_grid_shifts.push_back(Grid_Shift(x,y,z));
                }
            }
        }
    }
    // this->tee << "grid_shifts size -> " << grid_shifts.size() << std::endl;
    //initialize grid map
    // this->tee << "process rna atom grid_map..." << std::endl;
    std::map<std::string,tmd::Grid> lj_grid_map;
    std::set<int> lj_grid_x_set;
    std::set<int> lj_grid_y_set;
    std::set<int> lj_grid_z_set;
    int lj_rna_atom_index = 0;
    assert(lj_grid_width>tmd::k_epsilon);
    for(const atom_info& ra : P.rna.noH) {
        const int atom_x_grid = static_cast<int>(ra.x/lj_grid_width);
        const int atom_y_grid = static_cast<int>(ra.y/lj_grid_width);
        const int atom_z_grid = static_cast<int>(ra.z/lj_grid_width);
        for(const Grid_Shift& gs : lj_grid_shifts) {
            //get absolute grid index of this RNA atom
            const int actual_x_grid = atom_x_grid + gs.x;
            const int actual_y_grid = atom_y_grid + gs.y;
            const int actual_z_grid = atom_z_grid + gs.z;

            lj_grid_x_set.insert(actual_x_grid);
            lj_grid_y_set.insert(actual_y_grid);
            lj_grid_z_set.insert(actual_z_grid);

            const std::string& grid_name = std::to_string(actual_x_grid)+"-"+std::to_string(actual_y_grid)+"-"+std::to_string(actual_z_grid);
            if(lj_grid_map.find(grid_name) != lj_grid_map.end()) {
                lj_grid_map.at(grid_name).rigid_atoms.push_back(lj_rna_atom_index);
            } else {
                tmd::Grid grid;
                grid.flag = true;
                grid.rigid_atoms.push_back(lj_rna_atom_index);
                lj_grid_map.insert({grid_name,grid});
            }
        }
        lj_rna_atom_index++;
    }
    const int lj_min_grid_x = (*lj_grid_x_set.begin());
    const int lj_min_grid_y = (*lj_grid_y_set.begin());
    const int lj_min_grid_z = (*lj_grid_z_set.begin());
    const int lj_max_grid_x = (*lj_grid_x_set.end());
    const int lj_max_grid_y = (*lj_grid_y_set.end());
    const int lj_max_grid_z = (*lj_grid_z_set.end());
    const int lj_size_x = lj_max_grid_x - lj_min_grid_x + 1;
    const int lj_size_y = lj_max_grid_y - lj_min_grid_y + 1;
    const int lj_size_z = lj_max_grid_z - lj_min_grid_z + 1;
    P.lj_grids = tmd::Grids(lj_size_x,lj_size_y,lj_size_z,lj_min_grid_x,lj_min_grid_y,lj_min_grid_z,lj_grid_width,tmd::Grid(false));

    for(int i = lj_min_grid_x; i <= lj_max_grid_x; ++i) {
        for(int j = lj_min_grid_y; j <= lj_max_grid_y; ++j) {
            for(int k = lj_min_grid_z; k <= lj_max_grid_z; ++k) {
                const std::string& grid_name = std::to_string(i)+"-"+std::to_string(j)+"-"+std::to_string(k);
                if(lj_grid_map.find(grid_name) != lj_grid_map.end()) {
                    P.lj_grids.assign(i,j,k,lj_grid_map.at(grid_name));
                }
            }
        }
    }
    // this->tee << "finish initializing grids... size: " << this->grids.dim_1() << " " << this->grids.dim_2() << " " << this->grids.dim_3() << std::endl;

}

inline const double RL_Score::evaluate(const tmd::RNA& yz_rna, const tmd::Ligand& yz_lig)
{
    // clock_t start,mid1,mid2,end;
    // start=clock();
    // mid1=clock();
    // end = clock();

    //update lig's conformation and complex
    int index_noH=0,index_wH=0,index_HH=0;
    const tmd::Atoms& yz_lig_atoms_ref = yz_lig.get_atoms_reference();
    // std::cout << yz_lig_atoms_ref.size() << " " << P.lig.wH.size() << std::endl;
    assert(yz_lig_atoms_ref.size()==P.lig.wH.size());
    for(const tmd::Atom& atom : yz_lig.get_atoms_reference()) {
        const tmd::Vec3d& xyz = atom.get_coord();
        if(P.lig.wH[index_wH].element != "H" && P.lig.wH[index_wH].element != "h")
        {
            P.lig.noH[index_noH].x = xyz[0];
            P.lig.noH[index_noH].y = xyz[1];
            P.lig.noH[index_noH].z = xyz[2];
            P.complex.noH[index_noH+P.rna_heavy_size].x = xyz[0];
            P.complex.noH[index_noH+P.rna_heavy_size].y = xyz[1];
            P.complex.noH[index_noH+P.rna_heavy_size].z = xyz[2];
            index_noH++;
        }
        else
        {
            P.lig.HH[index_HH].x = xyz[0];
            P.lig.HH[index_HH].y = xyz[1];
            P.lig.HH[index_HH].z = xyz[2];
            P.complex.HH[index_HH+P.rna_HH_size].x = xyz[0];
            P.complex.HH[index_HH+P.rna_HH_size].y = xyz[1];
            P.complex.HH[index_HH+P.rna_HH_size].z = xyz[2];
            index_HH++;
        }
        P.lig.wH[index_wH].x = xyz[0];
        P.lig.wH[index_wH].y = xyz[1];
        P.lig.wH[index_wH].z = xyz[2];
        P.complex.wH[index_wH+P.rna_wH_size].x = xyz[0];
        P.complex.wH[index_wH+P.rna_wH_size].y = xyz[1];
        P.complex.wH[index_wH+P.rna_wH_size].z = xyz[2];
        index_wH++;
    }

    // std::cout << "complete dis mat and update dis mat" << std::endl;

    //complete dis mat and update dis mat
    for (int i = 0; i < P.complex_heavy_size; i++)
    {
        for (int j = P.rna_heavy_size; j < P.complex_heavy_size; j++)
        {
            if (i < j)
            {
                double dis_square =
                    (P.complex.noH[j].x - P.complex.noH[i].x) * (P.complex.noH[j].x - P.complex.noH[i].x) +
                    (P.complex.noH[j].y - P.complex.noH[i].y) * (P.complex.noH[j].y - P.complex.noH[i].y) +
                    (P.complex.noH[j].z - P.complex.noH[i].z) * (P.complex.noH[j].z - P.complex.noH[i].z);

                if(tmd::eq(dis_square,0)) {
                    dis_square = 0.0000000000000001;
                }
                P.complex_dis_square_mat(i,j) = dis_square;
                P.complex_dis_mat(i,j) = sqrt(dis_square);
            }
        }
    }

    // if(this->RL_Score_Mode == RL_ALL) {
        // std::cout << "update born raidus for ligand every time" << std::endl;
        //update born raidus for ligand every time
        ///////////////////////////////////////////////
        /// ligand //////////////////////////////////////
        P.lig_born_radius_table = vector<double>(P.lig_heavy_size,0.0);
        for (int i = 0; i < P.lig_heavy_size; i++)
        { // start to calculate rb0
            double Lij, Uij, dis;
            double intb = 0;
            double radi = P.lig.noH[i].radius;
            for (int j = 0; j < P.lig_heavy_size; j++)
            {
                if (j!= i)
                {
                    if(i < j) {
                        dis = P.complex_dis_mat(i+P.rna_heavy_size,j+P.rna_heavy_size);
                    } else {
                        dis = P.complex_dis_mat(j+P.rna_heavy_size,i+P.rna_heavy_size);
                    }
                    double radj = P.lig.noH[j].radius;
                    double radc = dis + P.lig.noH[j].born_scale * radj;
                    if (radi >= radc)
                    {
                        Lij = 1.;
                        Uij = 1.;
                    }
                    else
                    {
                        if (radi > (dis - P.lig.noH[j].born_scale * radj))
                        {
                            Lij = radi;
                        }
                        else
                        {
                            Lij = dis - P.lig.noH[j].born_scale * radj;
                        }
                        Uij = dis + P.lig.noH[j].born_scale * radj;
                    }
                    intb = intb + 0.5 * ((1. / Lij - 1. / Uij) + (P.lig.noH[j].born_scale * P.lig.noH[j].born_scale * radj * radj / (4. * dis) - dis / 4.) * (1. / (Lij * Lij) - 1. / (Uij * Uij)) + (1. / (2. * dis)) * log(Lij / Uij));
                }
            }
            P.lig_born_radius_table[i] = 1. / (1. / radi - intb);
            //  printf("%d %lf %lf\n",i,rb0_ligand[i],rrr_ligand[i]);
        }

        //update complex born radius every time
        // std::cout << "start Born_Radius_Complex" << std::endl;
        P.complex_born_radius_table  = vector<double>(P.complex_heavy_size,0.0);
        for (int i = 0; i < P.complex_heavy_size; i++)
        { // start to calculate rb0
            double Lij, Uij, dis;
            double intb = 0;
            double radi = P.complex.noH[i].radius;
            for (int j = 0; j < P.complex_heavy_size; j++)
            {
                if(i < P.rna_heavy_size && j < P.rna_heavy_size) {
                    continue;
                }
                if (j!= i)
                {
                    if(i < j) {
                        dis = P.complex_dis_mat(i,j);
                    } else {
                        dis = P.complex_dis_mat(j,i);
                    }
                    double radj = P.complex.noH[j].radius;
                    double radc = dis + P.complex.noH[j].born_scale * radj;
                    if (radi >= radc)
                    {
                        Lij = 1.;
                        Uij = 1.;
                    }
                    else
                    {
                        if (radi > (dis - P.complex.noH[j].born_scale * radj))
                        {
                            Lij = radi;
                        }
                        else
                        {
                            Lij = dis - P.complex.noH[j].born_scale * radj;
                        }
                        Uij = dis + P.complex.noH[j].born_scale * radj;
                    }
                    intb = intb + 0.5 * ((1. / Lij - 1. / Uij) + (P.complex.noH[j].born_scale * P.complex.noH[j].born_scale * radj * radj / (4. * dis) - dis / 4.) * (1. / (Lij * Lij) - 1. / (Uij * Uij)) + (1. / (2. * dis)) * log(Lij / Uij));
                }
                if(isnan(intb)) {
                    std::cout << i << " " << j << std::endl;
                    std::cout << P.complex.noH[i].born_scale << std::endl;
                    std::cout << P.complex.noH[i].radius << std::endl;
                    std::cout << P.complex.noH[j].born_scale << std::endl;
                    std::cout << P.complex.noH[j].radius << std::endl;
                    if(i < j) {
                        std::cout << P.complex_dis_mat(i,j) << std::endl;
                    } else {
                        std::cout << P.complex_dis_mat(j,i) << std::endl;
                    }
                    std::cout << Lij << std::endl;
                    std::cout << Uij << std::endl;
                    assert(false);
                }
            }
            if(i < P.rna_heavy_size) {
                P.complex_born_radius_table[i] = 1. / (1. / radi - intb - P.rna_fixed_born_part_table[i]);
            } else {
                P.complex_born_radius_table[i] = 1. / (1. / radi - intb);
            }
            //  printf("%d %lf %lf\n",i,rb0_complex[i],rrr_complex[i]);
        }
    // }
    // P.rna_after_born_radius_table = vector<double>(P.rna_heavy_size,0.0);
    // P.lig_after_born_radius_table = vector<double>(P.lig_heavy_size,0.0);
    // for (int i = 0; i < P.complex_heavy_size; i++)
    // {
    //     if (i < P.rna_heavy_size)
    //         P.rna_after_born_radius_table[i] = P.complex_born_radius_table[i];
    //     if (i >= P.rna_heavy_size)
    //         P.lig_after_born_radius_table[i - P.rna_heavy_size] = P.complex_born_radius_table[i];
    // }

    // std::cout << "cal LJ" << std::endl;
    // RNA-ligand LJ step 2
    P.LJenergy = 0.0;
    for(int i = 0; i < P.lig_heavy_size; i++) {
        const double& lig_x = P.lig.noH[i].x;
        const double& lig_y = P.lig.noH[i].y;
        const double& lig_z = P.lig.noH[i].z;
        const int atom_x_grid = static_cast<int>(lig_x/P.lj_grids.width);
        const int atom_y_grid = static_cast<int>(lig_y/P.lj_grids.width);
        const int atom_z_grid = static_cast<int>(lig_z/P.lj_grids.width);

        const tmd::Grid& grid = P.lj_grids.at(atom_x_grid,atom_y_grid,atom_z_grid);
        if(grid.flag) {
            for(const int& r_j : grid.rigid_atoms) {
                // const double& rna_x = P.rna.noH[r_j].x;
                // const double& rna_y = P.rna.noH[r_j].y;
                // const double& rna_z = P.rna.noH[r_j].z;
                // const double& dis = sqrt((rna_x-lig_x)*(rna_x-lig_x)+(rna_y-lig_y)*(rna_y-lig_y)+(rna_z-lig_z)*(rna_z-lig_z));
                const double& dis = P.complex_dis_mat(r_j, i+P.rna_heavy_size);
                const double equ_dis = (P.rna.noH[r_j].radius + P.lig.noH[i].radius);
                if (dis < P.LJ_cut * equ_dis) {
                    const int iLJ_dis = (dis / (P.LJ_cut * equ_dis)) * P.LJstep + 1;
                    P.LJenergy = P.LJenergy + P.LJ_table[P.rna.noH[r_j].type_idx][P.lig.noH[i].type_idx][iLJ_dis];
                }
            }
        }
    }
    // for (int i = 0; i < P.lig_heavy_size; i++)
    // {
    //     for (int iatom_RNA = 0; iatom_RNA < P.rna_heavy_size; iatom_RNA++)
    //     {
    //         const double equ_dis = (P.rna.noH[iatom_RNA].radius + P.lig.noH[i].radius);
    //         const double cal_dis = sqrt(
    //             (P.lig.noH[i].x - P.rna.noH[iatom_RNA].x) * (P.lig.noH[i].x - P.rna.noH[iatom_RNA].x)
    //           + (P.lig.noH[i].y - P.rna.noH[iatom_RNA].y) * (P.lig.noH[i].y - P.rna.noH[iatom_RNA].y)
    //           + (P.lig.noH[i].z - P.rna.noH[iatom_RNA].z) * (P.lig.noH[i].z - P.rna.noH[iatom_RNA].z)
    //         );
    //         P.complex_dis_mat(iatom_RNA,i+P.rna_heavy_size) = cal_dis;
    //         // P.complex_dis_mat[i+P.rna_heavy_size][iatom_RNA] = cal_dis;
    //         if (cal_dis < P.LJ_cut * equ_dis)
    //         {
    //             const int iLJ_dis = (cal_dis / (P.LJ_cut * equ_dis)) * P.LJstep + 1;
    //             P.LJenergy = P.LJenergy + P.LJ_table[P.rna.noH[iatom_RNA].type_idx][P.lig.noH[i].type_idx][iLJ_dis];
    //         }
    //     }
    // }


    // P.ligand_ELEenergy = 0.0;
    P.ligand_LJenergy = 0.0;
    P.ligand_POLenergy = 0.0;
    P.ligand_SELFenergy = 0.0;
    for (int i = 0; i < P.lig_heavy_size; i++)
    {
        // if(this->RL_Score_Mode == RL_ALL) {
            P.ligand_SELFenergy = P.ligand_SELFenergy + 0.5 * (1. / P.e2 - 1. / P.e1) * P.lB0 * P.lig.noH[i].charge * P.lig.noH[i].charge * (1. / P.complex_born_radius_table[i+P.rna_heavy_size] - 1. / P.lig_born_radius_table[i]);
        // }
        for (int j = 0; j < P.lig_heavy_size; j++)
        {
            if(i < j) {
                const double& cal_dis = P.complex_dis_mat(i+P.rna_heavy_size,j+P.rna_heavy_size);
                const double& cal_dis_square = P.complex_dis_square_mat(i+P.rna_heavy_size,j+P.rna_heavy_size);
                // if(this->RL_Score_Mode == RL_ALL) {
                    P.ligand_POLenergy = P.ligand_POLenergy + P.lB0 * (1. / P.e2 - 1. / P.e1) * P.lig.noH[i].charge * P.lig.noH[j].charge / sqrt(cal_dis_square + P.lig_born_radius_table[i] * P.lig_born_radius_table[j] * exp(-cal_dis_square / (4. * P.lig_born_radius_table[i] * P.lig_born_radius_table[j])));
                // }

                if(P.bond_type_mat(i+P.rna_heavy_size,j+P.rna_heavy_size) == tmd::DUMMY_BOND) {
                    continue;
                }
                // P.ligand_ELEenergy = P.ligand_ELEenergy + P.lB0 / P.e1 * P.lig.noH[i].charge * P.lig.noH[j].charge / cal_dis;

                const double& equ_dis = (P.lig.noH[j].radius + P.lig.noH[i].radius);
                // P.complex_dis_mat(i+P.rna_heavy_size,j+P.rna_heavy_size) = cal_dis;
                // P.complex_dis_mat[i+P.rna_heavy_size][j+P.rna_heavy_size] = cal_dis;
                if (cal_dis < P.LJ_cut * equ_dis)
                {
                    const int& iLJ_dis = (cal_dis / (P.LJ_cut * equ_dis)) * P.LJstep + 1;
                    P.ligand_LJenergy = P.ligand_LJenergy + P.LJ_table[P.lig.noH[j].type_idx][P.lig.noH[i].type_idx][iLJ_dis];
                }
            }
        }
    }

    // if(this->RL_Score_Mode == RL_ALL) {
        P.RNA_SELFenergy = 0.0;
        P.ELEenergy = 0.0;
        P.HBenergy = 0.0;
        for (int iatom_RNA = 0; iatom_RNA < P.rna_heavy_size; iatom_RNA++)
        {
            P.RNA_SELFenergy = P.RNA_SELFenergy + 0.5 * (1. / P.e2 - 1. / P.e1) * P.lB0 * P.rna.noH[iatom_RNA].charge * P.rna.noH[iatom_RNA].charge * (1. / P.complex_born_radius_table[iatom_RNA] - 1. / P.rna_born_radius_table[iatom_RNA]);
            for (int i = 0; i < P.lig_heavy_size; i++)
            {
                P.ELEenergy = P.ELEenergy + P.lB0 / P.e1 * P.lig.noH[i].charge * P.rna.noH[iatom_RNA].charge / P.complex_dis_mat(iatom_RNA,i+P.rna_heavy_size);


                ///////////////////////////////////////////////////////////////////////
                // if(P.ELEenergy < -400) {
                //     std::cout << P.ELEenergy << " " << iatom_RNA << " " << i << " " << P.lig.noH[i].charge << " " << P.rna.noH[iatom_RNA].charge << " " << P.complex_dis_mat(iatom_RNA,i+P.rna_heavy_size) << std::endl;
                //     yz_lig.write(std::cout);

                //     const double& dis = P.complex_dis_mat(iatom_RNA, i+P.rna_heavy_size);
                //     const double equ_dis = (P.rna.noH[iatom_RNA].radius + P.lig.noH[i].radius);
                //     if (dis < P.LJ_cut * equ_dis) {
                //         const int iLJ_dis = (dis / (P.LJ_cut * equ_dis)) * P.LJstep + 1;
                //         // P.LJenergy = P.LJenergy + P.LJ_table[P.rna.noH[iatom_RNA].type_idx][P.lig.noH[i].type_idx][iLJ_dis];
                //         std::cout << dis << " " << equ_dis << " " << iLJ_dis << " " << P.rna.noH[iatom_RNA].type_idx << " " << P.lig.noH[i].type_idx << " " << P.LJ_table[P.rna.noH[iatom_RNA].type_idx][P.lig.noH[i].type_idx][iLJ_dis] << std::endl;
                //     }
                //     std::cout << P.LJenergy << std::endl;

                //     const double& lig_x = P.lig.noH[i].x;
                //     const double& lig_y = P.lig.noH[i].y;
                //     const double& lig_z = P.lig.noH[i].z;
                //     std::cout << "atomtype: " << P.lig.noH[i].name << " " << P.lig.noH[i].yz_atom_index << " lig coord: " << lig_x << " " << lig_y << " " << lig_z << std::endl;
                //     const int atom_x_grid = static_cast<int>(lig_x/P.lj_grids.width);
                //     const int atom_y_grid = static_cast<int>(lig_y/P.lj_grids.width);
                //     const int atom_z_grid = static_cast<int>(lig_z/P.lj_grids.width);

                //     const tmd::Grid& grid = P.lj_grids.at(atom_x_grid,atom_y_grid,atom_z_grid);
                //     std::cout << "lj grids: " << P.lj_grids.dim_1() << " " << P.lj_grids.dim_2() << " " << P.lj_grids.dim_3() << " " << P.lj_grids.width << " " << P.lj_grids.min_i << " " << P.lj_grids.min_j << " " << P.lj_grids.min_k << " " << std::endl;
                //     std::cout << grid.flag << std::endl;
                //     if(grid.flag) {
                //         std::cout << "rigid atoms size: " << grid.rigid_atoms.size() << std::endl;
                //         for(const int& r_j : grid.rigid_atoms) {
                //             const double& rna_x = P.rna.noH[r_j].x;
                //             const double& rna_y = P.rna.noH[r_j].y;
                //             const double& rna_z = P.rna.noH[r_j].z;
                //             std::cout << r_j << std::endl;
                //             std::cout << "rna atomtype: " << P.rna.noH[r_j].name << " " << P.rna.noH[r_j].yz_atom_index << " rna coord: " << rna_x << " " << rna_y << " " << rna_z << std::endl;
                //             // // const double& dis = sqrt((rna_x-lig_x)*(rna_x-lig_x)+(rna_y-lig_y)*(rna_y-lig_y)+(rna_z-lig_z)*(rna_z-lig_z));
                //             // const double& dis = P.complex_dis_mat(r_j, i+P.rna_heavy_size);
                //             // const double equ_dis = (P.rna.noH[r_j].radius + P.lig.noH[i].radius);
                //             // if (dis < P.LJ_cut * equ_dis) {
                //             //     const int iLJ_dis = (dis / (P.LJ_cut * equ_dis)) * P.LJstep + 1;
                //             //     P.LJenergy = P.LJenergy + P.LJ_table[P.rna.noH[r_j].type_idx][P.lig.noH[i].type_idx][iLJ_dis];
                //             // }
                //         }
                //     }


                //     assert(false);
                // }
                ////////////////////////////////////////////////////////////////////////
                // POLenergy = POLenergy + lB0 * (1. / e2 - 1. / e1) * charge_ligand[i] * charge_RNA[iatom_RNA] / sqrt(dis_ok[i][iatom_RNA] * dis_ok[i][iatom_RNA] + Rborn[isite][atomtype_ligand[i]] * rb0[iatom_RNA] * exp(-dis_ok[i][iatom_RNA] * dis_ok[i][iatom_RNA] / (4. * Rborn[isite][atomtype_ligand[i]] * rb0[iatom_RNA])));
                // double equ_dis = P.radius_with_water[P.lig.noH[i].type_idx] + P.radius_with_water[P.rna.noH[iatom_RNA].type_idx];
                // if (dis_ok[i][iatom_RNA] < equ_dis)
                // {
                //     iLJ_dis = (dis_ok[i][iatom_RNA] / (equ_dis)) * LJstep + 1;
                //     SASAenergy = SASAenergy + SASA[atomtype_RNA[iatom_RNA]][atomtype_ligand[i]][iLJ_dis];
                // }
                if (P.lig.noH[i].numH > 0 || P.rna.noH[iatom_RNA].numH > 0)
                {
                    const double& equ_dis = P.Hbond_max * (P.radius_type[P.lig.noH[i].type_idx] + P.radius_type[P.rna.noH[iatom_RNA].type_idx]);
                    if (P.complex_dis_mat(iatom_RNA,i+P.rna_heavy_size) < equ_dis)
                    {
                        const int& iLJ_dis = (P.complex_dis_mat(iatom_RNA,i+P.rna_heavy_size) / (equ_dis)) * P.LJstep + 1;
                        P.HBenergy = P.HBenergy - P.Hbond_table[P.rna.noH[iatom_RNA].type_idx][P.lig.noH[i].type_idx][iLJ_dis];
                    }
                }
            }
        }

        //complex pol
        P.complex_POLenergy = 0.0;
        for (int i = 0; i < P.complex_heavy_size; i++)
        {
            for (int j = 0; j < P.complex_heavy_size; j++)
            {
                if (i < j)
                {
                    const double& dis = P.complex_dis_mat(i,j);
                    const double& dis_square = P.complex_dis_square_mat(i,j);
                    P.complex_POLenergy = P.complex_POLenergy + P.lB0 * (1. / P.e2 - 1. / P.e1) * P.complex.noH[i].charge * P.complex.noH[j].charge / sqrt(dis_square + P.complex_born_radius_table[i] * P.complex_born_radius_table[j] * exp(-dis_square / (4. * P.complex_born_radius_table[i] * P.complex_born_radius_table[j])));

                    // if(std::isnan(P.complex_POLenergy)) {
                    //     std::cout << "complex polenergy is nan" << std::endl;
                    //     std::cout << dis_square << " " << P.rna_heavy_size << " " << i << " " << j << " " << P.complex_born_radius_table[i] << " " << P.complex_born_radius_table[j] << std::endl;
                    //     assert(false);
                    // }
                }
            }
        }
        P.POLenergy = P.complex_POLenergy - P.ligand_POLenergy - P.RNA_POLenergy;


        //cal the change in sasa after sample a new pose
        // std::cout << "start SASA_After" << std::endl;
        // SASA_After(P, P.complex, P.rna, P.lig);
        // P.SumSASA_complex = SumSASAtmp;
        // P.SASAenergy = P.gamma_SASA * (P.SumSASA_RNA + P.SumSASA_ligand) / P.rate_kcal_to_kt;

        P.sasa = 0.0;
        // std::cout << "cal SASA" << std::endl;
        std::map<int,std::vector<int>> sasa_rna_near_atom_map;
        for(int i = 0; i < P.lig_heavy_size; i++) {
            const double& lig_x = P.lig.noH[i].x;
            const double& lig_y = P.lig.noH[i].y;
            const double& lig_z = P.lig.noH[i].z;
            const int& atom_x_grid = static_cast<int>(lig_x/P.sasa_grids.width);
            const int& atom_y_grid = static_cast<int>(lig_y/P.sasa_grids.width);
            const int& atom_z_grid = static_cast<int>(lig_z/P.sasa_grids.width);
            const double& r1 = P.lig.noH[i].radius + P.radius_water;
            const double& r1_square = r1*r1;

            const tmd::Grid& grid = P.sasa_grids.at(atom_x_grid,atom_y_grid,atom_z_grid);
            if(grid.flag) {
                std::vector<int> sasa_rna_atom_indices;
                for(const int& r_j : grid.rigid_atoms) {
                    // const double& rna_x = P.rna.noH[r_j].x;
                    // const double& rna_y = P.rna.noH[r_j].y;
                    // const double& rna_z = P.rna.noH[r_j].z;
                    // const double& r2 = P.rna.noH[r_j].radius + P.radius_water;
                    // const double& dis = sqrt((rna_x-lig_x)*(rna_x-lig_x)+(rna_y-lig_y)*(rna_y-lig_y)+(rna_z-lig_z)*(rna_z-lig_z));
                    const double& dis = P.complex_dis_mat(r_j,i+P.rna_heavy_size);
                    if (dis < (P.rna.noH[r_j].radius + P.lig.noH[i].radius + 2.0 * P.radius_water)) {
                        sasa_rna_atom_indices.push_back(r_j);
                        if(sasa_rna_near_atom_map.find(r_j)!=sasa_rna_near_atom_map.end()) {
                            sasa_rna_near_atom_map[r_j].push_back(i);
                        } else {
                            sasa_rna_near_atom_map.insert({r_j,{i}});
                        }
                    }
                }
                for(const sasa_polar& sap : P.sasa_polar_table[P.lig.noH[i].type_idx]) {
                    const double& Rx = sap.Rx + lig_x;
                    const double& Ry = sap.Ry + lig_y;
                    const double& Rz = sap.Rz + lig_z;
                    for(const int& r_j : sasa_rna_atom_indices) {
                        const double& rna_x = P.rna.noH[r_j].x;
                        const double& rna_y = P.rna.noH[r_j].y;
                        const double& rna_z = P.rna.noH[r_j].z;
                        const double& r2 = P.rna.noH[r_j].radius + P.radius_water;
                        if ( (Rx-rna_x)*(Rx-rna_x)+(Ry-rna_y)*(Ry-rna_y)+(Rz-rna_z)*(Rz-rna_z) <= r2*r2 ) {
                            P.sasa += r1_square * sap.sin_sita * sap.dsita * sap.dphi;
                            break;
                        }
                    }
                }
            }
        }
        for(const auto& rna_near : sasa_rna_near_atom_map) {
            const int& ri = rna_near.first;
            const double& rna_x = P.rna.noH[ri].x;
            const double& rna_y = P.rna.noH[ri].y;
            const double& rna_z = P.rna.noH[ri].z;
            const double& r2 = P.rna.noH[ri].radius + P.radius_water;
            const double& r2_square = r2*r2;

            for(const sasa_polar& sap : P.sasa_polar_table[P.rna.noH[ri].type_idx]) {
                const double& Rx = sap.Rx + rna_x;
                const double& Ry = sap.Ry + rna_y;
                const double& Rz = sap.Rz + rna_z;
                for(const int& lj : rna_near.second) {
                    const double& lig_x = P.lig.noH[lj].x;
                    const double& lig_y = P.lig.noH[lj].y;
                    const double& lig_z = P.lig.noH[lj].z;
                    const double& r1 = P.lig.noH[lj].radius + P.radius_water;
                    if ( (Rx-lig_x)*(Rx-lig_x)+(Ry-lig_y)*(Ry-lig_y)+(Rz-lig_z)*(Rz-lig_z) <= r1*r1 ) {
                        P.sasa += r2_square * sap.sin_sita * sap.dsita * sap.dphi;
                        break;
                    }
                }
            }
        }
        // std::cout << "sasa: " <<  P.gamma_SASA * P.sasa / P.rate_kcal_to_kt << std::endl;
        P.SASAenergy = P.gamma_SASA * P.sasa / P.rate_kcal_to_kt;
    // }

    ////////////////////////////////////////////////////////////
    // std::cout << "cal whole energy" << std::endl;
    const double energy =
            P.re_weight[0] * P.LJenergy +
            P.re_weight[1] * P.ELEenergy +
            P.re_weight[2] * P.POLenergy +
            P.re_weight[3] * P.ligand_SELFenergy +
            P.re_weight[4] * P.SASAenergy +
            P.re_weight[5] * P.HBenergy +
            P.re_weight[6] * P.RNA_SELFenergy
            // +P.ligand_ELEenergy
            +P.ligand_LJenergy;

    // if(energy < 23.18 && yz_lig.rmsd_with_respect_to_ref_atoms() < 20) {
    //     std::cout << "score mode: " << this->RL_Score_Mode << std::endl;
    //     yz_lig.write(std::cout);
    //     std::cout << "LJenergy: " << P.LJenergy << std::endl;
    //     std::cout << "ELEenergy: " << P.ELEenergy << std::endl;
    //     std::cout << "POLenergy: " << P.POLenergy  << " --> complex: " << P.complex_POLenergy << " rna: " << P.RNA_POLenergy << " lig: " << P.ligand_POLenergy << std::endl;
    //     std::cout << "ligand_SELFenergy: " << P.ligand_SELFenergy << std::endl;
    //     std::cout << "SASAenergy: " << P.SASAenergy << std::endl;
    //     std::cout << "HBenergy: " << P.HBenergy << std::endl;
    //     std::cout << "RNA_SELFenergy: " << P.RNA_SELFenergy << std::endl;
    //     std::cout << "ligand_LJenergy: " << P.ligand_LJenergy << std::endl;
    //     // std::cout << "ligand_ELEenergy: " << P.ligand_ELEenergy << std::endl;
    //     std::cout << "total energy: " << energy << std::endl;
    //     // P.print_lig_born_radius_table();
    //     // P.print_rna_born_radius_table();
    //     // P.print_lig_after_born_radius_table();
    //     // P.print_rna_after_born_radius_table();
    //     // P.print_complex_born_radius_table();
    //     std::cout << "----------------------------------------------" << std::endl;
    //     // if(energy<-1000) {
    //     //     exit(2);
    //     // }
    //     for(int i = 0; i < P.lig_heavy_size; ++i) {
    //         std::cout << P.lig.noH[i].index << " " << P.lig.noH[i].name << " " << P.lig.noH[i].x << " " << P.lig.noH[i].y << " " << P.lig.noH[i].z << " " << P.lig.noH[i].charge << " " << P.lig_born_radius_table[i] << " " << P.lig.noH[i].radius << " " << P.lig.noH[i].born_scale << std::endl;
    //     }
    //     for(int i = P.rna_heavy_size; i < P.complex.noH.size(); ++i) {
    //         std::cout << P.complex.noH[i].index << " " << P.complex.noH[i].name << " " << P.complex.noH[i].x << " " << P.complex.noH[i].y << " " << P.complex.noH[i].z << " " << P.complex.noH[i].charge << " " << P.complex_born_radius_table[i] << " " << P.complex.noH[i].radius << " " << P.complex.noH[i].born_scale << std::endl;
    //     }
    //     for(int i = 0; i < P.rna_heavy_size; ++i) {
    //         std::cout << P.rna.noH[i].index << " " << P.rna.noH[i].name << " " << P.rna.noH[i].x << " " << P.rna.noH[i].y << " " << P.rna.noH[i].z << " " << P.rna.noH[i].charge << " " << P.rna_born_radius_table[i] << " " << P.rna.noH[i].radius << " " << P.rna.noH[i].born_scale << std::endl;
    //     }
    //     for(int i = 0; i < P.rna_heavy_size; ++i) {
    //         std::cout << P.complex.noH[i].index << " " << P.complex.noH[i].name << " " << P.complex.noH[i].x << " " << P.complex.noH[i].y << " " << P.complex.noH[i].z << " " << P.complex.noH[i].charge << " " << P.complex_born_radius_table[i] << " " << P.complex.noH[i].radius << " " << P.complex.noH[i].born_scale << std::endl;
    //     }

    //     for (int i = 0; i < P.rna_heavy_size; i++)
    //     {
    //         for (int j = P.rna_heavy_size; j < P.complex_heavy_size; j++)
    //         {
    //             if (i < j)
    //             {
    //                 const double& dis_square =
    //                     (P.complex.noH[j].x - P.complex.noH[i].x) * (P.complex.noH[j].x - P.complex.noH[i].x) +
    //                     (P.complex.noH[j].y - P.complex.noH[i].y) * (P.complex.noH[j].y - P.complex.noH[i].y) +
    //                     (P.complex.noH[j].z - P.complex.noH[i].z) * (P.complex.noH[j].z - P.complex.noH[i].z);
    //                 if(sqrt(dis_square) < 5) {
    //                     std::cout << i << " " << j << " " << sqrt(dis_square) << " " << P.complex_dis_mat(i,j) << " " << dis_square << " " << P.complex_dis_square_mat(i,j) << std::endl;
    //                 }
    //             }
    //         }
    //     }
    //     exit(2);
    // }

    return energy/static_cast<double>(P.lig_heavy_size);

}

}