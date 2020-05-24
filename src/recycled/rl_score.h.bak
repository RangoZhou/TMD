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

void print_vector(const std::string s, const std::vector<double>& v) {
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
    int reseq;
    double charge;
    string element;
    double radius;
    double born_scale;
    int type_idx;
    int numH = 0;
    tmd::Atom_Index yz_atom_index;
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
                string::size_type position;
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
            wH.push_back(a);
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


class parameter
{
public:
    // int num_threads = 1;
    const double pi = 3.1415926;
    // int num_atom_type = 9;
    const double radius_small_sphere = 1.5;
    const double radius_large_sphere = 6.0;
    const int LJstep = 1000;
    const double LJ_cut = 2.50;
    double equ_LJ = 0.8;
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

    std::vector<double> re_weight = {3.300,1.320,0.780,2.520,0.920,0.120,4.960};


    double LJenergy = 0.;
    double ligand_LJenergy = 0.0;
    double ligand_ELEenergy = 0.0;
    double ELEenergy = 0.;
    double HBenergy = 0.;
    double POLenergy = 0.0;
    double RNA_POLenergy = 0.;
    double ligand_POLenergy = 0.;
    double complex_POLenergy = 0.;
    double RNA_SELFenergy = 0.;
    double ligand_SELFenergy = 0.;
    double SASAenergy = 0.0;

    // tmd::Strictly_Triangular_Matrix<tmd::DISTANCE_TYPE> dis_type_mat;
    tmd::Strictly_Triangular_Matrix<tmd::BOND_TYPE> bond_type_mat;

    vector<vector<vector<double>>> LJ_table;
    vector<vector<vector<double>>> SASA_table;
    vector<vector<vector<double>>> Hbond_table;
    vector<double> rna_born_radius_table;
    vector<double> lig_born_radius_table;
    vector<double> complex_born_radius_table;
    vector<double> rna_after_born_radius_table;
    vector<double> lig_after_born_radius_table;
    // tmd::Strictly_Triangular_Matrix<double> rna_dis_mat;
    // tmd::Strictly_Triangular_Matrix<double> lig_dis_mat;
    tmd::Strictly_Triangular_Matrix<double> complex_dis_mat;
    double SumSASA_ligand;
    double SumSASA_RNA;
    double SumSASA_complex;

    void print_lig_after_born_radius_table() {
        print_vector("lig_after_born_radius_table",this->lig_after_born_radius_table);
    }
    void print_rna_after_born_radius_table() {
        print_vector("rna_after_born_radius_table",this->rna_after_born_radius_table);
    }
    void print_lig_born_radius_table() {
        print_vector("lig_born_radius_table",this->lig_born_radius_table);
    }void print_rna_born_radius_table() {
        print_vector("rna_born_radius_table",this->rna_born_radius_table);
    }
    void print_complex_born_radius_table() {
        print_vector("complex_born_radius_table",this->complex_born_radius_table);
    }
};


void LJ(parameter& P)
{
    vector<double> tmp_lj(P.LJstep+1,0.0);
    vector<vector<double>> itmp_lj(P.radius_type.size(),tmp_lj);
    P.LJ_table = vector<vector<vector<double>>>(P.radius_type.size(),itmp_lj);
    ////////LJ index start from 1 in calculation 0 can not be divined.
    for(int i=0;i!=P.radius_type.size();i++)
    {
        for(int j=0;j!=P.radius_type.size();j++)
        {
            const double rdis = P.equ_LJ*(P.radius_type[i] + P.radius_type[j]);
            const double dr = P.LJ_cut*(P.radius_type[i] + P.radius_type[j])/P.LJstep;
            for(int istep=1; istep<=P.LJstep; istep++)
            {
                double dis = istep*dr;
                double dx = rdis/dis;
                double x6 = dx*dx*dx*dx*dx*dx;
                P.LJ_table[i][j][istep] = x6*x6 - x6;
            }
        }
    }
}

const double max(const double a, const double b)
{
    if (a > b)
        return (a);
    if (a < b)
        return (b);
    if (a == b)
        return (a);
}
const double min(const double a, const double b)
{
    if (a > b)
        return (b);
    if (a < b)
        return (a);
    if (a == b)
        return (b);
}

void SASA(parameter& P) {
    double dismax, dismax2, dismax3;
    double dr, dis, Rlarge, Rsmall, Hlarge, Hsmall;
    double Vlarge, Vsmall, Vall1, Vall2, Vall3;
    double anglelarge, anglesmall;

    vector<double> tmp_sasa(P.LJstep+1,0.0);
    vector<vector<double>> itmp_sasa(P.radius_type.size(),tmp_sasa);
    P.SASA_table = vector<vector<vector<double>>>(P.radius_type.size(),itmp_sasa);

    for (int i=0;i!=P.radius_type.size();i++) {
        for (int j=0;j!=P.radius_type.size();j++) {
            dismax = P.radius_with_water[i] + P.radius_with_water[j];
            dismax2 = P.radius_type[i] + P.radius_type[j];
            dismax3 = P.radius_with_water[i] + P.radius_type[j];
            dr = dismax / P.LJstep;
            for (int istep = 1; istep <= P.LJstep; istep++) {
                dis = dr * istep;
                Rlarge = max(P.radius_with_water[i], P.radius_with_water[j]);
                Rsmall = min(P.radius_with_water[i], P.radius_with_water[j]);
                Vlarge = 0.;
                Vsmall = 0.;
                Vall1 = 0.;
                if (dis < Rlarge + Rsmall)
                {
                    if (Rlarge > dis + Rsmall)
                    {
                        Vall1 = 4. * P.pi * Rsmall * Rsmall;
                    }
                    else
                    {
                        anglelarge = (Rlarge * Rlarge + dis * dis - Rsmall * Rsmall) / (2. * Rlarge * dis);
                        Hlarge = Rlarge * (1. - anglelarge);
                        anglesmall = (Rsmall * Rsmall + dis * dis - Rlarge * Rlarge) / (2. * Rsmall * dis);
                        Hsmall = Rsmall * (1. - anglesmall);
                        Vlarge = 2. * P.pi * Rlarge * Hlarge;
                        Vsmall = 2. * P.pi * Rsmall * Hsmall;
                        Vall1 = Vlarge + Vsmall;
                    }
                }
                Rlarge = max(P.radius_type[i], P.radius_with_water[j]);
                Rsmall = min(P.radius_type[i], P.radius_with_water[j]);
                Vlarge = 0.;
                Vsmall = 0.;
                Vall2 = 0.;
                if (dis < Rlarge + Rsmall)
                {
                    if (Rlarge > dis + Rsmall)
                    {
                        Vall2 = 4. * P.pi * Rsmall * Rsmall;
                    }
                    else
                    {
                        anglelarge = (Rlarge * Rlarge + dis * dis - Rsmall * Rsmall) / (2. * Rlarge * dis);
                        Hlarge = Rlarge * (1. - anglelarge);
                        anglesmall = (Rsmall * Rsmall + dis * dis - Rlarge * Rlarge) / (2. * Rsmall * dis);
                        Hsmall = Rsmall * (1. - anglesmall);
                        Vlarge = 2. * P.pi * Rlarge * Hlarge;
                        Vsmall = 2. * P.pi * Rsmall * Hsmall;
                        Vall2 = Vlarge + Vsmall;
                    }
                }
                Rlarge = max(P.radius_with_water[i], P.radius_type[j]);
                Rsmall = min(P.radius_with_water[i], P.radius_type[j]);
                Vlarge = 0.;
                Vsmall = 0.;
                Vall3 = 0.;
                if (dis < Rlarge + Rsmall)
                {
                    if (Rlarge > dis + Rsmall)
                    {
                        Vall3 = 4. * P.pi * Rsmall * Rsmall;
                    }
                    else
                    {
                        anglelarge = (Rlarge * Rlarge + dis * dis - Rsmall * Rsmall) / (2. * Rlarge * dis);
                        Hlarge = Rlarge * (1. - anglelarge);
                        anglesmall = (Rsmall * Rsmall + dis * dis - Rlarge * Rlarge) / (2. * Rsmall * dis);
                        Hsmall = Rsmall * (1. - anglesmall);
                        Vlarge = 2. * P.pi * Rlarge * Hlarge;
                        Vsmall = 2. * P.pi * Rsmall * Hsmall;
                        Vall3 = Vlarge + Vsmall;
                    }
                }
                if (i != 1 && j != 1)
                    P.SASA_table[i][j][istep] = P.gamma_SASA * (Vall1 - Vall2 + Vall1 - Vall3) / P.rate_kcal_to_kt;
                else
                    P.SASA_table[i][j][istep] = 0.;
                //   if(i==2&&j==1)printf("%d %d %d %lf %lf %lf %lf %lf \n",i,j,istep,dis,Vall1,Vall2,Vall3,SASA[i][j][istep]);
            }
        }
    }
}

void Hbond(parameter& P)
{
    double dismax, dr, dis, dismin;

    vector<double> tmp_hbond(P.LJstep+1,0.0);
    vector<vector<double>> itmp_bond(P.radius_type.size(),tmp_hbond);
    P.Hbond_table = vector<vector<vector<double>>>(P.radius_type.size(),itmp_bond);

    for (int i = 0; i < P.radius_type.size(); i++) {
        for (int j = 0; j < P.radius_type.size(); j++) {
            dismax = P.Hbond_max * (P.radius_type[i] + P.radius_type[j]);
            dismin = P.Hbond_min * (P.radius_type[i] + P.radius_type[j]);
            dr = dismax / P.LJstep;
            for(int istep = 1; istep <= P.LJstep; istep++)
            {
                dis = dr * istep;
                if (dis > dismin)
                {
                    P.Hbond_table[i][j][istep] = 1. - (dis - dismin) / (dismax - dismin);
                }
                else
                {
                    P.Hbond_table[i][j][istep] = 1.;
                }
                //   if(i==2&&j==3)printf("%d %d %d %lf %lf \n",i,j,istep,dis, Hbond[i][j][istep]);
            }
        }
    }
}

void Born_Radius_RNA(parameter& P, const case_info& rna) {

    double intb, Lij, Uij;
    double radi, radj, radc, dcp, dis;
    ////////////////////////////////////////
    //// pocket ////////////////////////////
    // for (i = 1; i <= num_pocket_site; i++)
    // {
    //     for (k = 1; k <= num_atom_type - 2; k++)
    //     {
    //         intb = 0;
    //         radi = radius_type[k];
    //         for (j = 1; j <= Nat_RNA; j++)
    //         {
    //             dis = (xc0_RNA[j] - xx_pocket_site[i]) * (xc0_RNA[j] - xx_pocket_site[i]) +
    //                   (yc0_RNA[j] - yy_pocket_site[i]) * (yc0_RNA[j] - yy_pocket_site[i]) +
    //                   (zc0_RNA[j] - zz_pocket_site[i]) * (zc0_RNA[j] - zz_pocket_site[i]);
    //             dis = sqrt(dis);
    //             radj = rrr_RNA[j];
    //             radc = dis + sb_RNA[j] * radj;
    //             if (radi >= radc)
    //             {
    //                 Lij = 1.;
    //                 Uij = 1.;
    //             }
    //             else
    //             {
    //                 if (radi > (dis - sb_RNA[j] * radj))
    //                 {
    //                     Lij = radi;
    //                 }
    //                 else
    //                 {
    //                     Lij = dis - sb_RNA[j] * radj;
    //                 }
    //                 Uij = dis + sb_RNA[j] * radj;
    //             }
    //             intb = intb + 0.5 * ((1. / Lij - 1. / Uij) + (sb_RNA[j] * sb_RNA[j] * radj * radj / (4. * dis) - dis / 4.) * (1. / (Lij * Lij) - 1. / (Uij * Uij)) + (1. / (2. * dis)) * log(Lij / Uij));
    //         }
    //         Rborn[i][k] = 1. / (1. / radi - intb);
    //     } // k
    // }     // i
    ///////////////////////////////////////////////
    /// RNA //////////////////////////////////////
    P.rna_born_radius_table = vector<double>(P.rna_heavy_size,0.0);
    for (int i = 0; i < P.rna_heavy_size; i++)
    { // start to calculate rb0
        intb = 0;
        radi = rna.noH[i].radius;
        for (int j = 0; j < P.rna_heavy_size; j++)
        {
            if (j!= i)
            {
                dis = (rna.noH[j].x - rna.noH[i].x) * (rna.noH[j].x - rna.noH[i].x) +
                      (rna.noH[j].y - rna.noH[i].y) * (rna.noH[j].y - rna.noH[i].y) +
                      (rna.noH[j].z - rna.noH[i].z) * (rna.noH[j].z - rna.noH[i].z);
                dis = sqrt(dis);
                radj = rna.noH[j].radius;
                radc = dis + rna.noH[j].born_scale * radj;
                if (radi >= radc)
                {
                    Lij = 1.;
                    Uij = 1.;
                }
                else
                {
                    if (radi > (dis - rna.noH[j].born_scale * radj))
                    {
                        Lij = radi;
                    }
                    else
                    {
                        Lij = dis - rna.noH[j].born_scale * radj;
                    }
                    Uij = dis + rna.noH[j].born_scale * radj;
                }
                intb = intb + 0.5 * ((1. / Lij - 1. / Uij) + (rna.noH[j].born_scale * rna.noH[j].born_scale * radj * radj / (4. * dis) - dis / 4.) * (1. / (Lij * Lij) - 1. / (Uij * Uij)) + (1. / (2. * dis)) * log(Lij / Uij));
            }
        }
        P.rna_born_radius_table[i] = 1. / (1. / radi - intb);
        //  printf("%d %lf %lf\n",i,rb0[i],rrr_RNA[i]);
    }
}

void Born_Radius_ligand(parameter& P, const case_info& lig) {

    double intb, Lij, Uij;
    double radi, radj, radc, dcp, dis;
    ///////////////////////////////////////////////
    /// ligand //////////////////////////////////////
    P.lig_born_radius_table = vector<double>(P.lig_heavy_size,0.0);
    for (int i = 0; i < P.lig_heavy_size; i++)
    { // start to calculate rb0
        intb = 0;
        radi = lig.noH[i].radius;
        for (int j = 0; j < P.lig_heavy_size; j++)
        {
            if (j!= i)
            {
                dis = (lig.noH[j].x - lig.noH[i].x) * (lig.noH[j].x - lig.noH[i].x) +
                      (lig.noH[j].y - lig.noH[i].y) * (lig.noH[j].y - lig.noH[i].y) +
                      (lig.noH[j].z - lig.noH[i].z) * (lig.noH[j].z - lig.noH[i].z);
                dis = sqrt(dis);
                radj = lig.noH[j].radius;
                radc = dis + lig.noH[j].born_scale * radj;
                if (radi >= radc)
                {
                    Lij = 1.;
                    Uij = 1.;
                }
                else
                {
                    if (radi > (dis - lig.noH[j].born_scale * radj))
                    {
                        Lij = radi;
                    }
                    else
                    {
                        Lij = dis - lig.noH[j].born_scale * radj;
                    }
                    Uij = dis + lig.noH[j].born_scale * radj;
                }
                intb = intb + 0.5 * ((1. / Lij - 1. / Uij) + (lig.noH[j].born_scale * lig.noH[j].born_scale * radj * radj / (4. * dis) - dis / 4.) * (1. / (Lij * Lij) - 1. / (Uij * Uij)) + (1. / (2. * dis)) * log(Lij / Uij));
            }
        }
        P.lig_born_radius_table[i] = 1. / (1. / radi - intb);
        //  printf("%d %lf %lf\n",i,rb0_ligand[i],rrr_ligand[i]);
    }
}

void Born_Radius_Complex(parameter& P, const case_info& complex, const case_info& rna, const case_info& lig) {

    double intb, Lij, Uij;
    double radi, radj, radc, dcp, dis;

    P.complex_born_radius_table  = vector<double>(P.complex_heavy_size,0.0);
    for (int i = 0; i < P.complex_heavy_size; i++)
    { // start to calculate rb0
        intb = 0;
        radi = complex.noH[i].radius;
        for (int j = 0; j < P.complex_heavy_size; j++)
        {
            if (j!= i)
            {
                if(i < j) {
                    dis = P.complex_dis_mat(i,j);
                } else {
                    dis = P.complex_dis_mat(j,i);
                }
                radj = complex.noH[j].radius;
                radc = dis + complex.noH[j].born_scale * radj;
                if (radi >= radc)
                {
                    Lij = 1.;
                    Uij = 1.;
                }
                else
                {
                    if (radi > (dis - complex.noH[j].born_scale * radj))
                    {
                        Lij = radi;
                    }
                    else
                    {
                        Lij = dis - complex.noH[j].born_scale * radj;
                    }
                    Uij = dis + complex.noH[j].born_scale * radj;
                }
                intb = intb + 0.5 * ((1. / Lij - 1. / Uij) + (complex.noH[j].born_scale * complex.noH[j].born_scale * radj * radj / (4. * dis) - dis / 4.) * (1. / (Lij * Lij) - 1. / (Uij * Uij)) + (1. / (2. * dis)) * log(Lij / Uij));
            }
            if(isnan(intb)) {
                std::cout << i << " " << j << std::endl;
                std::cout << complex.noH[i].born_scale << std::endl;
                std::cout << complex.noH[i].radius << std::endl;
                std::cout << complex.noH[j].born_scale << std::endl;
                std::cout << complex.noH[j].radius << std::endl;
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
        P.complex_born_radius_table[i] = 1. / (1. / radi - intb);
        //  printf("%d %lf %lf\n",i,rb0_complex[i],rrr_complex[i]);
    }
    P.rna_after_born_radius_table = vector<double>(P.rna_heavy_size,0.0);
    P.lig_after_born_radius_table = vector<double>(P.lig_heavy_size,0.0);
    for (int i = 0; i < P.complex_heavy_size; i++)
    {
        if (i < P.rna_heavy_size)
            P.rna_after_born_radius_table[i] = P.complex_born_radius_table[i];
        if (i >= P.rna_heavy_size)
            P.lig_after_born_radius_table[i - P.rna_heavy_size] = P.complex_born_radius_table[i];
    }
}

void SASA_After(parameter& P, const case_info& complex, const case_info& rna, const case_info& lig)
{
    double Rwater, step;
    double Rx, Ry, Rz;
    double radij;
    int Near[P.complex_heavy_size], LabNear[P.complex_heavy_size][100], NgSASA;
    double dsita, dphi, sita, phi;
    double SubSumSASA;
    int l, ll;

    Rwater = P.radius_water;
    step = P.step_SASA;

    double SumSASAtmp = 0.;
    for (int i = 0; i < P.rna_heavy_size; i++)
    {
        Near[i] = 0;
        for (int j = P.rna_heavy_size; j < P.complex_heavy_size; j++)
        {
            if(i < j) {
                radij = P.complex_dis_mat(i,j);
            } else {
                radij = P.complex_dis_mat(j,i);
            }
            if (i!= j)
            {
                if (radij < (complex.noH[i].radius + complex.noH[j].radius + 2 * Rwater))
                {
                    Near[i] = Near[i] + 1;
                    LabNear[i][Near[i]] = j;
                }
            }
        }
    }
    for (int i = 0; i < P.rna_heavy_size; i++)
    {
        SubSumSASA = 0.;
        dsita = step / (complex.noH[i].radius + Rwater);
        for (int j = 1; j <= P.pi / dsita; j++)
        {
            sita = j * dsita;
            dphi = step / (complex.noH[i].radius + Rwater) * sin(sita);
            for (int k = 1; k <= 2. * P.pi / dphi + 1; k++)
            {
                phi = k * dphi;
                Rx = (complex.noH[i].radius + Rwater) * sin(sita) * cos(phi) + complex.noH[i].x;
                Ry = (complex.noH[i].radius + Rwater) * sin(sita) * sin(phi) + complex.noH[i].y;
                Rz = (complex.noH[i].radius + Rwater) * cos(sita) + complex.noH[i].z;
                NgSASA = 0;
                for (ll = 1; ll <= Near[i]; ll++)
                {
                    l = LabNear[i][ll];
                    if ((pow(Rx - complex.noH[l].x, 2.) + pow(Ry - complex.noH[l].y, 2.) + pow(Rz - complex.noH[l].z, 2.)) <= pow(complex.noH[l].radius + Rwater, 2.))
                    {
                        NgSASA = 1;
                        break;
                    }
                }
                if (NgSASA > 0)
                    SubSumSASA = SubSumSASA + pow(complex.noH[i].radius + Rwater, 2.) * sin(sita) * dsita * dphi;
            }
        }
        SumSASAtmp = SumSASAtmp + SubSumSASA;
        // printf("%d %d %lf\n",i,Near[i],SubSumSASA);
    }
    P.SumSASA_RNA = SumSASAtmp;

    SumSASAtmp = 0.;
    for (int i = P.rna_heavy_size; i < P.complex_heavy_size; i++)
    {
        Near[i] = 0;
        for (int j = 0; j < P.rna_heavy_size; j++)
        {
            if(i < j) {
                radij = P.complex_dis_mat(i,j);
            } else {
                radij = P.complex_dis_mat(j,i);
            }
            if (i!= j)
            {
                if (radij < (complex.noH[i].radius + complex.noH[j].radius + 2 * Rwater))
                {
                    Near[i] = Near[i] + 1;
                    LabNear[i][Near[i]] = j;
                }
            }
        }
    }
    for (int i = P.rna_heavy_size; i < P.complex_heavy_size; i++)
    {
        SubSumSASA = 0.;
        dsita = step / (complex.noH[i].radius + Rwater);
        for (int j = 1; j <= P.pi / dsita; j++)
        {
            sita = j * dsita;
            dphi = step / (complex.noH[i].radius + Rwater) * sin(sita);
            for (int k = 1; k <= 2. * P.pi / dphi + 1; k++)
            {
                phi = k * dphi;
                Rx = (complex.noH[i].radius + Rwater) * sin(sita) * cos(phi) + complex.noH[i].x;
                Ry = (complex.noH[i].radius + Rwater) * sin(sita) * sin(phi) + complex.noH[i].y;
                Rz = (complex.noH[i].radius + Rwater) * cos(sita) + complex.noH[i].z;
                NgSASA = 0;
                for (ll = 1; ll <= Near[i]; ll++)
                {
                    l = LabNear[i][ll];
                    if ((pow(Rx - complex.noH[l].x, 2.) + pow(Ry - complex.noH[l].y, 2.) + pow(Rz - complex.noH[l].z, 2.)) <= pow(complex.noH[l].radius + Rwater, 2.))
                    {
                        NgSASA = 1;
                        break;
                    }
                }
                if (NgSASA > 0)
                    SubSumSASA = SubSumSASA + pow(complex.noH[i].radius + Rwater, 2.) * sin(sita) * dsita * dphi;
            }
        }
        SumSASAtmp = SumSASAtmp + SubSumSASA;
        // printf("%d %d %lf\n",i,Near[i],SubSumSASA);
    }
    P.SumSASA_ligand = SumSASAtmp;
}


void rl_score_init(parameter& P, const tmd::RNA& yz_rna, const tmd::Ligand& yz_lig) {
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


    // std::cout << "start init LJ" << std::endl;
    LJ(P);
    // std::cout << "start init SASA" << std::endl;
    SASA(P);
    // std::cout << "start init Hbond" << std::endl;
    Hbond(P);

    // P.dis_type_mat.resize(P.complex_heavy_size,tmd::DISTANCE_FLEXIBLE);
    P.bond_type_mat.resize(P.complex_heavy_size,tmd::NONE_BOND);
    P.complex_dis_mat.resize(P.complex_heavy_size,0.0);

    for (int i = 0; i < P.lig_heavy_size; i++)
    {
        for (int j = 0; j < P.lig_heavy_size; j++)
        {
            int real_i = i + P.rna_heavy_size;
            int real_j = j + P.rna_heavy_size;
            if(real_i < real_j) {
                const tmd::Atom& ai = yz_lig.get_atoms_reference()[P.complex.noH[real_i].yz_atom_index];
                const tmd::Atom_Index& aj = P.complex.noH[real_j].yz_atom_index;
                for(const tmd::Bond& b : ai.get_bonds()) {
                    if(b.get_bonded_atom_index() == aj && b.get_bond_type() != tmd::NONE_BOND && b.get_bond_type() != tmd::NOT_CONNECTED_BOND) {
                        P.bond_type_mat(real_i,real_j) = tmd::DUMMY_BOND;
                        break;
                    }
                }
            }
        }
    }

    Born_Radius_RNA(P, P.rna);

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
            }
        }
    }

}

const double rl_score_evaluate(parameter& P, const tmd::RNA& yz_rna, const tmd::Ligand& yz_lig)
{
    clock_t start,mid1,mid2,end;
    start=clock();
    mid1=clock();
    end = clock();

    //update lig's conformation
    int index_noH=0,index_wH=0,index_HH=0;
    const tmd::Atoms& yz_lig_atoms_ref = yz_lig.get_atoms_reference();
    assert(yz_lig_atoms_ref.size()==P.lig.wH.size());
    for(const tmd::Atom& atom : yz_lig.get_atoms_reference()) {
        const tmd::Vec3d xyz = atom.get_coord();
        if(P.lig.wH[index_wH].element != "H" && P.lig.wH[index_wH].element != "h")
        {
            P.lig.noH[index_noH].x = xyz[0];
            P.lig.noH[index_noH].y = xyz[1];
            P.lig.noH[index_noH].z = xyz[2];
            index_noH++;
        }
        else
        {
            P.lig.HH[index_HH].x = xyz[0];
            P.lig.HH[index_HH].y = xyz[1];
            P.lig.HH[index_HH].z = xyz[2];
            index_HH++;
        }
        P.lig.wH[index_wH].x = xyz[0];
        P.lig.wH[index_wH].y = xyz[1];
        P.lig.wH[index_wH].z = xyz[2];
        index_wH++;
    }

    Born_Radius_ligand(P, P.lig);

    // std::cout << "cal LJ" << std::endl;
    // RNA-ligand LJ step 2
    P.LJenergy = 0.0;
    for (int i = 0; i < P.lig_heavy_size; i++)
    {
        for (int iatom_RNA = 0; iatom_RNA < P.rna_heavy_size; iatom_RNA++)
        {
            const double equ_dis = (P.rna.noH[iatom_RNA].radius + P.lig.noH[i].radius);
            const double cal_dis = sqrt(
                (P.lig.noH[i].x - P.rna.noH[iatom_RNA].x) * (P.lig.noH[i].x - P.rna.noH[iatom_RNA].x)
              + (P.lig.noH[i].y - P.rna.noH[iatom_RNA].y) * (P.lig.noH[i].y - P.rna.noH[iatom_RNA].y)
              + (P.lig.noH[i].z - P.rna.noH[iatom_RNA].z) * (P.lig.noH[i].z - P.rna.noH[iatom_RNA].z)
            );
            P.complex_dis_mat(iatom_RNA,i+P.rna_heavy_size) = cal_dis;
            // P.complex_dis_mat[i+P.rna_heavy_size][iatom_RNA] = cal_dis;
            if (cal_dis < P.LJ_cut * equ_dis)
            {
                const int iLJ_dis = (cal_dis / (P.LJ_cut * equ_dis)) * P.LJstep + 1;
                P.LJenergy = P.LJenergy + P.LJ_table[P.rna.noH[iatom_RNA].type_idx][P.lig.noH[i].type_idx][iLJ_dis];
            }
        }
    }

    // ligand LJ
    P.ligand_LJenergy = 0.0;
    for (int i = 0; i < P.lig_heavy_size; i++)
    {
        for (int j = 0; j < P.lig_heavy_size; j++)
        {
            if(i < j) {
                if(P.bond_type_mat(i+P.rna_heavy_size,j+P.rna_heavy_size) == tmd::DUMMY_BOND) {
                    continue;
                }
                const double equ_dis = (P.lig.noH[j].radius + P.lig.noH[i].radius);
                const double cal_dis = sqrt(
                    (P.lig.noH[i].x - P.lig.noH[j].x) * (P.lig.noH[i].x - P.lig.noH[j].x)
                + (P.lig.noH[i].y - P.lig.noH[j].y) * (P.lig.noH[i].y - P.lig.noH[j].y)
                + (P.lig.noH[i].z - P.lig.noH[j].z) * (P.lig.noH[i].z - P.lig.noH[j].z)
                );
                P.complex_dis_mat(i+P.rna_heavy_size,j+P.rna_heavy_size) = cal_dis;
                // P.complex_dis_mat[i+P.rna_heavy_size][j+P.rna_heavy_size] = cal_dis;
                if (cal_dis < P.LJ_cut * equ_dis)
                {
                    const int iLJ_dis = (cal_dis / (P.LJ_cut * equ_dis)) * P.LJstep + 1;
                    P.ligand_LJenergy = P.ligand_LJenergy + P.LJ_table[P.lig.noH[j].type_idx][P.lig.noH[i].type_idx][iLJ_dis];
                }
            }
        }
    }

    //lig ele
    P.ligand_ELEenergy = 0.0;
    for (int i = 0; i < P.lig_heavy_size; i++)
    {
        for (int j = 0; j < P.lig_heavy_size; j++)
        {
            if(i < j) {
                if(P.bond_type_mat(i+P.rna_heavy_size,j+P.rna_heavy_size) == tmd::DUMMY_BOND) {
                    continue;
                }
                P.ligand_ELEenergy = P.ligand_ELEenergy + P.lB0 / P.e1 * P.lig.noH[i].charge * P.lig.noH[j].charge / P.complex_dis_mat(i+P.rna_heavy_size,j+P.rna_heavy_size);
            }
        }
    }

    //other step2
    //RNA-ligand ele and hbond
    P.ELEenergy = 0.0;
    P.HBenergy = 0.0;
    for (int i = 0; i < P.lig_heavy_size; i++)
    {
        for (int iatom_RNA = 0; iatom_RNA < P.rna_heavy_size; iatom_RNA++)
        {
            P.ELEenergy = P.ELEenergy + P.lB0 / P.e1 * P.lig.noH[i].charge * P.rna.noH[iatom_RNA].charge / P.complex_dis_mat(iatom_RNA,i+P.rna_heavy_size);
            // POLenergy = POLenergy + lB0 * (1. / e2 - 1. / e1) * charge_ligand[i] * charge_RNA[iatom_RNA] / sqrt(dis_ok[i][iatom_RNA] * dis_ok[i][iatom_RNA] + Rborn[isite][atomtype_ligand[i]] * rb0[iatom_RNA] * exp(-dis_ok[i][iatom_RNA] * dis_ok[i][iatom_RNA] / (4. * Rborn[isite][atomtype_ligand[i]] * rb0[iatom_RNA])));
            // double equ_dis = P.radius_with_water[P.lig.noH[i].type_idx] + P.radius_with_water[P.rna.noH[iatom_RNA].type_idx];
            // if (dis_ok[i][iatom_RNA] < equ_dis)
            // {
            //     iLJ_dis = (dis_ok[i][iatom_RNA] / (equ_dis)) * LJstep + 1;
            //     SASAenergy = SASAenergy + SASA[atomtype_RNA[iatom_RNA]][atomtype_ligand[i]][iLJ_dis];
            // }
            if (P.lig.noH[i].numH > 0 || P.rna.noH[iatom_RNA].numH > 0)
            {
                const double equ_dis = P.Hbond_max * (P.radius_type[P.lig.noH[i].type_idx] + P.radius_type[P.rna.noH[iatom_RNA].type_idx]);
                if (P.complex_dis_mat(iatom_RNA,i+P.rna_heavy_size) < equ_dis)
                {
                    const int iLJ_dis = (P.complex_dis_mat(iatom_RNA,i+P.rna_heavy_size) / (equ_dis)) * P.LJstep + 1;
                    P.HBenergy = P.HBenergy - P.Hbond_table[P.rna.noH[iatom_RNA].type_idx][P.lig.noH[i].type_idx][iLJ_dis];
                }
            }
        }
        // SELFenergy = SELFenergy + 0.5 * (1. / e2 - 1. / e1) * lB0 * charge_ligand[i] * charge_ligand[i] * (1. / Rborn[isite][atomtype_ligand[i]] - 1. / rrr_ligand[i]);
    }

    /////////////////////////////////////////////////////
    /////// ligand alone ////////////////////////////////
    /////////////////////////////////////////////////////
    // std::cout << "cal ligand_POLenergy" << std::endl;
    P.ligand_POLenergy = 0.0;
    for (int i = 0; i < P.lig_heavy_size; i++)
    {
        for (int j = 0; j < P.lig_heavy_size; j++)
        {
            if (i < j)
            {
                double dis = (P.lig.noH[i].x - P.lig.noH[j].x) * (P.lig.noH[i].x - P.lig.noH[j].x) +
                      (P.lig.noH[i].y - P.lig.noH[j].y) * (P.lig.noH[i].y - P.lig.noH[j].y) +
                      (P.lig.noH[i].z - P.lig.noH[j].z) * (P.lig.noH[i].z - P.lig.noH[j].z);
                dis = sqrt(dis);
                P.ligand_POLenergy = P.ligand_POLenergy + P.lB0 * (1. / P.e2 - 1. / P.e1) * P.lig.noH[i].charge * P.lig.noH[j].charge / sqrt(dis * dis + P.lig_born_radius_table[i] * P.lig_born_radius_table[j] * exp(-dis * dis / (4. * P.lig_born_radius_table[i] * P.lig_born_radius_table[j])));
                // P.lig_dis_mat[i][j] = dis;
                // P.lig_dis_mat[j][i] = dis;
                P.complex_dis_mat(i + P.rna_heavy_size,j + P.rna_heavy_size) = dis;
                // P.complex_dis_mat[j + P.rna_heavy_size][i + P.rna_heavy_size] = dis;
            }
        }
    }
    /////////////////////////////////////////////////////
    ///////// complex ///////////////////////////////////
    /////////////////////////////////////////////////////

    // for (int i = 0; i < P.rna_heavy_size; i++) {
    //     for (int j = 0 + P.rna_heavy_size; j < P.complex_heavy_size; j++)
    //     {
    //         double dis =
    //                 (complex.noH[i].x - complex.noH[j].x) * (complex.noH[i].x - complex.noH[j].x) +
    //                 (complex.noH[i].y - complex.noH[j].y) * (complex.noH[i].y - complex.noH[j].y) +
    //                 (complex.noH[i].z - complex.noH[j].z) * (complex.noH[i].z - complex.noH[j].z);
    //         dis = sqrt(dis);
    //         // Rij_complex[i][j] = dis;
    //         // Rij_complex[j][i] = dis;
    //     }
    // }

    // std::cout << "start Born_Radius_Complex" << std::endl;
    Born_Radius_Complex(P, P.complex, P.rna, P.lig);
    // std::cout << "cal complex_POLenergy" << std::endl;
    //complex pol
    P.complex_POLenergy = 0.0;
    for (int i = 0; i < P.complex_heavy_size; i++)
    {
        for (int j = 0; j < P.complex_heavy_size; j++)
        {
            if (i < j)
            {
                double dis = P.complex_dis_mat(i,j);
                P.complex_POLenergy = P.complex_POLenergy + P.lB0 * (1. / P.e2 - 1. / P.e1) * P.complex.noH[i].charge * P.complex.noH[j].charge / sqrt(dis * dis + P.complex_born_radius_table[i] * P.complex_born_radius_table[j] * exp(-dis * dis / (4. * P.complex_born_radius_table[i] * P.complex_born_radius_table[j])));
            }
        }
    }
    P.POLenergy = P.complex_POLenergy - P.ligand_POLenergy - P.RNA_POLenergy;

    // std::cout << "cal RNA_SELFenergy" << std::endl;
    P.RNA_SELFenergy = 0.0;
    for (int i = 0; i < P.rna_heavy_size; i++)
    {
        P.RNA_SELFenergy = P.RNA_SELFenergy + 0.5 * (1. / P.e2 - 1. / P.e1) * P.lB0 * P.rna.noH[i].charge * P.rna.noH[i].charge * (1. / P.rna_after_born_radius_table[i] - 1. / P.rna_born_radius_table[i]);
    }

    // std::cout << "cal ligand_SELFenergy" << std::endl;
    P.ligand_SELFenergy = 0.0;
    for (int i = 0; i < P.lig_heavy_size; i++)
    {
        P.ligand_SELFenergy = P.ligand_SELFenergy + 0.5 * (1. / P.e2 - 1. / P.e1) * P.lB0 * P.lig.noH[i].charge * P.lig.noH[i].charge * (1. / P.lig_after_born_radius_table[i] - 1. / P.lig_born_radius_table[i]);
    }

    // std::cout << "start SASA_After" << std::endl;
    SASA_After(P, P.complex, P.rna, P.lig);
    // P.SumSASA_complex = SumSASAtmp;
    P.SASAenergy = P.gamma_SASA * (P.SumSASA_RNA + P.SumSASA_ligand) / P.rate_kcal_to_kt;

    ////////////////////////////////////////////////////////////
    // std::cout << "cal whole energy" << std::endl;
    const double energy =
            P.re_weight[0] * P.LJenergy +
            P.re_weight[1] * P.ELEenergy +
            P.re_weight[2] * P.POLenergy +
            P.re_weight[3] * P.ligand_SELFenergy +
            P.re_weight[4] * P.SASAenergy +
            P.re_weight[5] * P.HBenergy +
            P.re_weight[6] * P.RNA_SELFenergy;

    std::cout << "LJenergy: " << P.LJenergy << std::endl;
    std::cout << "ELEenergy: " << P.ELEenergy << std::endl;
    std::cout << "POLenergy: " << P.POLenergy << std::endl;
    std::cout << "ligand_SELFenergy: " << P.ligand_SELFenergy << std::endl;
    std::cout << "SASAenergy: " << P.SASAenergy << std::endl;
    std::cout << "HBenergy: " << P.HBenergy << std::endl;
    std::cout << "RNA_SELFenergy: " << P.RNA_SELFenergy << std::endl;
    std::cout << "ligand_LJenergy: " << P.ligand_LJenergy << std::endl;
    std::cout << "ligand_ELEenergy: " << P.ligand_ELEenergy << std::endl;
    // P.print_lig_born_radius_table();
    // P.print_rna_born_radius_table();
    // P.print_lig_after_born_radius_table();
    // P.print_rna_after_born_radius_table();
    // P.print_complex_born_radius_table();
    // std::cout << "----------------------------------------------" << std::endl;
    exit(2);

    return energy;

}

}