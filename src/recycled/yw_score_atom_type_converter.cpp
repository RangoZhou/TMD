
//我自己写的mol2 转新原子类型的文件，所有转换规则都在 New_Atom_Type_Lig() 里面
//一般涉及到的是原子电荷，成键数，以及成键原子中H的数目或者非H的数目。lig是去H的，H的电荷赋值给成键原子的
//实际过程中，一些应该是带正点的基团，比如-NH3,合并电荷后，N是带-0.001这样，所以在需要判断正负的时候>-0.1都算+，<0.1算-
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
#include <unistd.h>
// #include "list/nr3.h"
// #include "list/ran.h"

using namespace std;

enum enum_type1{RNA,LIG};

struct struct_basic_parameter

{

    int contact_number = 70;

    vector<string> contact_type;

    string input_dir = "0FOOD/";

    string output_dir = "1Transition/";



    double rmax = 10.0;

    double dr = 0.2;

    int usize = int(rmax/dr);

    double Pi = 3.1415926;

    double ds = 0.0001;

    double spd = 1.0;

    double umax =1.68787*3.0;//in the unit of kBT = 3 kcal/mol ;1kcal/mol = 1.68787kBT

    double gmin = exp(-umax);

    double kcal2kbT = 1.68787;

    int num_iterative;

    int num_thread = 7;

    double rmsd_max = 4.0;

   

    bool clear_flag =false;

    int ensemble_size =4000;

    int vi_ensemble_size =2000;

    double non_near_native_ratio = 0.95;

    int RL_size = 10;

    // vector< vector<int> > decoy_category= {{0,1},{2,3},{4,5}};

    std::map<char,double> LJ_rmin = {{'O',1.6612},{'C',1.9080},{'N',1.8240},{'S',2.0},{'P',2.1},{'H',2.0}};

    std::map<char,double> LJ_eps= {{'O',0.2100},{'C',0.086},{'N',1.7},{'S',0.25},{'P',0.2},{'H',0.265}};

    vector<string> special_case = {"4P5J_SPM_1"};

    bool input_u = false;

    double Vmax = 4.0/3.0*Pi*rmax*rmax*rmax;

    vector<double> dv;



    set<string> mol2_ad;

    set<string> zou_ad;

    //   O           1.6612  0.2100             OPLS

    //   C           1.9080  0.0860             Spellmeyer

    //   N           1.8240  0.1700             OPLS

    //   S           2.0000  0.2500             W. Cornell CH3SH and CH3SCH3 FEP's

    //   P           2.1000  0.2000             JCC ,7,(1986),230;

    //   H           2.0     0.265              Halgon 

    //   F           1.75    0.061              Gough et al. JCC 13,(1992),963.

    //   Cl          1.948   0.265              Fox, JPCB,102,8070,(98),flex.mdl CHCl3

    //   Br          2.22    0.320              Junmei(?)

    //   I           2.35    0.40               JCC,7,(1986),230;  

    // double dangle = 20.0;

    // double angle2rad = Pi/180.0;

    // double rad_bin = dangle*angle2rad;

    // int angle_bin_num = 180.0/dangle;

    // std::map<char,vector< vector<string> > > torsion = {

    //     {'A',{{"N1","C2"},{"C2","N3"},{"N3","C4"},{"C4","C5"},{"C5","C6"},{"C6","N1"},{"C4","N9"},{"N9","C8"},{"C8","N7"},{"N7","C5"},{"C6","N6"},{"N9","C1'"},{"C1'","C2'"},{"C2'","C3'"},{"C3'","C4'"},{"C4'","O4'"},{"O4'","C1'"},{"C3'","O3'"},{"C4'","C5'"},{"C5'","O5'"},{"C2'","O2'"},{"O5'","P"},{"P","OP1"},{"P","OP2"}}},

    //     {'C',{{"N1","C2"},{"C2","N3"},{"N3","C4"},{"C4","C5"},{"C5","C6"},{"C6","N1"},{"C2","O2"},{"C4","N4"},{"N1","C1'"},{"C1'","C2'"},{"C2'","C3'"},{"C3'","C4'"},{"C4'","O4'"},{"O4'","C1'"},{"C3'","O3'"},{"C4'","C5'"},{"C5'","O5'"},{"C2'","O2'"},{"O5'","P"},{"P","OP1"},{"P","OP2"}}},

    //     {'G',{{"N1","C2"},{"C2","N3"},{"N3","C4"},{"C4","C5"},{"C5","C6"},{"C6","N1"},{"C4","N9"},{"N9","C8"},{"C8","N7"},{"N7","C5"},{"C6","O6"},{"N9","C1'"}, {"C1'","C2'"},{"C2'","C3'"},{"C3'","C4'"},{"C4'","O4'"},{"O4'","C1'"},{"C3'","O3'"},{"C4'","C5'"},{"C5'","O5'"},{"C2'","O2'"},{"O5'","P"},{"P","OP1"},{"P","OP2"}}},

    //     {'U',{{"N1","C2"},{"C2","N3"},{"N3","C4"},{"C4","C5"},{"C5","C6"},{"C6","N1"},{"C2","O2"},{"C4","O4"},{"N1","C1'"},{"C1'","C2'"},{"C2'","C3'"},{"C3'","C4'"},{"C4'","O4'"},{"O4'","C1'"},{"C3'","O3'"},{"C4'","C5'"},{"C5'","O5'"},{"C2'","O2'"},{"O5'","P"},{"P","OP1"},{"P","OP2"}}},

    // };



    std::vector<string> lig_atom_type={"C.2", "C.3", "C.ar", "C.cat", "N.2", "N.3", "N.4", "N.ar", "N.am", "N.pl3", "O.2", "O.3", "O.co2", "S.2", "S.3", "P.3", "Ha"};

    std::vector<double> pair_weight;

    std::vector<double> ori_weight;



    double emax;

    double emin;

    double eother;

    

    double Tc = 25.0;

    double eps_water =(87.740-0.4008*Tc+9.398*1e-4*Tc*Tc-1.41*1e-6*Tc*Tc*Tc); //e2

    double eps_rna = 20.0;//e1

    double lB = 7.15*(273+25)*eps_water/((273+Tc)*eps_water);

    double lB0=eps_water*lB;

    double rO = 1.5;//1.66;  //1.5

    double rP = 1.9;//2.10;  // 1.9

    double rH = 1.0;//2.0;   //1.0

    double rC = 1.7;//1.908;  //1.7

    double rN = 1.65;//1.8240;  //1.65

    double rS = 1.8;//2.00;   //1.80

    // double rO = 1.66;//1.5;//1.66;  //1.5

    // double rP = 2.10;//1.9;//2.10;  // 1.9

    // double rH = 2.0;//1.0;//2.0;   //1.0

    // double rC = 1.908;//1.7;//1.908;  //1.7

    // double rN = 1.824;//1.65;//1.8240;  //1.65

    // double rS = 2.0;//1.8;//2.00;   //1.80

    double rrest = 1.9;//2.0

    double sO = 0.85;

    double sP = 0.86;

    double sH = 0.85;

    double sC = 0.72;

    double sN = 0.79;

    double sS = 0.80;

    double srest = 0.85;



    //Tc    25.000    C 

    //E1    20.000 

    //        O    P    H   C     N     S

    //radius 1.5  1.9  1.0  1.7  1.65  1.80

    //Born scaling 0.85  0.86  0.85  0.72  0.79  0.80 



    double step_sasa=0.25;

    double gamma_sasa=0.005;

    double rwater = 1.4;

    double rate_kcal_to_kt=0.593*(273+Tc)/(273+25);



   





};



struct bond_info

{

    int num_bond= 0;

    int bond_index;

    std::vector<int> bonded_atom;

    std::vector<string>bonded_mol2;

    std::vector<string>bonded_type;

    std::vector<double>bonded_charge; 

};

struct atom_info

{

    int index;

    string atom_type;

    string atom_type_pdbqt;

    string residue_name;

    int res_index;

    double x;

    double y;

    double z;

    int mol2index;

    double charge=0.;

    string element_type;

    string ad;//donor or acceptor



    int index_mol2;

    int index_pdbqt;



    string atom_name_mol2;

    string atom_type_mol2;

    string residue_name_mol2;

    bond_info bond;



    string new_type;

    string don_acc;

};





struct obs_pdbqt

{   

    vector<atom_info> pdbqt_original;

    vector<atom_info> pdbqt_heavy_atom;

};

struct contact_parameter

{

    // double rou_ave = 0.;

    vector< int > cnt;

    vector<double >f;

    // vector<double >fmax;

    // vector<double >frest;

    // vector< bool > firstflag;

    vector<double >g;

    double gb = 0.;//g_below

    vector<double >u;

    int max_index;

    // double dV;

    // double F;

    // double Fmax;

    // double Frest;

    bool F_first_flag = true;

    int cnt_total;

};



struct statis_info

{

    vector<int> cnt;

    int cnt_total;

};



struct exp_pdbqt

{

    double rmsd;

    double score;

    double eng;

    double eng_pair;

    // vector<atom_info> pdbqt_original;

    vector<atom_info> exp;

    std::map<string, statis_info >exp_cnt;



    double lig_pol;

    double complex_pol;

    double lig_sasa;

    double complex_sasa;

    double delta_sasa;



};

struct residue_info

{

    char residue_name;

    std::map<string,atom_info> atoms;

};

struct complex_info

{

    string name;



    bool have_ba = false;

    double ba_eng = 0.;

    

    double eng_max;

    double eng_rest;

    vector<atom_info> obs_lig;

    vector<atom_info> obs_lig_mol2;

    vector<atom_info> rec;

    vector<atom_info> rec_mol2;

    std::map<string, statis_info >obs_cnt;

    vector<exp_pdbqt> exp_pack;

    double eng_max_pair;

    double eng_rest_pair;



    



};

struct input_flag

{

    int near = 0;//rmsd 0-3

    int moderate = 0;//rmsd 3-6

    int far = 0; //rmsd 6-9

};

struct thread_parameter

{

    int num_begin;

    int num_end;

};










void New_Atom_Type_Lig(atom_info& a)

{

    if(a.atom_type_mol2=="C.1"||a.atom_type_mol2=="C.2"||a.atom_type_mol2=="C.ar"||a.atom_type_mol2=="C.cat")

    {

        bool carbon_amide_flag=false;

        bool carbon_positive_flag=false;

        bool carbon_negative_flag=false;

        bool carbon_O2_flag=false;

        for(int i=0;i!=a.bond.num_bond;i++)

        {

            if(a.bond.bonded_mol2[i]=="N.am"){carbon_amide_flag=true;break;}

            char element = a.bond.bonded_mol2[i][0];

            double charge = a.bond.bonded_charge[i];

            if(element=='N'&&charge>-0.1)carbon_positive_flag=true;

            if(element=='O'&&charge<(-0.1))carbon_negative_flag=true;

            if(a.bond.bonded_mol2[i]=="O.2"){carbon_O2_flag=true;}

          



        }

          if(carbon_amide_flag){a.new_type="C2N";}

            else if(carbon_negative_flag){a.new_type="C2M";}

            else if(carbon_positive_flag){a.new_type="C2P";}

            else if(carbon_O2_flag){a.new_type="C2O";}

            else {a.new_type="C2X";}



    }

    else if(a.atom_type_mol2=="C.3")

    {

        bool hydrogen_carbon_only =true;

        for(int i=0;i!=a.bond.num_bond;i++)

        {

            char element = a.bond.bonded_mol2[i][0];

            if(element!='H'&&element!='C'){hydrogen_carbon_only=false;break;}

        }

        if(hydrogen_carbon_only){a.new_type="C3F";}

        else{a.new_type="C3X";}

    }

    else if(a.atom_type_mol2=="N.4"){a.new_type="NC";}

    else if(a.atom_type_mol2=="N.3")

    {

        int num_hydrogen = 0;

        for(int i=0;i!=a.bond.num_bond;i++)

        {

            if(a.bond.bonded_mol2[i]=="H")num_hydrogen++;

        }

        if(num_hydrogen!=0){a.new_type="NC";}

        else{a.new_type="N3X";}

    }

    else if(a.atom_type_mol2=="N.1"){a.new_type="N1";}

    else if(a.atom_type_mol2=="N.am"){a.new_type="N2N";}

    else if(a.atom_type_mol2=="N.2"||a.atom_type_mol2=="N.ar"||a.atom_type_mol2=="N.pl3")

    {

        // if(a.charge>0.1){a.new_type="NC";}

        // else{

            int non_hydrogen = 0;

            int num_hydrogen = 0;

            for(int i=0;i!=a.bond.num_bond;i++)

            {

                if(a.bond.bonded_mol2[i]=="H"){num_hydrogen++;}

                else{non_hydrogen++;}

            }

            if(non_hydrogen==0){a.new_type="NC";}

            else if(non_hydrogen==1){a.new_type="N21";}

            else if(non_hydrogen==2){a.new_type="N22";}

            else{a.new_type="N2X";}



        // }

    }

    else if(a.atom_type_mol2=="O.co2"){a.new_type="OC";}

    else if(a.atom_type_mol2=="O.2"){a.new_type="O2";}

    else if(a.atom_type_mol2=="O.3")

    {

        int non_hydrogen = 0;

        for(int i=0;i!=a.bond.num_bond;i++)

        {

           if(a.bond.bonded_mol2[i]!="H")non_hydrogen++;

        }

       if(non_hydrogen>=2){a.new_type="O32";}

       else if(non_hydrogen<=1){a.new_type="O31";}

       else {a.new_type="O33";cout<<a.bond.num_bond<<" "<<endl;

       }

        //a.new_type="O32";

       

    }

    else if(a.atom_type_mol2=="S.2"||a.atom_type_mol2=="S.3"||a.atom_type_mol2=="S.o"||a.atom_type_mol2=="S.o2"||a.atom_type_mol2=="S.O"||a.atom_type_mol2=="S.O2")

    {

        bool o2_flag = false;

        int non_hydrogen = 0;

        for(int i=0;i!=a.bond.num_bond;i++)

        {

            if(a.bond.bonded_mol2[i]=="O.2")o2_flag=true;

            if(a.bond.bonded_type[i]=="1" && a.bond.bonded_mol2[i]!="H")non_hydrogen++;

        }

        if(o2_flag){a.new_type="SX";}

        // else if(non_hydrogen==1){a.new_type="S1";}

        else{a.new_type="SX";}

    }

    else if(a.atom_type_mol2=="Ha"){a.new_type="Ha";}

    else if(a.atom_type_mol2=="P.3"){a.new_type="P";}

    else {

        a.new_type="MET";

        // cout<<"WRONG "<<a.index<<" "<<a.atom_type_mol2<<" "<<a.x<<" "<<a.charge<<endl;

    }

    







}



void Input_decoy(struct_basic_parameter& P, vector<complex_info>& all_sub)

{

   

    for(auto& aa:all_sub)

    {

        string inputstr = aa.name;

        complex_info cmp_info_temp;

        cmp_info_temp.name = inputstr;

        string pdbqtstr;

        int index,number_atoms,num_bonds;

        ////////lig input

        ifstream output_pdbqt;

  

        // system(("rm "+P.output_dir+inputstr+"/"+inputstr+"_rec.txt").c_str());

        // ifstream rec_dat_file((P.output_dir+inputstr+"/"+inputstr+"_rec_q.txt").c_str());

       

            output_pdbqt.open((P.input_dir+inputstr+"/"+inputstr+"_rec_wH.pdbqt").c_str());

            index=0;

            cout<<cmp_info_temp.name<<endl;

            std::vector<atom_info> donor_rec;

            while(getline(output_pdbqt,pdbqtstr))

            {

                // cout<<pdbqtstr<<endl;

                if(pdbqtstr.find("ATOM  ")!=std::string::npos || pdbqtstr.find("HETATM")!=std::string::npos)

                {

                    atom_info a;

                    a.index = index;

                    a.atom_type_pdbqt = pdbqtstr.substr(12,4);

                    a.residue_name = pdbqtstr.substr(17,3);

                    // a.res_index = std::stod(pdbqtstr.substr(22,4));

                    a.x = std::stod(pdbqtstr.substr(30,8));

                    a.y = std::stod(pdbqtstr.substr(38,8));

                    a.z = std::stod(pdbqtstr.substr(46,8));

                    a.charge = std::stod(pdbqtstr.substr(70,6));

                    a.element_type = pdbqtstr.substr(77,1);

                    string str1 = pdbqtstr.substr(77,1);

                    string str2 = pdbqtstr.substr(78,1);

                    if(str1 == "A" || str2 == "A")a.don_acc="acc";



                    if(a.element_type =="A")a.element_type ="C";

                    if(a.element_type =="H"){donor_rec.push_back(a);}

                    else{

                        cmp_info_temp.rec.push_back(a);

                        index++;

                    }

                    

                }

            }

            for(auto&b:donor_rec)

            {

                double dis_min=1000000.0;int index_min;

                for(int i=0;i!=cmp_info_temp.rec.size();i++)

                {

                    atom_info a = cmp_info_temp.rec[i];

                    double dis_2 = (b.x-a.x)*(b.x-a.x)+(b.y-a.y)*(b.y-a.y)+(b.z-a.z)*(b.z-a.z);

                    if(dis_2<dis_min){dis_min = dis_2;index_min = i;}

                }

                cmp_info_temp.rec[index_min].don_acc = "don";

                cmp_info_temp.rec[index_min].charge +=b.charge;



            }

          



            output_pdbqt.close();

            cout<<cmp_info_temp.name<<" "<<cmp_info_temp.rec.size()<<endl;





            int index_mol2=0;

            output_pdbqt.open((P.input_dir+inputstr+"/"+inputstr+"_rec_wH.mol2").c_str());

            while(getline(output_pdbqt,pdbqtstr))

            {

                //  cout<<pdbqtstr<<endl;

                if(pdbqtstr.find("@<TRIPOS>MOLECULE")!=std::string::npos)

                {

                     

                    std::getline(output_pdbqt,pdbqtstr);

                    // cout<<pdbqtstr<<endl;

                    std::getline(output_pdbqt,pdbqtstr);

                    // cout<<pdbqtstr<<endl;

                    std::stringstream ss(pdbqtstr);

                    ss >> number_atoms;

                    ss>>num_bonds;//cout

                }

                // cout<<pdbqtstr<<endl;

                if(pdbqtstr.find("@<TRIPOS>ATOM")!=std::string::npos)

                {

                    for(int i = 1; i <= number_atoms; ++i)

                    {

                        std::getline(output_pdbqt,pdbqtstr);

                        std::istringstream ss(pdbqtstr);

                        std::string buf;

                        std::vector<std::string> token;

                        while(ss >> buf) token.push_back(buf);

                        // cout<<token[0]<<" "<<token[2]<<" "<<endl;



                        int idx_record = std::stod(token[0]);

                        while(index_mol2<idx_record){

                            atom_info fake_mol2;

                            fake_mol2.index = -1;

                            cmp_info_temp.rec_mol2.push_back(fake_mol2);

                            index_mol2++;

                        }

                        

                        double xx = std::stod(token[2]);

                        double yy = std::stod(token[3]);

                        double zz = std::stod(token[4]);

                        double q = std::stod(token[8]);

                        int res_index = std::stoi(token[6]);

                        string mol2type = (token[5]);

                        if(mol2type=="O"){cout<<inputstr<<" "<<xx<<" "<<i<<" "<<endl;break;}

                        string mol2_atom_name = (token[1]);

                        string mol2_residue_name = (token[7]);

                        if(mol2type=="F"||mol2type=="Cl"||mol2type=="Br"||mol2type=="I"||mol2type=="BR"){mol2type="Ha";}

                        

                        int index_pdbqt = -1; 

                        string don_acc;

                        if(mol2type=="H")

                        {

                        }

                        else

                        {

                            for(auto&a:cmp_info_temp.rec)

                            {

                                if(abs(xx-a.x)<0.001&&abs(yy-a.y)<0.001&&abs(zz-a.z)<0.001)

                                {

                                    a.residue_name_mol2 = mol2_residue_name;

                                    a.atom_name_mol2 = mol2_atom_name;

                                    a.atom_type_mol2 = mol2type;

                                    a.res_index = res_index;

                                    q=a.charge ;

                                    a.index_mol2 = index_mol2;

                                    index_pdbqt = a.index;

                                    don_acc=a.don_acc;

                                    

                                    break;

                                    // a.atom_type_pdbqt = token[1];

                                    // a.res_index = std::stoi(token[6]);

                                    // a.residue_name = token[7];

                                }

                            }

                        }



                       

                        atom_info atom_mol2;

                        atom_mol2.index=index_mol2;

                        atom_mol2.x=xx;

                        atom_mol2.y=yy;

                        atom_mol2.z=zz;

                        atom_mol2.charge=q;

                        atom_mol2.residue_name_mol2 = mol2_residue_name;

                        atom_mol2.atom_name_mol2 = mol2_atom_name;

                        atom_mol2.atom_type_mol2 = mol2type;

                        atom_mol2.index_pdbqt=index_pdbqt;

                        atom_mol2.don_acc=don_acc;

                        cmp_info_temp.rec_mol2.push_back(atom_mol2);

                        index_mol2++;

                    }

                }

                if(pdbqtstr.find("@<TRIPOS>BOND")!=std::string::npos)

                {

                    for(int i = 1; i <= num_bonds; ++i)

                    {

                        std::getline(output_pdbqt,pdbqtstr);

                        std::istringstream ss(pdbqtstr);

                        std::string buf;

                        std::vector<std::string> token;

                        while(ss >> buf) token.push_back(buf);

                        int bondleft = std::stoi(token[1]);

                        int bondright = std::stoi(token[2]);

                        string bond_type = token[3];



                        // if((bondleft==353||bondright==353)&&inputstr=="1CS7-S02#1")cout<<bondleft<<" "<<bondright<<endl;

                        if(cmp_info_temp.rec_mol2[bondleft].index_pdbqt>=0)

                        {

                            int aa=cmp_info_temp.rec_mol2[bondright].index_pdbqt;

                            string aa_mol2=cmp_info_temp.rec_mol2[bondright].atom_type_mol2;

                            double q=cmp_info_temp.rec_mol2[bondright].charge;

                            cmp_info_temp.rec_mol2[bondleft].bond.num_bond++;

                            cmp_info_temp.rec_mol2[bondleft].bond.bonded_atom.push_back(aa);

                            cmp_info_temp.rec_mol2[bondleft].bond.bonded_type.push_back(bond_type);

                            cmp_info_temp.rec_mol2[bondleft].bond.bonded_mol2.push_back(aa_mol2);

                            cmp_info_temp.rec_mol2[bondleft].bond.bonded_charge.push_back(q);

                        }

                        if(cmp_info_temp.rec_mol2[bondright].index_pdbqt>=0)

                        {

                            int aa=cmp_info_temp.rec_mol2[bondleft].index_pdbqt;

                            string aa_mol2=cmp_info_temp.rec_mol2[bondleft].atom_type_mol2;

                            double q=cmp_info_temp.rec_mol2[bondleft].charge;



                            cmp_info_temp.rec_mol2[bondright].bond.num_bond++;

                            cmp_info_temp.rec_mol2[bondright].bond.bonded_atom.push_back(aa);

                            cmp_info_temp.rec_mol2[bondright].bond.bonded_type.push_back(bond_type);

                            cmp_info_temp.rec_mol2[bondright].bond.bonded_mol2.push_back(aa_mol2);

                            cmp_info_temp.rec_mol2[bondright].bond.bonded_charge.push_back(q);

                        }









                    }

                }

                // 

            }

            for(auto&a: cmp_info_temp.rec)

                {

                    int mol2_index = a.index_mol2;

                    // cout<<a.index<<" "<<a.x<<" "<<mol2_index<<" "<<a.index_mol2<<endl;

                    a.bond = cmp_info_temp.rec_mol2[mol2_index].bond;

                }

                // exit(3);

            output_pdbqt.close();

            cout<<cmp_info_temp.name<<endl;

            // exit(3);

            for(auto&a:cmp_info_temp.rec)

            {

                // a.new_type=New_Atom_Type_Rec(a.atom_name_mol2,a.residue_name_mol2);

                New_Atom_Type_Lig(a);

                

                if(a.new_type=="MET")cout<<cmp_info_temp.name<<" rec "<<a.index<<" "<<a.x<<" "<<a.new_type<<" "<<a.atom_type_mol2<<endl;



                // cout<<a.index<<" "<<a.atom_type_mol2<<" "<<a.new_type<<endl;



            }

            std::ofstream rec_output_dat_file((P.output_dir+inputstr+"/"+inputstr+"_rec0.txt").c_str());

            for(auto&a:cmp_info_temp.rec)

            {

                rec_output_dat_file<<a.index<<" "<<a.atom_type_mol2<<" "<<a.new_type<<" "<<a.x<<" "<<a.y<<" "<<a.z<<" "<<a.atom_type_pdbqt<<" "<<a.residue_name_mol2<<" "<<a.res_index<<" "<<a.element_type<<" "<<a.charge<<endl;

            }

            rec_output_dat_file.close();





            // if(inputstr=="1CS7-S02#1")exit(3);

            ///////////////////////////lig



            output_pdbqt.open((P.input_dir+inputstr+"/"+inputstr+"_lig_wH.pdbqt").c_str());

            index=0;

            donor_rec.clear();

            while(getline(output_pdbqt,pdbqtstr))

            {

                if(pdbqtstr.find("ATOM  ")!=std::string::npos || pdbqtstr.find("HETATM")!=std::string::npos)

                {

                    atom_info a;

                    a.index = index;

                    a.atom_type_pdbqt = pdbqtstr.substr(12,4);

                    a.residue_name = pdbqtstr.substr(17,3);

                    a.x = std::stod(pdbqtstr.substr(30,8));

                    a.y = std::stod(pdbqtstr.substr(38,8));

                    a.z = std::stod(pdbqtstr.substr(46,8));

                    a.charge = std::stod(pdbqtstr.substr(70,6));

                    a.element_type = pdbqtstr.substr(77,1);

                    string str1 = pdbqtstr.substr(77,1);

                    string str2 = pdbqtstr.substr(78,1);

                    if(str1 == "A" || str2 == "A")a.don_acc="acc";

                    if(a.element_type =="A")a.element_type ="C";

                    if(a.element_type =="H"){donor_rec.push_back(a);}

                    else{

                        cmp_info_temp.obs_lig.push_back(a);

                        index++;

                    }

                    // native_pdbqt.push_back(a);

                }



            }

            for(auto&b:donor_rec)

            {

                double dis_min=1000000.0;int index_min;

                for(int i=0;i!=cmp_info_temp.obs_lig.size();i++)

                {

                    atom_info a = cmp_info_temp.obs_lig[i];

                    double dis_2 = (b.x-a.x)*(b.x-a.x)+(b.y-a.y)*(b.y-a.y)+(b.z-a.z)*(b.z-a.z);

                    if(dis_2<dis_min){dis_min = dis_2;index_min = i;}

                }

                cmp_info_temp.obs_lig[index_min].don_acc = "don";

                cmp_info_temp.obs_lig[index_min].charge +=b.charge;



            }

            output_pdbqt.close();

            // cout<<"lig_pdbqt"<<endl;

            /////mol2

            index_mol2=0;

            output_pdbqt.open((P.input_dir+inputstr+"/"+inputstr+"_lig_wH.mol2").c_str());

            while(getline(output_pdbqt,pdbqtstr))

            {

                // cout<<pdbqtstr<<endl;

                if(pdbqtstr.find("@<TRIPOS>MOLECULE")!=std::string::npos)

                {

                    std::getline(output_pdbqt,pdbqtstr);

                    std::getline(output_pdbqt,pdbqtstr);

                    std::stringstream ss(pdbqtstr);

                    ss >> number_atoms;

                    ss>>num_bonds;

                }

                if(pdbqtstr.find("@<TRIPOS>ATOM")!=std::string::npos)

                {

                    for(int i = 1; i <= number_atoms; ++i)

                    {

                         std::getline(output_pdbqt,pdbqtstr);

                        std::istringstream ss(pdbqtstr);

                        std::string buf;

                        std::vector<std::string> token;

                        while(ss >> buf) token.push_back(buf);

                        

                        int idx_record = std::stod(token[0]);

                        while(index_mol2<idx_record){

                            atom_info fake_mol2;

                            fake_mol2.index = -1;

                            cmp_info_temp.obs_lig_mol2.push_back(fake_mol2);

                            index_mol2++;

                        }

                        

                        double xx = std::stod(token[2]);

                        double yy = std::stod(token[3]);

                        double zz = std::stod(token[4]);

                        double q = std::stod(token[8]);

                        int res_index = std::stoi(token[6]);

                        string mol2type = (token[5]);

                        string mol2_atom_name = (token[1]);

                        string mol2_residue_name = (token[7]);



                        if(mol2type=="O"){cout<<inputstr<<" "<<xx<<" "<<i<<" "<<endl;break;}

                       if(mol2type=="F"||mol2type=="Cl"||mol2type=="Br"||mol2type=="I"||mol2type=="BR")mol2type="Ha";



                        int index_pdbqt = -1; 

                        if(mol2type=="H")

                        {}

                        else

                        {

                            for(auto&a:cmp_info_temp.obs_lig)

                            {

                                // cout<<a.index<<" "<<a.x<<" "<<xx<<" "<<a.y<<" "<<yy<<" "<<a.z<<" "<<zz<<endl;

                                if(abs(xx-a.x)<0.001&&abs(yy-a.y)<0.001&&abs(zz-a.z)<0.001)

                                {

                                    a.residue_name_mol2 = mol2_residue_name;

                                    a.atom_name_mol2 = mol2_atom_name;

                                    a.atom_type_mol2 = mol2type;

                                    a.res_index = res_index;

                                    q = a.charge;

                                    a.index_mol2 = index_mol2;

                                    index_pdbqt = a.index;

                                    break;

                                }





                            }

                        }

                         

                        atom_info atom_mol2;

                        atom_mol2.index=index_mol2;

                        atom_mol2.x=xx;

                        atom_mol2.y=yy;

                        atom_mol2.z=zz;

                        atom_mol2.charge=q;

                        atom_mol2.residue_name_mol2 = mol2_residue_name;

                        atom_mol2.atom_name_mol2 = mol2_atom_name;

                        atom_mol2.atom_type_mol2 = mol2type;

                        atom_mol2.index_pdbqt=index_pdbqt;

                        cmp_info_temp.obs_lig_mol2.push_back(atom_mol2);

                        index_mol2++;

                        // if(cmp_info_temp.name=="1LC4_TOY_1")cout<<atom_mol2.index<<" "<<index_pdbqt<<" "<<xx<<" "<<yy<<endl;

                    }

                }

                if(pdbqtstr.find("@<TRIPOS>BOND")!=std::string::npos)

                {

                    for(int i = 1; i <= num_bonds; ++i)

                    {

                        std::getline(output_pdbqt,pdbqtstr);

                        std::istringstream ss(pdbqtstr);

                        std::string buf;

                        std::vector<std::string> token;

                        while(ss >> buf) token.push_back(buf);

                        int bondleft = std::stoi(token[1]);

                        int bondright = std::stoi(token[2]);

                        string bond_type = token[3];

 

                        if(cmp_info_temp.obs_lig_mol2[bondleft].index_pdbqt>=0)

                        {

                            int aa=cmp_info_temp.obs_lig_mol2[bondright].index_pdbqt;

                            string aa_mol2=cmp_info_temp.obs_lig_mol2[bondright].atom_type_mol2;

                            double q=cmp_info_temp.obs_lig_mol2[bondright].charge;



                            cmp_info_temp.obs_lig_mol2[bondleft].bond.num_bond++;

                            cmp_info_temp.obs_lig_mol2[bondleft].bond.bonded_atom.push_back(aa);

                            cmp_info_temp.obs_lig_mol2[bondleft].bond.bonded_type.push_back(bond_type);

                            cmp_info_temp.obs_lig_mol2[bondleft].bond.bonded_mol2.push_back(aa_mol2);

                            cmp_info_temp.obs_lig_mol2[bondleft].bond.bonded_charge.push_back(q);

                        }

                        if(cmp_info_temp.obs_lig_mol2[bondright].index_pdbqt>=0)

                        {

                            int aa=cmp_info_temp.obs_lig_mol2[bondleft].index_pdbqt;

                            string aa_mol2=cmp_info_temp.obs_lig_mol2[bondleft].atom_type_mol2;

                            double q=cmp_info_temp.obs_lig_mol2[bondleft].charge;



                            cmp_info_temp.obs_lig_mol2[bondright].bond.num_bond++;

                            cmp_info_temp.obs_lig_mol2[bondright].bond.bonded_atom.push_back(aa);

                            cmp_info_temp.obs_lig_mol2[bondright].bond.bonded_type.push_back(bond_type);

                            cmp_info_temp.obs_lig_mol2[bondright].bond.bonded_mol2.push_back(aa_mol2);

                            cmp_info_temp.obs_lig_mol2[bondright].bond.bonded_charge.push_back(q);

                        }

                    }

                }

                

            }

            for(auto&a: cmp_info_temp.obs_lig)

            {

                int mol2_index = a.index_mol2;

                a.bond = cmp_info_temp.obs_lig_mol2[mol2_index].bond;

            }

            output_pdbqt.close();



            for(auto&a:cmp_info_temp.obs_lig)

            {

                New_Atom_Type_Lig(a);

                // cout<<cmp_info_temp.name<<" lig "<<a.index<<" "<<a.new_type<<" "<<a.atom_type_mol2<<endl;

                if(a.new_type=="MET")cout<<cmp_info_temp.name<<" lig "<<a.index<<" "<<a.new_type<<" "<<a.atom_type_mol2<<endl;

            }

            // exit(3);

            //    cout<<"mol2_pdbqt"<<endl;

            //////

            exp_pdbqt obs_version;

            obs_version.rmsd = 0.;

            obs_version.exp = cmp_info_temp.obs_lig;

            obs_version.exp_cnt = cmp_info_temp.obs_cnt;

            cmp_info_temp.exp_pack.push_back(obs_version);

            cout<<"RL "<<cmp_info_temp.name<<" "<<cmp_info_temp.exp_pack.size()<<endl;

            if(access((P.output_dir+inputstr).c_str(),F_OK)!=0){system(("mkdir "+P.output_dir+inputstr).c_str());}

            ofstream lig_output_dat_file((P.output_dir+inputstr+"/"+inputstr+"_lig0.txt").c_str());

            index=0;

            lig_output_dat_file<<cmp_info_temp.exp_pack.size()<<" "<<cmp_info_temp.exp_pack[0].exp.size()<<endl;

            for(auto&a:cmp_info_temp.exp_pack)

            {

                lig_output_dat_file<<"MODEL "<<index<<" 0.0 "<<a.rmsd<<endl;index++;

                for(auto&b:a.exp)

                {

                    lig_output_dat_file<<b.index<<" "<<b.atom_type_mol2<<" "<<b.new_type<<" "<<b.x<<" "<<b.y<<" "<<b.z<<" "<<b.element_type<<" "<<b.charge<<endl;

                }

                lig_output_dat_file<<"ENDMDL"<<endl;

            }

            lig_output_dat_file.close();

                   

        aa = cmp_info_temp;

        // cout<<aa.name<<" finish"<<endl;

      

    } 

}









int main(int argc, char** argv)

{

    std::cout << "Have " << argc << " arguments:" << std::endl;

    for (int i = 0; i < argc; ++i){

        std::cout << argv[i] << std::endl;

    }

    string inputstr ;

    



    struct_basic_parameter P;

    string sline;

    ifstream contactlist("list/contact_new.list");

    

    

    int count=0;

    while( getline(contactlist,sline))

    {

        P.contact_type.push_back(sline);

        count++;

    }

    P.contact_number=count;

    cout<<"number of contact "<<P.contact_number<<endl;

    contactlist.close();

    





    vector<complex_info> all;

    ifstream all_list_file;

    string inputfile = "list/ribo.list";//180andBAandRNA

    all_list_file.open((inputfile).c_str());

    while(getline(all_list_file,sline))

    {

        

            complex_info cmp_info_temp; 

            cmp_info_temp.name = sline;

            // system(("mdkir "+P.output_dir+sline).c_str()); 

            all.push_back(cmp_info_temp);



    }

    all_list_file.close();



   



    if(all.size()<P.num_thread)P.num_thread = all.size();

    cout<<all.size()<<endl;

    // exit(3);

    vector<vector<complex_info> > all_input_sub(P.num_thread);

    for(int i=0;i!=all.size();i++)

    {

        int ii =i%P.num_thread;

        all_input_sub[ii].push_back(all[i]);

    }



    std::thread t_input[P.num_thread];

    for(unsigned int i = 0; i != P.num_thread; ++i)

    {

        //  cout<<i<<endl;

        t_input[i] = std::thread(Input_decoy,std::ref(P), std::ref(all_input_sub[i]));

    }

    

    

    for(unsigned int i = 0; i != P.num_thread; ++i)

    {

        t_input[i].join();

        for(auto&a:all_input_sub[i])

        {

            for(auto&b:all)

            {

                

                if(a.name==b.name){

                   b=a;

                }

            }

        }

    }

    cout<<"input finish "<<all.size()<<endl;



    // ofstream ff;

    // ff.open("mol2_ad.txt");

    // for(auto&a:P.mol2_ad){

    //     ff<<a<<endl;

    // }

    // ff.close();

    // ff.open("zou_ad.txt");

    // for(auto&a:P.zou_ad){

    //     ff<<a<<endl;

    // }

    // ff.close();







    // exit(3);





   







   



   

    



    return 0;

}