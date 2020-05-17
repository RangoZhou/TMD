
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

using namespace std;



class parameter 

{

    public:

    int num_threads = 1;



    double pi = 3.1415926;

    int num_atom_type = 9;

    double radius_small_sphere = 1.5;

    double radius_large_sphere = 6.0;

    int LJstep = 1000;

    double LJ_cut = 2.50;

    double equ_LJ = 0.8;

    double excludLJ = 0.6;

    double Num_in_sphere = 1000;

    int Nats_RNA = 5500;

    int Nats_ligand =200;

    int Nats_RNA_addH = 11000;

    int Nats_ligand_addH = 400;

    double gamma_SASA = -0.0054;

    double r_water = 1.4;

    double step_pocket = 0.5;

    int MaxN = 50;

    int Maxpose = 20000;

    double gamma_step1 = 10.0;

    double Max_pocket_site = 10000;

    double rmsd_cut_step1 = 0.5;

    double ligand_cut_step1 = 2.0;

    int TOP_step1 = 2;



    double LJ_cut_new =8.0;

    double lj_dr = LJ_cut_new/LJstep;

    double iseed1 = 25478553;

    double iseed2 = 25212112;

    double Tc = 25.0;

    double e1 = 20.0;

    double radius_O = 1.5;

    double radius_P = 1.9;

    double radius_H = 1.0;

    double radius_C = 1.7;

    double radius_N = 1.65;

    double radius_S = 1.80;

    // vector<double> atom_radius = {radius_H,radius_C,radius_N,radius_O,radius_P,radius_S,radius_C,radius_small_sphere,radius_large_sphere};

    vector<double> atom_radius = {radius_H,radius_C,radius_N,radius_O,radius_P,radius_S,radius_C};

    double born_scale_O = 0.85;

    double born_scale_P = 0.86;

    double born_scale_H = 0.85;

    double born_scale_C = 0.72;

    double born_scale_N = 0.79; 

    double born_scale_S = 0.80;

    double e2 = (87.740-0.4008*Tc+9.398*1e-4*Tc*Tc-1.41*1e-6*Tc*Tc*Tc);   // the dielctric consant of water in temperature T 

    double e25 = 87.740-0.4008*25+9.398*1e-4*25*25-1.41*1e-6*25*25*25;    // the dielctric consant of water in temperature 25

    double lB = 7.15*(273+25)*e25/((273+Tc)*e2);   // e^2/(ebs*kB*Tc) 

    double lB0=e2*lB;        // in A lb0=e^2/(4*pi*e0*kB*T) 

    double rate_kcal_to_kt=0.593*(273+Tc)/(273+25);



    vector<double> ux;

    vector<double> vy;

    vector<double> wz;

    void SetSphere(string filename)

    {

        ifstream spherefile(filename.c_str());

        string sline;

        while(getline(spherefile,sline))

        {

            std::istringstream ss(sline);

            std::string buf;

            std::vector<std::string> token;

            while(ss >> buf) token.push_back(buf);

            ux.push_back(stod(token[1]));

            vy.push_back(stod(token[2]));

            wz.push_back(stod(token[3]));

        }

    }



    vector<double> aa;

    vector<double> bb;

    vector<double> cc;

    vector<double> dd;

    vector<double> ee;

    vector<double> ff;

    vector<double> gg;

    vector<double> hh;

    vector<double> ii;



    vector< vector<double> > poses;



    void SetAngle()

    {

        int Ngamma = 360.0/gamma_step1+0.01;

        int index = 0;

        for(int ig =0; ig != Ngamma;ig++)

        {

            double dg = ig*(2*pi/Ngamma);

            double sing = sin(dg);

            double cosg = cos(dg);

            for(int i=0;i!=ux.size();i++)

            {

                aa.push_back(ux[i]*ux[i]+(vy[i]*vy[i]+wz[i]*wz[i])*cosg);

                bb.push_back(ux[i]*vy[i]*(1.-cosg)-wz[i]*sing);

                cc.push_back(ux[i]*wz[i]*(1.-cosg)+vy[i]*sing);

                dd.push_back(ux[i]*vy[i]*(1.-cosg)+wz[i]*sing);

                ee.push_back(vy[i]*vy[i]+(ux[i]*ux[i]+wz[i]*wz[i])*cosg);

                ff.push_back(vy[i]*wz[i]*(1.-cosg)-ux[i]*sing);

                gg.push_back(ux[i]*wz[i]*(1.-cosg)-vy[i]*sing);

                hh.push_back(vy[i]*wz[i]*(1.-cosg)+ux[i]*sing);

                ii.push_back(wz[i]*wz[i]+(ux[i]*ux[i]+vy[i]*vy[i])*cosg);

                

                index++;

            }

        }

    }

};

struct thread_parameter

{

    int ibegin = 0;

    int iend;

};

struct struct_MC_thread_info

{

    int seed1 = 25245822;

    int seed2 = 25311542;

    int sample_num = 0;

};



struct atom_info

{

    int index;

    string name;

    double x;

    double y;

    double z;

    string type;

    string resi_name;

    double q;



    string element;

    double rrr;

    double sb;

    int type_idx;

    

};

struct cordinate_info

{

    int index;

    double x;

    double y;

    double z;    

};





class Grids_info

{

    public:

    double igd;

    int Nnx;

    int Nny;

    int Nnz;

    double x0;

    double y0;

    double z0;

    vector< vector< vector<int> > > grids_pocket_flag;

    vector< vector< vector< vector<double> > > >grids_energy;    

    vector< vector< vector< vector<int> > > >lj_flag;

    void Set_grid_int(){

        grids_pocket_flag=vector<vector<vector< int > > >(Nnx, vector<vector< int > >(Nny, vector< int >(Nnz,0)));

    }

    void initial_grid_lj(const parameter& P){

        for(int i=0;i!=P.atom_radius.size();i++){

        grids_energy.push_back(vector<vector<vector< double > > >(Nnx, vector<vector< double > >(Nny, vector< double >(Nnz,0))) );

        lj_flag.push_back(vector<vector<vector< int > > >(Nnx, vector<vector< int > >(Nny, vector< int >(Nnz,0))) );

        }

    }

    void Find_RNA_grid(const parameter& P, int xi,int xa, int yi, int ya, int zi, int za, const atom_info& rec_atom){

        for(int i=xi; i<=xa; i++)

        for(int j=yi; j<=ya; j++)

        for(int k=zi; k<=za; k++)

        {

            if(grids_pocket_flag[i][j][k]==0)

            {

                double xtmp = P.step_pocket*i+x0 - rec_atom.x;

                dis  = xtmp*xtmp;

                double ytmp = P.step_pocket*j+y0 - rec_atom.y;

                dis += ytmp*ytmp;

                double ztmp = P.step_pocket*k+z0 - rec_atom.z;

                dis += ztmp*ztmp;

                dis = sqrt(dis);

                if( dis < rec_atom.rrr + P.radius_small_sphere) grids_pocket_flag[i][j][k]=1;

            }

        }

    }

    void Set_LJ(const parameter& P, const double& select_radius, const atom_info& rec_atom, const vector< vector< vector<double> > > lj)

    {

        xi = (rec_atom.x - x0 - select_radius)/igd;   //if(xi<0)xi=0;

        xa = (rec_atom.x - x0 + select_radius)/igd+1; //if(xa>=Nnx)xa=Nnx-1;



        yi = (rec_atom.y - y0 - select_radius)/igd;   //if(yi<0)yi=0;

        ya = (rec_atom.y - y0 + select_radius)/igd+1; //if(ya>=Nny)ya=Nny-1;



        zi = (rec_atom.z - z0 - select_radius)/igd;   //if(zi<0)zi=0;

        za = (rec_atom.z - z0 + select_radius)/igd+1; //if(za>=Nnz)za=Nnz-1;



        for(int i=xi; i<=xa; i++)

        for(int j=yi; j<=ya; j++)

        for(int k=zi; k<=za; k++)

        {

            double xtmp = igd*i+x0 - rec_atom.x;

            dis  = xtmp * xtmp;

            double ytmp = igd*j+y0 - rec_atom.y;

            dis += ytmp * ytmp;

            double ztmp = igd*k+z0 - rec_atom.z;

            dis += ztmp*ztmp;

            dis = sqrt(dis);

            int dis_idx = dis/P.lj_dr;

            for(int l=0;l!=P.atom_radius.size();l++)

            {

                if(dis < P.excludLJ*(P.atom_radius[l]+rec_atom.rrr))lj_flag[l][i][j][k]=1;

                grids_energy[l][i][j][k] += lj[rec_atom.type_idx][l][dis_idx];

            }

        }





    }

    private:

    double dis;

    int xi;

    int xa;

    int yi;

    int ya;

    int zi;

    int za;

    

    







};

struct label_idx

{

    int lig_idx;

    int initial_idx;

    double lj;

};

bool less_label_lj (const label_idx& l1, const label_idx& l2)

{

    return l1.lj < l2.lj;

}

class pocket_info

{

    public:

    double x;

    double y;

    double z;

    double lj_site_min=1000.0;

    vector< vector<label_idx> > LJ_min;

    

    // bool operator < (const pocket_info& p)

    // {

    //     return lj_site_min < p.lj_site_min;

    // }

};

bool less_pocket_lj (const pocket_info& p1, const pocket_info& p2)

{

    return p1.lj_site_min < p2.lj_site_min;

}

struct grid_energy_map

{

    vector<double> lj;

}

;



class case_info

{

    public:

    std::vector< atom_info > noH;

    std::vector< atom_info > wH;

    std::vector< atom_info > HH;

    void SetValue (const parameter P, string filename)

    {

        ifstream ff(filename.c_str());

        string sline;

        while (getline(ff,sline))

        {

            if(sline.find("@<TRIPOS>ATOM")!=std::string::npos)

            {

                int index=0,index_H=0;

                while (getline(ff,sline))

                {

                    if(sline.find("@<TRIPOS>BOND")!=std::string::npos)break;

                    std::istringstream ss(sline);

                    std::string buf;

                    std::vector<std::string> token;

                    while(ss >> buf) token.push_back(buf);



                    atom_info a;

                    a.index = index_H;index_H++;

                    a.name = token[1];

                    a.x = stod(token[2]);

                    a.y = stod(token[3]);

                    a.z = stod(token[4]);

                    a.type = token[5];

                    a.resi_name = token[7];

                    a.q = stod(token[8]);



                    if(a.type.find(".")!=std::string::npos)

                    {

                        string::size_type position;

                        position = a.type.find(".");

                        a.element = a.type.substr(0,position);

                    }

                    else

                    {

                        a.element = a.type;

                    }



                    if(a.element == "H" || a.element == "h")

                    {

                        a.rrr = P.radius_H;

                        a.sb = P.born_scale_H;

                        a.type_idx = 0;

                    }

                    else if(a.element == "C" || a.element == "c")

                    {

                        a.rrr = P.radius_C;

                        a.sb = P.born_scale_C;

                        a.type_idx = 1;



                    }

                    else if(a.element == "N" || a.element == "n")

                    {

                        a.rrr = P.radius_N;

                        a.sb = P.born_scale_N;

                        a.type_idx = 2;



                    }

                    else if(a.element == "O" || a.element == "o")

                    {

                        a.rrr = P.radius_O;

                        a.sb = P.born_scale_O;

                        a.type_idx = 3;

                    }

                    else if(a.element == "P" || a.element == "p")

                    {

                        a.rrr = P.radius_P;

                        a.sb = P.born_scale_P;

                        a.type_idx = 4;

                    }

                    else if(a.element == "S" || a.element == "s")

                    {

                        a.rrr = P.radius_S;

                        a.sb = P.born_scale_S;

                        a.type_idx = 5;

                    }

                    else

                    {

                        a.rrr = P.radius_C;

                        a.sb = P.born_scale_C;

                        a.type_idx = 6;

                    }



                    if(a.element != "H" && a.element != "h")

                    {

                        a.index = index;

                        index++;

                        noH.push_back(a);

                        xcen+=a.x;

                        ycen+=a.y;

                        zcen+=a.z;

                    }

                    else

                    {

                        HH.push_back(a);

                    }

                    wH.push_back(a);

                }

            }

        }

        xcen/=(1.0*noH.size());

        ycen/=(1.0*noH.size());

        zcen/=(1.0*noH.size());



    }

    void H_charge_transfer()

    {

        for(auto& a:HH)

        {

            int idx;

            double dmin = 10000000;

            for(auto& b:noH)

            {

                double dis = (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z);

                if(dis<dmin){dmin=dis;idx=b.index;}

            }

            noH[idx].q += a.q;    

        }

    };

    

    std::vector< atom_info > center0;

    



    std::vector< atom_info > center;

    void Move_to_Center_Easy()

    {

        center = noH;

        xcen=0,ycen=0,zcen=0;

        for(auto&a:noH){xcen+=a.x;ycen+=a.y;zcen+=a.z;}

        xcen/=1.0*noH.size();

        ycen/=1.0*noH.size();

        zcen/=1.0*noH.size();

        for(auto&a:center){a.x-=xcen;a.y-=ycen;a.z-=zcen;}

    }



    

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

    double x_max(){return xmax;}

    double y_max(){return ymax;}

    double z_max(){return zmax;}

    double x_min(){return xmin;}

    double y_min(){return ymin;}

    double z_min(){return zmin;}

     

    private:

    double xcen=0;

    double ycen=0;

    double zcen=0;

    double xmax;

    double ymax;

    double zmax;

    double xmin;

    double ymin;

    double zmin;



};



void Input_Mol2(const parameter& P, case_info& cc, string filename)

{

    ifstream ff(filename.c_str());

    string sline;

    while (getline(ff,sline))

    {

        if(sline.find("@<TRIPOS>ATOM")!=std::string::npos)

        {

            int index=0,index_H=0;

            while (getline(ff,sline))

            {

                if(sline.find("@<TRIPOS>BOND")!=std::string::npos)break;

                std::istringstream ss(sline);

                std::string buf;

                std::vector<std::string> token;

                while(ss >> buf) token.push_back(buf);



                atom_info a;

                a.index = index_H;index_H++;

                a.name = token[1];

                a.x = stod(token[2]);

                a.y = stod(token[3]);

                a.z = stod(token[4]);

                a.type = token[5];

                a.resi_name = token[7];

                a.q = stod(token[8]);



                if(a.type.find(".")!=std::string::npos)

                {

                    string::size_type position;

                    position = a.type.find(".");

                    a.element = a.type.substr(0,position);

                }

                else

                {

                    a.element = a.type;

                }



                if(a.element == "H" || a.element == "h")

                {

                    a.rrr = P.radius_H;

                    a.sb = P.born_scale_H;

                    a.type_idx = 0;

                }

                else if(a.element == "C" || a.element == "c")

                {

                    a.rrr = P.radius_C;

                    a.sb = P.born_scale_C;

                    a.type_idx = 1;



                }

                else if(a.element == "N" || a.element == "n")

                {

                    a.rrr = P.radius_N;

                    a.sb = P.born_scale_N;

                    a.type_idx = 2;



                }

                else if(a.element == "O" || a.element == "o")

                {

                    a.rrr = P.radius_O;

                    a.sb = P.born_scale_O;

                    a.type_idx = 3;

                }

                else if(a.element == "P" || a.element == "p")

                {

                    a.rrr = P.radius_P;

                    a.sb = P.born_scale_P;

                    a.type_idx = 4;

                }

                else if(a.element == "S" || a.element == "s")

                {

                    a.rrr = P.radius_S;

                    a.sb = P.born_scale_S;

                    a.type_idx = 5;

                }

                else

                {

                    a.rrr = P.radius_C;

                    a.sb = P.born_scale_C;

                    a.type_idx = 6;

                }



                if(a.element != "H" && a.element != "h")

                {

                    a.index = index;

                    index++;

                    cc.noH.push_back(a);

                }

                else

                {



                }

                cc.wH.push_back(a);



            }

        }



    }



}



void Find_pocket(const parameter& P, const Grids_info& gd, const vector<atom_info>& rec_atom,const thread_parameter& thread_p,vector<pocket_info>& output){

    for (int i = thread_p.ibegin; i<=thread_p.iend;i++)

    for (int j = 0; j != gd.Nny; j++)

    for (int k = 0; k != gd.Nnz; k++)

    {

        double xx = i*P.step_pocket+gd.x0;

        double yy = j*P.step_pocket+gd.y0;

        double zz = k*P.step_pocket+gd.z0;

        int totalblock,blockrx=0,blocklx=0,blockry=0,blockly=0,blockrz=0,blocklz=0;

        if(gd.grids_pocket_flag[i][j][k]==0){

            for(auto&a:rec_atom){

                double disxx = (xx - a.x)*(xx - a.x);

                double disyy = (yy - a.y)*(yy - a.y);

                double diszz = (zz - a.z)*(zz - a.z);



                if(blockrx==0){

                    if(a.x < xx && xx< a.x + P.radius_large_sphere){

                        double dd = sqrt(disyy+diszz);

                        if(dd < a.rrr + P.radius_small_sphere)blockrx=1;

                    }

                }



                if(blocklx==0){

                    if(a.x > xx && xx> a.x - P.radius_large_sphere){

                        double dd = sqrt(disyy+diszz);

                        if(dd < a.rrr + P.radius_small_sphere)blocklx=1;

                    }

                }



                if(blockry==0){

                    if(a.y < yy && yy< a.y + P.radius_large_sphere){

                        double dd = sqrt(disxx+diszz);

                        if(dd < a.rrr + P.radius_small_sphere)blockry=1;

                    }

                }



                if(blockly==0){

                    if(a.y > yy && yy> a.y - P.radius_large_sphere){

                        double dd = sqrt(disxx+diszz);

                        if(dd < a.rrr + P.radius_small_sphere)blockly=1;

                    }

                }



                if(blockrz==0){

                    if(a.z < zz && zz< a.z + P.radius_large_sphere){

                        double dd = sqrt(disyy+disxx);

                        if(dd < a.rrr + P.radius_small_sphere)blockrz=1;

                    }

                }



                if(blocklz==0){

                    if(a.z > zz && zz> a.z - P.radius_large_sphere){

                        double dd = sqrt(disyy+disxx);

                        if(dd < a.rrr + P.radius_small_sphere)blocklz=1;

                    }

                }

                totalblock=blocklx+blockrx+blockly+blockry+blocklz+blockrz;

                if(totalblock==6)break;

            }

             if(totalblock==6){

                pocket_info p_tmp;

                p_tmp.x = xx;

                p_tmp.y = yy;

                p_tmp.z = zz;

                output.push_back(p_tmp);

                if(output.size()>P.Max_pocket_site/P.num_threads){cout<<"Too much pocket sites"<<endl;exit(2);}

             }



        }



    }

}



void Find_Site(const parameter& P,const vector<atom_info>centers, const vector< vector<atom_info> >& initial_poses, const thread_parameter& thread_p,const Grids_info& gd, vector<pocket_info>& pockets)

{

    // cout<<thread_p.ibegin<<" "<<thread_p.iend<<endl;

    

    for(int ic = thread_p.ibegin; ic <= thread_p.iend;ic++)

    {

        int idx_initial = ic/initial_poses[0].size();

        int label = ic % initial_poses[0].size(); 

        vector< vector<atom_info> > use_poses;

        int num_pose = P.aa.size();

        int num_ligand = initial_poses[idx_initial].size();

        // cout<<"aaaa "<<P.aa.size()<<endl;

        clock_t aaa,bbb;

        aaa=clock();

        cordinate_info void_cordinate;

        void_cordinate.index=0;

        void_cordinate.x=0;

        void_cordinate.y=0;

        void_cordinate.z=0;



        vector< vector<cordinate_info> >pp(P.aa.size(),vector<cordinate_info>(initial_poses[idx_initial].size(),void_cordinate));

        double xpose[num_pose][num_ligand],ypose[num_pose][num_ligand],zpose[num_pose][num_ligand];

        double xposeuse[num_pose][num_ligand],yposeuse[num_pose][num_ligand],zposeuse[num_pose][num_ligand];

        for(int i=0;i!=P.aa.size();i++)

        {

            vector<atom_info> tmp_pose = initial_poses[idx_initial];

            vector<cordinate_info>temp_pose( initial_poses[idx_initial].size());

            for (int j=0;j!=initial_poses[idx_initial].size();j++)

            {

                double xtmp = initial_poses[idx_initial][j].x - centers[ic].x; 

                double ytmp = initial_poses[idx_initial][j].y - centers[ic].y;

                double ztmp = initial_poses[idx_initial][j].z - centers[ic].z;

                // tmp_pose[j].x = P.aa[i]*xtmp + P.bb[i]*ytmp + P.cc[i]*ztmp;

                // tmp_pose[j].y = P.dd[i]*xtmp + P.ee[i]*ytmp + P.ff[i]*ztmp;

                // tmp_pose[j].z = P.gg[i]*xtmp + P.hh[i]*ytmp + P.ii[i]*ztmp;

                pp[i][j].x = P.aa[i]*xtmp + P.bb[i]*ytmp + P.cc[i]*ztmp;

                pp[i][j].y = P.dd[i]*xtmp + P.ee[i]*ytmp + P.ff[i]*ztmp;

                pp[i][j].z = P.gg[i]*xtmp + P.hh[i]*ytmp + P.ii[i]*ztmp;

                xpose[i][j] = P.aa[i]*xtmp + P.bb[i]*ytmp + P.cc[i]*ztmp;

                ypose[i][j] = P.dd[i]*xtmp + P.ee[i]*ytmp + P.ff[i]*ztmp;

                zpose[i][j] = P.gg[i]*xtmp + P.hh[i]*ytmp + P.ii[i]*ztmp;





            }

            bool use=true;

            

            // for(int iuse=0;iuse!=use_poses.size();iuse++)

            // {

            //     double rmsd = 0;

            //     for(int j=0;j!=tmp_pose.size();j++)

            //     {

            //         // if(i==1)cout<<"i==1 "<<j<<" "<<tmp_pose[j].x <<" old "<<use_poses[iuse][j].x<<endl;

            //         rmsd += (tmp_pose[j].x - use_poses[iuse][j].x) * (tmp_pose[j].x - use_poses[iuse][j].x);

            //         rmsd += (tmp_pose[j].y - use_poses[iuse][j].y) * (tmp_pose[j].y - use_poses[iuse][j].y);

            //         rmsd += (tmp_pose[j].z - use_poses[iuse][j].z) * (tmp_pose[j].z - use_poses[iuse][j].z);

            //     }

            //     rmsd =sqrt(rmsd/tmp_pose.size());

            //     if(rmsd<P.rmsd_cut_step1){use = false;break;}

            // }

        

            if(use)use_poses.push_back(tmp_pose);

        } 

        vector<int>pose_flag(P.aa.size(),0);

        int nn=0,lig_num = initial_poses[idx_initial].size();

        for(int i=0;i!=pp.size();i++)

        {

            bool select=true;

            for(int j=0;j!=i;j++)

            {

                if(pose_flag[j]==1)

                {

                    double rmsd=0;

                    for(int k=0;k!=lig_num;k++)

                    {

                        // rmsd+=(pp[i][k].x-pp[j][k].x)*(pp[i][k].x-pp[j][k].x);

                        // rmsd+=(pp[i][k].y-pp[j][k].y)*(pp[i][k].y-pp[j][k].y);

                        // rmsd+=(pp[i][k].z-pp[j][k].z)*(pp[i][k].z-pp[j][k].z); 

                        rmsd+=(xpose[i][k]-xpose[j][k])*(xpose[i][k]-xpose[j][k]);

                        rmsd+=(ypose[i][k]-ypose[j][k])*(ypose[i][k]-ypose[j][k]);

                        rmsd+=(zpose[i][k]-zpose[j][k])*(zpose[i][k]-zpose[j][k]); 

                        

                    }

                    rmsd=sqrt(rmsd/lig_num);

                    if(rmsd<P.rmsd_cut_step1){select=false;break;}

                   

                }

            }

            if(select){pose_flag[i]=1;nn++;}



        }

        bbb=clock();

        cout<<pp.size()<<" "<<nn<<" "<<(bbb-aaa)/CLOCKS_PER_SEC<<endl;

        // cout<<"use "<<use_poses.size()<<endl;

        for(int isite = 0;isite!=pockets.size();isite++)

        {

            

            for(auto& pose:use_poses)

            {

                bool record=true;

                vector<atom_info>tmp_pose(pose.size());

                double ljenergy = 0.;

                for(int ilig =0; ilig!=pose.size();ilig++)

                {

                    tmp_pose[ilig].x = pockets[isite].x + pose[ilig].x;

                    tmp_pose[ilig].y = pockets[isite].y + pose[ilig].y;

                    tmp_pose[ilig].z = pockets[isite].z + pose[ilig].z;



                    int xidx = (tmp_pose[ilig].x + gd.igd*0.5)/gd.igd;

                    int yidx = (tmp_pose[ilig].y + gd.igd*0.5)/gd.igd;

                    int zidx = (tmp_pose[ilig].z + gd.igd*0.5)/gd.igd;



                    // cout<<xidx<<" "<<yidx<<" "<<zidx<<endl;

                    // cout<<gd.lj_flag[pose[ilig].type_idx][xidx][yidx][zidx]<<endl;



                    if(gd.lj_flag[pose[ilig].type_idx][xidx][yidx][zidx] != 0){record=false;break;}

                    ljenergy += gd.grids_energy[pose[ilig].type_idx][xidx][yidx][zidx];

                    if(ljenergy >50.0){record=false;break;}

                }exit(3);

                if(record){

                    if(ljenergy < pockets[isite].LJ_min[idx_initial][label].lj)

                    {

                        cout<<ljenergy<<endl;

                        pockets[isite].LJ_min[idx_initial][label].lj = ljenergy;

                        if(ljenergy < pockets[isite].lj_site_min)

                        {

                            pockets[isite].lj_site_min= ljenergy;

                        }

                    }   

                }                

            }

        }

    }



}



void LJ(const parameter& P, vector< vector< vector<double> > >& LJ_potential)

{

    ////////LJ index start from 1 in calculation 0 can not be divined.

    for(int i=0;i!=P.atom_radius.size();i++)

    {

        vector< vector<double> > itmp_lj;

        for(int j=0;j!=P.atom_radius.size();j++)

        {

            double rdis = P.equ_LJ*(P.atom_radius[i] + P.atom_radius[j]);

            // double dr = P.LJ_cut*(P.atom_radius[i] + P.atom_radius[j]);

            // if(dr<P.LJ_cut_new)dr = P.LJ_cut_new;

            // dr /= P.LJstep;



            vector<double> tmp_lj;

            for( int istep=1; istep<=P.LJstep;istep++)

            {

                double dis = istep*P.lj_dr;

                double dx = rdis/dis;

                double x6 = dx*dx*dx*dx*dx*dx;

                double var = x6*x6 - x6;

                tmp_lj.push_back(var);

            }

            itmp_lj.push_back(tmp_lj);

        }

        LJ_potential.push_back(itmp_lj);

    }

    



}





int main(int argc, char** argv )

{

    std::cout << "Have " << argc << " arguments:" << std::endl;

    vector<string>argv_str;

    for (int i = 0; i < argc; ++i){

        std::cout << argv[i] << std::endl;

        string argv_tmp = argv[i];

        argv_str.push_back(argv_tmp);

    }

    parameter P;

    P.SetSphere("../list/sphere.dat");

    P.SetAngle();

    cout<<"Parameter initialization"<<endl;

    //////////////////

    clock_t start,mid1,mid2,end;

    start=clock();



    case_info lig;

    string filename = "input/"+argv_str[1]+"/"+argv_str[1]+"_ligand.mol2";

    lig.SetValue(P,filename);

    if(lig.HH.size()>0)lig.H_charge_transfer();

    lig.Move_to_Center_Easy();

    cout<<"native Ligand input finish "<<lig.wH.size()<<" "<<lig.noH.size()<<endl;



    case_info rec;

    filename = "input/"+argv_str[1]+"/"+argv_str[1]+"_RNA.mol2";

    rec.SetValue(P,filename);

    rec.Find_Max_and_Min();

    if(rec.HH.size()>0)rec.H_charge_transfer();

    cout<<"Receptor input finish "<<rec.wH.size()<<" "<<rec.noH.size()<<endl;



    //////////Grid for binding site

   

    Grids_info Grid;

    Grid.Nnx = int((rec.x_max()+2*P.r_water - rec.x_min() + 2*P.r_water)/P.step_pocket + 0.1);

    Grid.Nny = int((rec.y_max()+2*P.r_water - rec.y_min() + 2*P.r_water)/P.step_pocket + 0.1);

    Grid.Nnz = int((rec.z_max()+2*P.r_water - rec.z_min() + 2*P.r_water)/P.step_pocket + 0.1);

    Grid.x0 = rec.x_min() - 2*P.r_water + P.step_pocket;

    Grid.y0 = rec.y_min() - 2*P.r_water + P.step_pocket;

    Grid.z0 = rec.z_min() - 2*P.r_water + P.step_pocket;

    Grid.Set_grid_int();

    cout<<"Grid for binding site initialization "<<Grid.Nnx<<" "<<Grid.Nny<<" "<<Grid.Nnz<<" "<<Grid.x0<<" "<<endl;

    



    //////kick out girds on RNA

    for(auto & rr:rec.noH)

    {

        int xmin= (rr.x-Grid.x0-(rr.rrr+P.radius_small_sphere))/P.step_pocket;

        if(xmin<0)xmin=0;

        int xmax= (rr.x-Grid.x0+(rr.rrr+P.radius_small_sphere))/P.step_pocket+1;

        if(xmax>=Grid.Nnx)xmax=Grid.Nnx-1;



        int ymin= (rr.y-Grid.y0-(rr.rrr+P.radius_small_sphere))/P.step_pocket;

        if(ymin<0)ymin=0;

        int ymax= (rr.y-Grid.y0+(rr.rrr+P.radius_small_sphere))/P.step_pocket+1;

        if(ymax>=Grid.Nny)ymax=Grid.Nny-1;



        int zmin= (rr.z-Grid.z0-(rr.rrr+P.radius_small_sphere))/P.step_pocket;

        if(zmin<0)zmin=0;

        int zmax= (rr.z-Grid.z0+(rr.rrr+P.radius_small_sphere))/P.step_pocket+1;

        if(zmax>=Grid.Nnz)zmax=Grid.Nnz-1;



        Grid.Find_RNA_grid(P,xmin,xmax,ymin,ymax,zmin,zmax,rr);

    }

    cout<<"Grid for binding site Finish"<<endl;



    ////////////////// 上下左右前后贯通

    thread t[P.num_threads];

    int num_per_thread = Grid.Nnx / P.num_threads;

    int num_plus = Grid.Nnx % P.num_threads;

    vector<thread_parameter> thread_p(P.num_threads);



    for(unsigned int i = 1; i != P.num_threads; ++i)

    {

        thread_p[i-1].iend =  thread_p[i-1].ibegin + num_per_thread-1;

        if(i<num_plus) thread_p[i-1].iend ++;

        thread_p[i].ibegin = thread_p[i-1].iend+1;

    }

    thread_p[P.num_threads-1].iend = Grid.Nnx -1;

    vector< vector<pocket_info> > pocket_sub(P.num_threads);

    vector<pocket_info> pockets;



    

    for(unsigned int i = 0; i != P.num_threads; ++i)

    {

        t[i] = std::thread(Find_pocket,std::ref(P), std::ref(Grid),std::ref(rec.noH),std::ref(thread_p[i]),std::ref(pocket_sub[i]));

    }    

    for(unsigned int i = 0; i != P.num_threads; ++i)

    {

        t[i].join();

        pockets.insert(pockets.end(),pocket_sub[i].begin(),pocket_sub[i].end());

    }

    cout<<"Find "<<pockets.size()<<" pockets"<<endl;

