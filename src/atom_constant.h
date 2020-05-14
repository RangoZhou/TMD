#pragma once

#include <string>
#include <array>
#include <map>

#include "common.h"

namespace tmd {

    enum ELEMENT_TYPE {
        EL_H=0,  EL_HE,  EL_LI,  EL_BE,   EL_B,   EL_C,   EL_N,   EL_O,   EL_F,  EL_NE,
        EL_NA,   EL_MG,  EL_AL,  EL_SI,   EL_P,   EL_S,  EL_CL,  EL_AR,   EL_K,  EL_CA,
        EL_SC,   EL_TI,   EL_V,  EL_CR,  EL_MN,  EL_FE,  EL_CO,  EL_NI,  EL_CU,  EL_ZN,
        EL_GA,   EL_GE,  EL_AS,  EL_SE,  EL_BR,  EL_KR,  EL_RB,  EL_SR,   EL_Y,  EL_ZR,
        EL_NB,   EL_MO,  EL_TC,  EL_RU,  EL_RH,  EL_PD,  EL_AG,  EL_CD,  EL_IN,  EL_SN,
        EL_SB,   EL_TE,   EL_I,  EL_XE,  EL_CS,  EL_BA,  EL_LA,  EL_CE,  EL_PR,  EL_ND,
        EL_PM,   EL_SM,  EL_EU,  EL_GD,  EL_TB,  EL_DY,  EL_HO,  EL_ER,  EL_TM,  EL_YB,
        EL_LU,   EL_HF,  EL_TA,   EL_W,  EL_RE,  EL_OS,  EL_IR,  EL_PT,  EL_AU,  EL_HG,
        EL_TL,   EL_PB,  EL_BI,  EL_PO,  EL_AT,  EL_RN,  EL_FR,  EL_RA,  EL_AC,  EL_TH,
        EL_PA,    EL_U,  EL_NP,  EL_PU,  EL_AM,  EL_CM,  EL_BK,  EL_CF,  EL_ES,  EL_FM,
        EL_MD,   EL_NO,  EL_LR,  EL_RF,  EL_DB,  EL_SG,  EL_BH,  EL_HS,  EL_MT,  EL_DS,
        EL_RG,   EL_CN, EL_UNK, EL_SIZE};

    struct Element_Info {
        const std::string name;
        const Size_Type atomic_number;
        const Float covalent_radius;
        const Float atomic_radius;
        const Float vdw_radius;
        Element_Info(const std::string n, const Size_Type an, const Float cr, const Float ar, const Float vr) : name(n), atomic_number(an), covalent_radius(cr), atomic_radius(ar), vdw_radius(vr) {}
    };
    //All Info comes from this webpage --->  http://periodictable.com/index.html
    //information below is for vdw radius
    //The newly determined radii_ (in Å) are:
    //Be, 1.53; B, 1.92; Al, 1.84; Ca, 2.31; Ge, 2.11; Rb, 3.03; Sr, 2.50; Sb, 2.06; Cs, 3.43; Ba, 2.68; Bi, 2.07; Po, 1.97; At, 2.02; Rn, 2.20; Fr, 3.48; and Ra, 2.83.
    //J Phys Chem A. 2009 May 14; 113(19): 5806–5812. doi:10.1021/jp8111556.
    //Consistent van der Waals Radii for the Whole Main Group
    //Manjeera Mantina, Adam C. Chamberlin, Rosendo Valero, Christopher J. Cramer, and Donald G. Truhlar *
    //Department of Chemistry and Supercomputing Institute, University of Minnesota, Minneapolis,MN 55455-0431
    const Float unkw = 1.70;//unkown radius
    const std::array<Element_Info,EL_SIZE> element_infos = {{
        //name, atomic_number, covalent_radius, atomic_radius, vdw_radius
        { "H",   1, 0.31, 0.53, 1.20},// EL_H
        {"He",   2, 0.28, 0.31, 1.40},//EL_HE
        {"Li",   3, 1.28, 1.67, 1.82},//EL_LI
        {"Be",   4, 0.96, 1.12, 1.53},//EL_BE
        { "B",   5, 0.84, 0.87, 1.92},// EL_B
        { "C",   6, 0.76, 0.67, 1.70},// EL_C
        { "N",   7, 0.71, 0.56, 1.55},// EL_N
        { "O",   8, 0.66, 0.48, 1.52},// EL_O
        { "F",   9, 0.57, 0.42, 1.47},// EL_F
        {"Ne",  10, 0.58, 0.38, 1.54},//EL_NE
        {"Na",  11, 1.66, 1.90, 2.27},//EL_NA
        {"Mg",  12, 1.41, 1.45, 1.73},//EL_MG
        {"Al",  13, 1.21, 1.18, 1.84},//EL_AL
        {"Si",  14, 1.11, 1.11, 2.10},//EL_SI
        { "P",  15, 1.07, 0.98, 1.80},// EL_P
        { "S",  16, 1.05, 0.88, 1.80},// EL_S
        {"Cl",  17, 1.02, 0.79, 1.75},//EL_CL
        {"Ar",  18, 1.06, 0.71, 1.88},//EL_AR
        { "K",  19, 2.03, 2.43, 2.75},// EL_K
        {"Ca",  20, 1.76, 1.94, 2.31},//EL_CA
        {"Sc",  21, 1.70, 1.84, unkw},//EL_SC
        {"Ti",  22, 1.60, 1.76, unkw},//EL_TI
        { "V",  23, 1.53, 1.71, unkw},// EL_V
        {"Cr",  24, 1.39, 1.66, unkw},//EL_CR
        {"Mn",  25, 1.39, 1.61, unkw},//EL_MN
        {"Fe",  26, 1.32, 1.56, unkw},//EL_FE
        {"Co",  27, 1.26, 1.52, unkw},//EL_CO
        {"Ni",  28, 1.24, 1.49, 1.63},//EL_NI
        {"Cu",  29, 1.32, 1.45, 1.40},//EL_CU
        {"Zn",  30, 1.22, 1.42, 1.39},//EL_ZN
        {"Ga",  31, 1.22, 1.36, 1.87},//EL_GA
        {"Ge",  32, 1.20, 1.25, 2.11},//EL_GE
        {"As",  33, 1.19, 1.14, 1.85},//EL_AS
        {"Se",  34, 1.20, 1.03, 1.90},//EL_SE
        {"Br",  35, 1.20, 0.94, 1.85},//EL_BR
        {"Kr",  36, 1.16, 0.88, 2.02},//EL_KR
        {"Rb",  37, 2.20, 2.65, 3.03},//EL_RB
        {"Sr",  38, 1.95, 2.19, 2.49},//EL_SR
        { "Y",  39, 1.90, 2.12, unkw},// EL_Y
        {"Zr",  40, 1.75, 2.06, unkw},//EL_ZR
        {"Nb",  41, 1.64, 1.98, unkw},//EL_NB
        {"Mo",  42, 1.54, 1.90, unkw},//EL_MO
        {"Tc",  43, 1.47, 1.83, unkw},//EL_TC
        {"Ru",  44, 1.46, 1.78, unkw},//EL_RU
        {"Rh",  45, 1.42, 1.73, unkw},//EL_RH
        {"Pd",  46, 1.39, 1.69, 1.63},//EL_PD
        {"Ag",  47, 1.45, 1.65, 1.72},//EL_AG
        {"Cd",  48, 1.44, 1.61, 1.58},//EL_CD
        {"In",  49, 1.42, 1.56, 1.93},//EL_IN
        {"Sn",  50, 1.39, 1.45, 2.17},//EL_SN
        {"Sb",  51, 1.39, 1.33, 2.06},//EL_SB
        {"Te",  52, 1.38, 1.23, 2.06},//EL_TE
        { "I",  53, 1.39, 1.15, 1.98},// EL_I
        {"Xe",  54, 1.40, 1.08, 2.16},//EL_XE
        {"Cs",  55, 2.44, 2.98, 3.43},//EL_CS
        {"Ba",  56, 2.15, 2.53, 2.68},//EL_BA
        {"La",  57, 2.07, unkw, unkw},//EL_LA
        {"Ce",  58, 2.04, unkw, unkw},//EL_CE
        {"Pr",  59, 2.03, 2.47, unkw},//EL_PR
        {"Nd",  60, 2.01, 2.06, unkw},//EL_ND
        {"Pm",  61, 1.99, 2.05, unkw},//EL_PM
        {"Sm",  62, 1.98, 2.38, unkw},//EL_SM
        {"Eu",  63, 1.98, 2.31, unkw},//EL_EU
        {"Gd",  64, 1.96, 2.33, unkw},//EL_GD
        {"Tb",  65, 1.94, 2.25, unkw},//EL_TB
        {"Dy",  66, 1.92, 2.28, unkw},//EL_DY
        {"Ho",  67, 1.92, 2.26, unkw},//EL_HO
        {"Er",  68, 1.89, 2.26, unkw},//EL_ER
        {"Tm",  69, 1.90, 2.22, unkw},//EL_TM
        {"Yb",  70, 1.87, 2.22, unkw},//EL_YB
        {"Lu",  71, 1.87, 2.17, unkw},//EL_LU
        {"Hf",  72, 1.75, 2.08, unkw},//EL_HF
        {"Ta",  73, 1.70, 2.00, unkw},//EL_TA
        { "W",  74, 1.62, 1.93, unkw},// EL_W
        {"Re",  75, 1.51, 1.88, unkw},//EL_RE
        {"Os",  76, 1.44, 1.85, unkw},//EL_OS
        {"Ir",  77, 1.41, 1.80, unkw},//EL_IR
        {"Pt",  78, 1.36, 1.77, 1.75},//EL_PT
        {"Au",  79, 1.36, 1.74, 1.66},//EL_AU
        {"Hg",  80, 1.32, 1.71, 1.55},//EL_HG
        {"Tl",  81, 1.45, 1.56, 1.96},//EL_TL
        {"Pb",  82, 1.46, 1.54, 2.02},//EL_PB
        {"Bi",  83, 1.48, 1.43, 2.07},//EL_BI
        {"Po",  84, 1.40, 1.35, 1.97},//EL_PO
        {"At",  85, 1.50, 1.27, 2.02},//EL_AT
        {"Rn",  86, 1.50, 1.20, 2.20},//EL_RN
        {"Fr",  87, 2.60, unkw, 3.48},//EL_FR
        {"Ra",  88, 2.21, unkw, 2.83},//EL_RA
        {"Ac",  89, 2.15, unkw, unkw},//EL_AC
        {"Th",  90, 2.06, unkw, unkw},//EL_TH
        {"Pa",  91, 2.00, unkw, unkw},//EL_PA
        { "U",  92, 1.96, unkw, 1.86},// EL_U
        {"Np",  93, 1.90, unkw, unkw},//EL_NP
        {"Pu",  94, 1.87, unkw, unkw},//EL_PU
        {"Am",  95, 1.80, unkw, unkw},//EL_AM
        {"Cm",  96, 1.69, unkw, unkw},//EL_CM
        {"Bk",  97, unkw, unkw, unkw},//EL_BK
        {"Cf",  98, unkw, unkw, unkw},//EL_CF
        {"Es",  99, unkw, unkw, unkw},//EL_ES
        {"Fm", 100, unkw, unkw, unkw},//EL_FM
        {"Md", 101, unkw, unkw, unkw},//EL_MD
        {"No", 102, unkw, unkw, unkw},//EL_NO
        {"Lr", 103, unkw, unkw, unkw},//EL_LR
        {"Rf", 104, unkw, unkw, unkw},//EL_RF
        {"Db", 105, unkw, unkw, unkw},//EL_DB
        {"Sg", 106, unkw, unkw, unkw},//EL_SG
        {"Bh", 107, unkw, unkw, unkw},//EL_BH
        {"Hs", 108, unkw, unkw, unkw},//EL_HS
        {"Mt", 109, unkw, unkw, unkw},//EL_MT
        {"Ds", 110, unkw, unkw, unkw},//EL_DS
        {"Rg", 111, unkw, unkw, unkw},//EL_RG
        {"Cn", 112, unkw, unkw, unkw},//EL_CN
        {"UNK",113, 0.76, 0.67, 1.70} //EL_UNK, set unkown as Carbon
    }};

    // inline const Float element_type_covalent_radius(const ELEMENT_TYPE et) {
    //     return element_infos.at(et).covalent_radius;
    // }
    // inline const Float element_type_atomic_radius(const ELEMENT_TYPE et) {
    //     return element_infos.at(et).atomic_radius;
    // }
    // inline const Float element_type_vdw_radius(const ELEMENT_TYPE et) {
    //     return element_infos.at(et).vdw_radius;
    // }
    // inline const std::string element_type_to_string(const ELEMENT_TYPE et) {
    //     return element_infos.at(et).name;
    // }

    enum SYBYL_TYPE {
        SYBYL_H = 0, SYBYL_H_spc, SYBYL_H_t3p, SYBYL_C_3  , SYBYL_C_2  , SYBYL_C_1  , SYBYL_C_ar ,
        SYBYL_C_cat, SYBYL_N_3  , SYBYL_N_2  , SYBYL_N_1  , SYBYL_N_ar , SYBYL_N_am , SYBYL_N_pl3, SYBYL_N_4  ,
        SYBYL_O_3  , SYBYL_O_2  , SYBYL_O_co2, SYBYL_O_spc, SYBYL_O_t3p, SYBYL_S_3  , SYBYL_S_2  , SYBYL_S_O  ,
        SYBYL_S_O2 , SYBYL_P_3  , SYBYL_F    , SYBYL_Cl   , SYBYL_Br   , SYBYL_I    , SYBYL_Li   , SYBYL_Na   ,
        SYBYL_Mg   , SYBYL_Al   , SYBYL_Si   , SYBYL_K    , SYBYL_Ca   , SYBYL_Cr_th, SYBYL_Cr_oh, SYBYL_Mn   ,
        SYBYL_Fe   , SYBYL_Co_oh, SYBYL_Cu   , SYBYL_Zn   , SYBYL_Se   , SYBYL_Mo   , SYBYL_Sn   , SYBYL_LP   ,
        SYBYL_Du   , SYBYL_Du_C , SYBYL_Any  , SYBYL_Hal  , SYBYL_Het  , SYBYL_Hev  , SYBYL_UNK  , SYBYL_SIZE};

        //mark's sybyl type
        // TYPE_H=1,    TYPE_C_2,   TYPE_C_3,  TYPE_C_ar,  TYPE_C_cat,  TYPE_N_2,
        // TYPE_N_3,    TYPE_N_4,  TYPE_N_ar,  TYPE_N_am,  TYPE_N_pl3,  TYPE_O_2,
        // TYPE_O_3,  TYPE_O_co2,   TYPE_S_2,   TYPE_S_3,    TYPE_P_3,    TYPE_F,
        //  TYPE_Cl,     TYPE_Br,     TYPE_I,   TYPE_UNK

    struct Sybyl_Info {
        const std::string name;
        const ELEMENT_TYPE EL_TYPE;
        const Float well_depth;
        const Float vdw_radius;
        const Float solvation;
        const Float vdw_volume;
        const Float covalent_radius;
        Sybyl_Info(const std::string n, const ELEMENT_TYPE et, const Float wd, const Float vr, const Float s, const Float vv, const Float cr) : name(n), EL_TYPE(et), well_depth(wd), vdw_radius(vr), solvation(s), vdw_volume(vv), covalent_radius(cr) {}
    };

    // from autodock vina's atom constants
    // use carbon parameter for atom not showed in autodock and not metal
    const ELEMENT_TYPE SYBYL_EL = EL_C;
    const Float sybyl_iwd = 0.15000;//well depth for atom not showed in autodock vina
    const Float sybyl_ivr = 2.00000;//vdw radius for atom not showed in autodock vina
    const Float sybyl_isp = -0.00143;//solvation_parameter for atom not showed in autodock vina
    const Float sybyl_ivv = 3.1415926*4.0/3.0*sybyl_ivr*sybyl_ivr*sybyl_ivr;//vdw volume for atom not showed int autodock vina
    const Float sybyl_icr = 0.77;//covalent radius for atom not showed int autodock vina

    //use Mg parameter
    const Float sybyl_mwd = 0.87500;//metal_well depth for atom not showed in autodock vina
    const Float sybyl_mvr = 0.65000;//metal_vdw radius for atom not showed in autodock vina
    const Float sybyl_msp = -0.00110;//metal_solvation_parameter for metal not in autodock vina
    const Float sybyl_mvv = 1.56000;//metal_vdw volume for atom not showed int autodock vina, does not equal to 3.1415926*4.0/3.0*sybyl_mvr*sybyl_mvr*sybyl_mvr
    const Float sybyl_mcr = 1.75;//metal_covalent_radius for metal not in autodock vina, this is not comes from Mg but from autodock vina
    // Se is equivalent to general S
    // solvation parameter for C and A in autodock are different C -0.00143 and -0.00052, but here we use C
    const std::array<Sybyl_Info,SYBYL_SIZE> sybyl_infos = {{
        //name  element_type  well_depth  vdw_radius  solvation  vdw_volume  covalent_radius
        {"H"    ,     EL_H,   0.02000,   1.00000,   0.00051,   0.00000,      0.37},//hydrogen
        {"H.spc",     EL_H,   0.02000,   1.00000,   0.00051,   0.00000,      0.37},//hydrogen in Single Point Charge (SPC) water model
        {"H.t3p",     EL_H,   0.02000,   1.00000,   0.00051,   0.00000,      0.37},//hydrogen in Transferable intermolecular Potential (TIP3P) water model
        {"C.3"  ,     EL_C,   0.15000,   2.00000,  -0.00143,  33.51030,      0.77},//carbon sp3
        {"C.2"  ,     EL_C,   0.15000,   2.00000,  -0.00143,  33.51030,      0.77},//carbon sp2
        {"C.1"  ,     EL_C,   0.15000,   2.00000,  -0.00143,  33.51030,      0.77},//carbon sp
        {"C.ar" ,     EL_C,   0.15000,   2.00000,  -0.00143,  33.51030,      0.77},//carbon aromatic
        {"C.cat",     EL_C,   0.15000,   2.00000,  -0.00143,  33.51030,      0.77},//carbocation (C+) used only in a guadinium group
        {"N.3"  ,     EL_N,   0.16000,   1.75000,  -0.00162,  22.44930,      0.75},//nitrogen sp3
        {"N.2"  ,     EL_N,   0.16000,   1.75000,  -0.00162,  22.44930,      0.75},//nitrogen sp2
        {"N.1"  ,     EL_N,   0.16000,   1.75000,  -0.00162,  22.44930,      0.75},//nitrogen sp
        {"N.ar" ,     EL_N,   0.16000,   1.75000,  -0.00162,  22.44930,      0.75},//nitrogen aromatic
        {"N.am" ,     EL_N,   0.16000,   1.75000,  -0.00162,  22.44930,      0.75},//nitrogen amide
        {"N.pl3",     EL_N,   0.16000,   1.75000,  -0.00162,  22.44930,      0.75},//nitrogen trigonal planar
        {"N.4"  ,     EL_N,   0.16000,   1.75000,  -0.00162,  22.44930,      0.75},//nitrogen sp3 positively charged
        {"O.3"  ,     EL_O,   0.20000,   1.60000,  -0.00251,  17.15730,      0.73},//oxygen sp3
        {"O.2"  ,     EL_O,   0.20000,   1.60000,  -0.00251,  17.15730,      0.73},//oxygen sp2
        {"O.co2",     EL_O,   0.20000,   1.60000,  -0.00251,  17.15730,      0.73},//oxygen in carboxylate and phosphate groups
        {"O.spc",     EL_O,   0.20000,   1.60000,  -0.00251,  17.15730,      0.73},//oxygen in Single Point Charge (SPC) water model
        {"O.t3p",     EL_O,   0.20000,   1.60000,  -0.00251,  17.15730,      0.73},//oxygen in Transferable Intermolecular Potential (TIP3P) water model
        {"S.3"  ,     EL_S,   0.20000,   2.00000,  -0.00214,  33.51030,      1.02},//sulfur sp3
        {"S.2"  ,     EL_S,   0.20000,   2.00000,  -0.00214,  33.51030,      1.02},//sulfur sp2
        {"S.O"  ,     EL_S,   0.20000,   2.00000,  -0.00214,  33.51030,      1.02},//sulfoxide sulfur
        {"S.O2" ,     EL_S,   0.20000,   2.00000,  -0.00214,  33.51030,      1.02},//sulfone sulfur
        {"P.3"  ,     EL_P,   0.20000,   2.10000,  -0.00110,  38.79240,      1.06},//phosphorous sp3
        {"F"    ,     EL_F,   0.08000,   1.54500,  -0.00110,  15.44800,      0.71},//fluorine
        {"Cl"   ,    EL_CL,   0.27600,   2.04500,  -0.00110,  35.82350,      0.99},//chlorine
        {"Br"   ,    EL_BR,   0.38900,   2.16500,  -0.00110,  42.56610,      1.14},//bromine
        {"I"    ,     EL_I,   0.55000,   2.36000,  -0.00110,  55.05850,      1.33},//iodine
        {"Li"   ,    EL_LI, sybyl_mwd, sybyl_mvr, sybyl_msp, sybyl_mvv, sybyl_mcr},//lithium
        {"Na"   ,    EL_NA, sybyl_mwd, sybyl_mvr, sybyl_msp, sybyl_mvv, sybyl_mcr},//sodium
        {"Mg"   ,    EL_MG,   0.87500,   0.65000,  -0.00110,   1.56000,      1.30},//magnesium
        {"Al"   ,    EL_AL, sybyl_mwd, sybyl_mvr, sybyl_msp, sybyl_mvv, sybyl_mcr},//aluminum
        {"Si"   ,    EL_SI, sybyl_mwd, sybyl_mvr, sybyl_msp, sybyl_mvv, sybyl_mcr},//silicon
        {"K"    ,     EL_K, sybyl_mwd, sybyl_mvr, sybyl_msp, sybyl_mvv, sybyl_mcr},//potassium
        {"Ca"   ,    EL_CA,   0.55000,   0.99000,  -0.00110,   2.77000,      1.74},//calcium
        {"Cr.th",    EL_CR, sybyl_mwd, sybyl_mvr, sybyl_msp, sybyl_mvv, sybyl_mcr},//chromium (tetrahedral)
        {"Cr.oh",    EL_CR, sybyl_mwd, sybyl_mvr, sybyl_msp, sybyl_mvv, sybyl_mcr},//chromium (octahedral)
        {"Mn"   ,    EL_MN,   0.87500,   0.65000,  -0.00110,   2.14000,      1.39},//manganese
        {"Fe"   ,    EL_FE,   0.01000,   0.65000,  -0.00110,   1.84000,      1.25},//iron
        {"Co.oh",    EL_CO, sybyl_mwd, sybyl_mvr, sybyl_msp, sybyl_mvv, sybyl_mcr},//cobalt (octahedral)
        {"Cu"   ,    EL_CU, sybyl_mwd, sybyl_mvr, sybyl_msp, sybyl_mvv, sybyl_mcr},//copper
        {"Zn"   ,    EL_ZN,   0.55000,   0.74000,  -0.00110,   1.70000,      1.31},//zinc
        {"Se"   ,    EL_SE,   0.20000,   2.00000,  -0.00214,  33.51030,      1.02},//selenium
        {"Mo"   ,    EL_MO, sybyl_mwd, sybyl_mvr, sybyl_msp, sybyl_mvv, sybyl_mcr},//molybdenum
        {"Sn"   ,    EL_SN, sybyl_mwd, sybyl_mvr, sybyl_msp, sybyl_mvv, sybyl_mcr},//tin
        {"LP"   , SYBYL_EL, sybyl_iwd, sybyl_ivr, sybyl_isp, sybyl_ivv, sybyl_icr},//lone pair
        {"Du"   , SYBYL_EL, sybyl_iwd, sybyl_ivr, sybyl_isp, sybyl_ivv, sybyl_icr},//dummy atom
        {"Du.C" , SYBYL_EL, sybyl_iwd, sybyl_ivr, sybyl_isp, sybyl_ivv, sybyl_icr},//dummy carbon
        {"Any"  , SYBYL_EL, sybyl_iwd, sybyl_ivr, sybyl_isp, sybyl_ivv, sybyl_icr},//any atom
        {"Hal"  , SYBYL_EL, sybyl_iwd, sybyl_ivr, sybyl_isp, sybyl_ivv, sybyl_icr},//halogen
        {"Het"  , SYBYL_EL, sybyl_iwd, sybyl_ivr, sybyl_isp, sybyl_ivv, sybyl_icr},//heteroatom = N, O, S, P
        {"Hev"  , SYBYL_EL, sybyl_iwd, sybyl_ivr, sybyl_isp, sybyl_ivv, sybyl_icr},//heavy atom (non hydrogen)
        {"UNK"  ,   EL_UNK, sybyl_iwd, sybyl_ivr, sybyl_isp, sybyl_ivv, sybyl_icr}
    }};

    const Size_Type sybyl_infos_size =  sybyl_infos.size();

    const std::map<std::string,SYBYL_TYPE> sybyl_type_lookup = {
        // name, SYBYL_TYPE
        {"H"    , SYBYL_H    },
        {"H.spc", SYBYL_H_spc},
        {"H.t3p", SYBYL_H_t3p},
        {"C.3"  , SYBYL_C_3  },
        {"C.2"  , SYBYL_C_2  },
        {"C.1"  , SYBYL_C_1  },
        {"C.ar" , SYBYL_C_ar },
        {"C.cat", SYBYL_C_cat},
        {"N.3"  , SYBYL_N_3  },
        {"N.2"  , SYBYL_N_2  },
        {"N.1"  , SYBYL_N_1  },
        {"N.ar" , SYBYL_N_ar },
        {"N.am" , SYBYL_N_am },
        {"N.pl3", SYBYL_N_pl3},
        {"N.4"  , SYBYL_N_4  },
        {"O.3"  , SYBYL_O_3  },
        {"O.2"  , SYBYL_O_2  },
        {"O.co2", SYBYL_O_co2},
        {"O.spc", SYBYL_O_spc},
        {"O.t3p", SYBYL_O_t3p},
        {"S.3"  , SYBYL_S_3  },
        {"S.2"  , SYBYL_S_2  },
        {"S.O"  , SYBYL_S_O  },
        {"S.O2" , SYBYL_S_O2 },
        {"P.3"  , SYBYL_P_3  },
        {"F"    , SYBYL_F    },
        {"Cl"   , SYBYL_Cl   },
        {"Br"   , SYBYL_Br   },
        {"I"    , SYBYL_I    },
        {"Li"   , SYBYL_Li   },
        {"Na"   , SYBYL_Na   },
        {"Mg"   , SYBYL_Mg   },
        {"Al"   , SYBYL_Al   },
        {"Si"   , SYBYL_Si   },
        {"K"    , SYBYL_K    },
        {"Ca"   , SYBYL_Ca   },
        {"Cr.th", SYBYL_Cr_th},
        {"Cr.oh", SYBYL_Cr_oh},
        {"Mn"   , SYBYL_Mn   },
        {"Fe"   , SYBYL_Fe   },
        {"Co.oh", SYBYL_Co_oh},
        {"Cu"   , SYBYL_Cu   },
        {"Zn"   , SYBYL_Zn   },
        {"Se"   , SYBYL_Se   },
        {"Mo"   , SYBYL_Mo   },
        {"Sn"   , SYBYL_Sn   },
        {"LP"   , SYBYL_LP   },
        {"Du"   , SYBYL_Du   },
        {"Du.C" , SYBYL_Du_C },
        {"Any"  , SYBYL_Any  },
        {"Hal"  , SYBYL_Hal  },
        {"Het"  , SYBYL_Het  },
        {"Hev"  , SYBYL_Hev  },
        {"UNK"  , SYBYL_UNK  }
    };





    // AutoDock4
    enum AD4_TYPE {
        AD4_C=0,
        AD4_A  ,
        AD4_N  ,
        AD4_O  ,
        AD4_P  ,
        AD4_S  ,
        AD4_SE ,// equivalent to AD4_S
        AD4_H  ,// non-polar hydrogen
        AD4_F  ,
        AD4_I  ,
        AD4_NA ,
        AD4_OA ,
        AD4_SA ,
        AD4_HD ,
        AD4_Mg ,
        AD4_Mn ,
        AD4_Zn ,
        AD4_Ca ,
        AD4_Fe ,
        AD4_Cl ,
        AD4_Br ,
        AD4_UNK,
        AD4_SIZE};


struct AD4_Info {
        const std::string name;
        const ELEMENT_TYPE EL_TYPE;
        const Float well_depth;
        const Float vdw_radius;
        const Float solvation;
        const Float vdw_volume;
        const Float covalent_radius;
        AD4_Info(const std::string n, const ELEMENT_TYPE et, const Float wd, const Float vr, const Float s, const Float vv, const Float cr) : name(n), EL_TYPE(et), well_depth(wd), vdw_radius(vr), solvation(s), vdw_volume(vv), covalent_radius(cr) {}
    };

const std::array<AD4_Info,AD4_SIZE> ad4_infos = {{
    // name, depth, vdw radius, solvation parameter, volume, covalent radius
	{ "C",    EL_C,  0.15000,  2.00000,  -0.00143,  33.51030,  0.77}, //
	{ "A",    EL_C,  0.15000,  2.00000,  -0.00052,  33.51030,  0.77}, //
	{ "N",    EL_N,  0.16000,  1.75000,  -0.00162,  22.44930,  0.75}, //
	{ "O",    EL_O,  0.20000,  1.60000,  -0.00251,  17.15730,  0.73}, //
	{ "P",    EL_P,  0.20000,  2.10000,  -0.00110,  38.79240,  1.06}, //
	{ "S",    EL_S,  0.20000,  2.00000,  -0.00214,  33.51030,  1.02}, //
    {"Se",   EL_SE,  0.20000,  2.00000,  -0.00214,  33.51030,  1.02}, //
	{ "H",    EL_H,  0.02000,  1.00000,   0.00051,   0.00000,  0.37}, //
	{ "F",    EL_F,  0.08000,  1.54500,  -0.00110,  15.44800,  0.71}, //
	{ "I",    EL_I,  0.55000,  2.36000,  -0.00110,  55.05850,  1.33}, //
	{"NA",    EL_N,  0.16000,  1.75000,  -0.00162,  22.44930,  0.75}, //
	{"OA",    EL_O,  0.20000,  1.60000,  -0.00251,  17.15730,  0.73}, //
	{"SA",    EL_S,  0.20000,  2.00000,  -0.00214,  33.51030,  1.02}, //
	{"HD",    EL_H,  0.02000,  1.00000,   0.00051,   0.00000,  0.37}, //
	{"Mg",   EL_MG,  0.87500,  0.65000,  -0.00110,   1.56000,  1.30}, //
	{"Mn",   EL_MN,  0.87500,  0.65000,  -0.00110,   2.14000,  1.39}, //
	{"Zn",   EL_ZN,  0.55000,  0.74000,  -0.00110,   1.70000,  1.31}, //
	{"Ca",   EL_CA,  0.55000,  0.99000,  -0.00110,   2.77000,  1.74}, //
	{"Fe",   EL_FE,  0.01000,  0.65000,  -0.00110,   1.84000,  1.25}, //
	{"Cl",   EL_CL,  0.27600,  2.04500,  -0.00110,  35.82350,  0.99}, //
	{"Br",   EL_BR,  0.38900,  2.16500,  -0.00110,  42.56610,  1.14}, //
    {"UNK", EL_UNK,  0.15000,  2.00000,  -0.00143,  33.51030,  0.77}  //set as carbon
}};

const Float ad4_msp = -0.00110;//metal_solvation_parameter
const Float ad4_mcr = 1.75; //metal_covalent_radius for metals not on the list // FIXME this info should be moved to non_ad_metals

const Size_Type ad4_infos_size =  ad4_infos.size();

const std::map<std::string,AD4_TYPE> ad4_type_lookup = {
{ "C" , AD4_C  },
{ "A" , AD4_A  },
{ "N" , AD4_N  },
{ "O" , AD4_O  },
{ "P" , AD4_P  },
{ "S" , AD4_S  },
{"Se" , AD4_SE },
{ "H" , AD4_H  },
{ "F" , AD4_F  },
{ "I" , AD4_I  },
{"NA" , AD4_NA },
{"OA" , AD4_OA },
{"SA" , AD4_SA },
{"HD" , AD4_HD },
{"Mg" , AD4_Mg },
{"Mn" , AD4_Mn },
{"Zn" , AD4_Zn },
{"Ca" , AD4_Ca },
{"Fe" , AD4_Fe },
{"Cl" , AD4_Cl },
{"Br" , AD4_Br },
{"UNK", AD4_UNK}
};


// struct Acceptor_Info {
//     Float radius;
//     Float depth;
//     Acceptor_Info(const Float r, const Float d) : radius(r), depth(d) {}
// };

// const acceptor_kind acceptor_kind_data[] = { // ad_type, optimal length, depth
// 	{AD_TYPE_NA, 1.9, 5.0},
// 	{AD_TYPE_OA, 1.9, 5.0},
// 	{AD_TYPE_SA, 2.5, 1.0}
// };


// inline bool ad_is_hydrogen(Size_Type ad) {
// 	return ad == AD_TYPE_H || ad == AD_TYPE_HD;
// }

// inline bool ad_is_heteroatom(Size_Type ad) { // returns false for ad >= AD_TYPE_SIZE
// 	return ad != AD_TYPE_A && ad != AD_TYPE_C  &&
// 		   ad != AD_TYPE_H && ad != AD_TYPE_HD &&
// 		   ad < AD_TYPE_SIZE;
// }

// inline Size_Type ad_type_to_el_type(Size_Type t) {
// 	switch(t) {
// 		case AD_TYPE_C    : return EL_TYPE_C;
// 		case AD_TYPE_A    : return EL_TYPE_C;
// 		case AD_TYPE_N    : return EL_TYPE_N;
// 		case AD_TYPE_O    : return EL_TYPE_O;
// 		case AD_TYPE_P    : return EL_TYPE_P;
// 		case AD_TYPE_S    : return EL_TYPE_S;
// 		case AD_TYPE_H    : return EL_TYPE_H;
// 		case AD_TYPE_F    : return EL_TYPE_F;
// 		case AD_TYPE_I    : return EL_TYPE_I;
// 		case AD_TYPE_NA   : return EL_TYPE_N;
// 		case AD_TYPE_OA   : return EL_TYPE_O;
// 		case AD_TYPE_SA   : return EL_TYPE_S;
// 		case AD_TYPE_HD   : return EL_TYPE_H;
// 		case AD_TYPE_Mg   : return EL_TYPE_Met;
// 		case AD_TYPE_Mn   : return EL_TYPE_Met;
// 		case AD_TYPE_Zn   : return EL_TYPE_Met;
// 		case AD_TYPE_Ca   : return EL_TYPE_Met;
// 		case AD_TYPE_Fe   : return EL_TYPE_Met;
// 		case AD_TYPE_Cl   : return EL_TYPE_Cl;
// 		case AD_TYPE_Br   : return EL_TYPE_Br;
// 		case AD_TYPE_SIZE : return EL_TYPE_SIZE;
// 		default: TMD_CHECK(false);
// 	}
// 	return EL_TYPE_SIZE; // to placate the compiler in case of warnings - it should never get here though
// }


// const std::string non_ad_metal_names[] = { // expand as necessary
// 	"Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"
// };

// inline bool is_non_ad_metal_name(const std::string& name) {
// 	const Size_Type s = sizeof(non_ad_metal_names) / sizeof(const std::string);
// 	for(Size_Type i = 0; i < s; ++i){
// 		if(non_ad_metal_names[i] == name){
// 			return true;
// 		}
// 	}
// 	return false;
// }

// inline bool xs_is_hydrophobic(Size_Type xs) {
// 	return xs == XS_TYPE_C_H ||
// 		   xs == XS_TYPE_F_H ||
// 		   xs == XS_TYPE_Cl_H ||
// 		   xs == XS_TYPE_Br_H ||
// 		   xs == XS_TYPE_I_H;
// }

// inline bool xs_is_acceptor(Size_Type xs) {
// 	return xs == XS_TYPE_N_A ||
// 		   xs == XS_TYPE_N_DA ||
// 		   xs == XS_TYPE_O_A ||
// 		   xs == XS_TYPE_O_DA;
// }

// inline bool xs_is_donor(Size_Type xs) {
// 	return xs == XS_TYPE_N_D ||
// 		   xs == XS_TYPE_N_DA ||
// 		   xs == XS_TYPE_O_D ||
// 		   xs == XS_TYPE_O_DA ||
// 		   xs == XS_TYPE_Met_D;
// }

// inline bool xs_donor_acceptor(Size_Type t1, Size_Type t2) {
// 	return xs_is_donor(t1) && xs_is_acceptor(t2);
// }

// inline bool xs_h_bond_possible(Size_Type t1, Size_Type t2) {
// 	return xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1);
// }

// inline const atom_kind& ad_type_property(Size_Type i) {
// 	assert(AD_TYPE_SIZE == atom_kinds_size);
//     assert(i < atom_kinds_size);
//     return atom_kind_data[i];
// }

// inline Size_Type string_to_ad_type(const std::string& name) { // returns AD_TYPE_SIZE if not found (no exceptions thrown, because metals unknown to AD4 are not exceptional)
//     for(Size_Type i = 0; i < atom_kinds_size; ++i){
// 		if(atom_kind_data[i].name == name){
// 			return i;
// 		}
// 	}
// 	for(Size_Type i = 0; i < atom_equivalences_size; ++i){
// 		if(atom_equivalence_data[i].name == name){
// 			return string_to_ad_type(atom_equivalence_data[i].to);
// 		}
// 	}
//     return AD_TYPE_SIZE;
// }

// inline Float max_covalent_radius() {
// 	Float tmp = 0;
// 	for(Size_Type i = 0; i < atom_kinds_size; ++i){
// 		if(atom_kind_data[i].covalent_radius > tmp){
// 			tmp = atom_kind_data[i].covalent_radius;
// 		}
// 	}
// 	return tmp;
// }
























    // const Size_Type acceptor_infos_size =  acceptor_infos.size();

    // inline const ELEMENT_TYPE sybyl_type_to_element_type(const SYBYL_TYPE at) {
    //     return sybyl_infos.at(at).EL_TYPE;
    // }

    // inline const Float covalent_radius(const SYBYL_TYPE at) {
    //     return sybyl_infos.at(at).covalent_radius;
    // }
    // inline const Float atomic_radius(const SYBYL_TYPE at) {
    //     return sybyl_infos.at(at).atomic_radius;
    // }
    // inline const Float vdw_radius(const SYBYL_TYPE at) {
    //     return sybyl_infos.at(at).vdw_radius;
    // }

    // inline const std::string sybyl_type_to_string(const SYBYL_TYPE at) {
    //     return sybyl_infos.at(at).name;
    // }

    // inline const SYBYL_TYPE string_to_sybyl_type(const std::string s) {
    //     SYBYL_TYPE at = (sybyl_type_lookup.find(s)==sybyl_type_lookup.end()) ? TYPE_UNK : sybyl_type_lookup.at(s);
    //     return at;
    // }

}