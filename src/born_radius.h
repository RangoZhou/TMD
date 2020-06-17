#pragma once

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
#include <set>
#include <array>

#include "atom.h"
#include "molecule.h"
#include "grid.h"
#include "matrix.h"
#include "score.h"

namespace lz {

inline void print_vector(const std::string s, const std::vector<double>& v) {
    std::cout << s << ": [min,max]: [" << *std::min_element(v.begin(), v.end()) << "," << *std::max_element(v.begin(), v.end()) << "]" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    for(const auto& a : v) {
        std::cout << a << " ";
    }
    std::cout << std::endl;
    std::vector<double> inverse(v.size());
    for(int i = 0; i < v.size(); ++i) {
        inverse[i] = 1.0 / v[i];
    }
    // std::sort(inverse.begin(), inverse.end());
    std::cout << "inverse: [min,max]: [" << *std::min_element(inverse.begin(), inverse.end()) << "," << *std::max_element(inverse.begin(), inverse.end()) << "]" << std::endl;
    for(const auto& a : inverse) {
        std::cout << a << " ";
    }
    std::cout << std::endl;
    std::cout << std::defaultfloat << std::setprecision(6) << std::endl;
}

// order  O P H C N
enum ATOM_TYPE {Type_O = 0, Type_P, Type_H, Type_C, Type_N, Type_Size};

const double radius_water = 1.4;
// order  O P H C N
// const double radius_S = 1.80;
//old rl vdw radii
// const std::array<double,Type_Size> vdw_radii = {{1.5, 1.9, 1.0, 1.7, 1.65}};
// mbondi3 radii set need to add detailed radius
// const double intrinsic_offset = 0.195141;
const std::array<double,Type_Size> vdw_radii = {{1.5, 1.85, 1.2, 1.7, 1.55}};

struct Born {
    // order  O P H C N
    // const double born_scale_S = 0.80;
    //old rl born scales
    // const std::array<double,Type_Size> scales = {{0.85, 0.86, 0.85, 0.72, 0.79}};
    //new scales
    const std::array<double,Type_Size> scales = {{0.184, 1.545, 1.697, 1.269, 1.426}};
    const std::array<double,Type_Size> alphas = {{0.606, 0.418, 0.537, 0.332, 0.686}};
    const std::array<double,Type_Size> betas = {{0.463, 0.290, 0.363, 0.197, 0.463}};
    const std::array<double,Type_Size> gammas = {{0.142, 0.106, 0.117, 0.093, 0.139}};
    const double offset = 0.195;
    const double neck_scale = 0.827;
    //row order: O P H C N
    //column order: O P H C N
    std::array<std::array<double,Type_Size>,Type_Size> d0_table;
    std::array<std::array<double,Type_Size>,Type_Size> m0_table;

    const double min_table_vdw_radius = 1.00;
    const double max_table_vdw_radius = 2.00;
    const double table_atom_radius_step = 0.05;
    const double table_size = 21;
    // GBNECKCUT: 2.8d0 (diameter of water) is "correct" value but
    // larger values give smaller discontinuities at the cut
    const double neck_cut = 6.8;
    // // d0
    // //These are the values used ford0in eqs 6 and 7. Distances andatom radii are in angstroms, Assuming a Solvent Molecule (Rw) Radiusof 1.4 Å
    // std::array<std::array<double,13>,13> d0_total_table = {{
    //     {{2.6797,  2.7250,  2.7719,  2.8188,  2.8656,  2.9125,  2.9609,  3.0078,  3.0562,  3.1047,  3.1531,  3.2016,  3.2500}},
    //     {{2.7359,  2.7813,  2.8281,  2.8750,  2.9219,  2.9688,  3.0156,  3.0641,  3.1109,  3.1594,  3.2078,  3.2563,  3.3047}},
    //     {{2.7922,  2.8375,  2.8844,  2.9297,  2.9766,  3.0234,  3.0719,  3.1188,  3.1672,  3.2141,  3.2625,  3.3109,  3.3594}},
    //     {{2.8500,  2.8953,  2.9406,  2.9859,  3.0328,  3.0797,  3.1266,  3.1750,  3.2219,  3.2703,  3.3172,  3.3656,  3.4141}},
    //     {{2.9062,  2.9516,  2.9969,  3.0422,  3.0891,  3.1359,  3.1828,  3.2297,  3.2766,  3.3250,  3.3719,  3.4203,  3.4688}},
    //     {{2.9625,  3.0078,  3.0531,  3.0984,  3.1437,  3.1906,  3.2375,  3.2844,  3.3313,  3.3797,  3.4266,  3.4750,  3.5234}},
    //     {{3.0188,  3.0641,  3.1078,  3.1547,  3.2000,  3.2469,  3.2922,  3.3391,  3.3875,  3.4344,  3.4813,  3.5297,  3.5781}},
    //     {{3.0750,  3.1203,  3.1641,  3.2094,  3.2563,  3.3016,  3.3484,  3.3953,  3.4422,  3.4891,  3.5359,  3.5844,  3.6313}},
    //     {{3.1313,  3.1750,  3.2203,  3.2656,  3.3109,  3.3563,  3.4031,  3.4500,  3.4969,  3.5438,  3.5906,  3.6391,  3.6859}},
    //     {{3.1875,  3.2313,  3.2766,  3.3203,  3.3656,  3.4125,  3.4578,  3.5047,  3.5516,  3.5984,  3.6453,  3.6922,  3.7406}},
    //     {{3.2437,  3.2875,  3.3313,  3.3766,  3.4219,  3.4672,  3.5125,  3.5594,  3.6063,  3.6531,  3.7000,  3.7469,  3.7953}},
    //     {{3.3000,  3.3422,  3.3875,  3.4312,  3.4766,  3.5219,  3.5688,  3.6141,  3.6609,  3.7078,  3.7547,  3.8016,  3.8484}},
    //     {{3.3547,  3.3984,  3.4422,  3.4875,  3.5313,  3.5766,  3.6234,  3.6688,  3.7156,  3.7625,  3.8094,  3.8563,  3.9031}}
    // }};
    // //m0
    // //These are the values used form0in eqs 6 and 7. Distances and atom radii are in angstroms. Assuming a Solvent Molecule (Rw) Radius of 1.4 Å
    // std::array<std::array<double,13>,13> m0_total_table = {{
    //     {{0.35281,  0.36412,  0.37516,  0.38594,  0.39645,  0.40670,  0.41670,  0.42646,  0.43598,  0.44527,  0.45434,  0.46319,  0.47183}},
    //     {{0.31853,  0.32889,  0.33902,  0.34890,  0.35855,  0.36797,  0.37717,  0.38615,  0.39492,  0.40348,  0.41185,  0.42001,  0.42799}},
    //     {{0.28847,  0.29798,  0.30728,  0.31637,  0.32525,  0.33392,  0.34240,  0.35069,  0.35878,  0.36669,  0.37441,  0.38196,  0.38934}},
    //     {{0.26199,  0.27074,  0.27930,  0.28768,  0.29587,  0.30387,  0.31170,  0.31936,  0.32684,  0.33416,  0.34131,  0.34830,  0.35514}},
    //     {{0.23859,  0.24666,  0.25455,  0.26228,  0.26985,  0.27725,  0.28449,  0.29158,  0.29851,  0.30529,  0.31193,  0.31842,  0.32477}},
    //     {{0.21783,  0.22528,  0.23258,  0.23972,  0.24673,  0.25358,  0.26029,  0.26686,  0.27330,  0.27959,  0.28575,  0.29179,  0.29769}},
    //     {{0.19935,  0.20624,  0.21300,  0.21962,  0.22611,  0.23247,  0.23870,  0.24480,  0.25078,  0.25664,  0.26237,  0.26799,  0.27349}},
    //     {{0.18285,  0.18923,  0.19550,  0.20165,  0.20767,  0.21358,  0.21938,  0.22505,  0.23062,  0.23607,  0.24141,  0.24665,  0.25178}},
    //     {{0.16807,  0.17400,  0.17982,  0.18553,  0.19114,  0.19664,  0.20203,  0.20732,  0.21251,  0.21759,  0.22258,  0.22747,  0.23226}},
    //     {{0.15480,  0.16031,  0.16573,  0.17104,  0.17626,  0.18139,  0.18642,  0.19135,  0.19620,  0.20095,  0.20561,  0.21018,  0.21466}},
    //     {{0.14285,  0.14798,  0.15303,  0.15798,  0.16285,  0.16764,  0.17233,  0.17694,  0.18147,  0.18591,  0.19027,  0.19455,  0.19875}},
    //     {{0.13207,  0.13685,  0.14155,  0.14618,  0.15073,  0.15520,  0.15959,  0.16390,  0.16814,  0.17230,  0.17638,  0.18039,  0.18433}},
    //     {{0.12231,  0.12677,  0.13117,  0.13549,  0.13975,  0.14393,  0.14804,  0.15208,  0.15605,  0.15995,  0.16378,  0.16754,  0.17124}}
    // }};

    // Lookup tables for position (atom separation, r) and value of the maximum
    // of the neck function for given atomic radii ri and rj. Values of neck
    // maximum are already divided by 4 Pi to save time. Values are given
    // for each 0.05 angstrom  between 1.0 and 2.0 (inclusive), so map to index with
    // nint((r-1.0)*20)).  Values were numerically determined in Mathematica;
    // note FORTRAN column-major array storage, so the data below may be
    // transposed from how you might expect it.

    std::array<std::array<double,21>,21> d0_total_table = {{
        {{2.26685,2.32548,2.38397,2.44235,2.50057,2.55867,2.61663,2.67444,2.73212,2.78965,2.84705,2.9043,2.96141,3.0184,3.07524,3.13196,3.18854,3.24498,3.30132,3.35752,3.4136}},
        {{2.31191,2.37017,2.4283,2.48632,2.5442,2.60197,2.65961,2.71711,2.77449,2.83175,2.88887,2.94586,3.00273,3.05948,3.1161,3.1726,3.22897,3.28522,3.34136,3.39738,3.45072}},
        {{2.35759,2.41549,2.47329,2.53097,2.58854,2.646,2.70333,2.76056,2.81766,2.87465,2.93152,2.98827,3.0449,3.10142,3.15782,3.21411,3.27028,3.32634,3.3823,3.43813,3.49387}},
        {{2.4038,2.46138,2.51885,2.57623,2.63351,2.69067,2.74773,2.80469,2.86152,2.91826,2.97489,3.0314,3.08781,3.1441,3.20031,3.25638,3.31237,3.36825,3.42402,3.4797,3.53527}},
        {{2.45045,2.50773,2.56492,2.62201,2.679,2.7359,2.7927,2.8494,2.90599,2.9625,3.0189,3.07518,3.13138,3.18748,3.24347,3.29937,3.35515,3.41085,3.46646,3.52196,3.57738}},
        {{2.4975,2.5545,2.61143,2.66825,2.72499,2.78163,2.83818,2.89464,2.95101,3.00729,3.06346,3.11954,3.17554,3.23143,3.28723,3.34294,3.39856,3.45409,3.50952,3.56488,3.62014}},
        {{2.54489,2.60164,2.6583,2.71488,2.77134,2.8278,2.88412,2.94034,2.9965,3.05256,3.10853,3.16442,3.22021,3.27592,3.33154,3.38707,3.44253,3.49789,3.55316,3.60836,3.66348}},
        {{2.59259,2.6491,2.70553,2.76188,2.81815,2.87434,2.93044,2.98646,3.04241,3.09827,3.15404,3.20974,3.26536,3.32089,3.37633,3.4317,3.48699,3.54219,3.59731,3.65237,3.70734}},
        {{2.64054,2.69684,2.75305,2.80918,2.86523,2.92122,2.97712,3.03295,3.0887,3.14437,3.19996,3.25548,3.31091,3.36627,3.42156,3.47677,3.5319,3.58695,3.64193,3.69684,3.75167}},
        {{2.68873,2.74482,2.80083,2.85676,2.91262,2.96841,3.02412,3.07976,3.13533,3.19082,3.24623,3.30157,3.35685,3.41205,3.46718,3.52223,3.57721,3.63213,3.68696,3.74174,3.79644}},
        {{2.73713,2.79302,2.84884,2.90459,2.96027,3.01587,3.0714,3.12686,3.18225,3.23757,3.29282,3.34801,3.40313,3.45815,3.51315,3.56805,3.6229,3.67767,3.73237,3.78701,3.84159}},
        {{2.78572,2.84143,2.89707,2.95264,3.00813,3.06356,3.11892,3.17422,3.22946,3.28462,3.33971,3.39474,3.44971,3.5046,3.55944,3.61421,3.66891,3.72356,3.77814,3.83264,3.8871}},
        {{2.83446,2.89,2.94547,3.00088,3.05621,3.11147,3.16669,3.22183,3.27689,3.33191,3.38685,3.44174,3.49656,3.55132,3.60602,3.66066,3.71523,3.76975,3.82421,3.8786,3.93293}},
        {{2.88335,2.93873,2.99404,3.04929,3.10447,3.15959,3.21464,3.26963,3.32456,3.37943,3.43424,3.48898,3.54366,3.5983,3.65287,3.70737,3.76183,3.81622,3.87056,3.92484,3.97905}},
        {{2.93234,2.9876,3.04277,3.09786,3.15291,3.20787,3.26278,3.31764,3.37242,3.42716,3.48184,3.53662,3.591,3.64551,3.69995,3.75435,3.80867,3.86295,3.91718,3.97134,4.02545}},
        {{2.98151,3.0366,3.09163,3.14659,3.20149,3.25632,3.3111,3.36581,3.42047,3.47507,3.52963,3.58411,3.63855,3.69293,3.74725,3.80153,3.85575,3.90991,3.96403,4.01809,4.07211}},
        {{3.03074,3.08571,3.14061,3.19543,3.25021,3.30491,3.35956,3.41415,3.46869,3.52317,3.57759,3.63196,3.68628,3.74054,3.79476,3.84893,3.90303,3.95709,4.01111,4.06506,4.11897}},
        {{3.08008,3.13492,3.1897,3.2444,3.29905,3.35363,3.40815,3.46263,3.51704,3.57141,3.62572,3.67998,3.73418,3.78834,3.84244,3.8965,3.95051,4.00447,4.05837,4.11224,4.16605}},
        {{3.12949,3.18422,3.23888,3.29347,3.348,3.40247,3.45688,3.51124,3.56554,3.6198,3.674,3.72815,3.78225,3.83629,3.8903,3.94425,3.99816,4.05203,4.10583,4.15961,4.21333}},
        {{3.17899,3.23361,3.28815,3.34264,3.39706,3.45142,3.50571,3.55997,3.61416,3.66831,3.72241,3.77645,3.83046,3.8844,3.93831,3.99216,4.04598,4.09974,4.15347,4.20715,4.26078}},
        {{3.22855,3.28307,3.33751,3.39188,3.4462,3.50046,3.55466,3.6088,3.6629,3.71694,3.77095,3.82489,3.8788,3.93265,3.98646,4.04022,4.09395,4.14762,4.20126,4.25485,4.3084}}

    }};

    std::array<std::array<double,21>,21> m0_total_table = {{
        {{0.0381511,0.0338587,0.0301776,0.027003,0.0242506,0.0218529,0.0197547,0.0179109,0.0162844,0.0148442,0.0135647,0.0124243,0.0114047,0.0104906,0.00966876,0.008928,0.0082587,0.00765255,0.00710237,0.00660196,0.00614589}},
        {{0.0396198,0.0351837,0.0313767,0.0280911,0.0252409,0.0227563,0.0205808,0.0186681,0.0169799,0.0154843,0.014155,0.0129696,0.0119094,0.0109584,0.0101031,0.00933189,0.0086348,0.00800326,0.00742986,0.00690814,0.00643255}},
        {{0.041048,0.0364738,0.0325456,0.0291532,0.0262084,0.0236399,0.0213897,0.0194102,0.0176622,0.0161129,0.0147351,0.0135059,0.0124061,0.0114192,0.0105312,0.00973027,0.00900602,0.00834965,0.0077535,0.00721091,0.00671609}},
        {{0.0424365,0.0377295,0.0336846,0.0301893,0.0271533,0.0245038,0.0221813,0.0201371,0.018331,0.0167295,0.0153047,0.014033,0.0128946,0.0118727,0.0109529,0.0101229,0.00937212,0.00869147,0.00807306,0.00751003,0.00699641}},
        {{0.0437861,0.0389516,0.0347944,0.0311998,0.0280758,0.0253479,0.0229555,0.0208487,0.0189864,0.0173343,0.0158637,0.0145507,0.0133748,0.0123188,0.0113679,0.0105096,0.0097329,0.00902853,0.00838835,0.00780533,0.0072733}},
        {{0.0450979,0.0401406,0.0358753,0.0321851,0.0289761,0.0261726,0.0237125,0.0215451,0.0196282,0.017927,0.0164121,0.0150588,0.0138465,0.0127573,0.0117761,0.0108902,0.0100882,0.00936068,0.00869923,0.00809665,0.00754661}},
        {{0.0463729,0.0412976,0.0369281,0.0331456,0.0298547,0.026978,0.0244525,0.0222264,0.0202567,0.0185078,0.0169498,0.0155575,0.0143096,0.0131881,0.0121775,0.0112646,0.010438,0.00968781,0.00900559,0.00838388,0.00781622}},
        {{0.0476123,0.0424233,0.0379534,0.034082,0.0307118,0.0277645,0.0251757,0.0228927,0.0208718,0.0190767,0.0174768,0.0160466,0.0147642,0.0136112,0.0125719,0.0116328,0.0107821,0.0100099,0.00930735,0.00866695,0.00808206}},
        {{0.0488171,0.0435186,0.038952,0.0349947,0.0315481,0.0285324,0.0258824,0.0235443,0.0214738,0.0196339,0.0179934,0.0165262,0.0152103,0.0140267,0.0129595,0.0119947,0.0111206,0.0103268,0.00960445,0.00894579,0.00834405}},
        {{0.0499883,0.0445845,0.0399246,0.0358844,0.032364,0.0292822,0.0265729,0.0241815,0.0220629,0.0201794,0.0184994,0.0169964,0.0156479,0.0144345,0.0133401,0.0123504,0.0114534,0.0106386,0.00989687,0.00922037,0.00860216}},
        {{0.0511272,0.0456219,0.040872,0.0367518,0.0331599,0.0300142,0.0272475,0.0248045,0.0226392,0.0207135,0.0189952,0.0174574,0.0160771,0.0148348,0.0137138,0.0126998,0.0117805,0.0109452,0.0101846,0.00949067,0.00885636}},
        {{0.0522348,0.0466315,0.0417948,0.0375973,0.0339365,0.030729,0.0279067,0.0254136,0.023203,0.0212363,0.0194809,0.0179092,0.016498,0.0152275,0.0140807,0.013043,0.012102,0.0112466,0.0104676,0.00975668,0.00910664}},
        {{0.0533123,0.0476145,0.042694,0.0384218,0.0346942,0.0314268,0.0285507,0.026009,0.0237547,0.0217482,0.0199566,0.018352,0.0169108,0.0156128,0.0144408,0.0133801,0.0124179,0.011543,0.010746,0.0100184,0.00935302}},
        {{0.0543606,0.0485716,0.04357,0.0392257,0.0354335,0.0321082,0.02918,0.0265913,0.0242943,0.0222492,0.0204225,0.0187859,0.0173155,0.0159908,0.0147943,0.0137111,0.0127282,0.0118343,0.0110197,0.0102759,0.00959549}},
        {{0.0553807,0.0495037,0.0444239,0.0400097,0.0361551,0.0327736,0.0297949,0.0271605,0.0248222,0.0227396,0.0208788,0.0192111,0.0177122,0.0163615,0.0151413,0.0140361,0.013033,0.0121206,0.0112888,0.0105292,0.00983409}},
        {{0.0563738,0.0504116,0.0452562,0.0407745,0.0368593,0.0334235,0.0303958,0.0277171,0.0253387,0.0232197,0.0213257,0.0196277,0.0181013,0.0167252,0.0154817,0.0143552,0.0133325,0.0124019,0.0115534,0.0107783,0.0100688}},
        {{0.0573406,0.0512963,0.0460676,0.0415206,0.0375468,0.0340583,0.030983,0.0282614,0.0258441,0.0236896,0.0217634,0.020036,0.0184826,0.017082,0.0158158,0.0146685,0.0136266,0.0126783,0.0118135,0.0110232,0.0102998}},
        {{0.0582822,0.0521584,0.0468589,0.0422486,0.038218,0.0346784,0.0315571,0.0287938,0.0263386,0.0241497,0.0221922,0.0204362,0.0188566,0.0174319,0.0161437,0.0149761,0.0139154,0.0129499,0.0120691,0.0112641,0.0105269}},
        {{0.0591994,0.0529987,0.0476307,0.042959,0.0388734,0.0352843,0.0321182,0.0293144,0.0268225,0.0246002,0.0226121,0.0208283,0.0192232,0.0177751,0.0164654,0.015278,0.0141991,0.0132167,0.0123204,0.0115009,0.0107504}},
        {{0.0600932,0.053818,0.0483836,0.0436525,0.0395136,0.0358764,0.0326669,0.0298237,0.0272961,0.0250413,0.0230236,0.0212126,0.0195826,0.0181118,0.0167811,0.0155744,0.0144778,0.0134789,0.0125673,0.0117338,0.0109702}},
        {{0.0609642,0.0546169,0.0491183,0.0443295,0.0401388,0.036455,0.0332033,0.030322,0.0277596,0.0254732,0.0234266,0.0215892,0.0199351,0.018442,0.0170909,0.0158654,0.0147514,0.0137365,0.0128101,0.0119627,0.0111863}}
    }};

    const double& get_d0(const int& r_1, const int& r_2) const {
        //table is column major
        return d0_total_table[r_2][r_1];
    }
    const double& get_m0(const int& r_1, const int& r_2) const {
        //colum is atom 1, row is atom 2
        return m0_total_table[r_2][r_1];
    }
    Born() {
        for(int i = 0; i < Type_Size; ++i) {
            for(int j = 0; j < Type_Size; ++j) {
                const double& r_i = vdw_radii[i];
                const double& r_j = vdw_radii[j];
                int ri_index = static_cast<int>( (r_i - min_table_vdw_radius) / table_atom_radius_step );
                int rj_index = static_cast<int>( (r_j - min_table_vdw_radius) / table_atom_radius_step );
                //need to check after finish
                // if(ri_index < 0) ri_index = 0;
                // if(rj_index < 0) rj_index = 0;
                // if(ri_index >= 13) ri_index = 12;
                // if(rj_index >= 13) rj_index = 12;
                assert(ri_index >= 0 && ri_index < table_size);
                assert(rj_index >= 0 && rj_index < table_size);
                m0_table[i][j] = get_m0(ri_index, rj_index);
                d0_table[i][j] = get_d0(ri_index, rj_index);
            }
        }
    }
};

struct atom_info {
    int index;
    ATOM_TYPE type_idx;
    int yz_atom_index;
    std::string name;
    double x;
    double y;
    double z;
    std::string sybyl_type;
    std::string resname;
    // int reseq;
    double charge;
    std::string element;
    int numH = 0;
};

class case_info {
public:
    std::vector< atom_info > noH;
    std::vector< atom_info > wH;
    std::vector< atom_info > HH;
    int heavy_size;
    int wH_size;
    int HH_size;
    void SetValue (const tmd::Atoms& atoms) {
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

            if(a.sybyl_type.find(".")!=std::string::npos) {
                int position;
                position = a.sybyl_type.find(".");
                a.element = a.sybyl_type.substr(0,position);
            } else {
                a.element = a.sybyl_type;
            }

            if(a.element == "H" || a.element == "h") {
                a.type_idx = Type_H;
            } else if(a.element == "C" || a.element == "c") {
                a.type_idx = Type_C;
            } else if(a.element == "N" || a.element == "n") {
                a.type_idx = Type_N;
            } else if(a.element == "O" || a.element == "o") {
                a.type_idx = Type_O;
            } else if(a.element == "P" || a.element == "p") {
                a.type_idx = Type_P;
            } else if(a.element == "S" || a.element == "s") {
                a.type_idx = Type_C;
            } else {
                a.type_idx = Type_C;
            }
            wH.push_back(a);

            if(a.element != "H" && a.element != "h") {
                a.index = index;
                index++;
                noH.push_back(a);
            } else {
                HH.push_back(a);
            }
        }
    }

    void H_charge_transfer() {
        for(auto& a : HH) {
            int idx;
            double dmin = 10000000;
            for(auto& b : noH) {
                double dis = (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z);
                if(dis < dmin) {
                    dmin = dis;
                    idx = b.index;
                }
            }
            noH[idx].charge += a.charge;
            noH[idx].numH++;
        }
    };

    void add_case(const case_info& ci) {
        const int old_noH_size = this->noH.size();
        this->heavy_size = old_noH_size+ci.noH.size();
        this->noH.resize(this->heavy_size);
        for(int i = old_noH_size; i < old_noH_size+ci.noH.size(); ++i) {
            this->noH[i] = ci.noH[i-old_noH_size];
        }

        const int old_wH_size = this->wH.size();
        this->wH_size = old_wH_size+ci.wH.size();
        this->wH.resize(this->wH_size);
        for(int i = old_wH_size; i < old_wH_size+ci.wH.size(); ++i) {
            this->wH[i] = ci.wH[i-old_wH_size];
        }

        const int old_HH_size = this->HH.size();
        this->HH_size = old_HH_size+ci.HH.size();
        this->HH.resize(this->HH_size);
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

    void Find_Max_and_Min() {
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

struct I_MS {
    double vdw = 0;
    double neck = 0;
};

struct Fixed {
    // tmd::Strictly_Triangular_Matrix<tmd::DISTANCE_TYPE> dis_type_mat;
    // tmd::Strictly_Triangular_Matrix<tmd::BOND_TYPE> bond_type_mat;

    std::vector<double> rna_born_radius_table;
    std::vector<I_MS> rna_fixed_born_part_table;
    void print_rna_born_radius_table() {
        print_vector("rna_born_radius_table",this->rna_born_radius_table);
    }
};


struct Variable {
    std::vector<double> lig_born_radius_table;
    std::vector<double> complex_born_radius_table;
    tmd::Strictly_Triangular_Matrix<double> complex_dis_mat;
    tmd::Strictly_Triangular_Matrix<double> complex_dis_square_mat;
    void print_lig_born_radius_table() {
        print_vector("lig_born_radius_table",this->lig_born_radius_table);
    }
    void print_complex_born_radius_table() {
        print_vector("complex_born_radius_table",this->complex_born_radius_table);
    }
};



inline const I_MS born_I_vdw_and_neck(const int& i, const int& start_index, const int& end_index, const Born& born, const std::vector<atom_info>& atoms, const tmd::Strictly_Triangular_Matrix<double>& dis_mat) {
    I_MS I_ms;
    double Lij, Uij;
    I_ms.vdw = 0;
    I_ms.neck = 0;
    double radi = vdw_radii[atoms[i].type_idx]-born.offset;
    for (int j = start_index; j < end_index; j++) {
        if (j!= i) {
            //first cal I_vdw
            double dis;
            if(i < j) {
                dis = dis_mat(i,j);
            } else {
                dis = dis_mat(j,i);
            }

            double radj = vdw_radii[atoms[j].type_idx]-born.offset;
            double radc = dis + born.scales[atoms[j].type_idx] * radj;
            if (radi >= radc) {
                Lij = 1.;
                Uij = 1.;
            } else {
                if (radi >= (dis - born.scales[atoms[j].type_idx] * radj)) {
                    Lij = radi;
                } else {
                    Lij = dis - born.scales[atoms[j].type_idx] * radj;
                }
                Uij = dis + born.scales[atoms[j].type_idx] * radj;
            }
            I_ms.vdw = I_ms.vdw + 0.5 * ((1. / Lij - 1. / Uij) + (born.scales[atoms[j].type_idx] * born.scales[atoms[j].type_idx] * radj * radj / (4. * dis) - dis / 4.) * (1. / (Lij * Lij) - 1. / (Uij * Uij)) + (1. / (2. * dis)) * log(Lij / Uij));

            //second cal I_neck
            if(dis <= born.neck_cut) {
                const double d_d0 = dis - born.d0_table[atoms[i].type_idx][atoms[j].type_idx];
                const double d_d0_square = d_d0*d_d0;
                I_ms.neck += born.m0_table[atoms[i].type_idx][atoms[j].type_idx] / (1 + d_d0_square + 0.3*d_d0_square*d_d0_square*d_d0_square );
            }
            // std::cout << j << " " << radi << " " << radj << " " << born.scales[atoms[i].type_idx] << " " << born.scales[atoms[j].type_idx] << " " << dis << " " << I_ms.vdw << " " << I_ms.neck << std::endl;
        }
        // std::cout << j << " " << I_ms.vdw << std::endl;
    }
    return I_ms;
}


//return born radius and I_MS
inline const double born_radius(const ATOM_TYPE& at, const Born& born, const I_MS& I_ms) {
    const double rho = vdw_radii[at];
    const double rho_tilde = rho - born.offset;
    const double omega = rho_tilde * (I_ms.vdw + born.neck_scale * I_ms.neck);
    //old
    // return 1. / (1. / rho - I_ms.vdw);
    //new
    const double tanh_arg = (born.alphas[at] - born.betas[at]*omega + born.gammas[at]*omega*omega)*omega;
    const double b_r = 1.0/( (1.0/rho_tilde) - (1.0/rho) * std::tanh(tanh_arg) );

    // std::cout << rho << " " << rho_tilde << " " << omega << " " << tanh_arg << " " << std::tanh(tanh_arg) << " " << b_r << std::endl;
    // if(b_r > 9) {
    //     exit(2);
    // }
    return b_r;
}


class RL_Score : public tmd::Scoring_Function {
    const tmd::RNA& tmd_rna;
    const tmd::Ligand& tmd_lig;
    std::ostream& tee;

    case_info rna;
    case_info lig;
    case_info complex;

    Fixed fixed;
    Variable var;
    // LJ lj;
    // HBond hbond;
    // SASA sasa;
    Born born;
public:
    RL_Score(const tmd::RNA& r, const tmd::Ligand& l, std::ostream& lg) : tmd_rna(r), tmd_lig(l), tee(lg) {
        tee << "RL_Score construct" << std::endl;
    }
    void init();
    const tmd::Float evaluate();
};


inline void RL_Score::init() {
    tee << "RL_Score init" << std::endl;
    lig.SetValue(this->tmd_lig.get_atoms_reference());
    if(lig.HH.size()>0) lig.H_charge_transfer();

    rna.SetValue(this->tmd_rna.get_atoms_reference());
    rna.Find_Max_and_Min();
    if(rna.HH.size()>0) rna.H_charge_transfer();

    complex = rna;
    complex.add_case(lig);

    rna.heavy_size = rna.noH.size();
    lig.heavy_size = lig.noH.size();
    complex.heavy_size = complex.noH.size();
    rna.wH_size = rna.wH.size();
    rna.HH_size = rna.HH.size();
    lig.wH_size = lig.wH.size();
    lig.HH_size = lig.HH.size();
    complex.wH_size = complex.wH.size();
    complex.HH_size = complex.HH.size();


    //init complex_dis_mat to 0.0
    var.complex_dis_mat.resize(complex.wH_size,0.0);
    var.complex_dis_square_mat.resize(complex.wH_size,0.0);
    //cal dis for complex for init
    for(int i = 0; i < complex.wH_size; ++i) {
        for(int j = 0; j < complex.wH_size; ++j) {
            if(i < j) {
                const double dis_square =   (complex.wH[i].x - complex.wH[j].x) * (complex.wH[i].x - complex.wH[j].x) +
                                            (complex.wH[i].y - complex.wH[j].y) * (complex.wH[i].y - complex.wH[j].y) +
                                            (complex.wH[i].z - complex.wH[j].z) * (complex.wH[i].z - complex.wH[j].z);
                const double dis = std::sqrt(dis_square);
                var.complex_dis_mat(i,j) = dis;
                var.complex_dis_square_mat(i,j) = dis_square;
            }
        }
    }

    //init born radius for rna
    //////////////////////////////////////////////
    /// RNA //////////////////////////////////////
    fixed.rna_born_radius_table = std::vector<double>(rna.wH_size,0.0);
    fixed.rna_fixed_born_part_table = std::vector<I_MS>(rna.wH_size);
    for (int i = 0; i < rna.wH_size; i++) { // start to calculate rb0
        const I_MS I_ms = born_I_vdw_and_neck(i, 0, rna.wH_size, born, complex.wH, var.complex_dis_mat);
        const double born_r = born_radius(complex.wH[i].type_idx, born, I_ms);
        fixed.rna_fixed_born_part_table[i] = I_ms;
        fixed.rna_born_radius_table[i] = born_r;
        // std::cout << i << " " << I_ms.vdw << " " << born_r << std::endl;
    }
}

inline const tmd::Float RL_Score::evaluate() {
    // clock_t start,mid1,mid2,end;
    // start=clock();
    // mid1=clock();
    // end = clock();

    //update lig's conformation and complex
    int index_noH=0,index_wH=0,index_HH=0;
    const tmd::Atoms& yz_lig_atoms_ref = this->tmd_lig.get_atoms_reference();
    // std::cout << yz_lig_atoms_ref.size() << " " << P.lig.wH.size() << std::endl;
    assert(yz_lig_atoms_ref.size() == lig.wH.size());
    for(const tmd::Atom& atom : this->tmd_lig.get_atoms_reference()) {
        const tmd::Vec3d& xyz = atom.get_coord();
        if(lig.wH[index_wH].element != "H" && lig.wH[index_wH].element != "h") {
            lig.noH[index_noH].x = xyz[0];
            lig.noH[index_noH].y = xyz[1];
            lig.noH[index_noH].z = xyz[2];
            complex.noH[index_noH+rna.wH_size].x = xyz[0];
            complex.noH[index_noH+rna.wH_size].y = xyz[1];
            complex.noH[index_noH+rna.wH_size].z = xyz[2];
            index_noH++;
        } else {
            lig.HH[index_HH].x = xyz[0];
            lig.HH[index_HH].y = xyz[1];
            lig.HH[index_HH].z = xyz[2];
            complex.HH[index_HH+rna.HH_size].x = xyz[0];
            complex.HH[index_HH+rna.HH_size].y = xyz[1];
            complex.HH[index_HH+rna.HH_size].z = xyz[2];
            index_HH++;
        }
        lig.wH[index_wH].x = xyz[0];
        lig.wH[index_wH].y = xyz[1];
        lig.wH[index_wH].z = xyz[2];
        complex.wH[index_wH+rna.wH_size].x = xyz[0];
        complex.wH[index_wH+rna.wH_size].y = xyz[1];
        complex.wH[index_wH+rna.wH_size].z = xyz[2];
        index_wH++;
    }

    // std::cout << "complete dis mat and update dis mat" << std::endl;

    //complete dis mat and update dis mat
    for (int i = 0; i < complex.wH_size; i++) {
        for (int j = rna.wH_size; j < complex.wH_size; j++) {
            if (i < j) {
                double dis_square =
                    (complex.wH[j].x - complex.wH[i].x) * (complex.wH[j].x - complex.wH[i].x) +
                    (complex.wH[j].y - complex.wH[i].y) * (complex.wH[j].y - complex.wH[i].y) +
                    (complex.wH[j].z - complex.wH[i].z) * (complex.wH[j].z - complex.wH[i].z);

                if(tmd::eq(dis_square,0)) {
                    dis_square = 0.0000000000000001;
                }
                var.complex_dis_square_mat(i,j) = dis_square;
                var.complex_dis_mat(i,j) = sqrt(dis_square);
            }
        }
    }

    ///////////////////////////////////////////////
    /// ligand //////////////////////////////////////
    var.lig_born_radius_table = std::vector<double>(lig.wH_size,0.0);
    for (int i = 0; i < lig.wH_size; i++) { // start to calculate rb0
        const I_MS I_ms = born_I_vdw_and_neck(i+rna.wH_size, rna.wH_size, complex.wH_size, born, complex.wH, var.complex_dis_mat);
        const double born_r = born_radius(complex.wH[i+rna.wH_size].type_idx, born, I_ms);
        var.lig_born_radius_table[i] = born_r;
    }

    //update complex born radius every time
    // std::cout << "start Born_Radius_Complex" << std::endl;
    var.complex_born_radius_table  = std::vector<double>(complex.wH_size,0.0);
    for (int i = 0; i < complex.wH_size; i++) { // start to calculate rb0
        if(i < rna.wH_size) {
            I_MS I_ms = born_I_vdw_and_neck(i, rna.wH_size, complex.wH_size, born, complex.wH, var.complex_dis_mat);
            I_ms.neck += fixed.rna_fixed_born_part_table[i].neck;
            I_ms.vdw += fixed.rna_fixed_born_part_table[i].vdw;
            var.complex_born_radius_table[i] = born_radius(complex.wH[i].type_idx, born, I_ms);
        } else {
            const I_MS I_ms = born_I_vdw_and_neck(i, 0, complex.wH_size, born, complex.wH, var.complex_dis_mat);
            var.complex_born_radius_table[i] = born_radius(complex.wH[i].type_idx, born, I_ms);
        }
    }


    std::cout << "rna-------------------------------------------" << std::endl;
    fixed.print_rna_born_radius_table();
    std::cout << "lig-------------------------------------------" << std::endl;
    var.print_lig_born_radius_table();
    std::cout << "complex---------------------------------------" << std::endl;
    var.print_complex_born_radius_table();
    std::cout << "----------------------------------------------" << std::endl;
    exit(2);
    // return energy/static_cast<double>(lig.heavy_size);
    return 0;

}

}
