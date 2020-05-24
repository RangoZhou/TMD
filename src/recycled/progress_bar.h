#pragma once

template<typename T1, typename T2>
void Progress_Bar(const T1 &counted, const T2 &total);

#include <iostream>
#include <iomanip>
namespace tmd {

template<typename T1, typename T2>
void progress_bar(const T1 &counted, const T2 &total) {
    Float progress = Float(counted)/Float(total);
    unsigned barWidth = 50;
    std::cout << "[";
    unsigned pos = barWidth * progress;
    for(unsigned i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(2) << (progress * 100.0) << "%\r";
    std::cout.flush();
    if(counted == total){
        std::cout << std::defaultfloat << std::setprecision(6) << std::endl;
    }
}

}//namespace