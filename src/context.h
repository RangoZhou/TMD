#pragma once

#include <string>
#include <vector>
#include <iomanip>

#include "atom.h"

namespace tmd {

//Context stores the atom line and its position in the file content while reading
struct Context {
    int atom_index = -1;
    std::string content = "";
    Context() {}
    Context(const int& ai, const std::string& ct) : atom_index(ai), content(ct) {
        assert(ai >= -1);
    }
};

using Contexts = std::vector<Context>;

inline void write_contexts(const Contexts& cs, const Atoms& as, std::ostream& out = std::cout) {
    for(const Context& c : cs) {
        std::string tmp_content = c.content;
        std::ostringstream oss;
        oss.fill(' ');
        if(c.atom_index != -1) {
            const Atom& a = as[c.atom_index];
            oss << std::setw(8) << std::fixed << std::setprecision(3) << a.get_coord()[0] << std::setw(8) << a.get_coord()[1] << std::setw(8) << a.get_coord()[2];
            tmp_content.replace(30, 24, oss.str());//24 is the string size, default fill with space
        }
        out << tmp_content << std::endl;
    }
}


}