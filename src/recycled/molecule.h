#pragma once

#include "block.h"
#include "atom.h"

// #include <vector>
#include <string>
#include <map>
#include <set>
// #include <memory>
// #include <iostream>
// #include <algorithm>
// #include <cstring>//memcmp

namespace mol{
    typedef Block<Ref_Atom> Selection;

    //Molecule class
    class Molecule : public Block<Atom> {
        friend void ReadPDB(std::istream& in, Molecule& m);
        friend void ReadMol2(std::istream& in, Molecule& m);
    public:
        Molecule(){}
        // Molecule(const std::string id) : id_(id) {}
        Molecule(const std::string id, const std::string type) : id_(id), type_(type) {}
        Molecule(const std::string id, const std::string type, std::initializer_list<Atom> il) : id_(id), type_(type), Block<Atom>(il) {}
        Molecule(const Molecule& m) : Block<Atom>(m), id_(m.id_), type_(m.type_) {}
        Molecule(const Selection& sel) {
            for(const auto& s : sel)
            this->push_back(s);
        };

        void set_id(const std::string id) {id_ = id;}
        void set_type(const std::string type) {type_ = type;}
        std::string get_id() const {return id_;};
        std::string get_type() const {return type_;};

        // unsigned int get_model_num() const {
        //     std::set<unsigned int> model_set;
        //     for(const auto& a : *this){
        //         model_set.insert(a.model_serial);
        //     }
        //     return model_set.size();
        // };
        unsigned int get_atom_num() const {return this->size();}

        Selection select(const size_type begin, const size_type end) {
            Selection sel;
            for(size_type i = begin; i < end; ++i) {
                //use [] operator will return a reference
                //reference will be treat as value when it is inserted into the Selection<Atom>?
                sel.push_back((*this)[i]);
            }
            return sel;
        }


    private:
        std::string id_ = "";
        std::string type_ = "";

        // void preProcess();
    };

    //the selection part is here
    inline Selection AtomSelect(Molecule& m, const Molecule::size_type begin, const Molecule::size_type end) {
        return m.select(begin,end);
    }

    //     //SelectionParser class
    // class SelectionParser{
    //     friend Selection AtomSelect(Molecule& m, const std::string str);
    //     public:
    //         SelectionParser(const std::string str);
    //     private:
    //         bool parse(const Atom& atom);
    //         bool parse_flag_ = false;
    //         std::map<std::string,std::string> sel_map_{{"name","none"},{"res_serial","none"},{"res_name","none"},{"chain_name","none"},{"serial","none"},{"model_serial","none"},{"alt_loc","none"},{"element","none"}};
    // };

}//namespace bio