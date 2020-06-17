#include <iostream>
#include <boost/program_options.hpp>
#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?
#include "TMDConfig.h"
#include "common.h"
#include "molecule.h"
#include "mol2_reader.h"
#include "tmd_reader.h"
#include "print.h"
#include "score.h"
#include "yw_score.h"
// #include "rl_score.h"
// #include "rl_score_new_born.h"
#include "born_radius.h"

#include "random.h"
#include "sampling.h"
// #include "file.h"
// #include "tee.h"

#ifdef PROFILE
#include <gperftools/profiler.h>
#endif

namespace po = boost::program_options;

struct Usage_Error : public std::runtime_error {
	Usage_Error(const std::string& message) : std::runtime_error(message) {}
};

struct Options_Occurrence {
	bool some;
	bool all;
	Options_Occurrence() : some(false), all(true) {} // convenience
	Options_Occurrence& operator+=(const Options_Occurrence& x) {
		some = some || x.some;
		all  = all  && x.all;
		return *this;
	}
};

Options_Occurrence get_occurrence(po::variables_map& vm,po::options_description& d) {
	Options_Occurrence tmp;
	for(int i = 0; i < d.options().size(); ++i) {
		if(vm.count((*d.options()[i]).long_name()))
			tmp.some = true;
		else
			tmp.all = false;
    }
	return tmp;
}

void check_occurrence(po::variables_map& vm, po::options_description& d) {
	for(int i = 0; i < d.options().size(); ++i) {
		const std::string& str = (*d.options()[i]).long_name();
		if(!vm.count(str))
			std::cerr << "Required parameter --" << str << " is missing!\n";
	}
}

std::string default_output(const std::string& input_path) {
	std::string tmp = input_path;
	if(tmp.size() >= 4 && tmp.substr(tmp.size()-4, 4) == ".tmd")
		tmp.resize(tmp.size() - 4); // FIXME?
	return tmp + "_out.tmd";
}

std::string default_log(const std::string& input_path) {
	std::string tmp = input_path;
	if(tmp.size() >= 4 && tmp.substr(tmp.size()-4, 4) == ".tmd")
		tmp.resize(tmp.size() - 4); // FIXME?
	return tmp + "_out.log";
}


int main(int argc, char* argv[]) {
#ifdef PROFILE
ProfilerStart("/tmp/profile.out");
#endif
try {
    int seed;
    int cpu;
    std::string rigid_path, ligand_path, flex_path, config_path, out_path, log_path;
    bool score_only = false, local_only = false, randomize_only = false;
    tmd::Float box_center_x, box_center_y, box_center_z, box_size_x, box_size_y, box_size_z;
    bool help = false, help_advanced = false, version = false;

    std::string sample_type, score_type, optimizer_type;

    po::positional_options_description positional;//remains empty
    po::options_description inputs("Input");
    inputs.add_options()
        ("receptor", po::value<std::string>(&rigid_path), "rigid part of the receptor (TMD)")
        // ("flex", po::value<std::string>(&flex_path), "flexible side chains, if any (TMD)")
        ("ligand", po::value<std::string>(&ligand_path), "ligand (TMD)")
        ("sample_type", po::value<std::string>(&sample_type), "sampling_mode: flexible or rigid")
        ("score_type", po::value<std::string>(&score_type), "scoring_function: yw_score, vdw_score")
        ("optimizer_type", po::value<std::string>(&optimizer_type), "optimizer: simulated_annealing or quasi_newton")
    ;
    //options_description search_area("Search area (required, except with --score_only)");
    po::options_description search_area("Search space (required)");
    search_area.add_options()
        ("center_x", po::value<tmd::Float>(&box_center_x), "X coordinate of the center")
        ("center_y", po::value<tmd::Float>(&box_center_y), "Y coordinate of the center")
        ("center_z", po::value<tmd::Float>(&box_center_z), "Z coordinate of the center")
        ("size_x", po::value<tmd::Float>(&box_size_x), "size in the X dimension (Angstroms)")
        ("size_y", po::value<tmd::Float>(&box_size_y), "size in the Y dimension (Angstroms)")
        ("size_z", po::value<tmd::Float>(&box_size_z), "size in the Z dimension (Angstroms)")
    ;

    po::options_description outputs("Output (optional)");
    outputs.add_options()
        ("out", po::value<std::string>(&out_path), "output models (TMD), the default is chosen based on the ligand file name")
        ("log", po::value<std::string>(&log_path), "optionally, write log file")
    ;


    po::options_description advanced("Advanced options (see the manual)");
    advanced.add_options()
        ("score_only",     po::bool_switch(&score_only),     "score only - search space can be omitted")
        ("local_only",     po::bool_switch(&local_only),     "do local search only")
        ("randomize_only", po::bool_switch(&randomize_only), "randomize input, attempting to avoid clashes")
        // ("weight_gauss1", po::value<tmd::Float>(&weight_gauss1)->default_value(weight_gauss1),                "gauss_1 weight")
        // ("weight_gauss2", po::value<tmd::Float>(&weight_gauss2)->default_value(weight_gauss2),                "gauss_2 weight")
        // ("weight_repulsion", po::value<tmd::Float>(&weight_repulsion)->default_value(weight_repulsion),       "repulsion weight")
        // ("weight_hydrophobic", po::value<tmd::Float>(&weight_hydrophobic)->default_value(weight_hydrophobic), "hydrophobic weight")
        // ("weight_hydrogen", po::value<tmd::Float>(&weight_hydrogen)->default_value(weight_hydrogen),          "Hydrogen bond weight")
        // ("weight_rot", po::value<tmd::Float>(&weight_rot)->default_value(weight_rot),                         "N_rot weight")
    ;
    po::options_description misc("Misc (optional)");
    misc.add_options()
        ("cpu", po::value<int>(&cpu), "the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1)")
        ("seed", po::value<int>(&seed)->default_value(743243), "explicit random seed")
        // ("exhaustiveness", po::value<int>(&exhaustiveness)->default_value(8), "exhaustiveness of the global search (roughly proportional to time): 1+")
        // ("num_modes", po::value<int>(&num_modes)->default_value(9), "maximum number of binding modes to generate")
        // ("energy_range", po::value<tmd::Float>(&energy_range)->default_value(3.0), "maximum energy difference between the best binding mode and the worst one displayed (kcal/mol)")
    ;
    po::options_description config("Configuration file (optional)");
    config.add_options()
        ("config", po::value<std::string>(&config_path), "the above options can be put here")
    ;
    po::options_description info("Information (optional)");
    info.add_options()
        ("help",          po::bool_switch(&help), "display usage summary")
        ("help_advanced", po::bool_switch(&help_advanced), "display usage summary with advanced options")
        ("version",       po::bool_switch(&version), "display program version")
    ;
    po::options_description desc, desc_config, desc_simple;
    desc       .add(inputs).add(search_area).add(outputs).add(advanced).add(misc).add(config).add(info);
    desc_config.add(inputs).add(search_area).add(outputs).add(advanced).add(misc);
    desc_simple.add(inputs).add(search_area).add(outputs).add(misc).add(config).add(info);

    po::variables_map vm;
    try {
        //store(parse_command_line(argc, argv, desc, command_line_style::default_style ^ command_line_style::allow_guessing), vm);
        po::store(po::command_line_parser(argc, argv)
            .options(desc)
            .style(po::command_line_style::default_style ^ po::command_line_style::allow_guessing)
            .positional(positional)
            .run(),
            vm);
        notify(vm);
    }
    catch(po::error& e) {
        std::cerr << "Command line parse error: " << e.what() << std::endl << "\nCorrect usage:\n" << desc_simple << std::endl;
        return 1;
    }

    if(vm.count("config")) {
        try {
            std::ifstream config_stream(config_path);
            assert(config_stream);
            po::store(po::parse_config_file(config_stream, desc_config), vm);
            po::notify(vm);
        }
        catch(po::error& e) {
            std::cerr << "Configuration file parse error: " << e.what() << std::endl << "\nCorrect usage:\n" << desc_simple << std::endl;
            return 1;
        }
    }
    if(help) {
        std::cout << desc_simple << std::endl;
        return 0;
    }
    if(help_advanced) {
        std::cout << desc << std::endl;
        return 0;
    }
    // if(version) {
    //     std::cout << version_string << std::endl;
    //     return 0;
    // }

    //process arguments
    bool search_box_needed = !score_only; // randomize_only and local_only still need the search space
    bool output_produced   = !score_only;
    bool receptor_needed   = !randomize_only;


    if(receptor_needed) {
        if(vm.count("receptor") <= 0) {
            std::cerr << "Missing receptor.\n" << "\nCorrect usage:\n" << desc_simple << std::endl;
            return 1;
        }
    }
    if(vm.count("ligand") <= 0) {
        std::cerr << "Missing ligand.\n" << "\nCorrect usage:\n" << desc_simple << std::endl;
        return 1;
    }
    if(vm.count("sample_type") <= 0) {
        std::cerr << "Missing sampling_mode.\n" << "\nCorrect usage:\n" << desc_simple << std::endl;
        return 1;
    }
    if(vm.count("score_type") <= 0) {
        std::cerr << "Missing scoring_function.\n" << "\nCorrect usage:\n" << desc_simple << std::endl;
        return 1;
    }
    if(vm.count("optimizer_type") <= 0) {
        std::cerr << "Missing optimizer.\n" << "\nCorrect usage:\n" << desc_simple << std::endl;
        return 1;
    }
    // if(cpu < 1)
    //     cpu = 1;
    if(vm.count("seed") == 0)
        assert(false);
    //     seed = auto_seed();
    // if(exhaustiveness < 1)
    //     throw usage_error("exhaustiveness must be 1 or greater");
    // if(num_modes < 1)
    //     throw usage_error("num_modes must be 1 or greater");
    // sz max_modes_sz = static_cast<sz>(num_modes);

    // if(vm.count("flex") && !vm.count("receptor"))
    //     throw usage_error("Flexible side chains are not allowed without the rest of the receptor"); // that's the only way parsing works, actually

    std::ofstream tee = (vm.count("log") > 0) ? std::ofstream(log_path) : std::ofstream(default_log(ligand_path));

    if(search_box_needed) {
        Options_Occurrence oo = get_occurrence(vm, search_area);
        if(!oo.all) {
            check_occurrence(vm, search_area);
            std::cerr << "\nCorrect usage:\n" << desc_simple << std::endl;
            return 1;
        }
        if(box_size_x <= 0 || box_size_y <= 0 || box_size_z <= 0)
            throw Usage_Error("Search space dimensions should be positive");
    }

    tee << "cite_message" << std::endl;

    if(search_box_needed && box_size_x * box_size_y * box_size_z > 27e3) {
        tee << "WARNING: The search space volume > 27000 Angstrom^3 (See FAQ)\n";
    }

    if(output_produced) {
        //if out exist won't excute
        if(!vm.count("out")) {
            out_path = default_output(ligand_path);
            tee << "Output will be " << out_path << std::endl;
        }
    }

    if(vm.count("cpu") == 0) {
        int num_cpus = boost::thread::hardware_concurrency();
        if(num_cpus > 0)
            tee << "Detected " << num_cpus << " CPU" << ((num_cpus > 1) ? "s" : "") << std::endl;
        else
            tee << "Could not detect the number of CPUs, using 1\n";
        if(num_cpus > 0)
            cpu = num_cpus;
        else
            cpu = 1;
    }
    if(cpu < 1)
        cpu = 1;
    // if(verbosity > 1 && exhaustiveness < cpu)
    //     tee << "WARNING: at low exhaustiveness, it may be impossible to utilize all CPUs\n";

    // int seed = tmd::auto_seed();
    tee << "Using seed: " << seed << std::endl;
    tmd::RNGType generator(static_cast<tmd::RNGType::result_type>(seed));
    tee << "hello" << std::endl;
    // std::string rna_path = argv[1];
    // std::string ligand_path = argv[2];
    tee << "decalare rna..." << std::endl;
    tmd::RNA rna;
    tee << "start read rna information..." << std::endl;
    tmd::read_rna_mol2(rigid_path,rna,tee);
    tee << "finish reading rna!" << std::endl;
    // if(rna.atom_num()>=2000) {
    //     tee << "RNA contain more than " << 2000 << " atoms, skip for now!" << std::endl;
    //     return 0;
    // }

    tee << "decalare ligand..." << std::endl;
    tmd::Ligand ligand;
    tee << "start read ligand information..." << std::endl;
    tmd::read_ligand_pdbqt(ligand_path,ligand,tee);
    tee << "finish reading ligand!" << std::endl;
    tmd::print(ligand, tee);

    int num_of_dof = ligand.dof_num();
    tee << "ligand has " << num_of_dof << " dofs" << std::endl;
    // tee << "print initial ligand:" << std::endl;
    // ligand.write(tee);

    tee << "main: scoring function construct and init..." << std::endl;
    tmd::Scoring_Function* scoring_function_ptr = nullptr;
    tmd::YW_Score yw_score = tmd::YW_Score(rna,ligand,tee);
    lz::RL_Score rl_score = lz::RL_Score(rna,ligand,tee);
    // tmd::SCORE_MODE Score_Type;
    if(score_type=="yw_score") {
        // Score_Type = tmd::YW_SCORE;
        scoring_function_ptr = &yw_score;
        scoring_function_ptr->init();
    } else if(score_type=="rl_score") {
        // Score_Type = tmd::RL_SCORE;
        scoring_function_ptr = &rl_score;
        scoring_function_ptr->init();
    } else if(score_type=="vdw_ligand") {
        // Score_Type = tmd::VDW_LIGAND;
        assert(false);
    } else if(score_type=="vdw_rna_ligand") {
        // Score_Type = tmd::VDW_RNA_LIGAND;
        assert(false);
    } else if(score_type=="all") {
        // Score_Type = tmd::ALL;
        assert(false);
    } else {
        assert(false);
    }
    tee << "using " << score_type << " scoring function" << std::endl;
    tmd::SAMPLE_MODE Sample_Type;
    if(sample_type=="rigid") {
        Sample_Type = tmd::SAMPLE_RIGID;
    } else if(sample_type=="flexible") {
        Sample_Type = tmd::SAMPLE_FLEXIBLE;
    } else {
        assert(false);
    }
    tee << "using " << sample_type << " sampling mode" << std::endl;
    tmd::OPTIMIZATION_MODE Optimizer_Type;
    if(optimizer_type=="simulated_annealing") {
        Optimizer_Type = tmd::SIMULATED_ANNEALING;
    } else if(optimizer_type=="quasi_newton") {
        Optimizer_Type = tmd::QUASI_NEWTON;
    } else {
        assert(false);
    }
    tee << "using " << optimizer_type << " optimizer" << std::endl;

    tmd::Vec3d center;
    tmd::Vec3d corner1;
    tmd::Vec3d corner2;
    if(search_box_needed) {
        center = tmd::Vec3d(box_center_x, box_center_y, box_center_z);
        corner1 = tmd::Vec3d(-box_size_x/2.0, -box_size_y/2.0, -box_size_z/2.0) + center;
        corner2 = tmd::Vec3d(box_size_x/2.0, box_size_y/2.0, box_size_z/2.0) + center;
        tee << "box center: ";
        tmd::print(center, tee);
        tee << std::endl;
        tee << "box size: " << box_size_x << " " << box_size_y << " " << box_size_z << std::endl;
        tee << "box corner1 and corner2: ";
        tmd::print(corner1, tee);
        tee << " ";
        tmd::print(corner2, tee);
        tee << std::endl;
    }

    // corner1 = tmd::Vec3d(-8/2.0, -8/2.0, -8/2.0) + center;
    // corner2 = tmd::Vec3d(8/2.0, 8/2.0, 8/2.0) + center;
    // std::ofstream outligand(out_path);
    // tmd::Floats ff(ligand.dof_num(),0.0);
    // ff[6] = tmd::k_pi;
    // // ff[7] = tmd::k_pi;
    // ligand.sample(ff,tmd::LIGAND_SAMPLE_FLEXIBLE);
    // std::cout << ligand.rmsd_with_respect_to_ref_atoms() << std::endl;
    // ligand.write(outligand);
    // outligand.close();
    // exit(2);

    //get pocket
    tmd::Pocket pocket(tee);
    // pocket.find_pocket(rna);
    tee << "find " << pocket.sites.size() << " pocket sites" << std::endl;

    // std::ofstream pocket_file("./pocket.mol2");
    // //out pocket to file in pdb format
    // // tmd::Atom pocket_atom("Mg", tmd::SYBYL_AD4);
    // // pocket_atom.set_name("Mg");
    // // pocket_atom.set_element_type(tmd::EL_MG);
    // // pocket_atom.set_res_name("Mg");
    // // pocket_atom.set_chain_name("A");
    // pocket_file << "@<TRIPOS>MOLECULE" << std::endl;
    // pocket_file << "pocket_sites" << std::endl;
    // pocket_file << pocket.sites.size() << " 26 1 0 0" << std::endl;
    // pocket_file << "PROTEIN" << std::endl;
    // pocket_file << "AMBER ff14SB" << std::endl;
    // pocket_file << "@<TRIPOS>ATOM" << std::endl;
    // int pocket_site_count = 0;
    // for(const auto& ps : pocket.sites) {
    //     pocket_site_count++;
    //     // pocket_atom.set_coord(tmd::Vec3d(ps.x,ps.y,ps.z));
    //     // tmd::print(pocket_atom,pocket_file);
    //     // pocket_atom.set_res_serial(pocket_site_count);
    //     // pocket_atom.set_serial(pocket_site_count);

    //     pocket_file << "  " << pocket_site_count << " " << "Mg" << " " << std::fixed << std::setprecision(4) << ps.x << " " << ps.y << " " << ps.z << " " << "Mg" << " " << pocket_site_count << " " << "Mg" << " " << 2.00 << std::endl;
	// 	pocket_file << std::defaultfloat << std::setprecision(6);
    // }
    // pocket_file.close();


    tmd::Sampling sampling(rna,ligand,scoring_function_ptr,tmd::Box(center,corner1,corner2),pocket,Sample_Type,Optimizer_Type,generator,tee);

    tee << "main: scoring evaluate for input conf ------> " << scoring_function_ptr->evaluate() << std::endl;

    if(randomize_only) {
        tee << "finish randomize only!" << std::endl;
        return 0;
    }

    if(score_only) {
        tee << "scoring evaluate: " << scoring_function_ptr->evaluate() << std::endl;
        tee << "finish score only!" << std::endl;
        return 0;
    }

    sampling.docking();

    tee << "finish docking" << std::endl;

    tmd::Conformers conformers = sampling.get_conformers();
    conformers.sort_by_score();
    conformers.prune_by_cutoff(100);
    tee << "output ligand: " << out_path << std::endl;
    // const tmd::Sampled_Ligands& clustered_ligands = sampling.cluster();
    tee << "has " << conformers.conformer_num() << " conformers stored" << std::endl;
    std::ofstream outligand_final(out_path);
    assert(outligand_final);
    int model_index = 0;
    int conformer_num_step = conformers.conformer_num() < 500 ? 1 : conformers.conformer_num()/500.0;
    for(int ci = 0; ci < conformers.conformer_num(); ci += conformer_num_step) {
        if(ci >= conformers.conformer_num()) break;
        tee << "conformer index: " << ci << " model index: " << model_index << " rmsd: " << conformers.get_conformers()[ci].rmsd << " score: " << conformers.get_conformers()[ci].score << std::endl;
        assert(ci < conformers.conformer_num());
        conformers.output_conformer(ci,++model_index,outligand_final);
        if(model_index>=500) {
            break;
        }
    }
    outligand_final.close();
    tee << std::endl;

    return 0;
}
// catch(tmd::File_Error& e) {
//     std::cerr << "\n\nError: could not open \"" << e.name.string() << "\" for " << (e.in ? "reading" : "writing") << ".\n";
//     return 1;
// }
// catch(std::filesystem_error& e) {
//     std::cerr << "\n\nFile system error: " << e.what() << std::endl;
//     return 1;
// }
catch(Usage_Error& e) {
    std::cerr << "\n\nUsage error: " << e.what() << ".\n";
    return 1;
}
catch(std::bad_alloc&) {
    std::cerr << "\n\nError: insufficient memory!\n";
    return 1;
}

// Errors that shouldn't happen:

catch(std::exception& e) {
    std::cerr << "\n\nAn error occurred: " << e.what() << ".\n";
    return 1;
}
catch(...) {
    std::cerr << "\n\nAn unknown error occurred.\n";
    return 1;
}

#ifdef PROFILE
ProfilerStop();
#endif
}
