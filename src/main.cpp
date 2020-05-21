#include <iostream>
#include <boost/program_options.hpp>
#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?
#include "TMDConfig.h"
#include "common.h"
#include "molecule.h"
#include "mol2.h"
#include "pdbqt.h"
#include "print.h"
#include "score.h"
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
	for(tmd::Size_Type i = 0; i < d.options().size(); ++i) {
		if(vm.count((*d.options()[i]).long_name()))
			tmp.some = true;
		else
			tmp.all = false;
    }
	return tmp;
}

void check_occurrence(po::variables_map& vm, po::options_description& d) {
	for(tmd::Size_Type i = 0; i < d.options().size(); ++i) {
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
    int seed, cpu;
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
        ("optimizer_type", po::value<std::string>(&optimizer_type), "optimizer: stimulated_anealing")
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

    std::ofstream log = (vm.count("log") > 0) ? std::ofstream(log_path) : std::ofstream(default_log(ligand_path));

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

    log << "cite_message" << std::endl;

    if(search_box_needed && box_size_x * box_size_y * box_size_z > 27e3) {
        log << "WARNING: The search space volume > 27000 Angstrom^3 (See FAQ)\n";
    }

    if(output_produced) {
        //if out exist won't excute
        if(!vm.count("out")) {
            out_path = default_output(ligand_path);
            log << "Output will be " << out_path << std::endl;
        }
    }

    if(vm.count("cpu") == 0) {
        unsigned num_cpus = boost::thread::hardware_concurrency();
        if(num_cpus > 0)
            log << "Detected " << num_cpus << " CPU" << ((num_cpus > 1) ? "s" : "") << std::endl;
        else
            log << "Could not detect the number of CPUs, using 1\n";
        if(num_cpus > 0)
            cpu = num_cpus;
        else
            cpu = 1;
    }
    if(cpu < 1)
        cpu = 1;
    // if(verbosity > 1 && exhaustiveness < cpu)
    //     log << "WARNING: at low exhaustiveness, it may be impossible to utilize all CPUs\n";

    // int seed = tmd::auto_seed();
    log << "Using seed: " << seed << std::endl;
    tmd::RNGType generator(static_cast<tmd::RNGType::result_type>(seed));
    log << "hello" << std::endl;
    // std::string rna_path = argv[1];
    // std::string ligand_path = argv[2];
    log << "decalare rna..." << std::endl;
    tmd::RNA rna;
    log << "start read rna information..." << std::endl;
    tmd::read_rna_mol2(rigid_path,rna,log);
    log << "finish reading rna!" << std::endl;
    // if(rna.atom_num()>=2000) {
    //     log << "RNA contain more than " << 2000 << " atoms, skip for now!" << std::endl;
    //     return 0;
    // }

    log << "decalare ligand..." << std::endl;
    tmd::Ligand ligand;
    log << "start read ligand information..." << std::endl;
    tmd::read_ligand_pdbqt(ligand_path,ligand,log);
    log << "finish reading ligand!" << std::endl;
    tmd::print(ligand, log);

    tmd::Floats::size_type num_of_dof = ligand.dof_num();
    log << "ligand has " << num_of_dof << " dofs" << std::endl;
    // log << "print initial ligand:" << std::endl;
    // ligand.write(log);

    //
    tmd::SCORE_MODE Score_Type;
    if(score_type=="yw_score") {
        Score_Type = tmd::YW_SCORE;
    } else if(score_type=="rl_score") {
        Score_Type = tmd::RL_SCORE;
    } else if(score_type=="vdw_ligand") {
        Score_Type = tmd::VDW_LIGAND;
    } else if(score_type=="vdw_rna_ligand") {
        Score_Type = tmd::VDW_RNA_LIGAND;
    } else if(score_type=="all") {
        Score_Type = tmd::ALL;
    } else {
        assert(false);
    }
    log << "using " << score_type << " scoring function" << std::endl;
    tmd::SAMPLE_MODE Sample_Type;
    if(sample_type=="rigid") {
        Sample_Type = tmd::SAMPLE_RIGID;
    } else if(sample_type=="flexible") {
        Sample_Type = tmd::SAMPLE_FLEXIBLE;
    } else {
        assert(false);
    }
    log << "using " << sample_type << " sampling mode" << std::endl;
    tmd::OPTIMIZATION_MODE Optimizer_Type;
    if(optimizer_type=="stimulated_anealing") {
        Optimizer_Type = tmd::SIMULATED_ANEALING;
    } else {
        assert(false);
    }
    log << "using " << optimizer_type << " optimizer" << std::endl;

    log << "scoring function init..." << std::endl;
    tmd::Scoring_Function scoring_function(rna,ligand,Score_Type,log);
    //SCORE_TYPE is enum and it has {YW_SCORE,VDW_LIGAND,VDW_RNA_LIGAND,ALL};
    // scoring_function.set_score_type(Score_Type);

    tmd::Vec3d center;
    tmd::Vec3d corner1;
    tmd::Vec3d corner2;
    if(search_box_needed) {
        center = tmd::Vec3d(box_center_x, box_center_y, box_center_z);
        corner1 = tmd::Vec3d(-box_size_x/2.0, -box_size_y/2.0, -box_size_z/2.0) + center;
        corner2 = tmd::Vec3d(box_size_x/2.0, box_size_y/2.0, box_size_z/2.0) + center;
        log << "box center: ";
        tmd::print(center, log);
        log << std::endl;
        log << "box size: " << box_size_x << " " << box_size_y << " " << box_size_z << std::endl;
        log << "box corner1 and corner2: ";
        tmd::print(corner1, log);
        log << " ";
        tmd::print(corner2, log);
        log << std::endl;
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

    tmd::Sampling sampling(rna,ligand,scoring_function,tmd::Box(center,corner1,corner2),Sample_Type,Optimizer_Type,generator,log);

    if(randomize_only) {
        log << "finish randomize only!" << std::endl;
        return 0;
    }

    if(score_only) {
        log << "scoring evaluate: " << scoring_function.evaluate() << std::endl;
        log << "finish score only!" << std::endl;
        return 0;
    }

    sampling.docking();
    log << "output ligand: " << out_path << std::endl;
    // const tmd::Sampled_Ligands& clustered_ligands = sampling.cluster();
    std::ofstream outligand_final(out_path);
    assert(outligand_final);
    unsigned int model_index = 0;
    unsigned int conformer_num_step = sampling.num_conformer() < 500 ? 1 : sampling.num_conformer()/500.0;
    for(int ci = sampling.num_conformer()-1; ci >= 0;) {
        // std::cout << ci << " " << model_index << std::endl;
        sampling.output_conformer(ci,++model_index,outligand_final);
        if(model_index>=500) {
            break;
        }
        ci -= conformer_num_step;
    }
    outligand_final.close();
    log << std::endl;

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
