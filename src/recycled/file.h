#pragma once

#include <filesystem>
#include <fstream>

namespace tmd {

typedef std::filesystem::path Path;

struct File_Error {
	Path name;
	bool in;
	File_Error(const Path& name_, bool in_) : name(name_), in(in_) {}
};

struct IFile : public std::ifstream {
	IFile(const Path& name) : std::ifstream(name) {
		if(!(*this))
			throw File_Error(name, true);
	}
	IFile(const Path& name, std::ios_base::openmode mode) : std::ifstream(name, mode) {
		if(!(*this))
			throw File_Error(name, true);
	}
};

struct OFile : public std::ofstream {
	OFile(const Path& name) : std::ofstream(name) {
		if(!(*this))
			throw File_Error(name, false);
	}
	OFile(const Path& name, std::ios_base::openmode mode) : std::ofstream(name, mode) {
		if(!(*this))
			throw File_Error(name, false);
	}
};

}