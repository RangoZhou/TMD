#pragma once

#include <iostream>
#include <fstream>
#include "file.h"

namespace tmd {

struct Tee {
	OFile* of;
	Tee() : of(NULL) {}
	void init(const Path& name) {
		of = new OFile(name);
	}
	virtual ~Tee() { delete of; }
	void flush() {
		std::cout << std::flush;
		if(of)
			(*of) << std::flush;
	}
	void endl() {
		std::cout << std::endl;
		if(of)
			(*of) << std::endl;
	}
	void setf(std::ios::fmtflags a) {
		std::cout.setf(a);
		if(of)
			of->setf(a);
	}
	void setf(std::ios::fmtflags a, std::ios::fmtflags b) {
		std::cout.setf(a, b);
		if(of)
			of->setf(a, b);
	}
};

template<typename T>
Tee& operator<<(Tee& out, const T& x) {
	std::cout << x;
	if(out.of)
		(*out.of) << x;
	return out;
}

}