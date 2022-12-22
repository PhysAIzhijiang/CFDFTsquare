// This file is under public domain.
#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <iostream>
#include <chrono>
#include <string>

// print code time consumption
#define TIMEIT_LOG(SS, CODE) {						\
		using namespace std::chrono;				\
		auto t0 = high_resolution_clock::now();			\
		{ CODE; }						\
		auto t1 = high_resolution_clock::now();			\
		SS << ""#CODE":\t" << duration_cast<milliseconds>(t1-t0).count() << "msec"<<std::endl; \
	}

#define TIMEIT(CODE) TIMEIT_LOG(std::cout, CODE)

#define INFO_LOG(SS,VAR) {					\
		SS << ""#VAR": " << (VAR) << std::endl;	\
	}

#define INFO(VAR) INFO_LOG(std::cout, VAR)

#define NVSTRING(VAR) (""#VAR": ")<<(VAR)

#include "type.hpp"   //to check variable type and ptint info to file SS
#define TYPEINFO_LOG(SS,VAR){						\
		SS << "decltype(VAR) is " << type_name<decltype(VAR)>() << std::endl; \
	}


namespace TicToc { // similar to matlab
	using namespace std::chrono;
	auto t0 = high_resolution_clock::now();
	void tic() {
		t0 = high_resolution_clock::now();
	}
	double toc(std::ostream &out=std::cout, std::string msg = "") {
		auto t1 = high_resolution_clock::now();
		const double dt = duration_cast<milliseconds>(t1-t0).count();
		out << dt << "msec"
			  << '\t' << msg << std::endl;
		return dt;
	}
}

#endif
