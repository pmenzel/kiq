/*************************************************
	KIQ

	Author: Peter Menzel <pmenzel@gmail.com>

	Copyright (C) 2018 Peter Menzel

	See the file README.md for documentation.
**************************************************/

#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <deque>
#include <stdexcept>
#include <cstring>

#include "version.hpp"
#include "util.hpp"
#include "kindex.hpp"
#include "kquery.hpp"
#include "kdb.hpp"
#include "kdump.hpp"
#ifdef KIQ_SRA
#include "ksra.hpp"
#endif



int main(int argc, char** argv) {

	struct HeaderDbFile hf;
	struct HeaderDbKmers hk;
	struct HeaderDbMetadata hm;
	hf.dbVer = 1;
	hk.numKmer = 2;
	hm.numExp = 2;

	Kmer k1 = 222;
	Kmer k2 = 3333222;
	ExperimentCount c1 = 2;
	ExperimentCount c2 = 0;
	ExperimentId i1 = 1;
	KmerCount kc1 = 200;
	ExperimentId i2 = 2;
	KmerCount kc2 = 500;

	std::ofstream os("testdb.bin", std::ios::out | std::ios::binary);
	if(!os.is_open()) {  error("Could not open file "); exit(EXIT_FAILURE); }
	os.write(reinterpret_cast<const char *>(&hf.kiq),sizeof(hf.kiq));
	os.write(reinterpret_cast<const char *>(&hf.dbVer),sizeof(hf.dbVer));
	os.write(reinterpret_cast<const char *>(&hk.numKmer),sizeof(hk.numKmer));
	os.write(reinterpret_cast<const char *>(&k1),sizeof(k1));
	os.write(reinterpret_cast<const char *>(&c1),sizeof(c1));
	os.write(reinterpret_cast<const char *>(&i1),sizeof(i1));
	os.write(reinterpret_cast<const char *>(&kc1),sizeof(kc1));
	os.write(reinterpret_cast<const char *>(&i2),sizeof(i2));
	os.write(reinterpret_cast<const char *>(&kc2),sizeof(kc2));
	os.write(reinterpret_cast<const char *>(&k2),sizeof(k2));
	os.write(reinterpret_cast<const char *>(&c2),sizeof(c2));
	os.write(reinterpret_cast<const char *>(&hm.label),sizeof(hm.label));
	os.write(reinterpret_cast<const char *>(&hm.numExp),sizeof(hm.numExp));
	os.close();

}

