/*************************************************
	KIQ

	Author: Peter Menzel <pmenzel@gmail.com>

	Copyright (C) 2018,2019 Peter Menzel

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

	boophf_t * bphf = new boomphf::mphf<u_int64_t,hasher_t>();
	load_index("index_fly_BSJ.bin", bphf);

	KmerIndex n_elem = (KmerIndex)bphf->nbKeys();
	std::vector<Kmer> initial_kmers;
	initial_kmers.reserve(n_elem+1);
	pCountMap * kmer2count = new pCountMap[n_elem](); // init new array of size n_elem

	ExpId2Name exp_id2name;
	ExpName2Id exp_name2id;
	ExpId2ReadCount exp_id2readcount;

	read_database("testdb.bin", initial_kmers, kmer2count, bphf, true, exp_id2name, exp_name2id, exp_id2readcount);

	

}

