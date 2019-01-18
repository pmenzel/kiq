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


int main(int argc, char** argv) {

	std::string filename_index;
	std::string filename_db_in;
	std::string filename_meta_in;
	std::string filename_db_out;
	
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "i:k:m:o:")) != -1) {
		switch (c)  {
			case 'm':
				filename_meta_in = optarg; break;
			case 'k':
				filename_db_in = optarg; break;
			case 'i':
				filename_index = optarg; break;
			case 'o':
				filename_db_out = optarg; break;
		}
	}

	if(filename_index.length() == 0) { error("Please specify the name of the index file, using the -i option."); exit(EXIT_FAILURE);  }
	if(filename_db_in.length() == 0) { error("Please specify the name of the database input file, using the -k option."); exit(EXIT_FAILURE); }
	if(filename_db_out.length() == 0) { error("Please specify the name of the database output file, using the -o option."); exit(EXIT_FAILURE); }
	if(filename_meta_in.length() == 0) { error("Please specify the name of the metadata input file, using the -m option."); exit(EXIT_FAILURE); }

	boophf_t * bphf = new boomphf::mphf<u_int64_t,hasher_t>();
	load_index(filename_index, bphf);

	KmerIndex n_elem = (KmerIndex)bphf->nbKeys();
	std::vector<Kmer> initial_kmers;
	initial_kmers.reserve(n_elem+1);
	pCountMap * kmer2count = new pCountMap[n_elem](); // init new array of size n_elem

	ExpId2Name exp_id2name;
	ExpId2Desc exp_id2desc;
	ExpName2Id exp_name2id;
	ExpId2ReadCount exp_id2readcount;

	read_kmer_database(filename_db_in,initial_kmers,kmer2count,bphf,true);
	read_experiment_database(filename_meta_in, exp_id2name, exp_name2id, exp_id2readcount);

	write_database(filename_db_out, initial_kmers, kmer2count, bphf, exp_id2name, exp_id2desc, exp_id2readcount);
	

}

