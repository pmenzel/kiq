#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <string>
#include <array>
#include <string>
#include <atomic>
#include <deque>
#include <stdexcept>

#include "BooPHF/BooPHF.h"
#include "util.hpp"


void usage_kdump() {
	print_usage_header();
	fprintf(stderr, "Usage:\n   kiq dump -i <file> -k <file> -p <mode>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of index file\n");
	fprintf(stderr, "   -k FILENAME   Name of k-mer count database file\n");
	fprintf(stderr, "   -p STRING     Mode is either db, metadata, stats, long\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -v            Enable verbose output\n");
	fprintf(stderr, "   -d            Enable debug output.\n");
	exit(EXIT_FAILURE);
}

int main_kdump(int argc, char** argv) {

	bool debug = false;
	bool verbose = false;

	std::string filename_index;
	std::string filename_db;
	std::string mode;

	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hdvi:k:p:")) != -1) {
		switch (c)  {
			case 'h':
				usage_kdump();
			case 'd':
				debug = true; break;
			case 'v':
				verbose = true; break;
			case 'k':
				filename_db = optarg; break;
			case 'i':
				filename_index = optarg; break;
			case 'p':
				mode = optarg; break;
			default:
				usage_kdump();
		}
	}
	if(filename_index.length() == 0) { error("Please specify the name of the index file, using the -i option."); usage_kdump(); }
	if(filename_db.length() == 0) { error("Please specify the name of the database file, using the -k option."); usage_kdump(); }
	if(mode.length()==0) { error("Specify mode with option -p."); usage_kdump(); }

	boophf_t * bphf = new boomphf::mphf<u_int64_t,hasher_t>();
	load_index(filename_index, bphf);

	KmerIndex n_elem = (KmerIndex)bphf->nbKeys();
	std::vector<Kmer> initial_kmers;
	pCountMap * kmer2countmap = new pCountMap[n_elem](); // init new array of size n_elem
	ExpId2Name exp_id2name;
	ExpId2Desc exp_id2desc;
	ExpName2Id exp_name2id;
	ExpId2ReadCount exp_id2readcount;

	try {
		read_database(filename_db, initial_kmers, kmer2countmap, bphf, true, exp_id2name, exp_id2desc, exp_name2id, exp_id2readcount);
	}
	catch(std::runtime_error e) {
		std::cerr << "Error while reading database (" << e.what() << ")." << std::endl;
		exit(EXIT_FAILURE);
	}

	if(mode=="db") {
		for(auto it : initial_kmers) {
			KmerIndex i = bphf->lookup(it);
			assert(i < n_elem);

			// print k-mer index
			std::cout << int_to_str(it) << "\t";
			// print number of experiments having this k-mer
			uint32_t num_exp = kmer2countmap[i]==nullptr ? 0 : (uint32_t)(kmer2countmap[i]->size());
			std::cout << num_exp;
			if(num_exp > 0) {
				for(auto const & it : *kmer2countmap[i]) {
					std::cout << "\tid=" << it.first << " count=" << it.second;
				}
			}
			std::cout << "\n";
		}
	}
	else if(mode=="long") {
		for(auto it : initial_kmers) {
			KmerIndex i = bphf->lookup(it);
			assert(i < n_elem);
			uint32_t num_exp = kmer2countmap[i]==nullptr ? 0 : (uint32_t)(kmer2countmap[i]->size());
			if(num_exp > 0) {
				std::string kmer = int_to_str(it) ;
				for(auto const & it : *kmer2countmap[i]) {
					double rpm = (double)it.second / (double)exp_id2readcount.at(it.first) * 1e6;
					std::string exp_name = exp_id2name.at(it.first);
					std::cout << kmer << "\t" << exp_name << "\t" << it.second << "\t" << rpm << "\n";
				}
			}
		}
	}
	else if(mode=="metadata") {
		// save experiment id to name mapping
		for(auto const & it : exp_id2name) {
			ExperimentId exp_id = it.first;
			std::string exp_name = it.second;
			ReadCount readcount = exp_id2readcount.at(exp_id);
			std::string exp_desc = "NA";
			auto it_desc = exp_id2desc.find(exp_id);
			if(it_desc != exp_id2desc.end() && it_desc->second.length() > 0) exp_desc = it_desc->second; 
			std::cout << exp_id << "\t" << exp_name << "\t" << readcount << "\t" << exp_desc << "\n";
		}
	}
	else if(mode=="stats") {
		// number of experiments in db
		std::cout << "Number of experiments\t" << exp_id2name.size() <<"\n";
		std::cout << "Number of k-mers\t" << n_elem <<"\n";
		// count number of k-mers with at least one experiment
		int count = 0;
		for(auto it : initial_kmers) {
			KmerIndex i = bphf->lookup(it);
			assert(i < n_elem);
			// print k-mer index
			uint32_t num_exp = kmer2countmap[i]==nullptr ? 0 : (uint32_t)(kmer2countmap[i]->size());
			if(num_exp > 0) {
				count++;
			}
		}
		std::cout << "K-mers with experiments\t" << count << "\n";
	}


	for(KmerIndex i = 0; i < n_elem;i++) {
		if(kmer2countmap[i] != nullptr) {
			delete kmer2countmap[i];
		}
	}
	delete[] kmer2countmap;
	delete bphf;

	return 0;

}


