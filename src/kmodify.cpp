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


void usage_kmodify() {
	print_usage_header();
	fprintf(stderr, "Usage:\n   kiq modify -i <file> -k <file> -c <file>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of index file\n");
	fprintf(stderr, "   -k FILENAME   Name of k-mer count database file\n");
	fprintf(stderr, "   -c FILENAME   Name of TSV file with modification commands\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -v            Enable verbose output\n");
	fprintf(stderr, "   -d            Enable debug output.\n");
	exit(EXIT_FAILURE);
}

int main_kmodify(int argc, char** argv) {

	bool debug = false;
	bool verbose = false;

	std::string filename_index;
	std::string filename_db;
	std::string filename_commands;

	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hdvi:k:c:")) != -1) {
		switch (c)  {
			case 'h':
				usage_kmodify();
			case 'd':
				debug = true; break;
			case 'v':
				verbose = true; break;
			case 'k':
				filename_db = optarg; break;
			case 'i':
				filename_index = optarg; break;
			case 'c':
				filename_commands = optarg; break;
			default:
				usage_kmodify();
		}
	}
	if(filename_index.length() == 0) { error("Please specify the name of the index file, using the -i option."); usage_kmodify(); }
	if(filename_db.length() == 0) { error("Please specify the name of the database file, using the -k option."); usage_kmodify(); }
	if(filename_commands.length() == 0) { error("Please specify the name of the modification file, using the -l option."); usage_kmodify(); }

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


	std::ifstream ifs(filename_commands);
	if(!ifs) { std::cerr << "Cannot open file " << filename_commands << std::endl; exit(EXIT_FAILURE); }

	std::string line;
	while(getline(ifs, line)) {
		if(line.length() == 0) { continue; }
		size_t tab = line.find_first_of('\t');
		if(tab == std::string::npos || tab == 0) { std::cerr << "Warning: Skipping line due to bad format: " << line << "\n"; continue; }

		std::string experiment_name = line.substr(0,tab);
		std::string command;
		std::string argument;
		size_t tab2 = line.find_first_of('\t',tab+1);
		if(tab2 == std::string::npos) { // no second tab found in line, then filename is from first tab to the line end
			command = line.substr(tab+1);
		}
		else { //second tab found in line, filename is between both tabs and description is from 2nd tab to line end
			command = line.substr(tab+1,tab2-tab-1);
			argument = line.substr(tab2+1);
		}

		if(exp_name2id.find(experiment_name) == exp_name2id.end()) {
			std::cerr << "Experiment with name " << experiment_name << " is not contained in database.\n";
			continue;
		}
		if(!(command=="delete" || command=="update_desc")) {
			std::cerr << "Warning: Wrong command in line: " << line << "\n";
			continue;
		}

		if(command=="delete") {
			ExperimentId exp_id = exp_name2id.at(experiment_name);
			// delete experiment from all k-mers
			for(Kmer it : initial_kmers) {
				KmerIndex index = bphf->lookup(it);
				ExperimentCount num_exp = kmer2countmap[index]==nullptr ? 0 : static_cast<ExperimentCount>(kmer2countmap[index]->size());
				if(num_exp > 0) {
					pCountMap m = kmer2countmap[index];
					auto it = m->find(exp_id);
					if(it != m->end()) {
						m->erase(it);
					}
				}
			}
			exp_name2id.erase(experiment_name);
			exp_id2name.erase(exp_id);
			exp_id2desc.erase(exp_id);
			exp_id2readcount.erase(exp_id);
		}
		else if(command=="update_desc") {
			if(argument.length()==0) {
				std::cerr << "Warning: Missing column 3 in line: " << line << "\n";
				continue;
			}
			ExperimentId exp_id = exp_name2id.at(experiment_name);
			exp_id2desc[exp_id] = argument;
		}

	}


	write_database(filename_db, initial_kmers, kmer2countmap, bphf, exp_id2name, exp_id2desc, exp_id2readcount);

	for(KmerIndex i = 0; i < n_elem;i++) {
		if(kmer2countmap[i] != nullptr) {
			delete kmer2countmap[i];
		}
	}
	delete[] kmer2countmap;
	delete bphf;

	return 0;

}


