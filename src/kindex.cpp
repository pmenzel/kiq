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
#include <stdexcept>

#include "BooPHF/BooPHF.h"
#include "util.hpp"

void usage_kindex() {
	print_usage_header();
	fprintf(stderr, "Usage:\n   kiq index -i <file> -k <file> -l <file>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i <file>   Name of index file\n");
	fprintf(stderr, "   -k <file>   Name of k-mer count database file\n");
	fprintf(stderr, "   -l <file>   Name of file containing k-mers to be indexed\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -v          Enable verbose output\n");
	fprintf(stderr, "   -d          Enable debug output.\n");
	exit(EXIT_FAILURE);
}

int main_kindex(int argc, char** argv) {


	std::string filename_kmers;
	std::string filename_db;
	std::string filename_index;
	bool debug = false;
	bool verbose = false;

	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hdvai:k:l:z:")) != -1) {
		switch (c)  {
			case 'h':
				usage_kindex();
			case 'd':
				debug = true; break;
			case 'v':
				verbose = true; break;
			case 'k':
				filename_db = optarg; break;
			case 'i':
				filename_index = optarg; break;
			case 'l':
				filename_kmers = optarg; break;
			default:
				usage_kindex();
		}
	}
	if(filename_index.length() == 0) { error("Please specify the name of the index file, using the -i option."); usage_kindex(); }
	if(filename_kmers.length() == 0) { error("Please specify the name of the file with k-mers, using the -l option."); usage_kindex(); }
	if(filename_db.length() == 0) { error("Please specify the name of the database file, using the -k option."); usage_kindex(); }

	std::ifstream filestream_kmers;
	filestream_kmers.open(filename_kmers);
	if(!filestream_kmers.is_open()) { error("Could not open file " + filename_kmers); exit(EXIT_FAILURE); }

	std::cerr << getCurrentTime() << " Start reading k-mers from file " << filename_kmers << "\n";
	std::string line_from_file;
	line_from_file.reserve(KMER_K+1);

	// first read all lines into a set, to remove duplicates
	std::set<Kmer> input_kmers_set;
	while(getline(filestream_kmers,line_from_file)) {
		if(line_from_file.length() == 0) { continue; }
		if(line_from_file.length() < KMER_K){ std::cerr << "Warning: Line too short: " << line_from_file << "\n"; continue; }
		if(line_from_file.length() > KMER_K){ std::cerr << "Warning: Line too long: " << line_from_file << "\n"; continue; }
		Kmer kmer = str_to_int(line_from_file);
		input_kmers_set.emplace(kmer);
	}
	filestream_kmers.close();

	// copy k-mers from set into vector
	std::vector<Kmer> initial_kmers;
	initial_kmers.insert(initial_kmers.end(), input_kmers_set.cbegin(), input_kmers_set.cend());

	// sort k-mers, so that they will be saved in proper order
	std::sort(initial_kmers.begin(), initial_kmers.end());

	std::cerr << getCurrentTime() << " Calculating hash functions for " << initial_kmers.size() << " k-mers\n";
	boophf_t * bphf = new boomphf::mphf<u_int64_t,hasher_t>(initial_kmers.size(),initial_kmers,1);

	write_initial_database(filename_db, initial_kmers);

	std::cerr << getCurrentTime() << " Writing index to file " << filename_index << "\n";
	std::ofstream os(filename_index, std::ofstream::out);
	if(!os.is_open()) { error("Could not open file " + filename_index); exit(EXIT_FAILURE); }
	bphf->save(os);
	os.close();

	delete bphf;

	return 0;

}


