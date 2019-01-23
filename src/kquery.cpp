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

#include "BooPHF/BooPHF.h"
#include "util.hpp"

void usage_kquery();
void run_query(const std::string & query, bool json, uint32_t threshold, uint32_t rpm_threshold, const std::vector<Kmer> & initial_kmers, boophf_t * bphf, pCountMap * kmer2countmap, const ExpId2Name & exp_id2name, const ExpId2Desc & exp_id2desc, const ExpId2ReadCount & exp_id2readcount);
void run_query_all(const std::string & query, bool first, std::set<ExperimentId> &, uint32_t threshold, uint32_t rpm_threshold, const std::vector<Kmer> & initial_kmers, boophf_t * bphf, pCountMap * kmer2countmap, const ExpId2ReadCount & exp_id2readcount);
void print_exp_set(std::set<ExperimentId> & exp_set, ExpId2Name & exp_id2name, ExpId2Desc & exp_id2desc, bool json);

int main_kquery(int argc, char** argv) {

	std::string filename_index;
	std::string filename_db;
	std::string filename_query;
	std::string arg_query;
	bool debug = false;
	bool verbose = false;
	bool json = false;
	bool all_kmers = false;
	uint32_t threshold = 0;
	uint32_t rpm_threshold = 0;

	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hjadvr:t:i:k:Q:q:z:")) != -1) {
		switch (c)  {
			case 'h':
				usage_kquery();
			case 'd':
				debug = true; break;
			case 'a':
				all_kmers = true; break;
			case 'j':
				json = true; break;
			case 'v':
				verbose = true; break;
			case 'k':
				filename_db = optarg; break;
			case 'i':
				filename_index = optarg; break;
			case 'q':
				filename_query = optarg; break;
			case 'Q':
				arg_query = optarg; break;
			case 't': {
				try {
					threshold = std::stoi(optarg);
				}
				catch(const std::invalid_argument& ia) {
					std::cerr << "Invalid argument in -t " << optarg << std::endl;
				}
				catch (const std::out_of_range& oor) {
					std::cerr << "Invalid argument in -t " << optarg << std::endl;
				}
				break;
			}
			case 'r': {
				try {
					rpm_threshold = std::stoi(optarg);
				}
				catch(const std::invalid_argument& ia) {
					std::cerr << "Invalid argument in -r " << optarg << std::endl;
				}
				catch (const std::out_of_range& oor) {
					std::cerr << "Invalid argument in -r " << optarg << std::endl;
				}
				break;
			}
			default:
				usage_kquery();
		}
	}
	if(filename_index.length() == 0) { error("Please specify the name of the index file, using the -i option."); usage_kquery(); }
	if(filename_db.length() == 0) { error("Please specify the name of the database file, using the -k option."); usage_kquery(); }
	if(filename_query.length() == 0 && arg_query.length() == 0) { error("Please specify either a query file with -q or the query k-mer(s) directly with -Q."); usage_kquery(); }
	if(filename_query.length() > 0 && arg_query.length() > 0) { error("Please specify only one of the options -q and -Q."); usage_kquery(); }

	boophf_t * bphf = new boomphf::mphf<u_int64_t,hasher_t>();
	load_index(filename_index,bphf);

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

	std::cerr << getCurrentTime() << " Runninq query " << arg_query << "\n";

	if(arg_query.length() > 0) {
		size_t pos = std::string::npos;
		if((pos = arg_query.find(",")) !=  std::string::npos) { // multi k-mer query
			std::string q = arg_query.substr(0,pos);
			if(!all_kmers) {
				if(json) { std::cout << "{ \"results\" : [ "; }
				run_query(q, json, threshold, rpm_threshold, initial_kmers, bphf, kmer2countmap, exp_id2name, exp_id2desc, exp_id2readcount);
				size_t start = pos+1;
				while((pos = arg_query.find(",",start)) != std::string::npos) {
					q = arg_query.substr(start,pos - start);
					if(json) { std::cout << ", "; }
					run_query(q, json, threshold, rpm_threshold, initial_kmers, bphf, kmer2countmap, exp_id2name, exp_id2desc, exp_id2readcount);
					start = pos +1;
				}
				if(json) { std::cout << " ] }\n"; }
			}
			else { // all k-mers
				std::set<ExperimentId> exp_set;
				// process first k-mer
				run_query_all(q, true, exp_set, threshold, rpm_threshold, initial_kmers, bphf, kmer2countmap, exp_id2readcount);
				if(!exp_set.empty()) {
					// do remaining k-mers
					size_t start = pos+1;
					while((pos = arg_query.find(",",start)) != std::string::npos) {
						q = arg_query.substr(start,pos - start);
						run_query_all(q, false, exp_set, threshold, rpm_threshold, initial_kmers, bphf, kmer2countmap, exp_id2readcount);
						start = pos +1;
					}
					// last k-mer after last ,
					q = arg_query.substr(start);
					run_query_all(q, false, exp_set, threshold, rpm_threshold, initial_kmers, bphf, kmer2countmap, exp_id2readcount);
				}
				print_exp_set(exp_set, exp_id2name, exp_id2desc, json);
			}

		}
		else { // single k-mer query
			run_query(arg_query, json, threshold, rpm_threshold, initial_kmers, bphf, kmer2countmap, exp_id2name, exp_id2desc, exp_id2readcount);
		}
	}
	else if(filename_query.length() > 0) {

		std::set<ExperimentId> exp_set;
		bool first = true;

		std::ifstream filestream_kmers;
		filestream_kmers.open(filename_query);
		if(!filestream_kmers.is_open()) { error("Could not open file " + filename_query); exit(EXIT_FAILURE); }

		std::cerr << getCurrentTime() << " Start reading k-mers from file " << filename_query << "\n";
		std::string line_from_file;
		line_from_file.reserve(KMER_K + 1);
		// query k-mers to DB:
		if(!all_kmers && json) { std::cout << "{ \"results\" : [ "; }
		while(getline(filestream_kmers,line_from_file)) {
			if(line_from_file.length() == 0) { continue; }
			if(!all_kmers) {
				if(json) { if(!first) { std::cout << ", "; } else { first = false; } }
				run_query(line_from_file, json, threshold, rpm_threshold, initial_kmers, bphf, kmer2countmap, exp_id2name, exp_id2desc, exp_id2readcount);
			}
			else {
				if(!first) {
					if(!exp_set.empty()) {
						run_query_all(line_from_file, false, exp_set, threshold, rpm_threshold, initial_kmers, bphf, kmer2countmap, exp_id2readcount);
					}
					else {
						break;
					}
				}
				else {
					run_query_all(line_from_file, true, exp_set, threshold, rpm_threshold, initial_kmers, bphf, kmer2countmap, exp_id2readcount);
					first = false;
				}
			}
		}
		if(!all_kmers && json) { std::cout << " ] }\n"; }
		filestream_kmers.close();

		if(all_kmers) {
			//Experiments that have all k-mers from input above thresholds are now in exp_set
			print_exp_set(exp_set, exp_id2name, exp_id2desc, json);
		}


	}

	// finished search, cleanup

	for(KmerIndex i = 0; i < n_elem;i++) {
		if(kmer2countmap[i] != nullptr) {
			delete kmer2countmap[i];
		}
	}

	delete[] kmer2countmap;
	delete bphf;

	return 0;

}

void run_query(const std::string & query, bool json, uint32_t threshold, uint32_t rpm_threshold, const std::vector<Kmer> & initial_kmers, boophf_t * bphf, pCountMap * kmer2countmap, const ExpId2Name & exp_id2name, const ExpId2Desc & exp_id2desc, const ExpId2ReadCount & exp_id2readcount) {
		if(query.length() < KMER_K){ printf("Warning, query too short:%s\n",query.c_str()); return; }
		if(query.length() > KMER_K){ printf("Warning, query too long:%s\n",query.c_str()); return; }
		Kmer kmer = str_to_int(query);
		//std::cerr << getCurrentTime() << " Searching " << query << "\n";
		if(json) std::cout << "{ \"query\" : \"" <<query << "\", \"experiments\" : [ ";
		if(std::binary_search(initial_kmers.begin(), initial_kmers.end(), kmer)) {
			// kmer was found in initial set, then check if it has counts in database
			KmerIndex index = bphf->lookup(kmer);
			assert(index < bphf->nbKeys());
			if(kmer2countmap[index] != nullptr) {
				bool first = true;
				for(auto const & it : *kmer2countmap[index]) { // go through all experiments that have counts for this k-mer
					ExperimentId exp_id = it.first;
					KmerCount count = it.second;
					assert(exp_id2readcount.find(exp_id) != exp_id2readcount.end());
					double rpm = (double)count / (double)exp_id2readcount.at(exp_id) * 1e6;
					if(count > threshold && rpm > rpm_threshold) {
						assert(exp_id2name.count(exp_id)>0);
						std::string exp_desc = "NA";
						auto it_desc = exp_id2desc.find(exp_id);
						if(it_desc != exp_id2desc.end() && it_desc->second.length() > 0) exp_desc = it_desc->second; 
						if(json) {
							if(!first) { printf(", ");} else { first = false; }
							printf("{ \"name\": \"%s\", \"count\":%d, \"rpm\":%f, \"desc\":\"%s\" } ", exp_id2name.at(exp_id).c_str(), count, rpm, exp_desc.c_str());
						}
						else {
							printf("%s\t%s\t%d\t%f\t%s\n",query.c_str(), exp_id2name.at(exp_id).c_str(), count, rpm, exp_desc.c_str());
						}
					}
				}
			}
			else {
				std::cerr << "K-mer " << query << " was not found in any experiment.\n";
			}
		}
		else {
			std::cerr << "K-mer " << query << " was not found in the inital k-mer set.\n";
		}
		if(json) std::cout << "]}\n";
}

void run_query_all(const std::string & query, bool first, std::set<ExperimentId> & exp_set, uint32_t threshold, uint32_t rpm_threshold, const std::vector<Kmer> & initial_kmers, boophf_t * bphf, pCountMap * kmer2countmap, const ExpId2ReadCount & exp_id2readcount) {
		if(query.length() < KMER_K){ printf("Warning, query too short:%s\n",query.c_str()); return; }
		if(query.length() > KMER_K){ printf("Warning, query too long:%s\n",query.c_str()); return; }

		Kmer kmer = str_to_int(query);
		std::cerr << getCurrentTime() << " Searching " << query << "\n";
		if(std::binary_search(initial_kmers.begin(), initial_kmers.end(), kmer)) {
			// kmer was found in initial set, then check if it has counts in database
			KmerIndex index = bphf->lookup(kmer);
			assert(index < bphf->nbKeys());
			std::set<ExperimentId> curr_exp_set;
			if(kmer2countmap[index] != nullptr) {
				for(auto const & it : *kmer2countmap[index]) { // go through all experiments that have counts for this k-mer
					ExperimentId exp_id = it.first;
					KmerCount count = it.second;
					assert(exp_id2readcount.find(exp_id) != exp_id2readcount.end());
					double rpm = (double)count / (double)exp_id2readcount.at(exp_id) * 1e6;
					if(count > threshold && rpm > rpm_threshold) {
						// experiment exp_id is good for k-mer kmer
						curr_exp_set.emplace(exp_id);
					}
				}
			}
			else {
				std::cerr << "K-mer " << query << " was not found in any experiment.\n";
			}
			if(first) {
				exp_set = curr_exp_set;
			}
			else {
				//intersect old_exp_set with curr_exp_set
				std::set<ExperimentId> intersection;
				std::set_intersection(exp_set.begin(), exp_set.end(), curr_exp_set.begin(), curr_exp_set.end(), std::inserter(intersection,intersection.begin()));
				exp_set = intersection;
			}
		}
		else {
			std::cerr << "K-mer " << query << " was not found in the inital k-mer set.\n";
		}
}

void print_exp_set(std::set<ExperimentId> & exp_set, ExpId2Name & exp_id2name, ExpId2Desc & exp_id2desc, bool json) {

	if(json) {
		std::cout << "{  \"experiments\" : [ ";
		bool first = true;
		for(auto const it : exp_set) {
			assert(exp_id2name.count(it)>0);
			std::string exp_desc = "NA";
			auto it_desc = exp_id2desc.find(it);
			if(it_desc != exp_id2desc.end() && it_desc->second.length() > 0) exp_desc = it_desc->second;
			if(!first) { printf(", ");} else { first = false; }
			std::cout << "{ \"name\" : \"" << exp_id2name.at(it) << "\", \"desc\" : \"" << exp_desc << "\" }";
		}
		std::cout << "]}\n";
	}
	else {
		for(auto it : exp_set) {
			assert(exp_id2name.count(it)>0);
			std::string exp_desc = "NA";
			auto it_desc = exp_id2desc.find(it);
			if(it_desc != exp_id2desc.end() && it_desc->second.length() > 0) exp_desc = it_desc->second;
			std::cout << exp_id2name.at(it) << "\t" << exp_desc << "\n";
		}
	}

}

void usage_kquery() {
	print_usage_header();
	fprintf(stderr, "Usage:\n   kiq query -i <file> -k <file> [-q <file> | -Q KMER[,KMER]* ]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of index file\n");
	fprintf(stderr, "   -k FILENAME   Name of k-mer count database file\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -q FILENAME   Name of file with query k-mers\n");
	fprintf(stderr, "   -Q KMER(S)    Single query k-mer or comma-separated list of query k-mers\n");
	fprintf(stderr, "   -t INT        Read count threshold\n");
	fprintf(stderr, "   -r FLOAT      RPM threshold\n");
	fprintf(stderr, "   -a            Only output experiments that contain all query k-mers\n");
	fprintf(stderr, "   -j            Output in JSON format\n");
	fprintf(stderr, "   -v            Enable verbose output.\n");
	fprintf(stderr, "   -d            Enable debug output.\n");
	fprintf(stderr, "   -h            Print this help.\n");
	exit(EXIT_FAILURE);
}
