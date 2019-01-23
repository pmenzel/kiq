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
#include <memory>

#include "ProducerConsumerQueue/ProducerConsumerQueue.hpp"
#include "zstr/zstr.hpp"
#include "BooPHF/BooPHF.h"
#include "util.hpp"
#include "ReadItem.hpp"
#include "CountThread.hpp"

void usage_kdb() {
	print_usage_header();
	fprintf(stderr, "Usage:\n   kiq db -i <file> -k <file> -l <file> [-a] [-z <int>]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i <file>   Name of index file\n");
	fprintf(stderr, "   -k <file>   Name of k-mer count database file\n");
	fprintf(stderr, "   -l <file>   Name of file with sample list in TSV format \n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -a          Append mode\n");
	fprintf(stderr, "   -z INT      Number of parallel threads for counting (default: 5)\n");
	fprintf(stderr, "   -v          Enable verbose output\n");
	fprintf(stderr, "   -d          Enable debug output.\n");
	exit(EXIT_FAILURE);
}

int main_kdb(int argc, char** argv) {

	size_t max_num_threads = 5;
	size_t curr_num_threads = 5;
	size_t max_queue_size = 9999;
	const float ma_alpha = 0.7f;
	bool append = false;
	bool debug = false;
	bool verbose = false;

	std::string filename_index;
	std::string filename_db;
	std::string filename_inputlist;

	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hdvai:k:l:z:")) != -1) {
		switch (c)  {
			case 'h':
				usage_kdb();
			case 'd':
				debug = true; break;
			case 'v':
				verbose = true; break;
			case 'a':
				append = true; break;
			case 'k':
				filename_db = optarg; break;
			case 'i':
				filename_index = optarg; break;
			case 'l':
				filename_inputlist = optarg; break;
			case 'z':
				max_num_threads = atoi(optarg);
				curr_num_threads = max_num_threads;
				break;
			default:
				usage_kdb();
		}
	}
	if(filename_index.length() == 0) { error("Please specify the name of the index file, using the -i option."); usage_kdb(); }
	if(filename_db.length() == 0) { error("Please specify the name of the database file, using the -k option."); usage_kdb(); }
	if(filename_inputlist.length() == 0) { error("Please specify the name of the sample list file, using the -l option."); usage_kdb(); }

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
		read_database(filename_db, initial_kmers, kmer2countmap, bphf, append, exp_id2name, exp_id2desc, exp_name2id, exp_id2readcount);
	}
	catch(std::runtime_error e) {
		std::cerr << "Error while reading database (" << e.what() << ")." << std::endl;
		exit(EXIT_FAILURE);
	}

	std::unordered_set<Kmer> initial_kmers_set;
	initial_kmers_set.insert(initial_kmers.cbegin(),initial_kmers.cend());

	std::atomic<KmerCount> * tmp_counts_atomic = new std::atomic<KmerCount>[n_elem];

	std::ifstream ifs_inputlist(filename_inputlist);
	if(!ifs_inputlist) { std::cerr << "Cannot open file " << filename_inputlist << std::endl; exit(EXIT_FAILURE); }

	std::string line;
	while(getline(ifs_inputlist, line)) {
		if(line.length() == 0) { continue; }
		size_t tab = line.find_first_of('\t');
		if(tab == std::string::npos || tab == 0) {
			continue;
		}
		std::string experiment_stringid = line.substr(0,tab);
		std::string filename_seq;
		std::string experiment_desc;
		size_t tab2 = line.find_first_of('\t',tab+1);
		if(tab2 == std::string::npos) { // no second tab found in line, then filename is from first tab to the line end
			filename_seq = line.substr(tab+1);
		}
		else { //second tab found in line, filename is between both tabs and description is from 2nd tab to line end
			filename_seq = line.substr(tab+1,tab2-tab-1);
			experiment_desc = line.substr(tab2+1);
		}

		ExperimentId experiment_numericid = 0;

		if(exp_name2id.count(experiment_stringid) > 0) {
			if(append) {
				std::cerr << "Experiment " << experiment_stringid << " is already in database. Skipping...\n" ;
				continue; // don't do that experiment again in append mode
			}
			experiment_numericid = exp_name2id.at(experiment_stringid);
		} else {
			experiment_numericid = get_next_experiment_id(exp_id2name);
			exp_id2name.emplace(experiment_numericid,experiment_stringid);
			exp_name2id.emplace(experiment_stringid,experiment_numericid);
			exp_id2desc.emplace(experiment_numericid,experiment_desc);
		}

		memset(tmp_counts_atomic,0,n_elem*sizeof(std::atomic<KmerCount>));

		auto queue = new queue_t(max_queue_size);
		std::deque<std::thread> threads;
		std::deque<std::unique_ptr<CountThread>> threadpointers;
		for(int i=0; i < static_cast<int>(curr_num_threads); i++) {
			std::unique_ptr<CountThread> p(new CountThread(queue,tmp_counts_atomic,bphf,&initial_kmers_set));
			threads.push_back(std::thread(&CountThread::doWork,p.get()));
			threadpointers.push_back(std::move(p));
		}

		std::ifstream ifstr(filename_seq);
		if(!ifstr.good()) {  error("Could not open file " + filename_seq); exit(EXIT_FAILURE); }
		zstr::istream in1_file(ifstr);
		if(!in1_file.good()) {  error("Could not open file " + filename_seq); exit(EXIT_FAILURE); }

		bool firstline_file1 = true;
		bool isFastQ_file1 = false;
		std::string line_from_file;
		line_from_file.reserve(2000);
		std::string sequence1;
		sequence1.reserve(2000);

		/* vars for progress indicator */
		unsigned char x[4] = { '/','-','\\','|'};
		size_t count=0;
		float ma = 1; // running average of queue size
		int i = 0;

		std::cerr << getCurrentTime() << " Start counting k-mers from file "<< filename_seq <<  "\n";

		while(getline(in1_file,line_from_file)) {
			if(line_from_file.length() == 0) { continue; }
			if(firstline_file1) {
				char fileTypeIdentifier = line_from_file[0];
				if(fileTypeIdentifier == '@') {
					isFastQ_file1 = true;
				}
				else if(fileTypeIdentifier != '>') {
					error("Auto-detection of file type for file " + filename_seq + " failed.");
					exit(EXIT_FAILURE);
				}
				firstline_file1 = false;
			}
			if(isFastQ_file1) {
				// read sequence line
				getline(in1_file,line_from_file);
				sequence1 = line_from_file;
				// skip + lin
				in1_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				// skip quality score line
				in1_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}
			else { //FASTA
				// read lines until next entry starts or file terminates
				sequence1.clear();
				while(!(in1_file.peek()=='>' || in1_file.peek()==EOF)) {
					getline(in1_file,line_from_file);
					sequence1.append(line_from_file);
				}
			} // end FASTA

			strip(sequence1); // remove non-alphabet chars
			if(sequence1.length() >= KMER_K) {
				queue->push(new ReadItem(sequence1));
			}
			if(count++ % 100000 == 0) {
				size_t curr_queue_size = queue->size();
				ma = ma * (1-ma_alpha) + static_cast<float>(curr_queue_size) * ma_alpha;
				fprintf(stderr,"\r%c %lu/%lu %4i ", x[i++], curr_num_threads, max_num_threads, static_cast<int>(ma));
				if(i>3) i=0;
			}
		} // end main loop around file1
		fprintf(stderr,"\r                              \r");
		fflush(stderr);

		// finish reading file

		queue->pushedLast();
		while(!threads.empty()) {
			threads.front().join();
			threads.pop_front();
			threadpointers.pop_front();
		}
		delete queue;

		std::cerr << getCurrentTime() << " Processed " << count << " sequences\n";

		if(exp_id2readcount.count(experiment_numericid) > 0) {
			exp_id2readcount.at(experiment_numericid) =  exp_id2readcount.at(experiment_numericid) + count;
		}
		else {
			exp_id2readcount.emplace(experiment_numericid,count);
		}

		//then finally go through all k-mers and add to big hash map
		for(KmerIndex i=0; i < n_elem; i++) {
			if(tmp_counts_atomic[i]==0) {
				continue;
			}
			if(kmer2countmap[i]==nullptr) { // experiment id has not been seen for this k-mer before
				kmer2countmap[i] = new CountMap();
				kmer2countmap[i]->emplace(experiment_numericid, tmp_counts_atomic[i]);
			}
			else {
				auto pos = kmer2countmap[i]->find(experiment_numericid);
				if(pos != kmer2countmap[i]->end()) {
					pos->second += tmp_counts_atomic[i];
				}
				else {
					kmer2countmap[i]->emplace(experiment_numericid, tmp_counts_atomic[i]);
				}
			}
		}

		// save database to file
		write_database(filename_db, initial_kmers, kmer2countmap, bphf, exp_id2name, exp_id2desc, exp_id2readcount);

	} // end while list of all experiments to read from files


	delete[] tmp_counts_atomic;
	for(KmerIndex i = 0; i < n_elem;i++) {
		if(kmer2countmap[i] != nullptr) {
			delete kmer2countmap[i];
		}
	}
	delete[] kmer2countmap;
	delete bphf;

	return 0;

}


