#include "util.hpp"
#include <algorithm>

void error(const std::string e) {
	std::cerr << "Error: " << e << std::endl << std::endl;
}

inline bool isalpha(const char & c) {
	return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}


void strip(std::string & s) {
		for(auto it = s.begin(); it!=s.end(); ++it) {
			if(!isalpha(*it)) {
				s.erase(it);
				it--;
			}
		}
}

void read_database(const std::string & filename,
										std::vector<Kmer> & initial_kmers,
										pCountMap * kmer2countmap,
										boophf_t * bphf,
										bool append,
										ExpId2Name & exp_id2name,
										ExpName2Id & exp_name2id,
										ExpId2ReadCount & exp_id2readcount) {

	std::cerr << getCurrentTime() << " Reading database file " << filename << "\n";
	std::ifstream ifs(filename, std::ios::in | std::ios::binary);
	if(!ifs) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }

	//read header
	struct HeaderDbFile h;
	ifs.read(reinterpret_cast<char*>(&h.kiq), sizeof(h.kiq));
	if(!ifs.good()) {  error("Error reading from file " + filename); exit(EXIT_FAILURE); }
	if(h.kiq[0] != 'K' || h.kiq[1] != 'I' || h.kiq[2] != 'Q') { error("Wrong file type detected in file " + filename); exit(EXIT_FAILURE); }
	ifs.read(reinterpret_cast<char*>(&h.dbVer), sizeof(h.dbVer));
	if(!ifs.good()) {  error("Error reading from file " + filename); exit(EXIT_FAILURE); }
	std::cerr << h.kiq[0] << h.kiq[1] << h.kiq[2] << " ver=" << (int)h.dbVer << "\n";

	//if(h.numKmer != bphf->nbKeys()) { error("Error: Mismatching number of k-mers in hash index and k-mer database " + filename); exit(EXIT_FAILURE); }

	// read k-mer section
	struct HeaderDbKmers k;
	ifs.read(reinterpret_cast<char*>(&k.numKmer), sizeof(k.numKmer));
	if(!ifs.good()) {  error("Error reading from file " + filename); exit(EXIT_FAILURE); }
	std::cerr << "numKmer=" << k.numKmer << "\n";

	for(uint64_t n = 1; n <= k.numKmer; n++) {
		Kmer kmer;
		ifs.read(reinterpret_cast<char*>(&kmer), sizeof(Kmer));
		if(!ifs.good()) {  error("Error reading from file " + filename); exit(EXIT_FAILURE); }
		KmerIndex index = bphf->lookup(kmer);
		std::cerr << "Read k-mer " << kmer << "="<< int_to_str(kmer)  <<" with index " << index << "\n";
		assert(index < bphf->nbKeys());
		initial_kmers.emplace_back(kmer);
		ExperimentCount num_exp = 0;
		ifs.read(reinterpret_cast<char*>(&num_exp), sizeof(ExperimentCount));
		if(!ifs.good()) {  error("Error reading from file " + filename); exit(EXIT_FAILURE); }
		std::cerr << "num exp=" << num_exp <<"\n";
		if(num_exp>0) {
			if(append) kmer2countmap[index] = new CountMap();
			for(int i=0; i < static_cast<int>(num_exp); i++) {
				ExperimentId exp_id = 0;
				ifs.read(reinterpret_cast<char*>(&exp_id), sizeof(ExperimentId));
				if(!ifs.good()) {  error("Error reading from file " + filename); exit(EXIT_FAILURE); }
				KmerCount count = 0;
				ifs.read(reinterpret_cast<char*>(&count), sizeof(KmerCount));
				if(!ifs.good()) {  error("Error reading from file " + filename); exit(EXIT_FAILURE); }
				std::cerr << " expid=" << exp_id << " count=" << count << "\n";
				if(append) kmer2countmap[index]->emplace(exp_id,count);
			}
		}
	}

	// read metadata section
	struct HeaderDbMetadata m;

	ifs.read(reinterpret_cast<char*>(&m.label), sizeof(m.label));
	if(!ifs.good()) {  error("Error reading from file " + filename); exit(EXIT_FAILURE); }
	if(m.label[0] != 'M' || m.label[1] != 'E' || m.label[7] != 'A') { error("Error reading from file " + filename); exit(EXIT_FAILURE); }

	ifs.read(reinterpret_cast<char*>(&m.numExp), sizeof(m.numExp));
	if(!ifs.good()) {  error("Error reading from file " + filename); exit(EXIT_FAILURE); }
	std::cerr << "numMetaExp=" << m.numExp << "\n";

	// continue reading metadata, add field for sample description 

}

void print_usage_header() {
	fprintf(stderr, "KIQ %s\n",KIQ_VERSION_STRING);
	fprintf(stderr, "Copyright 2018,2019 Peter Menzel\n");
	fprintf(stderr, "\n");
}

std::string getCurrentTime() {
	time_t t = time(0);
	char buffer[9] = {0};
	strftime(buffer, 9, "%H:%M:%S", localtime(&t));
	return std::string(buffer);
}

/**
* Converts a string of "ATCG" to a uint64_t
* where each character is represented by using only two bits
* from https://github.com/splatlab/squeakr/blob/master/kmer.h
*/
Kmer str_to_int(const std::string & str) {
	uint64_t strint = 0;
	for(auto it : str) {
		strint <<= 2;
		uint8_t curr = DNA_MAP::A;
		switch(it) {
			//case 'A': { curr = DNA_MAP::A; break; }
			case 'T': { curr = DNA_MAP::T; break; }
			case 't': { curr = DNA_MAP::T; break; }
			case 'C': { curr = DNA_MAP::C; break; }
			case 'c': { curr = DNA_MAP::C; break; }
			case 'G': { curr = DNA_MAP::G; break; }
			case 'g': { curr = DNA_MAP::G; break; }
			// TODO right now all non-CGT chars are used as A
		}
		strint = strint | curr;
	}
	return (Kmer)strint;
}

/**
* Converts a uint64_t to a string of "ACTG"
* where each character is represented by using only two bits
* from https://github.com/splatlab/squeakr/blob/master/kmer.h
*/
std::string int_to_str(const Kmer kmer) {
	uint8_t base;
	std::string str;
	for(int i=KMER_K; i>0; i--) {
		base = (kmer >> (i*2-2)) & 3ULL;
		char chr = 'A';
		switch(base) {
			//	case DNA_MAP::A: { return 'A'; }
				case DNA_MAP::T: { chr = 'T'; break; }
				case DNA_MAP::C: { chr = 'C'; break; }
				case DNA_MAP::G: { chr = 'G'; break; }
		}
		str.push_back(chr);
	}
	return str;
}


void read_kmer_database(const std::string & filename, std::vector<Kmer> & initial_kmers, pCountMap * kmer2countmap, boophf_t * bphf, bool append) {
	std::cerr << getCurrentTime() << " Reading database from file " << filename << "\n";
	std::ifstream ifs(filename, std::ios::binary);
	if(!ifs.is_open()) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }
	uint64_t count = 0;
	while(ifs.peek()!=EOF) {
		Kmer kmer;
		ifs.read(reinterpret_cast<char*>(&kmer), sizeof(Kmer));
		KmerIndex index = bphf->lookup(kmer);
		//std::cerr << "Read k-mer " << kmer << "="<< int_to_str(kmer)  <<" with index " << index << "\n";
		assert(index < bphf->nbKeys());
		initial_kmers.emplace_back(kmer);
		ExperimentCount num_exp = 0;
		//std::cerr << "kmer index=" << index <<"\n";
		ifs.read(reinterpret_cast<char*>(&num_exp), sizeof(ExperimentCount));
		//std::cerr << "num exp=" << num_exp <<"\n";
		if(num_exp>0) {
			if(append) kmer2countmap[index] = new CountMap();
			for(int i=0; i < static_cast<int>(num_exp); i++) {
				ExperimentId exp_id = 0;
				ifs.read(reinterpret_cast<char*>(&exp_id), sizeof(ExperimentId));
				KmerCount count = 0;
				ifs.read(reinterpret_cast<char*>(&count), sizeof(KmerCount));
				//std::cerr << "kmer_index i=" << index << " expid=" << exp_id << " count=" << count << "\n";
				if(append) kmer2countmap[index]->emplace(exp_id,count);
			}
		}
		++count;
	}
	ifs.close();
	assert(count == bphf->nbKeys());
}

void write_initial_kmer_database(const std::string & filename, const std::vector<Kmer> & initial_kmers) {
	std::cerr << getCurrentTime() << " Writing k-mer database to file " << filename << "\n";
	std::ofstream os(filename, std::ios::out | std::ios::binary);
	if(!os.is_open()) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }
	KmerCount zero = 0;
	for(Kmer it : initial_kmers) {
		// write k-mer
		Kmer kmer = it;
		os.write(reinterpret_cast<const char *>(&kmer),sizeof(Kmer));
		// write 0
		os.write(reinterpret_cast<const char *>(&zero),sizeof(KmerCount));
	}
	os.close();
}

void write_kmer_database(const std::string & filename,  boophf_t * bphf,  pCountMap * kmer2countmap, const std::vector<Kmer> & initial_kmers) {

	std::cerr << getCurrentTime() << " Writing k-mer database to file " << filename << "\n";
	std::ofstream os(filename, std::ios::out | std::ios::binary);
	if(!os.is_open()) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }
	for(Kmer it : initial_kmers) {
		KmerIndex index = bphf->lookup(it);
		assert(index < bphf->nbKeys());
		Kmer kmer = it;
		// write k-mer
		os.write(reinterpret_cast<const char *>(&kmer),sizeof(Kmer));
		// write number of experiments having this k-mer
		ExperimentCount num_exp = kmer2countmap[index]==nullptr ? 0 : static_cast<ExperimentCount>(kmer2countmap[index]->size());
		os.write(reinterpret_cast<const char *>(&num_exp),sizeof(ExperimentCount));

		if(num_exp > 0) {
			for(auto const & it : *kmer2countmap[index]) {
				ExperimentId exp_id = it.first;
				KmerCount count = kmer2countmap[index]->at(it.first);
				//std::cerr << "kmer_index i=" << i << " expid=" << exp_id << " count=" << count << "\n";
				os.write(reinterpret_cast<const char *>(&exp_id),sizeof(ExperimentId));
				os.write(reinterpret_cast<const char *>(&count),sizeof(KmerCount));
			}
		}
	}
	os.close();
}


void read_experiment_database(const std::string & filename, ExpId2Name & exp_id2name, ExpName2Id & exp_name2id, ExpId2ReadCount & exp_id2readcount) {

	std::cerr << getCurrentTime() << " Reading metadata from file " << filename << "\n";
	std::ifstream ifs(filename, std::ios::binary);
	if(!ifs.is_open()) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }
	while(ifs.peek()!=EOF) {
		ExperimentId exp_id = 0;
		ifs.read(reinterpret_cast<char*>(&exp_id), sizeof(ExperimentId));
		assert(!ifs.eof());
		assert(exp_id > 0);
		ReadCount readcount = 0;
		ifs.read(reinterpret_cast<char*>(&readcount), sizeof(ReadCount));
		assert(!ifs.eof());
		std::string exp_name;
		getline(ifs, exp_name,'\0');
		assert(exp_name.length() > 0);
		exp_id2name.emplace(exp_id,exp_name);
		exp_id2readcount.emplace(exp_id,readcount);
		exp_name2id.emplace(exp_name,exp_id);
	}
	ifs.close();

}

void load_index(const std::string & filename_index,  boophf_t * bphf) {
	std::cerr << getCurrentTime() << " Reading index from file " << filename_index << "\n";
	std::ifstream ifs(filename_index, std::ios::binary);
	if(!ifs.is_open()) { std::cerr << "Cannot open file " << filename_index << std::endl; exit(EXIT_FAILURE); }
	bphf->load(ifs);
}

void write_experiment_database(const std::string & filename, ExpId2Name & exp_id2name, ExpId2ReadCount & exp_id2readcount) {

	std::cerr << getCurrentTime() << " Writing metadata to file " << filename << "\n";
	std::ofstream os(filename, std::ios::out | std::ios::binary);
	if(!os.is_open()) {  error("Could not open file " + filename); exit(EXIT_FAILURE); }
	for(auto const & it : exp_id2name) {
		ExperimentId exp_id = it.first;
		std::string exp_name = it.second;
		assert(exp_id2readcount.count(exp_id) > 0);
		ReadCount readcount = exp_id2readcount.at(exp_id);

		os.write(reinterpret_cast<const char *>(&exp_id),sizeof(ExperimentId));
		os.write(reinterpret_cast<const char *>(&readcount),sizeof(ReadCount));
		os.write(exp_name.c_str(),exp_name.length() + 1);
	}
	os.close();

}

