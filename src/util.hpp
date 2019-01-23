#pragma once

#include <stdlib.h>
#include <string>
#include <iostream>
#include <set>
#include <time.h>
#include <map>
#include <fstream>

#include "BooPHF/BooPHF.h"
#include "version.hpp"

using hasher_t = boomphf::SingleHashFunctor<u_int64_t>;
using boophf_t = boomphf::mphf<u_int64_t, hasher_t>;

using ExperimentId = uint32_t;
using ExperimentCount = uint32_t;
using KmerCount = uint32_t;
using KmerIndex = uint64_t;
using Kmer = uint64_t;
using ReadCount = uint64_t;

using CountMap = std::map<ExperimentId, KmerCount>;
using pCountMap = CountMap*;
using ExpId2Name = std::map<ExperimentId, std::string>;
using ExpId2Desc = std::map<ExperimentId, std::string>;
using ExpName2Id = std::map<std::string, ExperimentId>;
using ExpId2ReadCount = std::map<ExperimentId, ReadCount>;

enum DNA_MAP {C, A, T, G};  // A=1, C=0, T=2, G=3

struct HeaderDbFile {
    uint8_t magic[4] = {'K','I','Q',0x0A}; // == KIQ\n
    uint32_t dbVer = 2;
};

struct HeaderDbKmers {
    uint64_t numKmer = 0;
};

struct HeaderDbMetadata {
    uint8_t label[8] = {'M','E','T','A','D','A','T','A'};
    uint64_t numExp = 0;
};



#define KMER_K 32

void error(const std::string e);

void strip(std::string &s);

bool isalpha(const char & c);

void print_usage_header();

std::string getCurrentTime();

Kmer str_to_int(const std::string & str);
std::string int_to_str(Kmer kmer);

void load_index(const std::string & filename_index,  boophf_t * bphf);

void read_database(const std::string & filename,
										std::vector<Kmer> & initial_kmers,
										pCountMap * kmer2countmap,
										boophf_t * bphf,
										bool append,
										ExpId2Name & exp_id2name,
										ExpId2Desc & exp_id2desc,
										ExpName2Id & exp_name2id,
										ExpId2ReadCount & exp_id2readcount);

void write_database(const std::string & filename,
										const std::vector<Kmer> & initial_kmers,
										pCountMap * kmer2countmap,
										boophf_t * bphf,
										const ExpId2Name & exp_id2name,
										const ExpId2Desc & exp_id2desc,
										const ExpId2ReadCount & exp_id2readcount);

void write_initial_database(const std::string & filename, const std::vector<Kmer> & initial_kmers);

ExperimentId get_next_experiment_id(const ExpId2Name & exp_id2name);

