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
#include "kmodify.hpp"
#ifdef KIQ_SRA
#include "ksra.hpp"
#endif


void usage();

int main(int argc, char** argv) {

	if(argc < 2) {
		usage();
		return 1;
	}

	if(strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0) {
		usage();
		return 0;
	}

	int ret = 0;
	if(strcmp(argv[1], "index") == 0)
		ret = main_kindex(argc-1, argv+1);
	else if(strcmp(argv[1], "query") == 0)
		ret = main_kquery(argc-1, argv+1);
	else if(strcmp(argv[1], "db") == 0)
		ret = main_kdb(argc-1, argv+1);
#ifdef KIQ_SRA
	else if(strcmp(argv[1], "sra") == 0)
		ret = main_ksra(argc-1, argv+1);
#endif
	else if(strcmp(argv[1], "dump") == 0)
		ret = main_kdump(argc-1, argv+1);
	else if(strcmp(argv[1], "modify") == 0)
		ret = main_kmodify(argc-1, argv+1);
	else {
		ret = 1;
		usage();
	}

	return ret;
}

void usage() {
	print_usage_header();
	fprintf(stderr, "Usage:\n   kiq [ index | db | sra | query | dump | modify ] ...\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "     index    create index from initial list of k-mers\n");
	fprintf(stderr, "     db       create k-mer count database from FASTQ files\n");
	fprintf(stderr, "     sra      create k-mer count database from SRA files\n");
	fprintf(stderr, "     query    query k-mers against count database\n");
	fprintf(stderr, "     dump     print database content / stats\n");
	fprintf(stderr, "     modify   modify database content\n");

}
