// Copyright Vladimir Prus 2002-2004.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/* The simplest usage of the library.
 */

#include <boost/program_options.hpp>
#include "Fasta.h"
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

int main(int argc, char* argv[])
{
	int optind = 1;
    Fasta FA(argc, argv);
//	if (strcmp(argv[optind],"-fasta") == 0)
//	{
//		Fasta FA(argc, argv);
//	}
    
    return 0;
}
