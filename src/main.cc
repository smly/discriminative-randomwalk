#include "dwalk.h"

extern "C" {
#include <unistd.h>
}
#ifndef HAVE_GETOPT_H
#include <getopt.h>
#endif

#if defined HAVE_GETOPT_H && defined HAVE_GETOPT_LONG
#include <getopt.h>
#else
extern "C" {
#include "getopt.h"
}
#endif

#define PREC_DEFAULT 8

#define HELP " [-p prec] -i graph_filename "\
  "-l label_filename -o output_prefix"

int
main(int argc, char** argv)
{
  std::string input_fn;
  std::string label_fn;
  std::string pref_fn;
  unsigned int prec = PREC_DEFAULT;
  std::string result;

  int opt;
  extern char *optarg;
  while ((opt = getopt(argc, argv, "i:l:o:p:h")) != -1) {
    switch (opt) {
      case 'i': input_fn = std::string(optarg); break;
      case 'l': label_fn = std::string(optarg); break;
      case 'o': pref_fn  = std::string(optarg); break;
      case 'p': prec = atoi(optarg); break;
      case 'h':
        std::cout << argv[0] << HELP << std::endl;
        exit(EXIT_SUCCESS);
      default:
        std::cout << argv[0] << HELP << std::endl;
        exit(EXIT_FAILURE);
    }
  }

  Dwalk d;
  d.load(input_fn.c_str(), label_fn.c_str());
  d.show_info();
  d.go(4, pref_fn);
  d.go(10, pref_fn);
  return 0;
}
