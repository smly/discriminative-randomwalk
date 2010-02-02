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
  unsigned int wl = 5;
  bool map_pred = false;
  bool symm = false;
  bool fixed_algorithm = false;

  int opt;
  extern char *optarg;
  while ((opt = getopt(argc, argv, "i:l:o:p:w:mshf")) != -1) {
    switch (opt) {
      case 'i': input_fn = std::string(optarg); break;
      case 'l': label_fn = std::string(optarg); break;
      case 'o': pref_fn  = std::string(optarg); break;
      case 'p': prec = atoi(optarg); break;
      case 'w': wl   = atoi(optarg); break;
      case 'm': map_pred = true; break;
      case 's': symm = true; break;
      case 'f': fixed_algorithm = true; break;
      case 'h':
        std::cout << argv[0] << HELP << std::endl;
        exit(EXIT_SUCCESS);
      default:
        std::cout << argv[0] << HELP << std::endl;
        exit(EXIT_FAILURE);
    }
  }

  gll::Dwalk d;
  d.load(input_fn.c_str(), label_fn.c_str(), symm);
  d.show_info();
  if (fixed_algorithm) {
    d.go3(wl, pref_fn, map_pred);
  } else {
    d.go(wl, pref_fn, map_pred);
  }
  return 0;
}
