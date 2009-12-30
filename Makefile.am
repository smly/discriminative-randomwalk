# TESTS -- programs run automatically by "make check"
# check_PROGRAMS -- programs built by "make check" but not necessarily run
ACLOCAL_AMFLAGS = -I m4
TESTS =
check_PROGRAMS =

SUBDIRS = src

if DEVEL
SUBDIRS += tests
TESTS += tests/test_simple
endif

EXTRA_DIST = README.md dat LICENSE

memcheck:
	valgrind --tool=memcheck src/main
