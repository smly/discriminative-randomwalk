# for gtest
AC_DEFUN([GTEST_LIB_CHECK],
[
AC_ARG_WITH([gtest],
        [AS_HELP_STRING([--with-gtest],
                        [Specifies how to find the gtest package.])],
        [],
        [with_gtest=yes])
AC_ARG_ENABLE([external-gtest],
              [AS_HELP_STRING([--disable-external-gtest],
                              [Disables any detection or use of a system
                              installed or user provided gtest. Any option to
                              '--with-gtest' is ignored. (Default is enabled.)])
              ], [], [enable_external_gtest=yes])
AS_IF([test "x$with_gtest" == "xno"],
      [AC_MSG_ERROR([dnl
Support for GoogleTest was explicitly disabled. Currently GoogleMock has a hard
dependency upon GoogleTest to build, please provide a version, or allow
GoogleMock to use any installed version and fall back upon its internal
version.])])
AC_ARG_VAR([GTEST_CONFIG],
           [The exact path of Google Test's 'gtest-config' script.])
AC_ARG_VAR([GTEST_CPPFLAGS],
           [C-like preprocessor flags for Google Test.])
AC_ARG_VAR([GTEST_CXXFLAGS],
           [C++ compile flags for Google Test.])
AC_ARG_VAR([GTEST_LDFLAGS],
           [Linker path and option flags for Google Test.])
AC_ARG_VAR([GTEST_LIBS],
           [Library linking flags for Google Test.])
AC_ARG_VAR([GTEST_VERSION],
           [The version of Google Test available.])
HAVE_BUILT_GTEST="no"
GTEST_MIN_VERSION="1.4.0"

AS_IF([test "x${enable_external_gtest}" = "xyes"],
      [# Begin filling in variables as we are able.
      AS_IF([test "x${with_gtest}" != "xyes"],
            [AS_IF([test -x "${with_gtest}/scripts/gtest-config"],
                   [GTEST_CONFIG="${with_gtest}/scripts/gtest-config"],
                   [GTEST_CONFIG="${with_gtest}/bin/gtest-config"])
            AS_IF([test -x "${GTEST_CONFIG}"], [],
                  [AC_MSG_ERROR([dnl
Unable to locate either a built or installed Google Test at '${with_gtest}'.])
                  ])])

      AS_IF([test -x "${GTEST_CONFIG}"], [],
            [AC_PATH_PROG([GTEST_CONFIG], [gtest-config])])
      AS_IF([test -x "${GTEST_CONFIG}"],
            [AC_MSG_CHECKING([for Google Test version >= ${GTEST_MIN_VERSION}])
            AS_IF([${GTEST_CONFIG} --min-version=${GTEST_MIN_VERSION}],
                  [AC_MSG_RESULT([yes])
                  HAVE_BUILT_GTEST="yes"],
                  [AC_MSG_RESULT([no])])])])

AS_IF([test "x${HAVE_BUILT_GTEST}" = "xyes"],
      [GTEST_CPPFLAGS=`${GTEST_CONFIG} --cppflags`
      GTEST_CXXFLAGS=`${GTEST_CONFIG} --cxxflags`
      GTEST_LDFLAGS=`${GTEST_CONFIG} --ldflags`
      GTEST_LIBS=`${GTEST_CONFIG} --libs`
      GTEST_VERSION=`${GTEST_CONFIG} --version`],
      [AC_CONFIG_SUBDIRS([gtest])
      # GTEST_CONFIG needs to be executable both in a Makefile environmont and
      # in a shell script environment, so resolve an absolute path for it here.
      GTEST_CONFIG=`which gtest-config`
      GTEST_CPPFLAGS='-I$(top_srcdir)/gtest/include'
      GTEST_CXXFLAGS='-g'
      GTEST_LDFLAGS=''
      GTEST_LIBS='$(top_builddir)/gtest/lib/libgtest.la'
      GTEST_VERSION="${GTEST_MIN_VERSION}"])
])

