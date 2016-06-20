CORE_SOURCE_FILES = $(top_srcdir)/src/mmio.cpp $(top_srcdir)/src/Matrix.cpp  $(top_srcdir)/src/CLI.cpp $(top_srcdir)/src/Utility.cc $(top_srcdir)/src/Timer.cpp $(top_srcdir)/src/PartitionLoader.cc $(top_srcdir)/src/IMatrix.cc $(top_srcdir)/src/CSegMatrix.cpp $(top_srcdir)/src/Configuration.cc $(top_srcdir)/src/NNZTag.cc $(top_srcdir)/src/SimplePartitionedMatrix.cc $(top_srcdir)/src/IRowColumnDS.cc $(top_srcdir)/src/ProfitPartitionMatrix.cc

# Include Directory for compilation
INCLUDES 	= -I$(top_srcdir)/src
#LIB = $(top_srcdir)/libs/libcsparse.a
#INCLUDE_HEADERS = mmio.h Matrix.h CLI.h

FORWARDDIFFERENCE_SOURCE_FILES = ForwardDifference.cc  StructuredJacobian.cc

## Programs to Build
bin_PROGRAMS = $(FORWARDDIFFERENCE) $(VERIFY) $(GCOLOR)

## Mappings of the executables. / Dictionary-Like
FORWARDDIFFERENCE = ForwardDifference
VERIFY = dsjmverify

ForwardDifference_SOURCES =$(FORWARDDIFFERENCE_SOURCE_FILES) $(CORE_SOURCE_FILES)
ForwardDifference_CXXFLAGS = $(CPP_PROFILING_FLAGS)
#ForwardDifference_LDADD = $(LIB)

dsjmverify_SOURCES = dsjmverify.cpp

GCOLOR = gcolor
gcolor_SOURCES = gcolor.cpp driver.cc $(CORE_SOURCE_FILES)
gcolor_CXXFLAGS = $(CPP_PROFILING_FLAGS)
