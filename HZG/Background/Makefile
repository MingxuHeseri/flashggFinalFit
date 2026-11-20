# Run inside CMSSW, or give the include root and library directory of Combine.
COMBINE_INCLUDE_ROOT ?= $(CMSSW_BASE)/src
COMBINE_LIBDIR ?= $(CMSSW_BASE)/lib/$(SCRAM_ARCH)
ROOT_CONFIG ?= root-config
CPPFLAGS += -I$(COMBINE_INCLUDE_ROOT) -Iinterface
CXXFLAGS += -O2 -fPIC $(shell $(ROOT_CONFIG) --cflags)
LDLIBS += -L$(COMBINE_LIBDIR) -Wl,-rpath,$(COMBINE_LIBDIR) -lHiggsAnalysisCombinedLimit \
          $(shell $(ROOT_CONFIG) --glibs) -lMinuit2 -lMinuit -lRooFit -lRooFitCore -lMathMore -lRooStats \
          -lboost_program_options -lboost_regex
ifneq ($(CMSSW_BASE),)
BOOST_BASE := $(shell scram tool tag boost BASE)
CPPFLAGS += -I$(BOOST_BASE)/include
LDFLAGS += -L$(BOOST_BASE)/lib
endif
ifeq ($(shell uname),Darwin)
LIB := lib/libHZGBackground.dylib
SHARED := -dynamiclib -install_name @rpath/libHZGBackground.dylib
EXE_RPATH := -Wl,-rpath,@loader_path/../lib
else
LIB := lib/libHZGBackground.so
SHARED := -shared -Wl,-soname,libHZGBackground.so
EXE_RPATH := -Wl,-rpath,'$$ORIGIN/../lib'
endif
TDR := ../../tdrStyle
OBJECTS := $(patsubst src/%.cc,build/%.o,$(wildcard src/*.cc)) build/tdrstyle.o build/CMS_lumi.o
PROGRAMS := background_fit check_gof_integrals check_plot_normalization

.PHONY: all check
all: $(addprefix bin/,$(PROGRAMS))

build bin lib:
	mkdir -p $@
build/%.o: src/%.cc | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@
build/%.o: $(TDR)/%.C | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -include TPad.h -MMD -MP -c $< -o $@
build/test_%.o: test/%.cpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@
$(LIB): $(OBJECTS) | lib
	$(CXX) $(SHARED) $(LDFLAGS) $^ $(LDLIBS) -o $@
bin/%: build/test_%.o $(LIB) | bin
	$(CXX) $(LDFLAGS) $(EXE_RPATH) $< -Llib -lHZGBackground $(LDLIBS) -o $@
check: all
	bin/check_gof_integrals
	bin/check_plot_normalization

.SECONDARY:
-include $(wildcard build/*.d)
