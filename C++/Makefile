# Simple Makefile to build DAIDALUS library and example applications
NAME=DAIDALUS
MAJOR=2
MINOR=.0.4
VERSION=$(NAME)$(MAJOR)
RELEASE=$(NAME)v$(MAJOR)$(MINOR)

PRECISION ?= double
OBJDIR = build/$(PRECISION)/obj
LIBDIR = lib
BINDIR = .
ifeq ($(PRECISION),float)
LIBDIR = build/float/lib
BINDIR = build/float
PRECISION_FLAGS = -DDAIDALUS_SINGLE_PRECISION
PRECISION_WARNINGS = -Werror=double-promotion -Wfloat-conversion
else ifneq ($(PRECISION),double)
$(error PRECISION must be double or float)
endif

SRC = $(wildcard src/*.cpp)
OBJS = $(patsubst src/%.cpp,$(OBJDIR)/%.o,$(SRC))
PROGRAMS = DaidalusExample DaidalusAlerting DaidalusBatch DAAEcho
ARCHIVE = $(LIBDIR)/$(RELEASE).a
INCLUDEFLAGS = -Iinclude
CXXFLAGS = -Wall -O

all: lib examples

lib: $(LIBDIR)/$(VERSION).a

$(LIBDIR)/$(VERSION).a: $(ARCHIVE)
	ln -sf $(RELEASE).a $@

$(ARCHIVE): $(OBJS)
	@mkdir -p $(LIBDIR)
	$(AR) rcs $@ $(OBJS)

$(OBJDIR)/%.o: src/%.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(INCLUDEFLAGS) $(CXXFLAGS) $(PRECISION_FLAGS) $(PRECISION_WARNINGS) -MMD -MP -c $< -o $@

examples: $(addprefix $(BINDIR)/,$(PROGRAMS))

$(addprefix $(BINDIR)/,$(PROGRAMS)): $(BINDIR)/%: examples/%.cpp $(ARCHIVE)
	$(CXX) $(CPPFLAGS) $(INCLUDEFLAGS) $(CXXFLAGS) $(PRECISION_FLAGS) $(LDFLAGS) $< $(ARCHIVE) $(LDLIBS) -o $@

-include $(OBJS:.o=.d)

doc:
	doxygen 

clean:
	-rm -f DaidalusExample DaidalusAlerting DaidalusBatch src/*.o examples/*.o lib/*.a
	-rm -rf build/double build/float

check:
	cppcheck --enable=all --cppcheck-build-dir=.cppcheck-config --suppressions-list=.cppcheck-config/cppcheck-suppressions.txt $(INCLUDEFLAGS) -q $(SRC) examples/

.PHONY: all lib examples doc clean check
