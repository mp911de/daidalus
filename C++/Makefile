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
DAIDALUS_FILE_IO ?= 1
ifeq ($(DAIDALUS_FILE_IO),0)
OBJDIR = build/$(PRECISION)/no-io/obj
LIBDIR = build/$(PRECISION)/no-io/lib
FILE_IO_FLAGS = -DDAIDALUS_FILE_IO=0
SRC := $(filter-out src/DAABatchProcessor.cpp src/DaidalusFileIO.cpp src/DaidalusFileWalker.cpp \
	src/DaidalusParametersFileIO.cpp src/SeparatedInput.cpp src/SeparatedOutput.cpp \
	src/SequenceReader.cpp src/StateReader.cpp,$(SRC))
PROGRAMS =
else ifneq ($(DAIDALUS_FILE_IO),1)
$(error DAIDALUS_FILE_IO must be 0 or 1)
endif
DAIDALUS_BATCH_FLOAT ?= 0
ifeq ($(DAIDALUS_BATCH_FLOAT),1)
ifneq ($(PRECISION),float)
$(error DAIDALUS_BATCH_FLOAT requires PRECISION=float)
endif
OBJDIR := $(OBJDIR)/batched
LIBDIR := $(LIBDIR)/batched
BINDIR := $(BINDIR)/batched
PRECISION_FLAGS += -DDAIDALUS_BATCH_FLOAT=1
endif
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
	$(CXX) $(CPPFLAGS) $(INCLUDEFLAGS) $(CXXFLAGS) $(FILE_IO_FLAGS) $(PRECISION_FLAGS) $(PRECISION_WARNINGS) -MMD -MP -c $< -o $@

examples: $(addprefix $(BINDIR)/,$(PROGRAMS))

$(addprefix $(BINDIR)/,$(PROGRAMS)): $(BINDIR)/%: examples/%.cpp $(ARCHIVE)
	@mkdir -p $(@D)
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
