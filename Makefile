# Build and test the C++ AKS programs.
#   make            build devVR/LenstraZnx.out and devVR/LenstraZ.out
#   make test       binary segmentation correctness test + AKS sweep (n <= 1000)
#   make clean
#
# macOS + Homebrew: NTL/GMP are found under $(brew --prefix) automatically.

CXX      ?= c++
CXXFLAGS ?= -g -O2 -std=c++17 -pthread
LDLIBS   = -lntl -lgmp -lm

BREW_PREFIX := $(shell brew --prefix 2>/dev/null)
ifneq ($(BREW_PREFIX),)
  CPPFLAGS += -I$(BREW_PREFIX)/include
  LDFLAGS  += -L$(BREW_PREFIX)/lib
endif

PROGRAMS = devVR/LenstraZnx.out devVR/LenstraZ.out
TESTS    = test/biSegCorrectnessTest.out test/probPrime.out
HEADERS  = $(wildcard devVR/*.h)

.PHONY: all test clean
all: $(PROGRAMS)

devVR/%.out: devVR/%.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

test/%.out: test/%.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

test: devVR/LenstraZnx.out $(TESTS)
	./test/biSegCorrectnessTest.out
	./test/aksSweep.sh 1000

clean:
	rm -f $(PROGRAMS) $(TESTS) test/*.out
