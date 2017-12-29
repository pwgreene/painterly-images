#
# Modify HALIDE_DIR to the appropriate path on your machine.
#
# Special instructions for Mac users
# ==================================
# You need to or must have installed libpng through Macports or Homebrew.
# Assuming that the installation succeeded, you should be able to run
#
# The brew command for installing libpng is
# brew install libpng
#
# libpng-config --I_opts
# libpng-config --L_opts
#
# Please add the output of the above commands to the following variables:
# PNG_INC
# PNG_LIB

MKDIR	:= mkdir -p
RM		:= rm -f
CP		:= cp -f
CXX		:= g++ -Wall -g3 -ggdb -std=c++11 -I.

# TODO: changeme to your halide path
HALIDE_DIR ?= $(HOME)/halide
HALIDE_LIB := $(HALIDE_DIR)/lib/libHalide.a

INC  := $(wildcard  *.h)
SRC  := $(wildcard  tutorial*.cpp)

# TODO: (OSX users) add the configuration flags for libpng
PNG_INC := `libpng-config --I_opts`
PNG_LIB := `libpng-config --L_opts`

CXXFLAGS := $(PNG_INC) -I$(HALIDE_DIR)/include/ -I$(HALIDE_DIR)/tools/ -I. -g -Wall
LDFLAGS  := $(PNG_LIB) -L$(HALIDE_DIR)/bin/     -lz -lpthread -ldl -lncurses -lpng

%: %.cpp $(HALIDE_LIB)
	$(CXX) $(CXXFLAGS) $< $(HALIDE_LIB) $(LDFLAGS) -o $@

painterly: src/a10_main.cpp src/a10.cpp src/basicImageManipulation.cpp src/filtering.cpp src/Image.cpp src/lodepng.cpp src/timing.cpp src/a10.h $(HALIDE_LIB)
	mkdir -p Output
	$(CXX) $(CXXFLAGS) src/a10_main.cpp src/a10.cpp src/basicImageManipulation.cpp src/filtering.cpp src/Image.cpp src/lodepng.cpp src/timing.cpp $(HALIDE_LIB) $(LDFLAGS) -o $@

run: painterly
	./painterly

clean:
	$(RM) Output/* painterly
	$(RM) -rf *.dSYM
	rm -rf Output
