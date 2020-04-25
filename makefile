ROOTLIBS=$(shell root-config --libs) -lMinuit2

CXXFLAGS=-std=c++14 -Wall -Wextra -Wpedantic -Wshadow -O3
CXX=g++

INCS=-I$(shell root-config --incdir) -I./interface -I$(CHIB_CHIC_POLFW_DIR)/general/interface -I$(EIGEN_INCLUDE_DIR)

ADD_FLAGS=

BIN=./bin/
SRC=./src/
RESULTS=./results/

SRCS=$(wildcard $(SRC)*.cc)
PROG=$(SRCS:$(SRC)%.cc=%)
LIST=$(addprefix $(BIN), $(PROG))

DEPS=$(LIST:%=%.d)
-include $(DEPS)

all: setup $(LIST)

debug: CXXFLAGS+=-ggdb3 -O0
debug: CXXFLAGS:=$(filter-out -O3,$(CXXFLAGS))
debug: all

setup:
	mkdir -p $(BIN) ${RESULTS}

$(BIN)%: $(SRC)%.cc
	$(CXX) $(INCS) $(CXXFLAGS) -MMD -MP $< -o $@ $(ROOTLIBS) $(ADD_FLAGS)

clean:
	rm -rf $(BIN)

