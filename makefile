ROOTLIBS=$(shell root-config --libs) -lMinuit2

CXXFLAGS=-std=c++14 -Wall -Wextra -Wpedantic -Wshadow -ggdb3
CXX=g++

INCS=-I$(shell root-config --incdir) -I./interface

BIN=./bin/
SRC=./src/

SRCS=$(wildcard $(SRC)*.cc)
PROG=$(SRCS:$(SRC)%.cc=%)
LIST=$(addprefix $(BIN), $(PROG))

DEPS=$(LIST:%=%.d)
-include $(DEPS)

all: setup $(LIST)

setup:
	mkdir -p $(BIN)

$(BIN)%: $(SRC)%.cc
	$(CXX) $(INCS) $(CXXFLAGS) -MMD -MP $< -o $@ $(ROOTLIBS)

clean:
	rm -rf $(BIN)

