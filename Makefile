.POSIX:

CXX=g++
CXXFLAGS=-g -Wall -fconcepts -std=c++17 #-fopenmp

SRC=src
SRCS=$(wildcard $(SRC)/*.cpp)

OBJS=$(patsubst $(SRC)/%.cpp, $(SRC)/%.o, $(SRCS))

TESTDIR=test

TARGET=svd


all: $(TARGET)

debug: CXXFLAGS += -DDEBUG=1
debug: $(TARGET)

tests:
	cd $(TESTDIR) && $(MAKE)


$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

$(SRC)/%.o: $(SRC)/%.cpp $(SRC)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SRC)/%.o: $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


.PHONY: all tests

clean:
	rm -r src/*.o results/* svd
	cd $(TESTDIR) && $(MAKE) clean
