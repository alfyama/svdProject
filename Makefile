.POSIX:

CXX=g++
CXXFLAGS=-g -Wall -std=c++17 -fopenmp

SRC=src
SRCS=$(wildcard $(SRC)/*.cpp)

OBJS=$(patsubst $(SRC)/%.cpp, $(SRC)/%.o, $(SRCS))

TARGET=svd


all: $(TARGET)

debug: CXXFLAGS += -DDEBUG=1
debug: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

$(SRC)/%.o: $(SRC)/%.cpp # $(SRC)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SRC)/%.o: $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -r src/*.o svd
	cd report && find . ! -name 'main.tex' -type f -exec rm -f {} +
	cd ..
