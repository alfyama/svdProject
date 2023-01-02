.POSIX:

CXX=g++
CXXFLAGS=-g -Wall -std=c++17

SRC=src
SRCS=$(wildcard $(SRC)/*.cpp)

OBJS=$(patsubst $(SRC)/%.cpp, $(SRC)/%.o, $(SRCS))


TARGET=svd

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

$(SRC)/%.o: $(SRC)/%.cpp # $(SRC)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SRC)/%.o: $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@



install: $(TARGET)
	mkdir -p $(DESTDIR)$(PREFIX)/bin 
	mkdir -p $(DESTDIR)$(PREFIX)/share/man/man1 
	cp -f game $(DESTDIR)$(PREFIX)/bin
	gzip < $(TARGET).1 > $(DESTDIR)$(PREFIX)/share/man/man1/$(TARGET).1.gz

clean:
	rm -r src/*.o zkt
