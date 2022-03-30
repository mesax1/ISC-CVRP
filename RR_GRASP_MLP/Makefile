appname := grasp_mlp

CXX := icpc
CXXFLAGS :=  -g -O3 -qopenmp -std=c++17 -lboost_system -Wall -Iinclude

srcfiles := $(shell find . -maxdepth 7 -name "*.cpp")
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) -o $(appname) $(objects)

clean:
	rm -f *.o $(objects)
