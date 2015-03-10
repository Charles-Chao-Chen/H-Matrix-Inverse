CC = g++

#EIGEN   = /usr/include/eigen3/
EIGEN   = ~/Software/eigen-eigen-10219c95fe65/
INCLUDE = -I ../include -I $(EIGEN)
CCFLAGS = -Wall $(INCLUDE)

#HEADER  = ../include/hmap.hpp

EXE = solve
SRC = hmat.cpp node.cpp helper.cpp dim2.cpp
# add prefix
OBJ = $(join $(addsuffix ../src/, $(dir $(SRC))), $(notdir $(SRC:.cpp=.o)))
OBJ += solve.o

all: $(EXE)

# debugging mode
debug: CCFLAGS += -DDEBUG -g
debug: $(EXE)

# profiling mode
profile: CCFLAGS += -pg
profile: $(EXE)

%.o : %.cpp
	$(CC) $(CCFLAGS) -c $< -o $@

$(EXE) : $(OBJ)
	$(CC) $(CCFLAGS) $(OBJ) -o $@

clean :
	rm $(EXE) $(OBJ) *~ ../src/*~ ../include/*~

