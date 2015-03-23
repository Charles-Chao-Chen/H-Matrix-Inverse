CC = g++

#EIGEN   = /usr/include/eigen3/
# "eigen" on corn
#EIGEN   = ~/Software/eigen-eigen-10219c95fe65/
# "eigen" on sapling
EIGEN   = /home/cchen10/Software/eigen-eigen-10219c95fe65
INCLUDE = -I ../include/ -I ../include/helper/ \
	  -I ../include/iterSolver/ -I $(EIGEN)
CCFLAGS = -Wall $(INCLUDE)

#HEADER  = ../include/hmap.hpp

EXE = solve
SRC = hmat.cpp node.cpp utility.cpp dim2.cpp 
SRC_HELP = zperm.cpp timer.cpp
SRC_ITER = solverBase.cpp cg.cpp fixedPoint.cpp

# add prefix
OBJ = $(join $(addsuffix ../src/, $(dir $(SRC))), $(notdir $(SRC:.cpp=.o)))
OBJ_HELP = $(join $(addsuffix ../src/helper/, $(dir $(SRC_HELP))),\
	$(notdir $(SRC_HELP:.cpp=.o)))
OBJ_ITER = $(join $(addsuffix ../src/iterSolver/, $(dir $(SRC_ITER))),\
	$(notdir $(SRC_ITER:.cpp=.o)))
OBJ += $(OBJ_HELP) $(OBJ_ITER) solve.o

# enable compiler optimization
all: CCFLAGS += -O2
all: $(EXE)

# debugging mode
debug: CCFLAGS += -DDEBUG -g
debug: $(EXE)

# profiling mode
profile: CCFLAGS += -pg
profile: $(EXE)

# run tests
check: test
	./test

%.o : %.cpp
	$(CC) $(CCFLAGS) -c $< -o $@

$(EXE) : $(OBJ)
	$(CC) $(CCFLAGS) $(OBJ) -o $@

test : $(OBJ) test.cpp
	$(CC) $(CCFLAGS) $(OBJ) test.cpp -o $@ -lcppunit

clean :
	rm $(EXE) $(OBJ) *~ 	\
	../src/*~		../include/*~    	\
	../src/iterSolver/*~	../include/iterSolver/*~  \
	../include/helper/*~	../src/helper/*~

