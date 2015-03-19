CC = g++

#EIGEN   = /usr/include/eigen3/
# "eigen" on corn
#EIGEN   = ~/Software/eigen-eigen-10219c95fe65/
# "eigen" on sapling
EIGEN   = /home/cchen10/Software/eigen-eigen-10219c95fe65
INCLUDE = -I ../include -I ../include/iterSolver/ -I $(EIGEN)
CCFLAGS = -Wall $(INCLUDE)

#HEADER  = ../include/hmap.hpp

EXE = solve
SRC = hmat.cpp node.cpp helper.cpp dim2.cpp zperm.cpp timer.cpp
SRC_ITER = solverBase.cpp cg.cpp fixedPoint.cpp

# add prefix
OBJ = $(join $(addsuffix ../src/, $(dir $(SRC))), $(notdir $(SRC:.cpp=.o)))
# add prefix
OBJ_ITER = $(join $(addsuffix ../src/iterSolver/, $(dir $(SRC_ITER))),\
	$(notdir $(SRC_ITER:.cpp=.o)))
OBJ += $(OBJ_ITER) solve.o

# enable compiler optimization
all: CCFLAGS += -O2
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
	rm $(EXE) $(OBJ) *~ 	\
	../src/*~     ../src/iterSolver/*~	\
	../include/*~ ../include/iterSolver/*~

