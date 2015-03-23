CC = g++

#============= Eigen path =============

# stanford corn machine
EIGEN   = ~/Software/eigen-eigen-10219c95fe65/

# cs sapling
ifeq ($(shell uname -n), sapling)
EIGEN   = /home/cchen10/Software/eigen-eigen-10219c95fe65
endif
#======================================

INCLUDE = -I ../include/ -I ../include/helper/ \
	  -I ../include/iterSolver/ -I $(EIGEN)
CCFLAGS = -Wall $(INCLUDE)

EXE = solve
SRC = hmat.cpp node.cpp utility.cpp dim2.cpp 
SRC_HELP = zperm.cpp timer.cpp
SRC_ITER = solverBase.cpp cg.cpp fixedPoint.cpp

# add prefix
OBJ      = $(join $(addsuffix ../src/, $(dir $(SRC))), \
	$(notdir $(SRC:.cpp=.o)))
OBJ_HELP = $(join $(addsuffix ../src/helper/, $(dir $(SRC_HELP))),\
	$(notdir $(SRC_HELP:.cpp=.o)))
OBJ_ITER = $(join $(addsuffix ../src/iterSolver/, $(dir $(SRC_ITER))),\
	$(notdir $(SRC_ITER:.cpp=.o)))
OBJ += $(OBJ_HELP) $(OBJ_ITER)

# enable compiler optimization
all: CCFLAGS += -O2
all: $(EXE)

# debugging mode
debug: CCFLAGS += -DDEBUG -g
debug: $(EXE)
	@echo "==================================="
	@echo "          debugging mode"
	@echo "==================================="

# profiling mode
profile: CCFLAGS += -pg
profile: $(EXE)
	@echo "==================================="
	@echo "          profiling mode"
	@echo "==================================="

# run tests
check: test
	./test

test : $(OBJ) test.o
	$(CC) $(CCFLAGS) $(OBJ) test.o -o $@ -lcppunit

%.o : %.cpp
	$(CC) $(CCFLAGS) -c $< -o $@

$(EXE) : $(OBJ) solve.o
	$(CC) $(CCFLAGS) $(OBJ) solve.o -o $@

clean :
	rm $(EXE) $(OBJ) *~ 	\
	../src/*~		../include/*~    	\
	../src/iterSolver/*~	../include/iterSolver/*~  \
	../include/helper/*~	../src/helper/*~

