CC = g++

EIGEN   = /usr/include/eigen3/
INCLUDE = -I ../include -I $(EIGEN)
CCFLAGS = -Wall $(INCLUDE)

#HEADER  = ../include/hmap.hpp

EXE = solve

OBJS =  solve.o \
	../src/hmat.o ../src/node.o \
	../src/helper.o ../src/dim2.o

all: $(EXE)

# debugging mode
debug: CCFLAGS += -DDEBUG -g
debug: $(EXE)

%.o : %.cpp
	$(CC) $(CCFLAGS) -c $< -o $@

$(EXE) : $(OBJS)
	$(CC) $(CCFLAGS) $(OBJS) -o $@

clean :
	rm $(EXE) $(OBJS) *~ ../src/*~ ../include/*~

