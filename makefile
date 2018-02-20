# Commet lines
# Here we define compiler option, libraries and the target
FLAGS = -O2 -Wall

# Here we make the executable file
SRCS = MD.cpp init.cpp integrator.cpp output.cpp forcepotential.cpp
OBJS = $(subst .cpp,.o,$(SRCS))
all: md

# Whereas here we create the object file
md: $(OBJS) constants.h
	mpic++ ${FLAGS} -o md $(OBJS)
#	export OMP_NUM_THREADS=20
	mpirun -np 4 ./md

%.o: %.cpp constants.h 
	mpic++ $(FLAGS) -c $<

#	g++ ${FLAGS} -c init.cpp
#	g++ ${FLAGS} -c integrator.cpp
#	g++ ${FLAGS} -c MD.cpp

# Clean
clean:
	rm *.o ./md
