SYSTEM     = x86-64_linux
#SYSTEM     = x86-64_osx
LIBFORMAT  = static_pic

# ---------------------------------------------------------------------
# Compiler selection
# ---------------------------------------------------------------------

#CCC = g++ -O3 -arch x86_64
CCC = g++ -O3
# CCC = g++ -pg -g

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------

CCOPT = -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
CCFLAGS = $(CCOPT)

all: leblanc_solver

leblanc_solver: solve_leblanc.o selfish_routing.o fr_solver.o
	$(CCC) $(CCFLAGS) -o leblanc_solver solve_leblanc.o selfish_routing.o fr_solver.o

solve_leblanc.o: src/solve_leblanc.cpp
	$(CCC) -c $(CCFLAGS) src/solve_leblanc.cpp

selfish_routing.o: src/selfish_routing.cpp src/selfish_routing.h
	$(CCC) -c $(CCFLAGS) src/selfish_routing.cpp

fr_solver.o: src/fr_solver.cpp src/fr_solver.h
	$(CCC) -c $(CCFLAGS) src/fr_solver.cpp

clean:
	rm -f *.o leblanc_solver
