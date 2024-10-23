SYSTEM     = x86-64_linux
#SYSTEM     = x86-64_osx
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

#CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio129/cplex
#CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio129/concert
CPLEXDIR      = /Applications/CPLEX_Studio221/cplex
CONCERTDIR    = /Applications/CPLEX_Studio221/concert
GUROBIDIR = /Library/gurobi1000/macos_universal2/
FLINT2DIR = $(HOME)/Desktop/flint2
FLINT2FMPZDIR = $(HOME)/Desktop/flint2/fmpz
GMPDIR = /opt/homebrew/include

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

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
GUROBILIBDIR = $(GUROBIDIR)/lib/
FLINTLIBDIR = $(FLINT2DIR)/lib/
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

# For dynamic linking
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIB      = cplex$(dynamic:yes=1290)

#CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) $(dynamic:yes=-L$(CPLEXBINDIR))
CCLNDIRS  = -L$(GUROBILIBDIR)
CLNDIRS   = -L$(CPLEXLIBDIR) $(dynamic:yes=-L$(CPLEXBINDIR))
#CCLNFLAGS = -lconcert -lilocplex -l$(CPLEXLIB) -lm -lpthread -ldl -lflint
#CCLNFLAGS = -lm -lpthread -ldl
CLNFLAGS  = -l$(CPLEXLIB) -lm -lpthread -ldl

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
GUROBIINCDIR = $(GUROBIDIR)/include/

#CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I$(GUROBIINCDIR) -I$(FLINT2DIR) -I$(FLINT2FMPZDIR) -I$(GMPDIR)
CCFLAGS = $(CCOPT) -I$(GUROBIINCDIR) -I$(FLINT2DIR) -I$(FLINT2FMPZDIR) -I$(GMPDIR)

all: leblanc_solver

leblanc_solver: solve_leblanc.o selfish_routing.o fr_solver.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o leblanc_solver solve_leblanc.o selfish_routing.o fr_solver.o $(CCLNFLAGS)

solve_leblanc.o: solve_leblanc.cpp
	$(CCC) -c $(CCFLAGS) solve_leblanc.cpp

selfish_routing.o: selfish_routing.cpp selfish_routing.h
	$(CCC) -c $(CCFLAGS) selfish_routing.cpp  $(CCLNFLAGS)

fr_solver.o: fr_solver.cpp fr_solver.h
	$(CCC) -c $(CCFLAGS) fr_solver.cpp  $(CCLNFLAGS)

clean:
	rm -f *.o leblanc_solver
