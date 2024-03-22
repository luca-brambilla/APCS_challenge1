# Luca Brambilla
# 10510718 - 919812


CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -Wall -I. -Wno-conversion-null -Wno-deprecated-declarations


EXEC     = main
LDFLAGS ?= 
LDLIBS  ?= 


# to use muparser
# change PACS_ROOT!
PACS_ROOT		?= /home/luca/Documents/pacs-examples/Examples
LDFLAGS			+= -L $(PACS_ROOT)/lib -Wl,-rpath=$(PACS_ROOT)/lib
# LDLIBS		    += -lmuparser
# local not working by copying the object files... in ./lib folder
# LDFLAGS			+= -L./lib -Wl,-rpath=./lib
# LDLIBS 			+= -l muparser

CPPFLAGS += -I./include -I$(mkBoostInc)
LDFLAGS  += -L$(mkBoostLib)
LDLIBS	 += -l boost_iostreams -l boost_system -l boost_filesystem

all: $(EXEC)

# COMPILER
main.o: main.cpp optimization.hpp muparser_fun.hpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $<

# LINKER and execute
$(EXEC): main.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@
	./$(EXEC)

clean:
	$(RM) *.o
	$(RM) ./report/*.aux ./report/*.log ./report/*.synctex.gz

distclean: clean
	$(RM) $(EXEC)
	$(RM) *~
