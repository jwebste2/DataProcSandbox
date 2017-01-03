
#
# Check that DELPHESPATH environment variable is set
# This should point to the directory containing Delphes-3.3.0, which is needed
# for processing Delphes output files.
#
ifndef DELPHESPATH
   $(error Please export DELPHESPATH=/path/to/Delphes-X.Y.Z before running make)
endif


CPP := g++
CPPFLAGS := -g -std=c++11 # -Wall

MYINCS := -I$(PWD) -I$(PWD)/$(notdir $(shell pwd)) 

ROOTLIBS := `root-config --libs` -lRooFitCore -lRooFit -lMinuit -lRooStats -lCore -lHistFactory
ROOTLIBFLAGS := `root-config --ldflags`
ROOTINCS := `root-config --cflags` 

DELPHESLIBS := -L$(DELPHESPATH) -lDelphes
DELPHESINCS := -I$(DELPHESPATH) # -I/users/jwebster/MG5/MG5_aMC_v2_3_3/Template/NLO/MCatNLO/include

BUILD_OBJ := $(CPP) $(CPPFLAGS) $(MYINCS) $(ROOTINCS) $(DELPHESINCS)
LINK_OBJ := $(CPP) -g $(ROOTLIBS) $(DELPHESLIBS) $(MYINCS) $(ROOTINCS) $(DELPHESINCS)

SOURCES := $(wildcard src/*.cpp)

# Don't automatically delete the object files after compiling
# .PRECIOUS : bin/%.o

all : $(patsubst src/%.cpp, %, $(SOURCES))

# Linking
% : bin/%.o
	@rm -f $@
	$(LINK_OBJ) $^ $(TTHBBLEPTONICOBJ) -o run/$@
	@rm -f bin/*.o

# Objects
bin/%.o : src/%.cpp
	@rm -f $@
	$(BUILD_OBJ) -c $< -o $@

clean : 
	rm -f bin/* run/* *~ */*~ *.o */*.gch log*



