CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated -fpermissive

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

LIBS 		= $(ROOTGLIBS) $(WCSIMDIR)/libWCSimRoot.so -lMinuit

INC = $(WCSIMDIR)/include
SRC= $(WCSIMDIR)/src

CXXFLAGS += -I$(LIB) -I$(SRC) -I$(INC) 

TARGET= analysis_absorption

all: $(TARGET)
analysis_absorption: analysis_absorption.o


%: %.o
	@echo "Now make $@"
	@$(CPP) -o $@ $< $(CXXFLAGS) $(LIBS) 
	@echo "..Compile done! "

%.o: %.c
	@echo "$<"
	@echo "Start Compiling $<"
	@$(CPP) $(CXXFLAGS) -c $<
	@echo ".. Compiling Object Files $<   --> done"
	@echo ""

clean: 
	@echo "Now Clean Up"
	rm -f $(TARGET) *~ *.o *.o~ core
