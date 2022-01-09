
CXX 		=g++
LIBS		=`root-config --libs` -lTMVA 

CXXFLAGS 	=-g3 -std=c++11  `root-config --cflags` 
INCPATH     =.

TARGET		= mvaDataMaker
all:	$(TARGET) 

mvaDataMaker : mvaDataMaker.cc Util.h TreeMaker.h TreeMaker.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

clean :
	@rm *.exe
