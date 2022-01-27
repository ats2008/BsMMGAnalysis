
CXX 		=g++
LIBS		=`root-config --libs` -lTMVA 

CXXFLAGS 	=-g3 -std=c++11  `root-config --cflags` 
INCPATH     =.

TARGET		= mvaDataMaker
all:	$(TARGET) 

mvaTrainer2 : photonID_bdt_MLP_classification.C 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

mvaTester : testMVAModels.cc MVATrainer.h
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

mvaDataMaker : mvaDataMaker.cc Util.h TreeMaker.h TreeMaker.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

mvaTrainer: makeMVAModels.cc  MVATrainer.h MVATrainer.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

clean :
	@rm *.exe
