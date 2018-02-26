MED2 : GeneInfo.o GeneSeq.o MED_start.o MainProcess.o \
	OftenUsedOperatLib.o SequenceTransform.o med.o
	g++ -v -o MED2 GeneInfo.o GeneSeq.o MED_start.o MainProcess.o \
	OftenUsedOperatLib.o SequenceTransform.o med.o
	
GeneInfo.o : GeneInfo.cpp GeneInfo.h
	g++ GeneInfo.cpp -c

GeneSeq.o : GeneSeq.cpp GeneSeq.h
	g++ GeneSeq.cpp -c
	
MED_start.o : MED_start.cpp MED_start.h
	g++ MED_start.cpp -c
	
MainProcess.o : MainProcess.cpp
	g++ MainProcess.cpp -c

OftenUsedOperatLib.o : OftenUsedOperatLib.cpp OftenUsedOperatLib.h
	g++ OftenUsedOperatLib.cpp -c

SequenceTransform.o : SequenceTransform.cpp SequenceTransform.h
	g++ SequenceTransform.cpp -c

med.o : med.cpp
	g++ med.cpp -c

clean :
	rm -f MED2 GeneInfo.o GeneSeq.o MED_start.o MainProcess.o \
	OftenUsedOperatLib.o SequenceTransform.o med.o
