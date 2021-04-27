CC=g++
INCDIR=../include/
SRC=../src/
MAIN_TEST=../main_test/
CPPFLAG=-I${INCDIR} -std=c++11 -llapacke -lblas

default:run

run:Test_1 Test_2 main_1 main_2
	./Test_1
	./Test_2

main_1:${MAIN_TEST}main_1.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^} -O3 -o ${@} ${CPPFLAG}

main_2:${MAIN_TEST}main_2.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^} -O3 -o ${@} ${CPPFLAG}

Test_1:${MAIN_TEST}Test_1.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^} -O3 -o ${@} ${CPPFLAG}

Test_2:${MAIN_TEST}Test_2.cpp ${INCDIR}*.h ${SRC}*.cpp
	g++ ${^} -O3 -o ${@} ${CPPFLAG}

.PHONY:clean

clean:
	-rm -f *.out *.m
