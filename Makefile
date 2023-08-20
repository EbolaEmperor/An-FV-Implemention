all: test1 test2 test3

FV_MOL.o: include/FV_MOL.h src/FV_MOL.cpp include/idpair.h include/RKTable.h
	g++ -c src/FV_MOL.cpp -Iinclude -O2 -O3 -Ofast

INSE.o: include/INSE.h src/INSE.cpp include/idpair.h include/RKTable.h
	g++ -c src/INSE.cpp -Iinclude -O2 -O3 -Ofast

amgSolver.o: include/amgSolver.h src/amgSolver.cpp
	g++ -c src/amgSolver.cpp -Iinclude -O2 -O3 -Ofast

matrix.o: include/matrix.h src/matrix.cpp
	g++ -c src/matrix.cpp -Iinclude -O2 -O3 -Ofast

sparseMatrix.o: include/sparseMatrix.h src/sparseMatrix.cpp
	g++ -c src/sparseMatrix.cpp -Iinclude -O2 -O3 -Ofast

avl.o: include/avl.h src/avl.cpp
	g++ -c src/avl.cpp -Iinclude -O2 -O3 -Ofast

norm.o: include/norm.h src/norm.cpp
	g++ -c src/norm.cpp -Iinclude -O2 -O3 -Ofast

TimeFunction2D.o: include/TimeFunction2D.h src/TimeFunction2D.cpp
	g++ -c src/TimeFunction2D.cpp -Iinclude -O2 -O3 -Ofast

test1.o: test1.cpp
	g++ -c test1.cpp -Iinclude -O2 -O3 -Ofast

test1: test1.o FV_MOL.o amgSolver.o avl.o matrix.o sparseMatrix.o norm.o TimeFunction2D.o
	g++ test1.o FV_MOL.o amgSolver.o avl.o matrix.o sparseMatrix.o norm.o TimeFunction2D.o -o test1 -O2 -O3 -Ofast

test2.o: test2.cpp
	g++ -c test2.cpp -Iinclude -O2 -O3 -Ofast

test2: test2.o FV_MOL.o amgSolver.o avl.o matrix.o sparseMatrix.o norm.o TimeFunction2D.o
	g++ test2.o FV_MOL.o amgSolver.o avl.o matrix.o sparseMatrix.o norm.o TimeFunction2D.o -o test2 -O2 -O3 -Ofast

test3.o: test3.cpp
	g++ -c test3.cpp -Iinclude -O2 -O3 -Ofast

test3: test3.o INSE.o amgSolver.o avl.o matrix.o sparseMatrix.o norm.o TimeFunction2D.o
	g++ test3.o INSE.o amgSolver.o avl.o matrix.o sparseMatrix.o norm.o TimeFunction2D.o -o test3 -O2 -O3 -Ofast

test4: test3.o INSE.o amgSolver.o avl.o matrix.o sparseMatrix.o norm.o TimeFunction2D.o
	g++ test3.o INSE.o amgSolver.o avl.o matrix.o sparseMatrix.o norm.o TimeFunction2D.o -o test4 -O2 -O3 -Ofast

clean:
	rm *.o test1 test2 test3