erdosSum: erdosSum.cpp
	g++ -o erdosSum main.cpp numberTheory.cpp erdosSum.cpp mordellSum.cpp -std=c++11 -O2 -lmpfr -lgmp -lgmpxx
