erdosSum: ErdosSum.cpp MordellSum.cpp
	g++ -o erdosSum main.cpp numberTheory.cpp ErdosSum.cpp MordellSum.cpp -std=c++11 -O2 -lmpfr -lgmp -lgmpxx
