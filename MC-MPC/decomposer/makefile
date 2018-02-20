all:
	g++ -o decompose -I/usr/local/include MPC.cpp decomposition.cpp ../util/utils.cpp main.cpp

test:
	g++ -o decomposition_tests -I/usr/local/include ../util/utils.cpp MPC.cpp decomposition.cpp tests/main.cpp tests/test_decomposition.cpp tests/test_MPC.cpp