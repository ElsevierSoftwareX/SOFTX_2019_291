CXX = g++ -std=c++11 -O0 -g


default: do_test run_test

do_test:
	$(CXX) -o mt_test main.cpp

run_test:
	./mt_test

clean:
	rm -f mt_test

