all: main.o
	g++ -std=c++11 -Wall -O3 main.o -o main

main.o: main.cpp
	g++ -std=c++11 -c main.cpp

run:
	./main <maxcut.in>maxcut.out
	
std:
	./main <maxcut.in

test1:
	./main <proj1_instances/unweighted_50.txt

test2:
	./main <proj1_instances/unweighted_100.txt

test3:
	./main <proj1_instances/unweighted_500.txt

test4:
	./main <proj1_instances/weighted_500.txt

test5:
	./main <proj1_instances/weighted_chimera_297.txt
	
clean:
	rm *.o main