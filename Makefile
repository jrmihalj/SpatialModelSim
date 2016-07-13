ravdebug:
	gcc -I/usr/local/include/gsl -O3 -c -g Joe_CostSim.c
	gcc -L/usr/local/lib/ Joe_CostSim.o -lgsl -lgslcblas -lm -o Joe_CostSim


clean:
	rm -f *.o
	rm -f *.exe
	rm -f FixMCInt
