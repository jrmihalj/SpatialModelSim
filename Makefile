ravdebug:
	gcc -I/usr/local/include/gsl -O3 -c -g SpatialModel.c
	gcc -L/usr/local/lib/ SpatialModel.o -lgsl -lgslcblas -lm -o SpatialModel


clean:
	rm -f *.o
	rm -f *.exe
	rm -f FixMCInt
