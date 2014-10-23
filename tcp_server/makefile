gpm_server: main.o tcp_server.o
	c++ -o gpm_server main.o tcp_server.o -O2 \
	-I /disk/scratch/lib/boost_1_55_0 

main.o: main.cpp tcp_server.h
	c++ -o main.o -c main.cpp -c -O2 \
	-I /disk/scratch/lib/boost_1_55_0 

tcp_server.o: tcp_server.cpp tcp_server.h
	c++ -o tcp_server.o -c tcp_server.cpp -O2 \
	-I /disk/scratch/lib/boost_1_55_0 


clean:
	rm -f main.o tcp_server.o \
	gpm_server

