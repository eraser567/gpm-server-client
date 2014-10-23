gpm_client: main.o tcp_client.o
	c++ -o gpm_client main.o tcp_client.o -O2

main.o: main.cpp tcp_client.h
	c++ -o main.o -c main.cpp -O2

tcp_client.o: tcp_client.cpp tcp_client.h
	c++ -o tcp_client.o -c tcp_client.cpp -O2

clean:
	rm -f main.o tcp_client.o \
	gpm_client

