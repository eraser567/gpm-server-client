#include <unistd.h>
#include <iostream>
#include <sys/socket.h>
#include <arpa/inet.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <ctype.h>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <stdint.h>
#include <inttypes.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

using namespace boost;
using namespace std;

// Note: BITMAPSIZE = BITMAPSIZE_INT * 4
#define BITMAPSIZE 20
#define BITMAPSIZE_INT 5

typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertex_iterator vertex_iter;
typedef graph_traits<Graph> GraphTraits;
typedef unsigned long SigType;

#define MAXSIZE 32768

class tcp_server
{
private:
        int socket_fd,accept_fd;
        sockaddr_in myserver;
        sockaddr_in remote_addr;

	map<int, vector<int> > g;
	map<int, int > label;
	set<int> labelset;
	map<int, vector<int> > inverted_index;
	map< int, SigType > sigs;

	int base;
	int lnum;

public:
	char *fn;
        tcp_server(int listen_port);
        int recv_msg();
	void init_graph_index(char *fn);
	void clean_up();
	int check(char *query_fn, int mode);
};
