#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <vector>
#include <queue>
#include <map>
#include <set>

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


struct VertexCompare {
public:
	map<int, int> graph1_label;
	map<int, int> graph2_label;

	int operator()(Vertex v1, Vertex v2)
	{
		
		return (graph1_label[(int)(v1)] == graph2_label[(int)(v2)]);
	}
};

int g_answer_count = 0;
map<int, int> *g_g_id_map;
map<int, int> *g_q_id_map;
set<int> outputnodes;
int id_outputnode;
// Default print_callback
template <typename Graph1, typename Graph2>
struct MyCallback {
        
    MyCallback(const Graph1& graph1, const Graph2& graph2) 
      : graph1_(graph1), graph2_(graph2) {}
    
    template <typename CorrespondenceMap1To2,
              typename CorrespondenceMap2To1>
    bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) const {
      
      // Print (sub)graph isomorphism map
      BGL_FORALL_VERTICES_T(v, graph1_, Graph1)
	{
		int qid = (*g_q_id_map)[get(vertex_index_t(), graph1_, v)];
		int gid = (*g_g_id_map)[get(vertex_index_t(), graph2_, get(f, v))];
		//printf("(%d, %d) ", qid, gid);
		if (qid == id_outputnode) {
			outputnodes.insert(gid);
		}
	}
      //printf("\n");
      
      g_answer_count ++;
	if (g_answer_count % 1000 == 0) {
		printf("%d answers are found, with %d distinct matchings for query node %d...\n", g_answer_count, outputnodes.size(), id_outputnode);
	}

      return true;
    }
    
  private:
    const Graph1& graph1_;
    const Graph2& graph2_;
};


int myhash (int value, int base) {
	return value%base;
}

void setbitmap (int hash, char *b) {
	b[hash/8] |= 1 << hash%8;
}

void printbyte(char *c) {
	for (int i=0; i<8; i++) {
		char t = 1 << (7-i);
		if (c[0] & t) printf("1");
		else printf("0");
	}
}

int main_testbitmap() {

	printf("%d\n", 1 <<5);

	char b[BITMAPSIZE];
	memset(b, 0, BITMAPSIZE*sizeof(char));
	int base = BITMAPSIZE*8;
	int value = 293;
	
	printf("myhash=%d\n", myhash(293, base));
	setbitmap(myhash(293, base), b);
	for (int i=0; i < BITMAPSIZE; i ++) {
		printf("%2d. ", i);
		printbyte(b+i);
		printf("\n");
	}
}

void setsig (int from, int from_l, int to, int to_l, int base, int lnum, map< int, char* > &sigs) {
	map< int, char* >::iterator it_from = sigs.find(from);
	map< int, char* >::iterator it_to = sigs.find(to);
	char *from_b, *to_b;
	if (it_from == sigs.end()) {
		from_b = new char[BITMAPSIZE];
		memset(from_b, 0, BITMAPSIZE*sizeof(char));
		sigs.insert(pair<int, char *>(from, from_b));
	}
	else {
		from_b = it_from->second;
	}
	if (it_to == sigs.end()) {
		to_b = new char[BITMAPSIZE];
		memset(to_b, 0, BITMAPSIZE*sizeof(char));
		sigs.insert(pair<int, char *>(to, to_b));
	}
	else {
		to_b = it_to->second;
	}
	
	//int hash1 = myhash(to_l, base);
	//int hash2 = myhash(base-from_l, base);
	//printf("from=%d(%d), to=%d(%d), base=%d, lnum=%d, ", from, from_l, to, to_l, base, lnum);
	//printf("hash1=%d, hash2=%d\n", hash1, hash2);
	setbitmap(myhash(lnum+to_l, base), from_b);
	setbitmap(myhash(lnum-from_l, base), to_b);
}

void add_to_graph (int &id, int &l, char *p, map<int, vector<int> > &g, map<int, int> &label, set<int> &labelset) {
	if (id == -1) {
		id = atoi(p);
		g.insert(pair<int, vector<int> >(id, vector<int>()));
		//printf("id=%d, ", id);
	}
	else if (l == -1) {
		l = atoi(p);
		labelset.insert(l);
		label[id] = l;
		//printf("label=%d: ", l);
	}
	else {
		int nei;
		nei = atoi(p);
		g[id].push_back(nei);
		//printf("%d, ", nei);
	}
}

char *trimwhitespace(char *str) {
	char *end;

	// Trim leading space
	while(isspace(*str)) str++;

	if(*str == 0)  // All spaces?
		return str;

	// Trim trailing space
	end = str + strlen(str) - 1;
	while(end > str && isspace(*end)) end--;

	// Write new null terminator
	*(end+1) = 0;

	return str;
}

void parse_line (char *buffer, map<int, vector<int> > &g, map<int, int> &label, set<int> &labelset) {
	int size = strlen(buffer);
	char *p = buffer;
	int id = -1;
	int l = -1;
	for (int i = 0; i < size; i ++) {
		if (buffer[i] == '\t') {
			add_to_graph(id, l, p, g, label, labelset);
			buffer[i] = 0;
			p = buffer + 1 + i;
		}
	}
	add_to_graph(id, l, p, g, label, labelset);
	//printf("\n");
}

// return: 0 for cycle, 1 for DAG
int toposort (map<int, vector<int> > &q, vector<int> &qrank, set<int> &roots) {
	map<int, vector<int> > parents;
	queue<int> que;
	map<int, int> indegree;

	// calculate parents
	map<int, vector<int> >::iterator itq = q.begin();
	map<int, vector<int> >::iterator itq_end = q.end();
	itq = q.begin();
	for (; itq != itq_end; itq ++) {
		parents.insert(pair<int, vector<int> >(itq->first, vector<int>()));
	}	
	itq = q.begin();
	for (; itq != itq_end; itq ++) {
		int from = itq->first;
		vector<int>::iterator itv = itq->second.begin();
		vector<int>::iterator itv_end = itq->second.end();
		for (; itv != itv_end; itv ++) {
			int to = *itv;
			parents[to].push_back(from);
		}
	}	
	
	// get all nodes with in-degree=0
	itq = q.begin();
	for (; itq != itq_end; itq ++) {
		int od = indegree[itq->first] = parents[itq->first].size();
		if (od == 0) {
			que.push(itq->first);
			roots.insert(itq->first);
			//printf("push into stack: %d\n", itq->first);
		}
	}
		
	int size = q.size();
	while (qrank.size() < size) {
		if (que.size() == 0) {
			printf("Q is not a DAG!\n");
			return 0;
		}
		
		int next = que.front();
		que.pop();

		qrank.push_back(next);

		vector<int>::iterator iter_v = q[next].begin();
		vector<int>::iterator iter_v_end = q[next].end();
		for (; iter_v != iter_v_end; iter_v ++) {
			int child_id = *iter_v;
			indegree[child_id] --;
			if (indegree[child_id] == 0) {
				que.push(child_id);
				//printf("push into stack: %d\n", stack.back());
			}
		}
	}
	
	//printf("topological ranking:\n");
	//for (int i = 0; i < qrank.size(); i ++) {
	//	printf("%d, ", qrank[i]);
	//}
	//printf("\n");

	return 1;
}

// return: -1 for invalid, 0 for cycle, 1 for DAG
int toposort_cycle (map<int, vector<int> > &q, vector<int> &qrank, set<int> &roots) {
	map<int, vector<int> > parents;
	queue<int> que;
	map<int, int> indegree;
	set<int> rest;

	int edge_count = 0;
	int ret = 1;

	// calculate parents
	map<int, vector<int> >::iterator itq = q.begin();
	map<int, vector<int> >::iterator itq_end = q.end();
	itq = q.begin();
	for (; itq != itq_end; itq ++) {
		parents.insert(pair<int, vector<int> >(itq->first, vector<int>()));
	}	
	itq = q.begin();
	for (; itq != itq_end; itq ++) {
		int from = itq->first;
		rest.insert(from);
		vector<int>::iterator itv = itq->second.begin();
		vector<int>::iterator itv_end = itq->second.end();
		for (; itv != itv_end; itv ++) {
			int to = *itv;
			parents[to].push_back(from);
			edge_count ++;
		}
	}	
	
	// get all nodes with in-degree=0
	itq = q.begin();
	for (; itq != itq_end; itq ++) {
		int od = indegree[itq->first] = parents[itq->first].size();
		if (od == 0) {
			que.push(itq->first);
			roots.insert(itq->first);
			//printf("push into stack: %d\n", itq->first);
		}
	}
		
	while (rest.size() > 0) {
		if (que.size() == 0 && rest.size() > 0) {
			printf("Q is not a DAG!\n");
			ret = 0;
			// choose a node with max out degree
			set<int>::iterator its = rest.begin();
			set<int>::iterator its_end = rest.end();
			int max = -1;
			int max_id = -1;
			for (; its != its_end; its ++) {
				//printf("out-degree of %d: %d, max=%d, max_id=%d\n", *its, q[*its].size(), max, max_id);
				if (q[*its].size() > max) {
					max = q[*its].size();
					max_id = *its;
				}
			}
			if (max_id > 0) {
				roots.insert(max_id);
				que.push(max_id);
			}
			else {
				max_id = *(rest.begin());
				roots.insert(max_id);
				que.push(max_id);
			}
			continue;
		}
		
		int next = que.front();
		que.pop();

		//printf("next=%d\n", next);

		qrank.push_back(next);
		rest.erase(next);
		
		vector<int>::iterator iter_v = q[next].begin();
		vector<int>::iterator iter_v_end = q[next].end();
		for (; iter_v != iter_v_end; iter_v ++) {
			int child_id = *iter_v;
			indegree[child_id] --;
			edge_count --;
			//printf("(%d, %d)\n", next, child_id);
			if (indegree[child_id] == 0 && rest.count(child_id) > 0) {
				que.push(child_id);
				//printf("add to queue: %d\n", child_id);
			}
		}
	}
	
	//printf("topological ranking:\n");
	//for (int i = 0; i < qrank.size(); i ++) {
	//	printf("%d, ", qrank[i]);
	//}
	//printf("\n");
	if (edge_count == 0) {
		return ret;
	}
	else {
		printf("Query has cycles, and there are %d edges cannot be covered by G_small.\n", edge_count);
		return -1;
	}
}

bool sigmatch (int *qs, int *gs) {
	for (int i=0; i < BITMAPSIZE_INT; i ++) {
		int check = qs[i]&gs[i];
		//printf("qs[%d] = %d, gs[%d] = %d, qs[%d]&gs[%d] = %d\n", i, qs[i], i, gs[i], i, i, check);
		if (check != qs[i]) {
			return false;
		}
	}
	return true;
}

int find_cand_in_g (int ql, char *qs, set<int> &cand, map<int, vector<int> > &g, map<int, int > &glabel, map< int, char* > &gsigs, map<int, vector<int> > &ivi, map<int, vector<int> > &gsmall) {
	//map<int, vector<int> >::iterator itg = g.begin();
	//map<int, vector<int> >::iterator itg_end = g.end();
	map<int, vector<int> >::iterator itivimap = ivi.find(ql);
	if (itivimap == ivi.end()) return 0;

	vector<int>::iterator itivi = itivimap->second.begin();
	vector<int>::iterator itivi_end = itivimap->second.end();
	map< int, char* >::iterator itsig_end = gsigs.end();
	int count = 0;
	//int no1 = 0, no2 = 0;
	//for (; itg != itg_end; itg ++) {
	for (; itivi != itivi_end; itivi ++) {
		//int gid = itg->first;
		//int gl = glabel[gid];
		int gid = *itivi;
		bool iscand = false;
		//if (ql == gl) {
			count ++;
			map< int, char* >::iterator itsig = gsigs.find(gid);
			if (itsig == itsig_end) {
				if (qs == NULL) {
					iscand = true;
				}
				//else no1 ++;
			}
			else {
				if (qs == NULL) {
					iscand = true;
				}
				else {
					char *gs = itsig->second;
					if (sigmatch((int *)qs, (int *)gs)) {
						iscand = true;
					}
					//else no2 ++;
					
					/*printf("signature:\n");	
					for (int i=0; i < BITMAPSIZE; i ++) {
						printf("%2d. ", i);
						printbyte(qs+i);
						printf("\t");
						printbyte(gs+i);
						printf("\n");
					}
					char fse;
					scanf("%c", &fse);*/
				}
			}
			if (iscand) {
				cand.insert(gid);
				if (gsmall.count(gid) == 0) gsmall.insert(pair<int, vector<int> >(gid, vector<int>())); 

			}
		//}
	}

	/*printf("signature:\n");	
	for (int i=0; i < BITMAPSIZE; i ++) {
		printf("%2d. ", i);
		printbyte(qs+i);
		printf("\n");
	}
	printf("same label(%d) count=%d, no1=%d, no2=%d\n", ql, count, no1, no2);*/
	//printf("%d out of %d", cand.size(), count);
	return cand.size();
}

int find_cand_in_neighors (set<int> &fromcand, int ql, char *qs, set<int> &cand, map<int, vector<int> > &g, map<int, int > &glabel, map< int, char* > &gsigs, map<int, vector<int> > &gsmall) {
	set<int>::iterator its = fromcand.begin();
	set<int>::iterator its_end = fromcand.end();
	map< int, char* >::iterator itsig_end = gsigs.end();
	int count = 0;
	//int no1=0, no2=0, no3=0;

	// for each candidate u of the "from node"	
	for (; its != its_end; its ++) {
		// find u's children, which satisfy ql and qs
		int uid = *its;
		vector<int>::iterator itv = g[uid].begin();
		vector<int>::iterator itv_end = g[uid].end();
		// for each child v
		for (; itv != itv_end; itv ++) {
			int gid = *itv;
			int gl = glabel[gid];
			bool iscand = false;
			if (ql == gl) {
				count ++;
				map< int, char* >::iterator itsig = gsigs.find(gid);
				if (itsig == itsig_end) {
					if (qs == NULL) {
						iscand = true;
					}
					//else no1 ++;
				}
				else {
					if (qs == NULL) {
						iscand = true;
					}
					else {
						char *gs = itsig->second;
						if (sigmatch((int *)qs, (int *)gs)) {
							iscand = true;
						}
						//else no2 ++;
					}
				}
				if (iscand) {
					cand.insert(gid);
					gsmall[uid].push_back(gid);
					if (gsmall.count(gid) == 0) gsmall.insert(pair<int, vector<int> >(gid, vector<int>())); 
				}
			}
		}
	}

	//printf("%d out of %d", cand.size(), count);
	//printf(", no1=%d, no2=%d", no1, no2);
	return cand.size();
}

int intersection (set<int> &small, set<int> &large, set<int> &result) {
	set<int>::iterator its = small.begin();
	set<int>::iterator its_end = small.end();
	set<int>::iterator large_end = large.end();
	set<int>::iterator found;
	for (; its != its_end; its ++) {
		found = large.find(*its);
		if (found != large_end) {
			result.insert(*its);
		}
	}
	return result.size();
}

void encode_newid (map<int, vector<int> > &g, map<int, vector<int> > &gnew, map<int, int > &id_new2old) {
	map<int, int> id_map;
	map<int, vector<int> >::iterator it = g.begin();
	map<int, vector<int> >::iterator it_end = g.end();
	int id = 0;
	for (; it != it_end; it ++) {
		id_map[it->first] = id;
		id_new2old[id] = it->first;
		gnew.insert(pair<int, vector<int> >(id, vector<int>()));
		id ++;
	}

	it = g.begin();
	for (; it != it_end; it ++) {
		id = id_map[it->first];
		vector<int>::iterator itv = it->second.begin();		
		vector<int>::iterator itv_end = it->second.end();		
		for (; itv != itv_end; itv ++) {
			gnew[id].push_back(id_map[*itv]);
		}
	}
}

void gen_graph_4_vf2 (map<int, vector<int> > &g, map<int, int> &glabel, map<int, int > &id_new2old, Graph &vf2graph, map<int, int> &newlabel) {
	map<int, int> id_map;
	map<int, vector<int> >::iterator it = g.begin();
	map<int, vector<int> >::iterator it_end = g.end();
	int id = 0;
	for (; it != it_end; it ++) {
		id_map[it->first] = id;
		id_new2old[id] = it->first;
		newlabel[id] = glabel[it->first];
		id ++;
	}

	it = g.begin();
	for (; it != it_end; it ++) {
		id = id_map[it->first];
		vector<int>::iterator itv = it->second.begin();		
		vector<int>::iterator itv_end = it->second.end();
		set<int> nei;
		for (; itv != itv_end; itv ++) {
			if (nei.count(*itv) == 0) {
				add_edge(id, id_map[*itv], vf2graph);
				nei.insert(*itv);
			}
		}
	}
}

void make_gsmall_smaller (map<int, vector<int> > &gsmall, map<int, set<int> > &cand) {
	set<int> candset;
	map<int, set<int> >::iterator itcand = cand.begin();
	map<int, set<int> >::iterator itcand_end = cand.end();
	for (; itcand != itcand_end; itcand ++) {
		candset.insert(itcand->second.begin(), itcand->second.end());
	}

	int r_edge_count = 0;
	int r_node_count = 0;

	map<int, vector<int> >::iterator itg = gsmall.begin();
	map<int, vector<int> >::iterator itg_end = gsmall.end();
	vector<int> remove;
	for (; itg != itg_end; itg ++) {
		if (candset.count(itg->first) == 0) {
			remove.push_back(itg->first);
			continue;
		}
		vector<int>::iterator itv = itg->second.begin();
		vector<int>::iterator itv_end = itg->second.end();
		vector<int> temp;
		for (; itv != itv_end; itv ++) {
			if (candset.count(*itv) > 0) {
				temp.push_back(*itv);
			}
			else r_edge_count ++;
		}
		itg->second.swap(temp);
	}

	vector<int>::iterator itv = remove.begin();
	vector<int>::iterator itv_end = remove.end();
	for (; itv != itv_end; itv ++) {
		gsmall.erase(*itv);
		r_node_count ++;
	}

	printf("make_gsmall_smaller: remove %d nodes and %d edges from g_small.\n", r_node_count, r_edge_count);
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		printf("argument is not enough!\n");
		printf("argument: 1.datagraph.filename\n");
		return -1;
	}
	
	map<int, vector<int> > g;
	map<int, int > label;
	set<int> labelset;

	// ******** load data graph ******** 
	char buffer[50000];
	FILE *fr = fopen(argv[1], "r");
	int line_count = 0;
	printf("\nload data graph...\n");
	while (fgets(buffer, 50000, fr) != NULL) {
		//printf("\n---------\n");
		char *p = buffer;
		p = trimwhitespace(p);
		if(strlen(buffer) == 0) continue;
		parse_line(buffer, g, label, labelset);
		if (line_count % 100000 == 0) {
			printf("load data graph progress = %d\n", line_count);
		}
		//if (line_count > 100000) break;
		line_count ++;
	}
	fclose(fr);

	// ******** build simple neighborhood signature ******** 
	printf("\nbuild signature...\n");
	int base = BITMAPSIZE*8;
	int lnum = labelset.size() + 3;
	line_count = 0;
	map< int, char* > sigs;
	map< int, vector<int> >::iterator itg = g.begin();
	map< int, vector<int> >::iterator itg_end = g.end();
	for (; itg != itg_end; itg ++) {
		int from = itg->first;
		int from_label = label[from];
		vector<int>::iterator itv = itg->second.begin();
		vector<int>::iterator itv_end = itg->second.end();
		for (; itv != itv_end; itv ++) {
			setsig(from, from_label, *itv, label[*itv], base, lnum, sigs);
		}

		if (line_count % 100000 == 0) {
			printf("build signature progress = %d\n", line_count);
		}
		//printf("%d\n", line_count);
		line_count ++;
	}

	// ******** build inverted index on labels for the datagraph ******** 
	printf("\nbuild inverted index on labels...\n");
	map<int, vector<int> > inverted_index;
	itg = g.begin();
	line_count = 0;
	for (; itg != itg_end; itg ++) {
		int id = itg->first;
		int l = label[id];
		
		if (inverted_index.count(l) == 0) {
			inverted_index.insert(pair<int, vector<int> >(l, vector<int>()));
		}
		inverted_index[l].push_back(id);

		if (line_count % 100000 == 0) {
			printf("build inverted index progress = %d\n", line_count);
		}
		//printf("%d\n", line_count);
		line_count ++;
	}

	while (true) {
		//int conti;
		//printf("continue? (-1 for EXIT, any other integers for CONTINUE.)\n");
		//scanf("%d", &conti);
		//if (conti == -1) break;
		
		// ******** input the query ******** 
		// ******** Note: the query should be (1) a DAG, and (2) connected. ******** 
		map<int, vector<int> > q;
		map<int, int > qlabel;
		set<int> qlabelset;
		map< int, char* > qsigs;
		int nq;
		printf("please 'paste' the query below (END with a line of a single integer -1, EXIT with -1):\n");	
		scanf("#%d\n", &nq);
		while (fgets(buffer, 50000, stdin) != NULL) {
			if(strlen(buffer) == 0 || buffer[0] == '\n' || buffer[0] == '#') continue;
			if(buffer[0] == '-' && buffer[1] == '1') break;
			//printf("LINE: %s\n", buffer);
			parse_line(buffer, q, qlabel, qlabelset);
		}
		if (q.size() == 0) break;

		printf("Please input the id of output node:\n");
		scanf("%d", &id_outputnode);

		bool has_answer = true;
		g_answer_count = 0;
		outputnodes.clear();

		// ******** build signature for the query ******** 
		map< int, vector<int> >::iterator itq = q.begin();
		map< int, vector<int> >::iterator itq_end = q.end();
		for (; itq != itq_end; itq ++) {
			int from = itq->first;
			int from_label = qlabel[from];
			//printf("%d(%d):", itq->first, qlabel[itq->first]);
			vector<int>::iterator itv = itq->second.begin();
			vector<int>::iterator itv_end = itq->second.end();
			for (; itv != itv_end; itv ++) {
				//printf("%d, ", *itv);
				setsig(from, from_label, *itv, qlabel[*itv], base, lnum, qsigs);
			}
			//printf("\n");
		}

		// ******** reduce the data graph to a smaller one, according to the query and signature ********
		map<int, vector<int> > gsmall;
		// 1. topo-sort for the query graph
		vector<int> qrank;
		set<int> roots;
		int isDAG = toposort_cycle(q, qrank, roots);
		if (isDAG == -1) {
			printf("The query edges are not covered!\n");
			continue;
		}
		/*printf("topo-sort: ");
		for (int i=0; i < qrank.size(); i ++) {
			printf("%d, ", qrank[i]);
		}
		printf("\n");*/

		map<int, set<int> > cand;
		vector<int>::iterator itv = qrank.begin();
		vector<int>::iterator itv_end = qrank.end();
		map< int, char* >::iterator itqsig_end = qsigs.end();
		for (; itv != itv_end; itv ++) {
			int qid = *itv;
			int ql = qlabel[qid];
			char *qs;
			map< int, char* >::iterator itqsig = qsigs.find(qid);
			if (itqsig == itqsig_end) qs = NULL;
			else qs = itqsig->second; 
		
			// 2. generate candidates for the roots	
			if (roots.count(qid) > 0) {
				set<int> ccand;
				//printf("candidates for [qid=%d], size=", qid);
				int ret = find_cand_in_g(ql, qs, ccand, g, label, sigs, inverted_index, gsmall);
				//printf("\n");
				if (ret == 0) {
					has_answer = false;
					break;
				}

				map<int, set<int> >::iterator itcand = cand.find(qid);
				if (itcand == cand.end()) {
					cand.insert(pair<int, set<int> >(qid, set<int>()));
					cand[qid].swap(ccand);
				}
				else {
					set<int> &old = itcand->second;
					set<int> result;
					if (ccand.size() >= old.size()) ret=intersection(old, ccand, result);
					else ret=intersection(ccand, old, result);
					
					if (ret == 0) {
						has_answer = false;
						break;
					}
					else itcand->second.swap(result);
				}
			}

			//printf("qid=%d\n", qid);

			// 3. generate candidates for non-root nodes
			vector<int>::iterator itc = q[qid].begin();
			vector<int>::iterator itc_end = q[qid].end();
			for (; itc != itc_end; itc ++) {
				int cid = *itc;
				int cl = qlabel[cid];
				char *cs;
				itqsig = qsigs.find(cid);
				if (itqsig == itqsig_end) cs = NULL;
				else cs = itqsig->second;
				
				//printf("q: %d --> %d\n", qid, cid);
				set<int> ccand;
 
				//printf("candidates for [cid=%d], size=", cid);
				int ret = find_cand_in_neighors(cand[qid], cl, cs, ccand, g, label, sigs, gsmall);
				//printf("\n");
				//printf("ccand.size=%d\n", ret);
				if (ret == 0) {
					has_answer = false;
					break;
				}


				map<int, set<int> >::iterator itcand = cand.find(cid);
				if (itcand == cand.end()) {
					cand.insert(pair<int, set<int> >(cid, set<int>()));
					cand[cid].swap(ccand);
				}
				else {
					set<int> &old = itcand->second;
					set<int> result;
					if (ccand.size() >= old.size()) ret=intersection(old, ccand, result);
					else ret=intersection(ccand, old, result);
					
					if (ret == 0) {
						has_answer = false;
						break;
					}
					else itcand->second.swap(result);
				}
			}
			if (!has_answer) {
				break;
			}
		}

		if (has_answer) {
			make_gsmall_smaller(gsmall, cand);
			// print candidates
			/*map<int, set<int> >::iterator itcand = cand.begin();
			map<int, set<int> >::iterator itcand_end = cand.end();
			for (; itcand != itcand_end; itcand ++) {
				printf("[%d] candidates (%d): ", itcand->first, itcand->second.size());
				set<int>::iterator itset = itcand->second.begin();
				set<int>::iterator itset_end = itcand->second.end();
				for (; itset != itset_end; itset ++) {
					printf("%d, ", *itset);
				}
				printf("\n");
			}*/
			// ******** generate a small data graph based on gsmall ******** 
			map<int, int > g_id_new2old;
			map<int, int > q_id_new2old;
			printf("G_small.size = %d\n", gsmall.size());
			Graph dgraph(gsmall.size());
			Graph qgraph(q.size());
			VertexCompare vcomp;
			gen_graph_4_vf2 (q, qlabel, q_id_new2old, qgraph, vcomp.graph1_label);
			gen_graph_4_vf2 (gsmall, label, g_id_new2old, dgraph, vcomp.graph2_label);
                        MyCallback<Graph, Graph> callback(qgraph, dgraph);

			g_q_id_map = &q_id_new2old;
			g_g_id_map = &g_id_new2old;

			// print gsmall, namely dgraph here
			/*printf("------G_small-----\n");
			GraphTraits::out_edge_iterator out_i, out_end;
      			GraphTraits::edge_descriptor e;
    			pair<vertex_iter, vertex_iter> vp;
    			for (vp = vertices(dgraph); vp.first != vp.second; ++vp.first) {
      				Vertex v = *vp.first;
				printf("[%d] <%d>: ", g_id_new2old[v], vcomp.graph2_label[v]);
				for (boost::tie(out_i, out_end) = out_edges(v, dgraph); out_i != out_end; ++out_i) {
					e = *out_i;
				        Vertex targ = target(e, dgraph);
				        printf("%d, ", g_id_new2old[targ]);
				}
				printf("\n");
    			}*/
			
			// print query
			/*printf("------Q-----\n");
    			for (vp = vertices(qgraph); vp.first != vp.second; ++vp.first) {
      				Vertex v = *vp.first;
				printf("[%d] <%d>: ", q_id_new2old[v], vcomp.graph1_label[v]);
				for (boost::tie(out_i, out_end) = out_edges(v, qgraph); out_i != out_end; ++out_i) {
					e = *out_i;
				        Vertex targ = target(e, qgraph);
				        printf("%d, ", q_id_new2old[targ]);
				}
				printf("\n\n");
    			}*/
			
			// ******** vf2 ******** 
			vf2_subgraph_iso(qgraph, 
				dgraph, 
				callback, 
				get(vertex_index, qgraph),
				get(vertex_index, dgraph),
				vertex_order_by_mult(qgraph),
				always_equivalent(), 
				vcomp
			);
			//printf("g_answer_count=%d\n", g_answer_count);
                        if (g_answer_count > 0) has_answer = true;
			else has_answer = false;
		}

		if (has_answer) {
			printf("HAS ANSWER!\n");
			printf("# of matchings: %d\n", g_answer_count);
			printf("# of matchings for query node %d: %d\n", id_outputnode, outputnodes.size());
		}
		else {
			printf("HAS NO ANSWER!\n");
		}

		// ******** release the memory for q ******** 
		map< int, char* >::iterator itsig = qsigs.begin();
		map< int, char* >::iterator itsig_end = qsigs.end();
		for (; itsig != itsig_end; itsig ++) {
			delete[] itsig->second;
		}

	}


	// ******** release the memory ******** 
	printf("\nrelease the memory...\n");
	map< int, char* >::iterator itsig = sigs.begin();
	map< int, char* >::iterator itsig_end = sigs.end();
	for (; itsig != itsig_end; itsig ++) {
		delete[] itsig->second;
	}

	return 0;
}
