
#ifndef DYNAMICTRIANGLELIST_H_
#define DYNAMICTRIANGLELIST_H_

#include<omp.h>
#include<string>
#include<vector>
#include<cstring>
#include<algorithm>
#include<cstdio>
#include<map>
#include<iostream>
#include<unordered_set>

using namespace std;

#define MAXN 2000000000

#define NONE 0
#define NEW_POS 1
#define OLD_POS 2
#define OLD_NEG 3

#define ALG_NAIVE 0
#define ALG_JOIN 1
#define ALG_AOT 2
#define ALG_AOT_DLT 3
#define ALG_NEW 4
#define ALG_HASH 5

#define MODE_INSERT 0
#define MODE_DELETE 1

class IntNode {
public:
	int val, id;
	bool operator < (const IntNode& v) const;
};

class Graph {
public:
	int n;
	long long m;
	vector<int> *con;
	vector<int> deg;
	vector<int> nid, oid;

	vector<int> *con_dlt;
	vector<int> deg_dlt;
	vector<int> list_v;

	vector<unordered_set<int> > hash;
	vector<unordered_set<int> > hash_dlt;

public:
	static void txt_to_bin(string path);
	static bool get_edge(char *line, int &a, int &b, int num_cnt);
	static int get_num_cnt(string path);
	static void get_order(vector<int> *con, int n, int *o);
	static int compare(int u, int du, int v, int dv);

public:
	Graph(string path);

	void run_batch(int mode, int n_thread, int alg, const vector<pair<int,int> > &l);
	void generate_edges(double percent, vector<pair<int,int> > &l);

	long long batch_insert_naive(const vector<pair<int,int> > &l);
	long long batch_delete_naive(const vector<pair<int,int> > &l);

	long long batch_insert_hash(const vector<pair<int,int> > &l);
	long long batch_delete_hash(const vector<pair<int,int> > &l);
	long long batch_hash(const vector<pair<int,int> > &l);
	void hash_join(unordered_set<int> *a, unordered_set<int> *b, vector<int> &v_join);

	long long batch_insert_join(const vector<pair<int,int> > &l);
	long long batch_delete_join(const vector<pair<int,int> > &l);
	long long batch_join(const vector<pair<int,int> > &l);

	long long batch_insert_new(const vector<pair<int,int> > &l);
	long long batch_delete_new(const vector<pair<int,int> > &l);
	long long batch_new(const vector<pair<int,int> > &l);

	long long batch_insert_aot(const vector<pair<int,int> > &l, bool is_dlt);
	long long batch_delete_aot(const vector<pair<int,int> > &l, bool is_dlt);
	void construct_graph_aot(vector<int> *con_new, vector<bool> *is_new, vector<int> &deg_new, vector<int> &deg_plus);
	void construct_graph_aot_dlt(vector<int> *con_new, vector<bool> *is_new, vector<int> &deg_new, vector<int> &deg_plus);
	long long batch_aot(vector<int> *con, vector<bool> *is_new, vector<int> &deg, vector<int> &deg_plus);

	void graph_batch_insert(bool use_hash = false);
	void graph_batch_delete(bool use_hash = false);

	void graph_single_insert(int u, int v);
	void graph_single_delete(int u, int v);

	void get_dlt(const vector<pair<int,int> > &l);
	void hash_l(const vector<pair<int,int> > &l);
	void clear_dlt();
	void print();
	int compare(int u, int v);

	void init_hash();

	~Graph();
};

bool IntNode::operator < (const IntNode& v) const {
	if(val == v.val) return id < v.id; else return val < v.val;
}

int Graph::get_num_cnt(string path) {
	FILE *fin = fopen( (path + "graph.txt").c_str(), "r" );
	char line[1000];
	int cnt = 0, min_cnt = 100;

	while( fgets( line, 1000, fin ) && cnt < 10 ) {
		if( !isdigit(line[0]) ) continue;
		vector<char*> v_num;
		int len = (int) strlen(line);
		for( int i = 0; i < len; ++i )
			if( !isdigit(line[i]) && line[i] != '.' ) line[i] = '\0';
			else if(i == 0 || !line[i-1]) v_num.push_back(line+i);
		if( (int) v_num.size() < 2 ) continue;
		min_cnt = min(min_cnt, (int)v_num.size());
		++cnt;
	}
	fclose( fin );
	return min_cnt;
}

bool Graph::get_edge(char *line, int &a, int &b, int num_cnt) {
	if( !isdigit(line[0]) ) return false;
	vector<char*> v_num;
	int len = (int) strlen(line);
	for( int i = 0; i < len; ++i )
		if( !isdigit(line[i]) && line[i] != '.' ) line[i] = '\0';
		else if(i == 0 || !line[i-1]) v_num.push_back(line+i);
	if( (int) v_num.size() != num_cnt ) return false;
	sscanf( v_num[0], "%d", &a );
	sscanf( v_num[1], "%d", &b );
	return true;
}

void Graph::get_order( vector<int> *con, int n, int *o) {
	IntNode *f = new IntNode[n];
	for( int i = 0; i < n; ++i )
		f[i].id = i, f[i].val = (int) con[i].size();
	sort(f, f + n);
	for(int i = 0; i < n; ++i) o[i] = f[n-1-i].id;
	delete[] f;
}

void Graph::txt_to_bin(string path) {
	FILE *fin = fopen( (path + "graph.txt").c_str(), "r" );
	char line[1024];
	int n = 0, a, b, num_cnt = get_num_cnt(path);
	long long cnt = 0, m = 0, p = 0;

	int *deg = new int[MAXN];
	memset(deg, 0, sizeof(int) * MAXN);
	printf( "Calculating size, num_cnt=%d...\n", num_cnt );
	while( fgets( line, 1024, fin ) ) {
			if( !get_edge(line, a, b, num_cnt) ) continue;
			if( a < 0 || b < 0 ) continue;
			n = max(n, a+1);
			n = max(n, b+1);
			if( a == b ) continue;
			++deg[a]; ++deg[b];
			m += 2;
			if( (++cnt) % (long long) 10000000 == 0 ) printf( "%lldM lines finished\n", cnt/1000000 );
		}
	fclose( fin );
	printf( "n=%d,m=%lld\n", n, m );

	int **con = new int*[n], *dat = new int[m];
	for( int i = 0; i < n; ++i ) {
		con[i] = dat + p;
		p += deg[i];
	}

	fin = fopen( (path + "graph.txt").c_str(), "r" );
	cnt = 0;
	memset(deg, 0, sizeof(int) * n);
	printf( "Loading text...\n" );
	while( fgets( line, 1024, fin ) ) {
		if( !get_edge(line, a, b, num_cnt) ) continue;
		if( a < 0 || b < 0 ) continue;
		if( a == b ) continue;
		con[a][deg[a]++] = b;
		con[b][deg[b]++] = a;
		if( (++cnt) % (long long) 10000000 == 0 ) printf( "%lldM lines finished\n", cnt/1000000 );
	}
	fclose( fin );


	int maxd = 0;
	m = 0;
	for( int i = 0; i < n; ++i )
		if( deg[i] > 0 ){
			sort( con[i], con[i]+deg[i] );
			int p = 0;
			for( int j = 0; j < deg[i]; ++j )
				if( (j == 0 || con[i][j-1] != con[i][j]) && con[i][j] != i ) con[i][p++] = con[i][j];
			deg[i] = p; m += p;
			maxd = max(maxd, p);
		}

	printf( "m=%lld, Reordering...\n", m );


	int *oid = new int[n], *nid = new int[n];
	//get_order(con, n, oid);

	IntNode *f = new IntNode[n];
	for( int i = 0; i < n; ++i ) f[i].id = i, f[i].val =  deg[i];
	sort(f, f + n);
	for(int i = 0; i < n; ++i) {oid[i] = f[n-1-i].id; deg[i] = f[n-1-i].val;}
	delete[] f;

	for( int i = 0; i < n; ++i ) nid[oid[i]] = i;


	printf( "Saving binary...\n" );

	FILE *fout = fopen( (path + "graph-sort.bin").c_str(), "wb" );
	fwrite( &n, sizeof(int), 1, fout );
	fwrite( &n, sizeof(int), 1, fout );
	fwrite( &m, sizeof(long long), 1, fout );
	fwrite( deg, sizeof(int), n, fout );

	int *nbr = new int[maxd];
	for( int i = 0; i < n; ++i ) {
		int u = oid[i], d = deg[i];
		for( int j = 0; j < d; ++j )
			nbr[j] = nid[con[u][j]];
		sort(nbr, nbr + d);
		fwrite(nbr, sizeof(int), d, fout);
	}

	fwrite( nid, sizeof(int), n, fout );
	fwrite( oid, sizeof(int), n, fout );
	fclose( fout );
	printf( "Created binary file, n = %d, m = %lld\n", n, m );
	delete[] dat;
	delete[] con; delete[] deg; delete[] oid; delete[] nid; delete[] nbr;
}

int Graph::compare(int u, int du, int v, int dv) {
	if(du > dv) return 1;
	if(du < dv) return -1;
	if(u > v) return 1;
	if(u < v) return -1;
	return 0;
}

Graph::Graph(string path) {
	printf( "path=%s\n", path.c_str() );
	FILE* fin = fopen( (path+"graph-sort.bin").c_str(), "rb" );
	fread( &n, sizeof(int), 1, fin );
	fread( &n, sizeof(int), 1, fin );
	fread( &m, sizeof(long long), 1, fin );
	int *dat = new int[m];
	deg.resize(n);
	deg_dlt.resize(n);
	nid.resize(n);
	oid.resize(n);
	con = new vector<int>[n];
	con_dlt = new vector<int>[n];

	printf( "Loading graph...\n" );
	fread( deg.data(), sizeof(int), n, fin );
	fread( dat, sizeof(int), m, fin );
	fread( nid.data(), sizeof(int), n, fin );
	fread( oid.data(), sizeof(int), n, fin );
	fclose(fin);

	long long p = 0;
	for( int i = 0; i < n; ++i ) {
		con[i].resize(deg[i]);
		memcpy(con[i].data(), dat+p, sizeof(int)*deg[i]);
		p+=deg[i];
	}
	int nown = 0;
	while(nown < n && deg[nown] > 1) ++nown;
	printf( "%s: n=%d,nown=%d,m=%lld\n", path.c_str(), n, nown, m );

	delete[] dat;
}

void Graph::init_hash() {
	hash.clear();
	hash.resize(n);
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < deg[i]; ++j)
			hash[i].insert(con[i][j]);
}

void Graph::print() {
	for(int i = 0; i < n; ++i) {
		printf( "%d:", i );
		for(int j = 0; j < deg[i]; ++j) printf( " %d", con[i][j] );
		printf( "\n" );
	}
}

void Graph::get_dlt(const vector<pair<int,int> > &l) {
	for(size_t i = 0; i < l.size(); ++i) {
		con_dlt[l[i].first].push_back(l[i].second);
		++deg_dlt[l[i].first];
		if(i == 0 || l[i].first != l[i-1].first) list_v.push_back(l[i].first);
	}
}

void Graph::hash_l(const vector<pair<int,int> > &l) {
	hash_dlt.clear();
	hash_dlt.resize(n);
	for(size_t i = 0; i < l.size(); ++i)
		hash_dlt[l[i].first].insert(l[i].second);
}


void Graph::clear_dlt() {
	for(size_t i = 0; i < list_v.size(); ++i) {
		con_dlt[list_v[i]].clear();
		deg_dlt[list_v[i]] = 0;
	}
	list_v.clear();
}

void Graph::generate_edges(double percent, vector<pair<int,int> > &l) {
	long long cnt = 0, t = 0;

	for(int u = 0; u < n; ++u)
		for(int i = 0; i < deg[u]; ++i) {
			int v = con[u][i];
			if(v <= u) continue;
			++cnt;
			if((t+1) <= percent * 0.01 * cnt) {
				++t;
				l.push_back(make_pair(u,v));
				l.push_back(make_pair(v,u));
			}
		}
	sort(l.begin(), l.end());
}

void Graph::run_batch(int mode, int n_thread, int alg, const vector<pair<int,int> > &l) {
	printf( "\n====Run Batch===\n" );
	if( n_thread > 0 ) omp_set_num_threads(n_thread);
	int np;
	#pragma omp parallel
	{
		if(omp_get_thread_num()==0) np = omp_get_num_threads();
	}
	printf( "Number of Threads = %d\n", alg==ALG_NAIVE?1:np );

	long long cnt = 0;
	double t = omp_get_wtime();
	if(mode == MODE_INSERT) {
		printf( "Batch Insertion, Algorithm: " );
		if(alg == ALG_NAIVE) { printf( "Naive.\n" ); cnt = batch_insert_naive(l); }
		else if(alg == ALG_JOIN ) { printf( "Join.\n" ); cnt = batch_insert_join(l); }
		else if( alg == ALG_AOT ) { printf( "AOT.\n"); cnt = batch_insert_aot(l, false); }
		else if( alg == ALG_AOT_DLT ) { printf( "AOT_DLT.\n"); cnt = batch_insert_aot(l, true); }
		else if( alg == ALG_HASH ) { printf( "Hash.\n"); cnt = batch_insert_hash(l); }
		else { printf( "New.\n" ); cnt = batch_insert_new(l); }
		printf( "%lld Triangles Inserted.\n", cnt );
	} else {
		printf( "Batch Deletion, Algorithm: " );
		if(alg == ALG_NAIVE) { printf( "Naive.\n" ); cnt = batch_delete_naive(l); }
		else if(alg == ALG_JOIN ) { printf( "Join.\n" ); cnt = batch_delete_join(l); }
		else if(alg == ALG_AOT ) { printf( "AOT.\n"); cnt = batch_delete_aot(l, false); }
		else if(alg == ALG_AOT_DLT ) { printf( "AOT_DLT.\n"); cnt = batch_delete_aot(l, true); }
		else if( alg == ALG_HASH ) { printf( "Hash.\n"); cnt = batch_delete_hash(l); }
		else { printf( "New.\n" ); cnt = batch_delete_new(l); }
		printf( "%lld Triangles Deleted.\n", cnt );
	}

	t = omp_get_wtime() - t;
	printf( "Total Processing Time = %0.5lf Seconds\n", t );
}

long long Graph::batch_insert_naive(const vector<pair<int,int> > &l) {
	long long cnt = 0;
	vector<int> v_join(n);
	for(size_t i = 0; i < l.size(); ++i) {
		int u = l[i].first, v = l[i].second;
		if(u >= v) continue;
		vector<int>::iterator it = set_intersection(con[u].begin(), con[u].begin() + deg[u], con[v].begin(), con[v].begin() + deg[v], v_join.begin());
		cnt += it - v_join.begin();
		graph_single_insert(u,v);
		graph_single_insert(v,u);
	}
	return cnt;
}

long long Graph::batch_delete_naive(const vector<pair<int,int> > &l) {
	long long cnt = 0;
	vector<int> v_join(n);
	for(size_t i = 0; i < l.size(); ++i) {
		int u = l[i].first, v = l[i].second;
		if(u >= v) continue;
		vector<int>::iterator it = set_intersection(con[u].begin(), con[u].begin() + deg[u], con[v].begin(), con[v].begin() + deg[v], v_join.begin());
		cnt += it - v_join.begin();
		graph_single_delete(u,v);
		graph_single_delete(v,u);
	}
	return cnt;
}

long long Graph::batch_insert_new(const vector<pair<int,int> > &l) {
	//get_dlt(l);
	long long cnt = batch_new(l);
	graph_batch_insert();
	//clear_dlt();
	return cnt;
}

long long Graph::batch_delete_new(const vector<pair<int,int> > &l) {
	//get_dlt(l);
	graph_batch_delete();
	long long cnt = batch_new(l);
	//clear_dlt();
	return cnt;
}

long long Graph::batch_new(const vector<pair<int,int> > &l) {
	long long cnt = 0;

	double t = omp_get_wtime();
	#pragma omp parallel
	{
		long long local_cnt = 0;
		vector<char> status(n,NONE);

		int np = omp_get_num_threads(), rank = omp_get_thread_num();

		//#pragma omp for schedule(dynamic)
		//for(size_t i = 0; i < list_v.size(); ++i) {
		for(size_t i = rank; i < list_v.size(); i+=np) {
			int u = list_v[i];
			for(int j = 0; j < deg[u]; ++j) status[con[u][j]] = compare(u,con[u][j]) > 0 ? OLD_POS : OLD_NEG;
			for(int j = 0; j < deg_dlt[u]; ++j) if(compare(u, con_dlt[u][j]) > 0 ) status[con_dlt[u][j]] = NEW_POS;

			for(int j = 0; j < deg_dlt[u]; ++j) {
				int v = con_dlt[u][j];
				if(compare(u,v) < 0) continue;
				for(int k = 0; k < deg[v]; ++k) {
					int w = con[v][k];
					if(status[w] == NONE) continue;
					if(status[w] != NEW_POS) {
						++local_cnt;
						//printf( "Type 1: (%d,%d,%d)\n", u, v, w );
					} else if(v > w) {
						++local_cnt;
						//printf( "Type 2.1: (%d,%d,%d)\n", u, v, w );
					}
				}
				for(int k = 0; k < deg_dlt[v]; ++k) {
					int w = con_dlt[v][k];
					if(status[w] == OLD_POS) {
						++local_cnt;
						//printf( "Type 2.2: (%d,%d,%d)\n", u, v, w );
					} else if(status[w] == NEW_POS && v > w) {
						++local_cnt;
						//printf( "Type 3: (%d,%d,%d)\n", u, v, w );
					}
				}
			}

			for(int j = 0; j < deg[u]; ++j) status[con[u][j]] = NONE;
			for(int j = 0; j < deg_dlt[u]; ++j) status[con_dlt[u][j]] = NONE;
		}

		#pragma omp critical
		{
			cnt += local_cnt;
		}
	}
	t = omp_get_wtime() - t;
	printf( "Triangle Processing Time = %0.6lf Seconds\n", t );
	return cnt;
}

void Graph::hash_join(unordered_set<int> *a, unordered_set<int> *b, vector<int> &v_join) {
	v_join.clear();
	if(a->size() > b->size()) {unordered_set<int> *c = a; a = b; b = c;}
	for(unordered_set<int>::iterator it = a->begin(); it != a->end(); it++) if(b->find(*it) != b->end()) v_join.push_back(*it);
}

long long Graph::batch_hash(const vector<pair<int,int> > &l) {
	long long cnt = 0;

	double t = omp_get_wtime();
	#pragma omp parallel
	{
		long long local_cnt = 0;
		vector<int> v_join(n);

		#pragma omp for schedule(dynamic)
		for(size_t i = 0; i < list_v.size(); ++i) {
			int u = list_v[i];
			for(int j = 0; j < deg_dlt[u]; ++j) {
				int v = con_dlt[u][j];
				hash_join(&hash_dlt[u], &hash[v], v_join);
				for(size_t k = 0; k < v_join.size(); ++k)
					if(v_join[k] < v) {
						++local_cnt;
						//printf( "Type 2: (%d,%d,%d)\n", u, v, v_join[k] );
					}
				if(v > u) continue;
				hash_join(&hash[u],&hash[v], v_join);
				local_cnt += (int)v_join.size();
				//for(int k = 0; k < it-v_join.begin(); ++k) printf( "Type 1: (%d,%d,%d)\n", u, v, v_join[k] );

				hash_join(&hash_dlt[u],&hash_dlt[v],v_join);
				for(size_t k = 0; k < v_join.size(); ++k)
					if(v_join[k] < v) {
						++local_cnt;
						//printf( "Type 3: (%d,%d,%d)\n", u, v, v_join[k] );
					}
			}
		}

		#pragma omp critical
		{
			cnt += local_cnt;
		}
	}
	t = omp_get_wtime() - t;
	printf( "Triangle Processing Time = %0.6lf Seconds\n", t );
	return cnt;
}

long long Graph::batch_insert_hash(const vector<pair<int,int> > &l) {
	//get_dlt(l, true);
	hash_l(l);
	long long cnt = batch_hash(l);
	graph_batch_insert(true);
	//clear_dlt();
	return cnt;
}

long long Graph::batch_delete_hash(const vector<pair<int,int> > &l) {
	//get_dlt(l, true);
	hash_l(l);
	graph_batch_delete(true);
	long long cnt = batch_hash(l);
	//clear_dlt();
	return cnt;
}


long long Graph::batch_join(const vector<pair<int,int> > &l) {
	long long cnt = 0;

	double t = omp_get_wtime();
	#pragma omp parallel
	{
		long long local_cnt = 0;
		vector<int> v_join(n);

		#pragma omp for schedule(dynamic)
		for(size_t i = 0; i < list_v.size(); ++i) {
			int u = list_v[i];
			for(int j = 0; j < deg_dlt[u]; ++j) {
				int v = con_dlt[u][j];
				vector<int>::iterator it = set_intersection(con_dlt[u].begin(), con_dlt[u].begin() + deg_dlt[u], con[v].begin(), con[v].begin() + deg[v], v_join.begin());
				int len = it - v_join.begin();
				for(int k = 0; k < len; ++k)
					if(v_join[k] < v) {
						++local_cnt;
						//printf( "Type 2: (%d,%d,%d)\n", u, v, v_join[k] );
					}
				if(v > u) continue;
				it = set_intersection(con[u].begin(), con[u].begin() + deg[u], con[v].begin(), con[v].begin() + deg[v], v_join.begin());
				local_cnt += it - v_join.begin();
				//for(int k = 0; k < it-v_join.begin(); ++k) printf( "Type 1: (%d,%d,%d)\n", u, v, v_join[k] );

				it = set_intersection(con_dlt[u].begin(), con_dlt[u].begin() + deg_dlt[u], con_dlt[v].begin(), con_dlt[v].begin() + deg_dlt[v], v_join.begin());
				len = it - v_join.begin();
				for(int k = 0; k < len; ++k)
					if(v_join[k] < v) {
						++local_cnt;
						//printf( "Type 3: (%d,%d,%d)\n", u, v, v_join[k] );
					}
			}
		}

		#pragma omp critical
		{
			cnt += local_cnt;
		}
	}
	t = omp_get_wtime() - t;
	printf( "Triangle Processing Time = %0.6lf Seconds\n", t );
	return cnt;
}

long long Graph::batch_insert_join(const vector<pair<int,int> > &l) {
	//get_dlt(l);
	long long cnt = batch_join(l);
	graph_batch_insert();
	//clear_dlt();
	return cnt;
}

long long Graph::batch_delete_join(const vector<pair<int,int> > &l) {
	//get_dlt(l);
	graph_batch_delete();
	long long cnt = batch_join(l);
	//clear_dlt();
	return cnt;
}


void Graph::graph_single_insert(int u, int v) {
	con[u].push_back(v);
	int p = deg[u]++;
	for(; p > 0 && con[u][p-1] > v; --p) con[u][p] = con[u][p-1];
	con[u][p] = v;
}

void Graph::graph_single_delete(int u, int v) {
	int p = 0;
	for(; con[u][p] != v; ++p);
	for(--deg[u]; p < deg[u]; ++p) con[u][p] = con[u][p+1];
	con[u].resize(deg[u]);
}

void Graph::graph_batch_insert(bool use_hash) {
	//double t = omp_get_wtime();
	//vector<char> used(n,0);
	//for(int i = 0; i < list_v.size(); ++i) used[list_v[i]] = 1;

	#pragma omp parallel
	{
		int np = omp_get_num_threads(), rank = omp_get_thread_num();
		//#pragma omp for schedule(dynamic)
		//for(size_t i = 0; i < list_v.size(); ++i) {
		for(size_t i = rank; i < list_v.size(); i+=np) {
			int v = list_v[i];
		//for(int v = rank; v < n; v+=np) {
		//	if(!used[v]) continue;
			int p0 = deg[v], p1 = deg_dlt[v];
			con[v].resize(p0+p1);
			for(int p = p0 + p1; p > 0; --p)
				if(p1 == 0 || (p0 > 0 && con[v][p0-1] > con_dlt[v][p1-1]))
					con[v][p-1] = con[v][--p0];
				else con[v][p-1] = con_dlt[v][--p1];
			deg[v] += deg_dlt[v];
			if(use_hash)
				for(int j = 0; j < deg_dlt[v]; ++j) hash[v].insert(con_dlt[v][j]);
		}
	}
	//printf( "Graph Update Time = %0.6lf Seconds\n", omp_get_wtime()-t );
}

void Graph::graph_batch_delete(bool use_hash) {
	//double t = omp_get_wtime();
	//vector<char> used(n,0);
	//for(int i = 0; i < list_v.size(); ++i) used[list_v[i]] = 1;
	#pragma omp parallel
	{
		int np = omp_get_num_threads(), rank = omp_get_thread_num();
		//#pragma omp for schedule(dynamic)
		//for(size_t i = 0; i < list_v.size(); ++i) {
		for(size_t i = rank; i < list_v.size(); i+=np) {
			int v = list_v[i];
		//for(int v = rank; v < n; v+=np) {
		//	if(!used[v]) continue;
			int p1 = 0, p = 0;
			for(int p0 = 0; p0 < deg[v]; ++p0 )
				if(p1 < deg_dlt[v] && con[v][p0] == con_dlt[v][p1])
					++p1;
				else con[v][p++] = con[v][p0];
			//set_difference(con[v].begin(),con[v].end(),con_dlt[v].begin(),con_dlt[v].end(),con[v].begin());
			deg[v] -= deg_dlt[v];
			con[v].resize(deg[v]);
			if(use_hash)
				for(int j = 0; j < deg_dlt[v]; ++j) hash[v].erase(con_dlt[v][j]);
		}
	}
	//printf( "Graph Update Time = %0.6lf Seconds\n", omp_get_wtime()-t );
}

int Graph::compare(int u, int v) {
	return compare(u, deg[u] + deg_dlt[u], v, deg[v] + deg_dlt[v]);
}

void Graph::construct_graph_aot(vector<int> *con_new, vector<bool> *is_new, vector<int> &deg_new, vector<int> &deg_plus) {
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for(int u = 0; u < n; ++u) {
			for(int i = 0; i < deg[u]; ++i) {
				int v = con[u][i];
				if(compare(u, v) < 0) {con_new[u].push_back(v);is_new[u].push_back(false);}
			}
			for(int i = 0; i < deg_dlt[u]; ++i) {
				int v = con_dlt[u][i];
				if(compare(u, v) < 0) {con_new[u].push_back(v);is_new[u].push_back(true);}
			}
			deg_plus[u] = (int) con_new[u].size();
			for(int i = 0; i < deg[u]; ++i) {
				int v = con[u][i];
				if(compare(u, v) > 0) {con_new[u].push_back(v);is_new[u].push_back(false);}
			}
			for(int i = 0; i < deg_dlt[u]; ++i) {
				int v = con_dlt[u][i];
				if(compare(u, v) > 0) {con_new[u].push_back(v);is_new[u].push_back(true);}
			}
			deg_new[u] = (int) con_new[u].size();
		}
	}
}

void Graph::construct_graph_aot_dlt(vector<int> *con_new, vector<bool> *is_new, vector<int> &deg_new, vector<int> &deg_plus) {
	vector<bool> flag(n,false);
	vector<int> deg_tmp(n,0);
	for(int u = 0; u < n; ++u)
		for(int i = 0; i < deg_dlt[u]; ++i) {
			int v = con_dlt[u][i];
			flag[u] = true; flag[v] = true;
		}
	for(int u = 0; u < n; ++u) {
		deg_tmp[u] = deg_dlt[u];
		for(int i = 0; i < deg[u]; ++i) {
			int v = con[u][i];
			if(flag[u] || flag[v]) ++deg_tmp[u];
		}
	}

	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for(int u = 0; u < n; ++u) {
			for(int i = 0; i < deg[u]; ++i) {
				int v = con[u][i];
				if((flag[u] || flag[v]) && compare(u, deg_tmp[u], v, deg_tmp[v]) < 0) {con_new[u].push_back(v);is_new[u].push_back(false);}
			}
			for(int i = 0; i < deg_dlt[u]; ++i) {
				int v = con_dlt[u][i];
				if(compare(u, deg_tmp[u], v, deg_tmp[v]) < 0) {con_new[u].push_back(v);is_new[u].push_back(true);}
			}
			deg_plus[u] = (int) con_new[u].size();
			for(int i = 0; i < deg[u]; ++i) {
				int v = con[u][i];
				if((flag[u] || flag[v]) && compare(u, deg_tmp[u], v, deg_tmp[v]) > 0) {con_new[u].push_back(v);is_new[u].push_back(false);}
			}
			for(int i = 0; i < deg_dlt[u]; ++i) {
				int v = con_dlt[u][i];
				if(compare(u, deg_tmp[u], v, deg_tmp[v]) > 0) {con_new[u].push_back(v);is_new[u].push_back(true);}
			}
			deg_new[u] = (int) con_new[u].size();
		}
	}
}

long long Graph::batch_insert_aot(const vector<pair<int,int> > &l, bool is_dlt) {
	//get_dlt(l);

	vector<int> *con_new = new vector<int>[n];
	vector<bool> *is_new = new vector<bool>[n];
	vector<int> deg_new(n,0);
	vector<int> deg_plus(n,0);

	if(is_dlt) construct_graph_aot_dlt(con_new, is_new, deg_new, deg_plus);
	else construct_graph_aot(con_new, is_new, deg_new, deg_plus);
	long long cnt = batch_aot(con_new, is_new, deg_new, deg_plus);

	delete[] con_new;
	delete[] is_new;

	graph_batch_insert();
	//clear_dlt();
	return cnt;
}

long long Graph::batch_delete_aot(const vector<pair<int,int> > &l, bool is_dlt) {
	//get_dlt(l);
	graph_batch_delete();

	vector<int> *con_new = new vector<int>[n];
	vector<bool> *is_new = new vector<bool>[n];
	vector<int> deg_new(n,0);
	vector<int> deg_plus(n,0);

	if(is_dlt) construct_graph_aot_dlt(con_new, is_new, deg_new, deg_plus);
	else construct_graph_aot(con_new, is_new, deg_new, deg_plus);
	long long cnt = batch_aot(con_new, is_new, deg_new, deg_plus);

	delete[] con_new;
	delete[] is_new;

	//clear_dlt();
	return cnt;
}


long long Graph::batch_aot(vector<int> *con, vector<bool> *is_new, vector<int> &deg, vector<int> &deg_plus) {
	printf( "Processing Triangles...\n" );
	long long cnt = 0;

	double t = omp_get_wtime();
	double alpha = 3.0;
	#pragma omp parallel
	{
		long long local_cnt = 0;
		int u, v, w;
		vector<char> last_use(n,0);

		#pragma omp for schedule(dynamic)
		for( u = 0; u < n; ++u ) {
			for( int i = 0; i < deg_plus[u]; ++i ) {
				v = con[u][i];
				last_use[v] = is_new[u][i] ? 2 : 1;
			}
			for( int i = 0; i < deg[u]; ++i ) {
				v = con[u][i];
				int c1 = is_new[u][i] ? 1 : 0;
				if( i >= deg_plus[u] && deg_plus[v] <= (int)(deg_plus[u] * alpha + 1e-8))
					for( int j = 0; j < deg_plus[v]; ++j ) {
						w = con[v][j]; if(compare(w,deg[w],u,deg[u]) <= 0) continue;
						int c2 = is_new[v][j] ? 1 : 0;
						if( last_use[w] ) {
							int c3 = last_use[w] == 2 ? 1 : 0;
							if(c1+c2+c3) ++local_cnt;
						}
					}
				else if( i < deg_plus[u] && deg_plus[u] > (int)(deg_plus[v] * alpha + 1e-8))
					for( int j = 0; j < deg_plus[v]; ++j ) {
						w = con[v][j];
						int c2 = is_new[v][j] ? 1 : 0;
						if( last_use[w] ) {
							int c3 = last_use[w] == 2 ? 1 : 0;
							if(c1+c2+c3) ++local_cnt;
						}
					}
			}
			for( int i = 0; i < deg_plus[u]; ++i ) {
				v = con[u][i];
				last_use[v] = 0;
			}
		}
		#pragma omp critical
		{
			cnt += local_cnt;
		}
	}
	t = omp_get_wtime() - t;
	printf( "Triangle Processing Time=%0.6lf Seconds.\n", t );
	return cnt;
}


Graph::~Graph() {
	delete[] con;
}
#endif /* DYNAMICTRIANGLELIST_H_ */
