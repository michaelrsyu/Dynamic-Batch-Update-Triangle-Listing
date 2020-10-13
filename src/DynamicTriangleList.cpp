#include "DynamicTriangleList.h"


void txt_to_bin(string path) {
	Graph::txt_to_bin(path);
}

void show_info(string path) {
	Graph *g = new Graph(path);
	printf( "n=%d,m=%lld\n", g->n, g->m );
	int maxd = 0;
	double sum_sqr = 0;
	for( int i = 0; i < g->n; ++i ) {
		maxd = max(maxd, g->deg[i]);
		sum_sqr += g->deg[i] * 1.0 * g->deg[i];
	}

	printf( "maxd = %d, sum_sqr = %0.3lfGB\n", maxd, sum_sqr/(1000000000.0) );

	delete g;
}

void show_raw(string path) {
	FILE *fin = fopen( (path + "graph.txt").c_str(), "r" );
	char line[1024], st[64];
	int n = 0, a, b, num_cnt = Graph::get_num_cnt(path);
	long long cnt = 0, m = 0;

	printf( "Loading text...\n" );
	while( fgets( line, 1024, fin ) ) {
		if( !Graph::get_edge(line, a, b, num_cnt) ) continue;
		n = max(max(n, a+1), b+1);
		if( (++cnt) % (long long) 10000000 == 0 ) printf( "%lld lines finished\n", cnt );
		++m;
	}
	fclose( fin );
	printf( "n=%d,m=%lld\n", n, m );
}

void run_batch(string path, double percent, int alg, int n_thread) {
	Graph g(path);
	vector<pair<int,int> > l;
	g.generate_edges(percent, l); g.get_dlt(l);
	printf( "%d (%0.3lf\%) edges generated\n", percent, (int)l.size() );
	if(alg == ALG_HASH) g.init_hash();
	g.run_batch(MODE_DELETE, n_thread, alg, l);
	g.run_batch(MODE_INSERT, n_thread, alg, l);
	g.clear_dlt();
}

int main(int argc, char *argv[]) {
	printf( "argc=%d\n", argc );
	for( int i = 0; i < argc; ++i )
		printf( "argv[%d]=%s\n", i, argv[i] );

	setvbuf(stdout, NULL, _IONBF, 0);
	setvbuf(stderr, NULL, _IONBF, 0);

	if( argc > 1 ) {
		if(strcmp( argv[1], "txt-to-bin" ) == 0)
			txt_to_bin( argv[2] );
		else if( strcmp(argv[1], "show-info" ) == 0 )
			show_info(argv[2]);
		else if( strcmp(argv[1], "show-raw" ) == 0 )
			show_raw(argv[2]);
		else {
			run_batch(argv[1], argc>2?atof(argv[2]):10.0, argc>3?atoi(argv[3]):ALG_NEW, argc>4?atoi(argv[4]):1);
		}
	}

	return 0;
}
