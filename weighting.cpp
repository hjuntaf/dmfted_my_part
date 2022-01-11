double *Wwn ;

#define ALPHA	1.
#define WWN	Wwn[i]

void allocWwn( int nmax ){
	Wwn = mkvectord( nmax );
}
void initWwn( double *matsu, int nmax ){
	allocWwn( nmax );
	double sum = 0.;
	for( int i=0; i<nmax; i++ ){
		WWN = pow( matsu[i] , -ALPHA );
		sum += WWN;
	}
	sum /= nmax;
	for( int i=0; i<nmax; i++ ) WWN /= sum;
	printf( "Weight-coeffcient : 1./%lf\n", sum );
}

int print_doubleArr( char fname[], double matsu[], double Green[], int start, int end){
        FILE *fw;
        fw = fopen( fname , "w" );
        for(int i=start; i<end; i++){
                fprintf( fw, "%.10lf\t", matsu[i] );
                fprintf( fw, "%.10lf\n", Green[i] );
        }
        fclose(fw);
        return 0;
}
