class Hamiltonianclass {
	public : 
		gsl_complex     **evecEcluster  ;
		double           *evalEcluster  ;
		Hamiltonianclass(void) {}
		void allocgsl( void ) {
			evecEcluster  = mkgscmatrixd(tNC,tNC) ;
			evalEcluster  = mkvectord(tNC) ;
		}
		void eigsystem( gsl_complex **Ecluster ) {
			eigen_lapack( Ecluster, evalEcluster, evecEcluster, tNC ,0) ; printf("Ecluster :: eigen_lapack completed.\n" ) ; //fflush(stdout) ; // Row-wise vectors, i.e. rdmevec[mu][] is the mu-th eigenvectors.
		}
		void freegsl( void ) {
			freegscmatrixd( evecEcluster , tNC ) ;
			free( evalEcluster ) ;
		}
		int dotransform_Ecluster( gsl_complex **transf, gsl_complex **Ecluster ) {
			gsl_complex **backup, **right, **pose, zero = gsl_complex_rect(0,0), **inv, **result ;
			int NN = tNC ;

			backup  = mkgscmatrixd(NN, NN);
			right   = mkgscmatrixd(NN, NN);
			pose    = mkgscmatrixd(NN, NN);
			inv     = mkgscmatrixd(NN, NN);
			result  = mkgscmatrixd(NN, NN);

			print_gscmatrixd("dotransform_Ecluster", transf, tNC, tNC);
			for(int mu=0; mu<NN; mu++) for(int nu=0; nu<NN; nu++){
				pose[mu][nu] = transf[nu][mu];
				backup[mu][nu] = Ecluster[mu][nu] ;
				result[mu][nu] = zero;
				right[mu][nu] = zero;
				inv[mu][nu] = gsl_complex_conjugate(transf[mu][nu]);
			}
			//print_gscmatrixd("excord", backup, DIMd4, DIMd4);
			//print_gscmatrixd("pose",   pose  , DIMd4, DIMd4);
			//print_gscmatrixd("inv",    inv   , DIMd4, DIMd4);
			for(int mu=0; mu<NN; mu++) for(int nu=0; nu<NN; nu++) for(int k=0; k<NN; k++)
				right[mu][nu]
					= gsl_complex_add(
							right[mu][nu],
							gsl_complex_mul(backup[mu][k], pose[k][nu])
							);
			for(int mu=0; mu<NN; mu++) for(int nu=0; nu<NN; nu++) for(int k=0; k<NN; k++)
				result[mu][nu] =
					gsl_complex_add(
							result[mu][nu],
							gsl_complex_mul( inv[mu][k], right[k][nu] )
						       );

			//print_gscmatrixd("original", backup, NN, NN);
			print_gscmatrixd("EclusterT", result, NN, NN);

			freegscmatrixd(right, NN);
			freegscmatrixd(pose, NN);
			freegscmatrixd(inv, NN);
			freegscmatrixd(backup, NN);
			freegscmatrixd(result, NN);

			return 0;
		}
} ;

class excitonclass {
	public : 
		gsl_complex **excorder  ;
		gsl_complex **excorderT ;
		excitonclass(void) {}
		
		void allocgsl( void ) {
			excorder	= mkgscmatrixd(tNC,tNC) ;
			excorderT	= mkgscmatrixd(tNC,tNC) ;
		}
		void freegsl( void ) {
			freegscmatrixd( excorder	, tNC ) ;
			freegscmatrixd( excorderT	, tNC ) ;
		}
		void obtain_exc( int gnd, int gndblock, typebasis **basis , int *table, gsl_complex *ground ) {
			typebasis bin[2];
			int phase , index ; 
			for(int mu=0; mu<tNC; mu++) for(int nu=0; nu<tNC; nu++) excorder[mu][nu] = zero ;
			for(int mu=0; mu<tNC; mu++) for(int nu=0; nu<tNC; nu++){
				for(int i = gnd; i<gnd+gndblock; i++){
					if( Zero( basis[i][mu%2], mu/2) && One( basis[i][nu%2], nu/2 ) ){
						bin[0] = basis[i][0];
						bin[1] = basis[i][1];
						phase  = permu( bin, nu/2, nu%2);       bin[nu%2] = Turnoff( bin[nu%2], nu/2 );
						phase *= permu( bin, mu/2, mu%2);       bin[mu%2] = Turnon( bin[mu%2], mu/2 );
						index = (bin[0]<<Ns) + bin[1];
						excorder[mu][nu]
							= gsl_complex_add(
									excorder[mu][nu],
									gsl_complex_mul( gsl_complex_conjugate(ground[table[index]-gnd]), gsl_complex_mul_real( ground[i-gnd], phase ) )
									);
					}
				}
			}
			for(int mu=0; mu<tNC; mu++) for(int i = gnd; i<gnd+gndblock; i++){
					if( One( basis[i][mu%2], mu/2 ) ){
						excorder[mu][mu]
							= gsl_complex_add(
									excorder[mu][mu],
									gcr( gsl_complex_abs2(ground[i-gnd])* ((int)One(basis[i][mu%2], mu/2) ) , 0.) 
									) ;
					}
			}
		}
		int dotransform_excorder( gsl_complex **transf ) {	// Transforming as O -> T^* O T^t 
			gsl_complex **backup, **right, **pose, zero = gsl_complex_rect(0,0), **inv;
			int NN = tNC ;

			backup  = mkgscmatrixd(NN, NN);
			right   = mkgscmatrixd(NN, NN);
			pose    = mkgscmatrixd(NN, NN);
			inv     = mkgscmatrixd(NN, NN);

			printf("dotransform_excorder\n") ;
			print_gscmatrixd("transf", transf, tNC, tNC);
			for(int mu=0; mu<NN; mu++) for(int nu=0; nu<NN; nu++){
				pose[mu][nu] = transf[nu][mu];
				backup[mu][nu] = excorder[mu][nu] ;
				excorderT[mu][nu] = zero;
				right[mu][nu] = zero;
				inv[mu][nu] = gsl_complex_conjugate(transf[mu][nu]);
			}
			//print_gscmatrixd("excorder", backup, DIMd4, DIMd4);
			//print_gscmatrixd("pose",   pose  , DIMd4, DIMd4);
			//print_gscmatrixd("inv",    inv   , DIMd4, DIMd4);
			for(int mu=0; mu<NN; mu++) for(int nu=0; nu<NN; nu++) for(int k=0; k<NN; k++)
				right[mu][nu]
					= gsl_complex_add(
							right[mu][nu],
							gsl_complex_mul(backup[mu][k], pose[k][nu])
							);
			for(int mu=0; mu<NN; mu++) for(int nu=0; nu<NN; nu++) for(int k=0; k<NN; k++)
				excorderT[mu][nu] =
					gsl_complex_add(
							excorderT[mu][nu],
							gsl_complex_mul( inv[mu][k], right[k][nu] )
						       );

			//print_gscmatrixd("original", backup, NN, NN);
			//print_gscmatrixd("excorderT", excorderT, NN, NN);

			freegscmatrixd(right, NN);
			freegscmatrixd(pose, NN);
			freegscmatrixd(inv, NN);
			freegscmatrixd(backup, NN);

			return 0;
		}
		int dotransforminv_excorder( gsl_complex **transf ) {	// Transforming as O -> T^t O T^*
			gsl_complex **backup, **right, **pose, zero = gsl_complex_rect(0,0), **inv;
			int NN = tNC ;

			backup  = mkgscmatrixd(NN, NN);
			right   = mkgscmatrixd(NN, NN);
			pose    = mkgscmatrixd(NN, NN);
			inv     = mkgscmatrixd(NN, NN);

			printf("dotransforminv_excorder\n") ;
			print_gscmatrixd("transf", transf, tNC, tNC);
			for(int mu=0; mu<NN; mu++) for(int nu=0; nu<NN; nu++){
				pose[mu][nu] = gsl_complex_conjugate(transf[mu][nu]);
				backup[mu][nu] = excorder[mu][nu] ;
				excorderT[mu][nu] = zero;
				right[mu][nu] = zero;
				inv[mu][nu] = transf[nu][mu];
			}
			//print_gscmatrixd("excorder", backup, DIMd4, DIMd4);
			//print_gscmatrixd("pose",   pose  , DIMd4, DIMd4);
			//print_gscmatrixd("inv",    inv   , DIMd4, DIMd4);
			for(int mu=0; mu<NN; mu++) for(int nu=0; nu<NN; nu++) for(int k=0; k<NN; k++)
				right[mu][nu]
					= gsl_complex_add(
							right[mu][nu],
							gsl_complex_mul(backup[mu][k], pose[k][nu])
							);
			for(int mu=0; mu<NN; mu++) for(int nu=0; nu<NN; nu++) for(int k=0; k<NN; k++)
				excorderT[mu][nu] =
					gsl_complex_add(
							excorderT[mu][nu],
							gsl_complex_mul( inv[mu][k], right[k][nu] )
						       );

			//print_gscmatrixd("original", backup, NN, NN);
			//print_gscmatrixd("excorderT", excorderT, NN, NN);

			freegscmatrixd(right, NN);
			freegscmatrixd(pose, NN);
			freegscmatrixd(inv, NN);
			freegscmatrixd(backup, NN);

			return 0;
		}
		int write_excdeg( char pararesult[] , int ic , char writemode[] ){
			FILE *fw ; 
			char fname[1025] ;
		       	sprintf( fname, "%s/excitonicdeg.dat", pararesult ) ;
			fw = fopen( fname , writemode );
			fprintf( fw , "%d\t", ic ) ;
			for(int mu=0; mu<tNC; mu++)  for(int nu=0; nu<tNC; nu++)  {
				fprintf( fw , "%19.16lf\t", GSL_REAL( excorder[mu][nu]) ) ;
				fprintf( fw , "%19.16lf\t", GSL_IMAG( excorder[mu][nu]) ) ;
			}
			fprintf( fw , "\n") ;
			fclose( fw ) ;
			print_gscmatrixd("excitonicdeg", excorder, tNC, tNC ) ;
			printf( "\n") ;
			return 0 ;
		}
		int write_excTdeg( char pararesult[] , int ic , char writemode[] ){
			FILE *fw ; 
			char fname[1025] ;
		       	sprintf( fname, "%s/excitonicTdeg.dat", pararesult ) ;
			fw = fopen( fname , writemode );
			fprintf( fw , "%d\t", ic ) ;
			for(int mu=0; mu<tNC; mu++)  for(int nu=0; nu<tNC; nu++)  {
				fprintf( fw , "%19.16lf\t", GSL_REAL( excorderT[mu][nu]) ) ;
				fprintf( fw , "%19.16lf\t", GSL_IMAG( excorderT[mu][nu]) ) ;
			}
			fprintf( fw , "\n") ;
			fclose( fw ) ;
			print_gscmatrixd("excitonicTdeg", excorderT, tNC, tNC ) ;
			printf( "\n") ;
			return 0 ;
		}
		int write_excTdeg_datname( char pararesult[] , int ic , char datname[] , char writemode[] ){
			FILE *fw ; 
			char fname[1025] ;
		       	sprintf( fname, "%s/%s", pararesult , datname ) ;
			fw = fopen( fname , writemode );
			fprintf( fw , "%d\t", ic ) ;
			for(int mu=0; mu<tNC; mu++)  for(int nu=0; nu<tNC; nu++)  {
				fprintf( fw , "%19.16lf\t", GSL_REAL( excorderT[mu][nu]) ) ;
				fprintf( fw , "%19.16lf\t", GSL_IMAG( excorderT[mu][nu]) ) ;
			}
			fprintf( fw , "\n") ;
			fclose( fw ) ;
			print_gscmatrixd( datname , excorderT, tNC, tNC ) ;
			printf( "\n") ;
			return 0 ;
		}
		int write_excTdegdiag_datname( char pararesult[] , int ic , char datname[] , char writemode[] ){
			FILE *fw ; 
			char fname[1025] ;
		       	sprintf( fname, "%s/%s", pararesult , datname ) ;
			fw = fopen( fname , writemode );
			fprintf( fw , "%d\t", ic ) ;
			for(int mu=0; mu<tNC; mu++)  
				fprintf( fw , "%19.16lf\t", GSL_REAL( excorderT[mu][mu]) ) ;
			fprintf( fw , "\n") ;
			fclose( fw ) ;
			print_gscmatrixd( datname , excorderT, tNC, tNC ) ;
			printf( "\n") ;
			return 0 ;
		}
		//void doall( int gnd, int gndblock, typebasis **basis , int *table, gsl_complex *ground, char pararesult[] , int ic , gsl_complex **transf ) {
		//	allocgsl() ; 
		//	obtain_exc( gnd, gndblock, basis, table, ground ) ;
		//	write_excdeg( pararesult, ic ) ;
		//	dotransform_excorder( transf ) ;
		//}
} ;

int write_excT( char pararesult[] , int ic , gsl_complex **excT , char writemode[] ){
	FILE *fw ; 
	char fname[1025] ;
	sprintf( fname, "%s/excitonicT.dat", pararesult ) ;
	fw = fopen( fname , writemode );
	fprintf( fw , "%d\t", ic ) ;
	for(int mu=0; mu<tNC; mu++)  for(int nu=0; nu<tNC; nu++)  {
		fprintf( fw , "%19.16lf\t", GSL_REAL( excT[mu][nu]) ) ;
		fprintf( fw , "%19.16lf\t", GSL_IMAG( excT[mu][nu]) ) ;
	}
	fprintf( fw , "\n") ;
	fclose( fw ) ;
	//print_gscmatrixd("excitonicT", excT, tNC, tNC ) ;
	printf( "\n") ;
	return 0 ;
}
int write_excname( char pararesult[] , int ic , gsl_complex **excT , char datname[] , char writemode[] ){
	FILE *fw ; 
	char fname[1025] ;
	sprintf( fname, "%s/%s.dat", pararesult , datname ) ;
	fw = fopen( fname , writemode );
	fprintf( fw , "%d\t", ic ) ;
	for(int mu=0; mu<tNC; mu++)  for(int nu=0; nu<tNC; nu++)  {
		fprintf( fw , "%19.16lf\t", GSL_REAL( excT[mu][nu]) ) ;
		fprintf( fw , "%19.16lf\t", GSL_IMAG( excT[mu][nu]) ) ;
	}
	fprintf( fw , "\n") ;
	fclose( fw ) ;
	//print_gscmatrixd("excitonicT", excT, tNC, tNC ) ;
	printf( "\n") ;
	return 0 ;
}
