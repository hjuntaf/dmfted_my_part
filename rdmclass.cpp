//#define DIMd 	64						//in lattice.h
//#define DIMd4	15						//in lattice.h
//#define gcr(zre, zim)	gsl_complex_rect( zre, zim )		//in lattice.h
//gsl_complex **transformMatd4LSJ = mkgscmatrixd(DIMd4, DIMd4);	//in static

class rdmclass {
	public : 
		gsl_complex **transformMatd4LSJ ;
		gsl_complex **transformMatd4jjJ ;
		gsl_complex **transformMatd4LSJcubic ;
		gsl_complex **transformMatd4t2g ;
		gsl_complex **transformMatd4Localdiag ;
		gsl_complex **rdmLSJ	   ;
		gsl_complex **rdmjjJ	   ;
		gsl_complex **rdmLSJcubic  ;
		gsl_complex **rdmt2g	   ;
		gsl_complex **rdmLocaldiag ;
		gsl_complex **jzopLSJ	   ;
		gsl_complex **rdmjeff	   ;
		gsl_complex **rdmorb0;
		gsl_complex **rdmorb1;
		gsl_complex **rdmorb2;
		rdmclass(void) {}
		void allocgsl( void ) {
			transformMatd4LSJ = mkgscmatrixd(DIMd4,DIMd4) ;
			transformMatd4jjJ = mkgscmatrixd(DIMd4,DIMd4) ;
			transformMatd4LSJcubic = mkgscmatrixd(DIMd4,DIMd4) ;
			transformMatd4t2g = mkgscmatrixd(DIMd4,DIMd4) ;
			transformMatd4Localdiag = mkgscmatrixd(DIMd4,DIMd4) ;
			rdmLSJ	          = mkgscmatrixd(DIMd4,DIMd4) ;
			rdmjjJ	          = mkgscmatrixd(DIMd4,DIMd4) ;
			rdmLSJcubic       = mkgscmatrixd(DIMd4,DIMd4) ;
			rdmt2g	          = mkgscmatrixd(DIMd4,DIMd4) ;
			rdmLocaldiag      = mkgscmatrixd(DIMd4,DIMd4) ;
			rdmorb0	          = mkgscmatrixd(4,4) ;
			rdmorb1	          = mkgscmatrixd(4,4) ;
			rdmorb2	          = mkgscmatrixd(4,4) ;
		}
		void freegsl( void ) {
			freegscmatrixd( transformMatd4LSJ , DIMd4 ) ;
			freegscmatrixd( transformMatd4jjJ , DIMd4 ) ;
			freegscmatrixd( transformMatd4LSJcubic , DIMd4 ) ;
			freegscmatrixd( transformMatd4t2g , DIMd4 ) ;
			freegscmatrixd( transformMatd4Localdiag , DIMd4 ) ;
			freegscmatrixd( rdmjjJ       , DIMd4 ) ;
			freegscmatrixd( rdmLSJ       , DIMd4 ) ;
			freegscmatrixd( rdmLSJcubic  , DIMd4 ) ;
			freegscmatrixd( rdmt2g       , DIMd4 ) ;
			freegscmatrixd( rdmLocaldiag , DIMd4 ) ;
			freegscmatrixd( rdmorb0      , 4 ) ;
			freegscmatrixd( rdmorb1      , 4 ) ;
			freegscmatrixd( rdmorb2      , 4 ) ;
		}
		int init_transformMatd4jjJ( void ){
			gsl_complex trans[DIMd4][DIMd4] = {
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(-1.,0), zero, zero, zero, zero, }, 
				{ gcr(0.866025,0), zero, zero, zero, zero, zero, zero, zero, gcr(-0.5,0), zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, gcr(-0.866025,0), zero, gcr(0.5,0), zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(-0.707107,0), zero, gcr(0.707107,0), }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(1.,0), zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(1.,0), zero, zero, zero, zero, zero, }, 
				{ gcr(0.5,0), zero, zero, zero, zero, zero, zero, zero, gcr(0.866025,0), zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, gcr(0.5,0), zero, gcr(0.866025,0), zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(0.707107,0), zero, gcr(0.707107,0), }, 
				{ zero, zero, gcr(-1.,0), zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, gcr(-1.,0), zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, gcr(-1.,0), zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, gcr(-1.,0), zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, gcr(-0.707107,0), zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(-0.707107,0), zero, }, 
				{ zero, gcr(-0.707107,0), zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(0.707107,0), zero, }
			};
			int mu, nu ;
			for(mu=0; mu<DIMd4; mu++) for(nu=0; nu<DIMd4; nu++){
				transformMatd4jjJ[mu][nu] = zero;
			}
			for(mu=0; mu<DIMd4; mu++) for(nu=0; nu<DIMd4; nu++){
				transformMatd4jjJ[mu][nu] = gsl_complex_rect( GSL_REAL(trans[mu][nu]) , GSL_IMAG(trans[mu][nu]) ) ;
			}
			print_gscmatrixd("init_transformMatd4jjJ", transformMatd4jjJ, DIMd4, DIMd4);
			return 0;
		}
		int init_transformMatd4LSJ( void ){
			gsl_complex trans[DIMd4][DIMd4] = {
				{ zero, gcr(-0.408248,0), zero, zero, zero, zero, zero, zero, zero, zero, gcr(0.816497,0), zero, zero, gcr(0.408248,0), zero, }, 
				{ gcr(0.866025,0), zero, zero, zero, zero, zero, zero, zero, gcr(-0.5,0), zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, gcr(-0.866025,0), zero, gcr(0.5,0), zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(-0.707107,0), zero, gcr(0.707107,0), }, 
				{ zero, zero, gcr(0.816497,0), zero, zero, zero, zero, zero, zero, zero, zero, gcr(0.57735,0), zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, gcr(0.816497,0), zero, zero, zero, gcr(0.57735,0), zero, zero, zero, zero, zero, }, 
				{ gcr(0.288675,0), zero, zero, zero, zero, zero, zero, gcr(0.816497,0), gcr(0.5,0), zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, gcr(0.816497,0), gcr(0.288675,0), zero, gcr(0.5,0), zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, gcr(0.57735,0), zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(0.408248,0), gcr(0.57735,0), gcr(0.408248,0), }, 
				{ zero, zero, gcr(-0.57735,0), zero, zero, zero, zero, zero, zero, zero, zero, gcr(0.816497,0), zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, gcr(-0.57735,0), zero, zero, zero, gcr(0.816497,0), zero, zero, zero, zero, zero, }, 
				{ gcr(0.408248,0), zero, zero, zero, zero, zero, zero, gcr(-0.57735,0), gcr(0.707107,0), zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, gcr(-0.57735,0), gcr(0.408248,0), zero, gcr(0.707107,0), zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, gcr(-0.408248,0), zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(0.57735,0), gcr(-0.408248,0), gcr(0.57735,0), }, 
				{ zero, gcr(-0.57735,0), zero, zero, zero, zero, zero, zero, zero, zero, gcr(-0.57735,0), zero, zero, gcr(0.57735,0), zero, }
		       	};
			int mu, nu ;
			for(mu=0; mu<DIMd4; mu++) for(nu=0; nu<DIMd4; nu++){
				transformMatd4LSJ[mu][nu] = zero;
			}
			for(mu=0; mu<DIMd4; mu++) for(nu=0; nu<DIMd4; nu++){
				transformMatd4LSJ[mu][nu] = gsl_complex_rect( GSL_REAL(trans[mu][nu]) , GSL_IMAG(trans[mu][nu]) ) ;
			}
			print_gscmatrixd("init_transformMatd4LSJ", transformMatd4LSJ, DIMd4, DIMd4);
			return 0;
		}
		int init_transformMatd4Localdiag( void ){
			gsl_complex trans[DIMd4][DIMd4] = {
				{ zero, gcr(0.19681,0), zero, zero, zero, zero, zero, zero, zero, zero, gcr(-0.873045,0), zero, gcr(-0.145068,0), gcr(-0.396192,0), gcr(-0.145068,0), }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(-0.707107,0), zero, gcr(0.707107,0), }, 
				{ gcr(-0.060308,0), zero, zero, gcr(-0.28797,0), gcr(0.620861,0), zero, gcr(-0.722674,0), gcr(0.0279722,0), gcr(0.0701976,0), zero, zero, zero, zero, zero, zero, }, 
				{ gcr(0.620861,0), zero, zero, gcr(-0.0279722,0), gcr(0.060308,0), zero, gcr(-0.0701976,0), gcr(-0.28797,0), gcr(-0.722674,0), zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, gcr(-0.691321,0), zero, zero, zero, zero, zero, zero, zero, zero, gcr(-0.722548,0), zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, gcr(0.691321,0), zero, zero, zero, gcr(0.722548,0), zero, zero, zero, zero, zero, }, 
				{ gcr(0.0481643,0), zero, zero, gcr(-0.521121,0), gcr(-0.736977,0), zero, gcr(-0.425494,0), gcr(0.0340573,0), gcr(0.0278077,0), zero, zero, zero, zero, zero, zero, }, 
				{ zero, gcr(-0.46291,0), zero, zero, zero, zero, zero, zero, zero, zero, gcr(0.205738,0), zero, gcr(-0.581914,0), gcr(-0.257172,0), gcr(-0.581914,0), }, 
				{ gcr(0.736977,0), zero, zero, gcr(0.0340573,0), gcr(0.0481643,0), zero, gcr(0.0278077,0), gcr(0.521121,0), gcr(0.425494,0), zero, zero, zero, zero, zero, zero, }, 
				{ zero, gcr(-0.492792,0), zero, zero, zero, zero, zero, zero, zero, zero, gcr(0.0835819,0), zero, gcr(0.364568,0), gcr(-0.695953,0), gcr(0.364568,0), }, 
				{ gcr(0.010692,0), zero, zero, gcr(0.801523,0), gcr(-0.255588,0), zero, gcr(-0.53897,0), gcr(-0.0335302,0), gcr(0.0225468,0), zero, zero, zero, zero, zero, zero, }, 
				{ gcr(0.255588,0), zero, zero, gcr(-0.0335302,0), gcr(0.010692,0), zero, gcr(0.0225468,0), gcr(-0.801523,0), gcr(0.53897,0), zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, gcr(0.722548,0), zero, zero, zero, zero, zero, zero, zero, zero, gcr(-0.691321,0), zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, gcr(-0.722548,0), zero, zero, zero, gcr(0.691321,0), zero, zero, zero, zero, zero, }, 
				{ zero, gcr(-0.710025,0), zero, zero, zero, zero, zero, zero, zero, zero, gcr(-0.43414,0), zero, gcr(0.0861478,0), gcr(0.540873,0), gcr(0.0861478,0), }
		       	};
			int mu, nu ;
			for(mu=0; mu<DIMd4; mu++) for(nu=0; nu<DIMd4; nu++){
				transformMatd4Localdiag[mu][nu] = zero;
			}
			for(mu=0; mu<DIMd4; mu++) for(nu=0; nu<DIMd4; nu++){
				transformMatd4Localdiag[mu][nu] = gsl_complex_rect( GSL_REAL(trans[mu][nu]) , GSL_IMAG(trans[mu][nu]) ) ;
			}
			print_gscmatrixd("init_transformMatd4Localdiag", transformMatd4Localdiag, DIMd4, DIMd4);
			return 0;
		}
		int init_transformMatd4LSJcubic( void ){
			gsl_complex trans[DIMd4][DIMd4] = {
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }	// Unsupported now
		       	};
			int mu, nu ;
			for(mu=0; mu<DIMd4; mu++) for(nu=0; nu<DIMd4; nu++){
				transformMatd4LSJcubic[mu][nu] = zero;
			}
			for(mu=0; mu<DIMd4; mu++) for(nu=0; nu<DIMd4; nu++){
				transformMatd4LSJcubic[mu][nu] = gsl_complex_rect( GSL_REAL(trans[mu][nu]) , GSL_IMAG(trans[mu][nu]) ) ;
			}
			print_gscmatrixd("init_transformMatd4LSJcubic", transformMatd4LSJcubic, DIMd4, DIMd4);
			return 0;
		}
		int init_transformMatd4t2g( void ){
			gsl_complex trans[DIMd4][DIMd4] = {
				{ zero, zero, zero, gcr(0.,0.57735), gcr(0.,-0.408248), zero, zero, zero, gcr(0.,-0.707107), zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, gcr(-0.57735,0), gcr(0.408248,0), zero, zero, zero, gcr(-0.707107,0), zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, gcr(-0.333333,0), zero, gcr(-0.471405,0), gcr(0.666667,0), gcr(-0.471405,0), }, 
				{ zero, zero, zero, zero, zero, gcr(-0.57735,0), zero, zero, zero, gcr(-0.408248,0), gcr(-0.333333,0), zero, gcr(-0.471405,0), gcr(-0.333333,0), gcr(0.235702,0), }, 
				{ zero, gcr(0.,0.5), gcr(0.,-0.288675), zero, zero, gcr(0.,0.288675), zero, zero, zero, gcr(0.,-0.408248), gcr(0.,-0.333333), gcr(0.,0.408248), gcr(0.,0.235702), gcr(0.,0.166667), gcr(0.,0.235702), }, 
				{ zero, gcr(-0.5,0), gcr(0.288675,0), zero, zero, gcr(0.288675,0), zero, zero, zero, gcr(-0.408248,0), gcr(-0.333333,0), gcr(-0.408248,0), gcr(0.235702,0), gcr(0.166667,0), gcr(0.235702,0), }, 
				{ zero, gcr(0.,0.5), gcr(0.,0.288675), zero, zero, gcr(0.,-0.288675), zero, zero, zero, gcr(0.,0.408248), gcr(0.,-0.333333), gcr(0.,-0.408248), gcr(0.,0.235702), gcr(0.,0.166667), gcr(0.,0.235702), }, 
				{ zero, zero, gcr(-0.57735,0), zero, zero, zero, zero, zero, zero, zero, gcr(-0.333333,0), gcr(-0.408248,0), gcr(0.235702,0), gcr(-0.333333,0), gcr(-0.471405,0), }, 
				{ gcr(0.,0.408248), zero, zero, zero, zero, zero, gcr(0.,0.707107), gcr(0.,-0.57735), zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, gcr(0.,0.57735), gcr(0.,0.816497), zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, gcr(-0.5,0), gcr(-0.288675,0), zero, zero, gcr(-0.288675,0), zero, zero, zero, gcr(0.408248,0), gcr(-0.333333,0), gcr(0.408248,0), gcr(0.235702,0), gcr(0.166667,0), gcr(0.235702,0), }, 
				{ zero, zero, gcr(0.,-0.57735), zero, zero, zero, zero, zero, zero, zero, gcr(0.,0.333333), gcr(0.,-0.408248), gcr(0.,-0.235702), gcr(0.,0.333333), gcr(0.,0.471405), }, 
				{ gcr(0.,0.816497), zero, zero, zero, zero, zero, zero, gcr(0.,0.57735), zero, zero, zero, zero, zero, zero, zero, }, 
				{ gcr(0.408248,0), zero, zero, zero, zero, zero, gcr(-0.707107,0), gcr(-0.57735,0), zero, zero, zero, zero, zero, zero, zero, }, 
				{ zero, zero, zero, zero, zero, gcr(0.,0.57735), zero, zero, zero, gcr(0.,0.408248), gcr(0.,-0.333333), zero, gcr(0.,-0.471405), gcr(0.,-0.333333), gcr(0.,0.235702), }
		       	};
			int mu, nu ;
			for(mu=0; mu<DIMd4; mu++) for(nu=0; nu<DIMd4; nu++){
				transformMatd4t2g[mu][nu] = zero;
			}
			for(mu=0; mu<DIMd4; mu++) for(nu=0; nu<DIMd4; nu++){
				transformMatd4t2g[mu][nu] = gsl_complex_rect( GSL_REAL(trans[mu][nu]) , GSL_IMAG(trans[mu][nu]) ) ;
			}
			print_gscmatrixd("init_transformMatd4t2g", transformMatd4t2g, DIMd4, DIMd4);
			return 0;
		}
		int obtainRDMorbt2g( void ) {
			for(int i=0; i<4; i++) for(int j=0; j<4; j++) {
				rdmorb0[i][j] = zero ;
				rdmorb1[i][j] = zero ;
				rdmorb2[i][j] = zero ;
			}
			int NN = DIMd4 ;
			for(int i=0; i<NN; i++) {
				unsigned int  ciup = (i>>NC) ,
					      cidn = i - ciup;
				int orbmu0 = One(ciup,0) + 2*One(cidn,0);
				int orbmu1 = One(ciup,1) + 2*One(cidn,1);
				int orbmu2 = One(ciup,2) + 2*One(cidn,2);
				for(int j=0; j<NN; j++){
					unsigned int  cjup = (j>>NC) ,
						      cjdn = j - cjup;
					int orbnu0 = One(cjup,0) + 2*One(cjdn,0);
					int orbnu1 = One(cjup,1) + 2*One(cjdn,1);
					int orbnu2 = One(cjup,2) + 2*One(cjdn,2);
					if( (orbmu1==orbnu1)&&(orbmu2==orbnu2) )
							rdmorb0[orbmu0][orbnu0] = gsl_complex_add(  rdmorb0[orbmu0][orbnu0] , rdmt2g[i][j]  ) ;
					if( (orbmu0==orbnu0)&&(orbmu2==orbnu2) )
							rdmorb1[orbmu1][orbnu1] = gsl_complex_add(  rdmorb1[orbmu1][orbnu1] , rdmt2g[i][j]  ) ;
					if( (orbmu0==orbnu0)&&(orbmu1==orbnu1) )
							rdmorb2[orbmu2][orbnu2] = gsl_complex_add(  rdmorb2[orbmu2][orbnu2] , rdmt2g[i][j]  ) ;
				}
			}
			return 0;
		}
		int dotransform_jjJd4( gsl_complex **original ) {
			int rdmvecIndd4[DIMd4] =  { 51, 53, 54, 60, 29, 45, 43, 39, 30, 57, 27, 23, 15, 46, 58 };
			int JzjjJfactor[DIMd4] =  {0, -1, 1, 0, -2, 2, -1, 1, 0, -2, 2, -1, 1, 0, 0};
			int Jzjefffactor[DIMd4]=  {-1, 0, -2, 1, 1, 2, 1, -1, -1, 2, 0, -2, 0, 0, 0};
			int mu, nu, k;
			gsl_complex **backup, **right, **pose, zero = gsl_complex_rect(0,0), **inv;
			int NN = DIMd4 ;
			double diagsum = 0 ;
			double Jzd4    = 0 ;

			backup  = mkgscmatrixd(NN, NN);
			right   = mkgscmatrixd(NN, NN);
			pose    = mkgscmatrixd(NN, NN);
			inv     = mkgscmatrixd(NN, NN);

			print_gscmatrixd("dotransformMatd4jjJ", transformMatd4jjJ, DIMd4, DIMd4);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
				pose[mu][nu] = transformMatd4jjJ[nu][mu];
				int selmu = rdmvecIndd4[mu] ;
				int selnu = rdmvecIndd4[nu] ;
				backup[mu][nu] = original[selmu][selnu] ;
				rdmjjJ[mu][nu] = zero;
				right[mu][nu] = zero;
				inv[mu][nu] = gsl_complex_conjugate(transformMatd4jjJ[mu][nu]);
			}
			print_gscmatrixd("rdmd4", backup, DIMd4, DIMd4);
			//print_gscmatrixd("pose",   pose  , DIMd4, DIMd4);
			//print_gscmatrixd("inv",    inv   , DIMd4, DIMd4);
			printf( "rdmd4diagsum : " ) ;
			diagsum = 0 ;
			Jzd4    = 0 ;
			for(int mu=0; mu<DIMd4; mu++) 
			{
				diagsum +=  GSL_REAL( backup[mu][mu]) ;
				Jzd4    +=  GSL_REAL( backup[mu][mu]) * Jzjefffactor[mu] ;
			}
			printf( "%19.16lf\n" , diagsum ) ;
			printf( "Jzd4jeff : %19.16lf\n", Jzd4 ) ;

			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				right[mu][nu]
					= gsl_complex_add(
							right[mu][nu],
							gsl_complex_mul(backup[mu][k], pose[k][nu])
							);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				rdmjjJ[mu][nu] =
					gsl_complex_add(
							rdmjjJ[mu][nu],
							gsl_complex_mul( inv[mu][k], right[k][nu] )
						       );

			//print_gscmatrixd("original", backup, NN, NN);
			print_gscmatrixd("rdmjjJ", rdmjjJ, NN, NN);

			printf( "rdmjjJd4diagsum : " ) ;
			diagsum = 0 ;
			Jzd4    = 0 ;
			for(int mu=0; mu<DIMd4; mu++)
			{
				diagsum +=  GSL_REAL( rdmjjJ[mu][mu]) ;
				Jzd4    +=  GSL_REAL( rdmjjJ[mu][mu]) * JzjjJfactor[mu] ;
			}
			printf( "%19.16lf\n" , diagsum ) ;
			printf( "Jzd4jjJ : %19.16lf\n", Jzd4 ) ;

			freegscmatrixd(right, NN);
			freegscmatrixd(pose, NN);
			freegscmatrixd(inv, NN);
			freegscmatrixd(backup, NN);

			return 0;
		}
		int dotransform_LSJd4( gsl_complex **original ) {
			int rdmvecIndd4[DIMd4] =  { 51, 53, 54, 60, 29, 45, 43, 39, 30, 57, 27, 23, 15, 46, 58 };
			int JzLSJfactor[DIMd4] =  {0, 1, -1, 0, 2, -2, 1, -1, 0, 2, -2, 1, -1, 0, 0};
			int Jzjefffactor[DIMd4]=  {-1, 0, -2, 1, 1, 2, 1, -1, -1, 2, 0, -2, 0, 0, 0};
			int mu, nu, k;
			gsl_complex **backup, **right, **pose, zero = gsl_complex_rect(0,0), **inv;
			int NN = DIMd4 ;
			double diagsum = 0 ;
			double Jzd4    = 0 ;

			backup  = mkgscmatrixd(NN, NN);
			right   = mkgscmatrixd(NN, NN);
			pose    = mkgscmatrixd(NN, NN);
			inv     = mkgscmatrixd(NN, NN);

			print_gscmatrixd("dotransformMatd4LSJ", transformMatd4LSJ, DIMd4, DIMd4);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
				pose[mu][nu] = transformMatd4LSJ[nu][mu];
				int selmu = rdmvecIndd4[mu] ;
				int selnu = rdmvecIndd4[nu] ;
				backup[mu][nu] = original[selmu][selnu] ;
				rdmLSJ[mu][nu] = zero;
				right[mu][nu] = zero;
				inv[mu][nu] = gsl_complex_conjugate(transformMatd4LSJ[mu][nu]);
			}
			print_gscmatrixd("backup", backup, DIMd4, DIMd4);
			//print_gscmatrixd("pose",   pose  , DIMd4, DIMd4);
			//print_gscmatrixd("inv",    inv   , DIMd4, DIMd4);
			printf( "rdmd4diagsum : " ) ;
			diagsum = 0 ;
			Jzd4    = 0 ;
			for(int mu=0; mu<DIMd4; mu++) 
			{
				diagsum +=  GSL_REAL( backup[mu][mu]) ;
				Jzd4    +=  GSL_REAL( backup[mu][mu]) * Jzjefffactor[mu] ;
			}
			printf( "%19.16lf\n" , diagsum ) ;
			printf( "Jzd4jeff : %19.16lf\n", Jzd4 ) ;

			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				right[mu][nu]
					= gsl_complex_add(
							right[mu][nu],
							gsl_complex_mul(backup[mu][k], pose[k][nu])
							);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				rdmLSJ[mu][nu] =
					gsl_complex_add(
							rdmLSJ[mu][nu],
							gsl_complex_mul( inv[mu][k], right[k][nu] )
						       );

			//print_gscmatrixd("original", backup, NN, NN);
			print_gscmatrixd("rdmLSJ", rdmLSJ, NN, NN);

			printf( "rdmLSJd4diagsum : " ) ;
			diagsum = 0 ;
			Jzd4    = 0 ;
			for(int mu=0; mu<DIMd4; mu++)
			{
				diagsum +=  GSL_REAL( rdmLSJ[mu][mu]) ;
				Jzd4    +=  GSL_REAL( rdmLSJ[mu][mu]) * JzLSJfactor[mu] ;
			}
			printf( "%19.16lf\n" , diagsum ) ;
			printf( "Jzd4LSJ : %19.16lf\n", Jzd4 ) ;

			freegscmatrixd(right, NN);
			freegscmatrixd(pose, NN);
			freegscmatrixd(inv, NN);
			freegscmatrixd(backup, NN);

			return 0;
		}
		int dotransform_LSJcubicd4( gsl_complex **original ) {
			int rdmvecIndd4[DIMd4] =  { 51, 53, 54, 60, 29, 45, 43, 39, 30, 57, 27, 23, 15, 46, 58 };
			int Jzjefffactor[DIMd4]=  {-1, 0, -2, 1, 1, 2, 1, -1, -1, 2, 0, -2, 0, 0, 0};
			int mu, nu, k;
			gsl_complex **backup, **right, **pose, zero = gsl_complex_rect(0,0), **inv;
			int NN = DIMd4 ;
			double diagsum = 0 ;

			backup  = mkgscmatrixd(NN, NN);
			right   = mkgscmatrixd(NN, NN);
			pose    = mkgscmatrixd(NN, NN);
			inv     = mkgscmatrixd(NN, NN);

			print_gscmatrixd("dotransformMatd4LSJcubic", transformMatd4LSJcubic, DIMd4, DIMd4);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
				pose[mu][nu] = transformMatd4LSJcubic[nu][mu];
				int selmu = rdmvecIndd4[mu] ;
				int selnu = rdmvecIndd4[nu] ;
				backup[mu][nu] = original[selmu][selnu] ;
				rdmLSJcubic[mu][nu] = zero;
				right[mu][nu] = zero;
				inv[mu][nu] = gsl_complex_conjugate(transformMatd4LSJcubic[mu][nu]);
			}
			print_gscmatrixd("backup", backup, DIMd4, DIMd4);
			//print_gscmatrixd("pose",   pose  , DIMd4, DIMd4);
			//print_gscmatrixd("inv",    inv   , DIMd4, DIMd4);
			printf( "rdmd4diagsum : " ) ;
			diagsum = 0 ;
			for(int mu=0; mu<DIMd4; mu++) 
			{
				diagsum +=  GSL_REAL( backup[mu][mu]) ;
			}
			printf( "%19.16lf\n" , diagsum ) ;

			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				right[mu][nu]
					= gsl_complex_add(
							right[mu][nu],
							gsl_complex_mul(backup[mu][k], pose[k][nu])
							);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				rdmLSJcubic[mu][nu] =
					gsl_complex_add(
							rdmLSJcubic[mu][nu],
							gsl_complex_mul( inv[mu][k], right[k][nu] )
						       );

			//print_gscmatrixd("original", backup, NN, NN);
			print_gscmatrixd("rdmLSJcubic", rdmLSJcubic, NN, NN);

			printf( "rdmLSJcubicd4diagsum : " ) ;
			diagsum = 0 ;
			for(int mu=0; mu<NN; mu++)
			{
				diagsum +=  GSL_REAL( rdmLSJcubic[mu][mu]) ;
			}
			printf( "%19.16lf\n" , diagsum ) ;

			freegscmatrixd(right, NN);
			freegscmatrixd(pose, NN);
			freegscmatrixd(inv, NN);
			freegscmatrixd(backup, NN);

			return 0;
		}
		int dotransform_t2gd4( gsl_complex **original ) {
			int rdmvecIndd4[DIMd4] =  { 51, 53, 54, 60, 29, 45, 43, 39, 30, 57, 27, 23, 15, 46, 58 };
			int mu, nu, k;
			gsl_complex **backup, **right, **pose, zero = gsl_complex_rect(0,0), **inv;
			int NN = DIMd4 ;
			double diagsum = 0 ;

			backup  = mkgscmatrixd(NN, NN);
			right   = mkgscmatrixd(NN, NN);
			pose    = mkgscmatrixd(NN, NN);
			inv     = mkgscmatrixd(NN, NN);

			print_gscmatrixd("dotransformMatd4t2g", transformMatd4t2g, DIMd4, DIMd4);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
				pose[mu][nu] = transformMatd4t2g[nu][mu];
				int selmu = rdmvecIndd4[mu] ;
				int selnu = rdmvecIndd4[nu] ;
				backup[mu][nu] = original[selmu][selnu] ;
				rdmt2g[mu][nu] = zero;
				right[mu][nu] = zero;
				inv[mu][nu] = gsl_complex_conjugate(transformMatd4t2g[mu][nu]);
			}
			print_gscmatrixd("backup", backup, DIMd4, DIMd4);
			//print_gscmatrixd("pose",   pose  , DIMd4, DIMd4);
			//print_gscmatrixd("inv",    inv   , DIMd4, DIMd4);
			printf( "rdmd4diagsum : " ) ;
			diagsum = 0 ;
			for(int mu=0; mu<DIMd4; mu++) 
			{
				diagsum +=  GSL_REAL( backup[mu][mu]) ;
			}
			printf( "%19.16lf\n" , diagsum ) ;

			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				right[mu][nu]
					= gsl_complex_add(
							right[mu][nu],
							gsl_complex_mul(backup[mu][k], pose[k][nu])
							);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				rdmt2g[mu][nu] =
					gsl_complex_add(
							rdmt2g[mu][nu],
							gsl_complex_mul( inv[mu][k], right[k][nu] )
						       );

			//print_gscmatrixd("original", backup, NN, NN);
			print_gscmatrixd("rdmt2g", rdmt2g, NN, NN);

			printf( "rdmt2gd4diagsum : " ) ;
			diagsum = 0 ;
			for(int mu=0; mu<DIMd4; mu++)
			{
				diagsum +=  GSL_REAL( rdmt2g[mu][mu]) ;
			}
			printf( "%19.16lf\n" , diagsum ) ;

			freegscmatrixd(right, NN);
			freegscmatrixd(pose, NN);
			freegscmatrixd(inv, NN);
			freegscmatrixd(backup, NN);

			return 0;
		}
		int dotransform_Localdiagd4( gsl_complex **original ) {
			int rdmvecIndd4[DIMd4] =  { 51, 53, 54, 60, 29, 45, 43, 39, 30, 57, 27, 23, 15, 46, 58 };
			int mu, nu, k;
			gsl_complex **backup, **right, **pose, zero = gsl_complex_rect(0,0), **inv;
			int NN = DIMd4 ;
			double diagsum = 0 ;
			double Jzd4    = 0 ;

			backup  = mkgscmatrixd(NN, NN);
			right   = mkgscmatrixd(NN, NN);
			pose    = mkgscmatrixd(NN, NN);
			inv     = mkgscmatrixd(NN, NN);

			print_gscmatrixd("dotransformMatd4Localdiag", transformMatd4Localdiag, DIMd4, DIMd4);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
				pose[mu][nu] = transformMatd4Localdiag[nu][mu];
				int selmu = rdmvecIndd4[mu] ;
				int selnu = rdmvecIndd4[nu] ;
				backup[mu][nu] = original[selmu][selnu] ;
				rdmLocaldiag[mu][nu] = zero;
				right[mu][nu] = zero;
				inv[mu][nu] = gsl_complex_conjugate(transformMatd4Localdiag[mu][nu]);
			}
			print_gscmatrixd("backup", backup, DIMd4, DIMd4);
			printf( "rdmd4diagsum in 'Localdiag'-basis : " ) ;
			diagsum = 0 ;
			for(int mu=0; mu<DIMd4; mu++) 
			{
				diagsum +=  GSL_REAL( backup[mu][mu]) ;
			}
			printf( "%19.16lf\n" , diagsum ) ;

			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				right[mu][nu]
					= gsl_complex_add(
							right[mu][nu],
							gsl_complex_mul(backup[mu][k], pose[k][nu])
							);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				rdmLocaldiag[mu][nu] =
					gsl_complex_add(
							rdmLocaldiag[mu][nu],
							gsl_complex_mul( inv[mu][k], right[k][nu] )
						       );

			//print_gscmatrixd("original", backup, NN, NN);
			print_gscmatrixd("rdmLocaldiag", rdmLocaldiag, NN, NN);

			printf( "rdmLocaldiag's d4diagsum : " ) ;
			diagsum = 0 ;
			for(int mu=0; mu<DIMd4; mu++)
			{
				diagsum +=  GSL_REAL( rdmLocaldiag[mu][mu]) ;
			}
			printf( "%19.16lf\n" , diagsum ) ;

			freegscmatrixd(right, NN);
			freegscmatrixd(pose, NN);
			freegscmatrixd(inv, NN);
			freegscmatrixd(backup, NN);

			return 0;
		}
		int dotransformevecd4( gsl_complex **original , gsl_complex **right, gsl_complex **transform ) { //original should be a matrix of all many-body states
			int mu, nu, k;
			gsl_complex **backup, **pose, zero = gsl_complex_rect(0,0) ;
			int NN = DIMd4 ;

			backup  = mkgscmatrixd(NN, NN);
			pose    = mkgscmatrixd(NN, NN);

			print_gscmatrixd("dotransformevecd4", transform, NN, NN);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
				pose[mu][nu] = gsl_complex_conjugate( transform[nu][mu] ) ;
				backup[mu][nu] = original[mu][nu] ;
				right[mu][nu] = zero;
			}
			print_gscmatrixd("originald4", backup, NN, NN);

			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				right[mu][nu]
					= gsl_complex_add(
							right[mu][nu],
							gsl_complex_mul(backup[mu][k], pose[k][nu])
							);

			print_gscmatrixd("rightd4", right, NN, NN);

			freegscmatrixd(pose, NN);
			freegscmatrixd(backup, NN);

			return 0;
		}
		int dotransform_testd4unit( gsl_complex **transform) {
			int mu, nu, k;
			int NN=DIMd4 ;
			gsl_complex **target, zero = gsl_complex_rect(0,0), **inv;
			inv	= mkgscmatrixd(NN, NN);
			target	= mkgscmatrixd(NN, NN);

			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
				target[mu][nu] = zero;
				inv[mu][nu] = gsl_complex_conjugate( transform[nu][mu] );
			}
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++){
				target[mu][nu]
					= gsl_complex_add(
							target[mu][nu],
							gsl_complex_mul(inv[mu][k], transform[k][nu])
							);
			}
			print_gscmatrixd("dotransform_testd4unit", target, NN, NN);
			for(mu=0; mu<NN; mu++){
				if( fabs(1-gsl_complex_abs2(target[mu][mu]))> 1e-5 ){
					printf("wrong diagonal part\n");
				}
				for(nu=mu+1; nu<NN; nu++){
					if( gsl_complex_abs2(target[mu][nu]) > 1e-5 ){
						printf("wrong off-diagonal part\n");
					}
				}
			}
			freegscmatrixd(inv, NN);
			freegscmatrixd(target, NN);

            return 0;
		}
		int dotransform_testd4( gsl_complex **original , gsl_complex **transform) {
			int mu, nu, k;
			gsl_complex **backup, **right, **pose, zero = gsl_complex_rect(0,0), **inv, **result;
			int NN = DIMd4 ;

			backup  = mkgscmatrixd(NN, NN);
			right   = mkgscmatrixd(NN, NN);
			pose    = mkgscmatrixd(NN, NN);
			inv     = mkgscmatrixd(NN, NN);
			result  = mkgscmatrixd(NN, NN);

			print_gscmatrixd("dotransform_testd4", transform, NN, NN);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
				pose[mu][nu] = transform[nu][mu];
				backup[mu][nu] = original[mu][nu] ;
				right[mu][nu] = zero;
				result[mu][nu] = zero;
				inv[mu][nu] = gsl_complex_conjugate(transform[mu][nu]);
			}
			print_gscmatrixd("backup", backup, NN, NN);

			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				right[mu][nu]
					= gsl_complex_add(
							right[mu][nu],
							gsl_complex_mul(backup[mu][k], pose[k][nu])
							);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				result[mu][nu] =
					gsl_complex_add(
							result[mu][nu],
							gsl_complex_mul( inv[mu][k], right[k][nu] )
						       );

			//print_gscmatrixd("original", backup, NN, NN);
			print_gscmatrixd("rdmd4_test", result, NN, NN);
			printf("rdmd4_test diag::\n") ;
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) printf("%19.16lf\t", GSL_REAL(result[mu][nu]) );
			printf("\n" ) ;

			freegscmatrixd(right, NN);
			freegscmatrixd(result, NN);
			freegscmatrixd(pose, NN);
			freegscmatrixd(inv, NN);
			freegscmatrixd(backup, NN);

			return 0;
		}
		int write_target( char pararesult[] , int ic , gsl_complex **target , char datname[]  , char writemode[] ){
		    return write_target( pararesult, ic , target , datname , writemode, "" );
        }
		int write_target( char pararesult[] , int ic , gsl_complex **target , char datname[]  , char writemode[], char prefix[] ){
			FILE *fw ; 
			char fname[1025] ;
		       	sprintf( fname, "%s/%s%s.dat", pararesult , prefix, datname ) ;
			fw = fopen( fname , writemode );
			fprintf( fw , "%d\t", ic ) ;
			int NN = DIMd4 ;
			for(int mu=0; mu<NN; mu++) for(int nu=0; nu<NN; nu++) fprintf( fw , "%19.16lf\t%19.16lf\t\t", GSL_REAL(target[mu][nu]), GSL_IMAG(target[mu][nu])   ) ;
			fprintf( fw , "\n") ;

			fclose(fw) ;
			return 0 ;
		}
		int write_targetd( char pararesult[] , int ic , double *target , char datname[]  , char writemode[] ){
		    return write_targetd( pararesult, ic , target , datname, writemode, "" );
        }
		int write_targetd( char pararesult[] , int ic , double *target , char datname[]  , char writemode[], char prefix[] ){
			FILE *fw ; 
			char fname[1025] ;
		       	sprintf( fname, "%s/%s%s.dat", pararesult , prefix, datname ) ;
			fw = fopen( fname , writemode );
			fprintf( fw , "%d\t", ic ) ;
			int NN = DIMd4 ;
			for(int mu=0; mu<NN; mu++) fprintf( fw , "%19.16lf\t", target[mu] ) ;
			fprintf( fw , "\n") ;

			fclose(fw) ;
			return 0 ;
		}

		double write_EE( char pararesult[] , int ic , double *rdmeval , char writemode[] ){
			return write_EE( pararesult, ic, rdmeval, writemode, "" );
		}
		double write_EE( char pararesult[] , int ic , double *rdmeval , char writemode[] , char prefix[] ){
			FILE *fw ; 
			char fname[1025] ;
		       	sprintf( fname, "%s/%sEEdeg.dat", pararesult, prefix ) ;
			fw = fopen( fname , writemode );
			fprintf( fw , "%d\t", ic ) ;
			double EE = 0 ; 
			for(int mu=0; mu<DIMd; mu++)  {
                double evalmu = rdmeval[mu];
                if( fabs(evalmu)<1e-15 ) evalmu += 1e-15;
				EE -= log( evalmu ) * evalmu ;
			}
			fprintf( fw , "%19.16lf\n", EE ) ;
			printf("EE : %19.16lf\n", EE  ) ;

			fclose(fw) ;
			return EE ;
		}
		int write_jjJd4( char pararesult[] , int ic , char writemode[] ){
			return write_jjJd4( pararesult, ic, writemode, "" );
		}
		int write_jjJd4( char pararesult[] , int ic , char writemode[], char prefix[] ){
			char fname[1025] ;
		       	sprintf( fname, "%s/%srdmjjJd4diag.dat", pararesult, prefix ) ;
			FILE *fw ; 
			fw = fopen( fname, writemode );
			fprintf( fw , "%d\t", ic ) ;
			for(int mu=0; mu<DIMd4; mu++)  {
				fprintf( fw , "%19.16lf\t", GSL_REAL( rdmjjJ[mu][mu]) ) ;
			}
			fprintf( fw , "\n") ;
			printf("rdmjjJd4diag : " ) ;
			for(int mu=0; mu<DIMd4; mu++)  {
				printf( "%19.16lf\t", GSL_REAL( rdmjjJ[mu][mu]) ) ;
			}
			printf( "\n") ;
			fclose(fw) ;
			return 0 ;
		}
		int write_LSJd4( char pararesult[] , int ic , char writemode[] ){
			return write_LSJd4( pararesult, ic, writemode, "" );
		}
		int write_LSJd4( char pararesult[] , int ic , char writemode[], char prefix[] ){
			char fname[1025] ;
		       	sprintf( fname, "%s/%srdmLSJd4diag.dat", pararesult, prefix ) ;
			FILE *fw ; 
			fw = fopen( fname, writemode );
			fprintf( fw , "%d\t", ic ) ;
			for(int mu=0; mu<DIMd4; mu++)  {
				fprintf( fw , "%19.16lf\t", GSL_REAL( rdmLSJ[mu][mu]) ) ;
			}
			fprintf( fw , "\n") ;
			printf("rdmLSJd4diag : " ) ;
			for(int mu=0; mu<DIMd4; mu++)  {
				printf( "%19.16lf\t", GSL_REAL( rdmLSJ[mu][mu]) ) ;
			}
			printf( "\n") ;
			fclose(fw) ;
			return 0 ;
		}
		int write_LSJcubicd4( char pararesult[] , int ic , char writemode[] ){
			return write_LSJcubicd4( pararesult, ic, writemode, "" );
		}
		int write_LSJcubicd4( char pararesult[] , int ic , char writemode[], char prefix[] ){
			FILE *fw ; 
			char fname[1025] ;
		       	sprintf( fname, "%s/%srdmLSJcubicd4diag.dat", pararesult, prefix ) ;
			fw = fopen( fname , writemode );
			fprintf( fw , "%d\t", ic ) ;
			for(int mu=0; mu<DIMd4; mu++)  {
				fprintf( fw , "%19.16lf\t", GSL_REAL( rdmLSJcubic[mu][mu]) ) ;
			}
			fprintf( fw , "\n") ;
			printf("rdmLSJcubicd4diag : " ) ;
			for(int mu=0; mu<DIMd4; mu++)  {
				printf( "%19.16lf\t", GSL_REAL( rdmLSJcubic[mu][mu]) ) ;
			}
			printf( "\n") ;
			fclose(fw) ;
			return 0 ;
		}
		int write_t2gd4( char pararesult[] , int ic , char writemode[] ){
			return write_t2gd4( pararesult, ic, writemode, "" );
		}
		int write_t2gd4( char pararesult[] , int ic , char writemode[], char prefix[] ){
			FILE *fw ; 
			char fname[1025] ;
		       	sprintf( fname, "%s/%srdmt2gd4diag.dat", pararesult, prefix ) ;
			fw = fopen( fname , writemode );
			fprintf( fw , "%d\t", ic ) ;
			for(int mu=0; mu<DIMd4; mu++)  {
				fprintf( fw , "%19.16lf\t", GSL_REAL( rdmt2g[mu][mu]) ) ;
			}
			fprintf( fw , "\n") ;
			printf("rdmt2gd4diag : " ) ;
			for(int mu=0; mu<DIMd4; mu++)  {
				printf( "%19.16lf\t", GSL_REAL( rdmt2g[mu][mu]) ) ;
			}
			printf( "\n") ;
			fclose(fw) ;
			return 0 ;
		}
		int write_Localdiagd4( char pararesult[] , int ic , char writemode[] ){
			return write_Localdiagd4( pararesult, ic, writemode, "" );
		}
		int write_Localdiagd4( char pararesult[] , int ic , char writemode[], char prefix[] ){
			FILE *fw ; 
			char fname[1025] ;
		       	sprintf( fname, "%s/%srdmLocaldiagd4diag.dat", pararesult, prefix ) ;
			fw = fopen( fname , writemode );
			fprintf( fw , "%d\t", ic ) ;
			for(int mu=0; mu<DIMd4; mu++)  {
				fprintf( fw , "%19.16lf\t", GSL_REAL( rdmLocaldiag[mu][mu]) ) ;
			}
			fprintf( fw , "\n") ;
			printf("rdmLocaldiagd4diag : " ) ;
			for(int mu=0; mu<DIMd4; mu++)  {
				printf( "%19.16lf\t", GSL_REAL( rdmLocaldiag[mu][mu]) ) ;
			}
			printf( "\n") ;
			fclose(fw) ;
			return 0 ;
		}
		int checkrdmHerimitian( gsl_complex **original ) {
			printf( "check rdm Herimitian : \n" ) ;
			int rdmvecIndd4[DIMd4] =  { 51, 53, 54, 60, 29, 45, 43, 39, 30, 57, 27, 23, 15, 46, 58 };
			for(int mu=0; mu<DIMd4; mu++) for(int nu=0; nu<DIMd4; nu++){
				int selmu = rdmvecIndd4[mu] ;
				int selnu = rdmvecIndd4[nu] ;
				if( gsl_complex_abs2( gsl_complex_sub( original[selmu][selnu] , gsl_complex_conjugate(original[selnu][selmu]) ) ) > 1e-5 ) 
					printf( "%d %d : %19.16lf +i %19.16lf , %19.16lf +i %19.16lf\n" , selmu,selnu,
							GSL_REAL( original[selmu][selnu] ) ,
							GSL_IMAG( original[selmu][selnu] ) ,
							GSL_REAL( original[selnu][selmu] ) ,
							GSL_IMAG( original[selnu][selmu] )  ) ;
			}
            return 0;
		}
		int obtain_jzmat_LSJd4( void ) {
			int JzLSJfactor[DIMd4] =  {0, 1, -1, 0, 2, -2, 1, -1, 0, 2, -2, 1, -1, 0, 0};
			int Jzjefffactor[DIMd4]=  {-1, 0, -2, 1, 1, 2, 1, -1, -1, 2, 0, -2, 0, 0, 0};
			int mu, nu, k;
			gsl_complex **backup, **right, **pose, zero = gsl_complex_rect(0,0), **inv;
			int NN = DIMd4 ;
			double diagsum = 0 ;
			double Jzd4    = 0 ;

			backup  = mkgscmatrixd(NN, NN);
			right   = mkgscmatrixd(NN, NN);
			pose    = mkgscmatrixd(NN, NN);
			inv     = mkgscmatrixd(NN, NN);
			jzopLSJ = mkgscmatrixd(NN, NN);

			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++){
				pose[mu][nu] = gsl_complex_conjugate(transformMatd4LSJ[nu][mu]);
				backup[mu][nu] = zero ;
				jzopLSJ[mu][nu] = zero;
				right[mu][nu] = zero;
				inv[mu][nu] = transformMatd4LSJ[mu][nu];
			}
			for(mu=0; mu<NN; mu++)  backup[mu][mu] = gsl_complex_rect(  Jzjefffactor[mu] , 0 ) ;
			print_gscmatrixd("backup", backup, DIMd4, DIMd4);
			//print_gscmatrixd("pose",   pose  , DIMd4, DIMd4);
			//print_gscmatrixd("inv",    inv   , DIMd4, DIMd4);
			//
			printf( "jzopjeffd4diag : " ) ;
			diagsum = 0 ;
			for(int mu=0; mu<DIMd4; mu++) 
			{
				double dum =  GSL_REAL( backup[mu][mu]) ;
				diagsum +=  dum ;
				printf( "%19.16lf " , dum ) ;
			}
			printf( "\n" ) ;
			printf( "jzopjeffd4diagsum : " ) ;
			printf( "%19.16lf\n" , diagsum ) ;

			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				right[mu][nu]
					= gsl_complex_add(
							right[mu][nu],
							gsl_complex_mul(backup[mu][k], pose[k][nu])
							);
			for(mu=0; mu<NN; mu++) for(nu=0; nu<NN; nu++) for(k=0; k<NN; k++)
				jzopLSJ[mu][nu] =
					gsl_complex_add(
							jzopLSJ[mu][nu],
							gsl_complex_mul( inv[mu][k], right[k][nu] )
						       );

			//print_gscmatrixd("original", backup, NN, NN);
			print_gscmatrixd("jzopLSJ", jzopLSJ, NN, NN);

			diagsum = 0 ;
			Jzd4    = 0 ;
			printf( "jzopLSJdiag : " ) ;
			for(int mu=0; mu<DIMd4; mu++)
			{
				double dum =  GSL_REAL( jzopLSJ[mu][mu]) ;
				diagsum +=  dum ;
				printf( "%19.16lf " , dum ) ;
			}
			printf( "\n" ) ;
			printf( "jzopLSJdiagsum : " ) ;
			printf( "%19.16lf\n" , diagsum ) ;

			freegscmatrixd(right, NN);
			freegscmatrixd(pose, NN);
			freegscmatrixd(inv, NN);
			freegscmatrixd(backup, NN);
			freegscmatrixd(jzopLSJ, NN);

			return 0;
		}
} ;
