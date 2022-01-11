double obtainMPS( gsl_complex **mps, gsl_complex *ground, typebasis **basis, int gnd, int gndblock ) {
	double mpssum = 0 ;
	#pragma omp parallel for default(shared) reduction(+:mpssum)
	for(int i = gnd; i<gnd+gndblock; i++){
		typebasis bin[2] ;
		bin[0] = basis[i][0];
		bin[1] = basis[i][1];
		//unsigned int  tstate = (bin[0]<<Ns) + bin[1]  ;
		unsigned int  bup = (bin[0]>>NC) ,
			      bdn = (bin[1]>>NC) ;
		unsigned int  bstate = (bup<<NB) + bdn ;
		unsigned int  cup = bin[0] - (bup<<NC) ,
			      cdn = bin[1] - (bdn<<NC) ;
		unsigned int  cstate = (cup<<NC) + cdn ;
		mps[bstate][cstate] = gsl_complex_add(  mps[bstate][cstate] , ground[i-gnd] ) ;
		mpssum += gsl_complex_abs2( mps[bstate][cstate] ) ;
	}
	return mpssum ;
}

double obtainMPS3orb( gsl_complex ****mps3orb, gsl_complex *ground, typebasis **basis, int gnd, int gndblock ) {
	double mpssum = 0 ;
	#pragma omp parallel for default(shared) reduction(+:mpssum)
	for(int i = gnd; i<gnd+gndblock; i++){
		typebasis bin[2] ;
		bin[0] = basis[i][0];
		bin[1] = basis[i][1];
		//unsigned int  tstate = (bin[0]<<Ns) + bin[1]  ;
		unsigned int  bup = (bin[0]>>NC) ,
			      bdn = (bin[1]>>NC) ;
		unsigned int  bstate = (bup<<NB) + bdn ;
		unsigned int  cup = bin[0] - (bup<<NC) ,
			      cdn = bin[1] - (bdn<<NC) ;
		unsigned int  cstate = (cup<<NC) + cdn ;
		int orbmu0 = One(cup,0) + 2*One(cdn,0);
		int orbmu1 = One(cup,1) + 2*One(cdn,1);
		int orbmu2 = One(cup,2) + 2*One(cdn,2);
		mps3orb[orbmu0][orbmu1][orbmu2][bstate] = gsl_complex_add(  mps3orb[orbmu0][orbmu1][orbmu2][bstate] , ground[i-gnd] ) ;
		mpssum += gsl_complex_abs2( mps3orb[orbmu0][orbmu1][orbmu2][bstate] ) ;
	}
	return mpssum ;
}

int obtainRDM( gsl_complex **rdm, gsl_complex **mps , int Powtnb, int Powtnc ) {
	#pragma omp parallel for shared(mps,rdm) collapse(2)
	for(int mu=0; mu<Powtnc; mu++) {
		for(int nu=0; nu<Powtnc; nu++){
			gsl_complex dum ;
			dum = gsl_complex_rect(0,0);
			for(int i=0; i<Powtnb ; i++) {
				dum = gsl_complex_add(
						dum ,
						gsl_complex_mul( mps[i][mu] ,
							gsl_complex_conjugate( mps[i][nu] ) 
							)
						) ;
			}
			GSL_REAL(rdm[mu][nu]) += GSL_REAL(dum) ;
			GSL_IMAG(rdm[mu][nu]) += GSL_IMAG(dum) ;
		}
	}
	return 0 ;
}

int obtainRDMorb( gsl_complex **rdm, gsl_complex ****mps3orb , int Powtnb , int orbind ) {
	if( orbind == 0 ) {
		#pragma omp parallel for default(shared) collapse(2)
		for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++) {
				gsl_complex dum ;
				dum = gsl_complex_rect(0,0);
				for(int orb1=0; orb1<4; orb1++) for(int orb2=0; orb2<4; orb2++) for(int i=0; i<Powtnb ; i++) {
					dum = gsl_complex_add(
							dum ,
							gsl_complex_mul( mps3orb[mu][orb1][orb2][i] ,
								gsl_complex_conjugate( mps3orb[mu][orb1][orb2][i] ) 
								)
							) ;
				}
				GSL_REAL(rdm[mu][nu]) += GSL_REAL(dum) ;
				GSL_IMAG(rdm[mu][nu]) += GSL_IMAG(dum) ;
			}
	}
	else if( orbind == 1 ) {
		#pragma omp parallel for default(shared) collapse(2)
		for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++) {
				gsl_complex dum ;
				dum = gsl_complex_rect(0,0);
				for(int orb1=0; orb1<4; orb1++) for(int orb2=0; orb2<4; orb2++) for(int i=0; i<Powtnb ; i++) {
					dum = gsl_complex_add(
							dum ,
							gsl_complex_mul( mps3orb[orb1][mu][orb2][i] ,
								gsl_complex_conjugate( mps3orb[orb1][mu][orb2][i] ) 
								)
							) ;
				}
				GSL_REAL(rdm[mu][nu]) += GSL_REAL(dum) ;
				GSL_IMAG(rdm[mu][nu]) += GSL_IMAG(dum) ;
			}
	}
	else if( orbind == 2 ) {
		#pragma omp parallel for default(shared) collapse(2)
		for(int mu=0; mu<4; mu++) for(int nu=0; nu<4; nu++) {
				gsl_complex dum ;
				dum = gsl_complex_rect(0,0);
				for(int orb1=0; orb1<4; orb1++) for(int orb2=0; orb2<4; orb2++) for(int i=0; i<Powtnb ; i++) {
					dum = gsl_complex_add(
							dum ,
							gsl_complex_mul( mps3orb[orb1][orb2][mu][i] ,
								gsl_complex_conjugate( mps3orb[orb1][orb2][mu][i] ) 
								)
							) ;
				}
				GSL_REAL(rdm[mu][nu]) += GSL_REAL(dum) ;
				GSL_IMAG(rdm[mu][nu]) += GSL_IMAG(dum) ;
			}
	}
	else {
		printf("ERROR : in obtainMPS3orb(), orbind is invalid.");
		exit(1);
	}
	return 0 ;
}

int calcNcRDMJzRDM( int ncrdm[], double jzrdm[] , int Powtnc ) {
	double jzjeff[tNC] = { 1.5, 0.5, -0.5, -1.5, 0.5, -0.5 } ;
	#pragma omp parallel for shared(ncrdm,jzrdm)
	for(int mu=0; mu<Powtnc; mu++) {
		unsigned int  cup = (mu>>NC) ,
			      cdn = mu - (cup<<NC) ;
		typebasis bin[2] ;
		bin[0] = cup ; 
		bin[1] = cdn ; 
		ncrdm[mu] = 0 ;
		jzrdm[mu] = 0 ;
		for(int nu=0 ; nu<tNC ; nu++ )  {
			if( One(bin[nu%2], nu/2) ) {
				ncrdm[mu] += 1 ;
				jzrdm[mu] += jzjeff[nu] ;
			}
		}
	}
	return 0 ;
}

gsl_complex calcOp( Op op1, gsl_complex *ground, typebasis **basis, int gnd, int gndblock , int *table) {
	gsl_complex opopval ;
	GSL_REAL(opopval) = 0;
	GSL_IMAG(opopval) = 0;
	double *dumre = new double[gndblock] ; 
	double *dumim = new double[gndblock] ; 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=1 ; 
		dumim[ig]=1 ; 
	}
	#pragma omp parallel for default(shared) 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=0 ; 
		dumim[ig]=0 ; 
		for( int bb=0 ; bb<op1.ndat ; bb++ )  {
			int mu1 = op1.cDat[bb].i ;
			int nu1 = op1.cDat[bb].j ;
			double op1val = op1.cDat[bb].dat ;
			if( One(  basis[i][nu1%2], nu1/2) ) {
				typebasis d1bin[2] ;
				d1bin[0]	= basis[i][0];
				d1bin[1]	= basis[i][1];
				int d1phase	= permu( d1bin, nu1/2, nu1%2);  d1bin[nu1%2] = Turnoff( d1bin[nu1%2], nu1/2 );
				if( Zero( d1bin[mu1%2], mu1/2) ) {
					typebasis op1bin[2] ;
					op1bin[0]	= d1bin[0];
					op1bin[1]	= d1bin[1];
					int op1phase	= d1phase * permu( op1bin, mu1/2, mu1%2);  op1bin[mu1%2] = Turnon(  op1bin[mu1%2], mu1/2 );
					int index = (op1bin[0]<<Ns) + op1bin[1];
					gsl_complex dum = gsl_complex_mul(
							gsl_complex_conjugate(ground[ig]),
							gsl_complex_mul_real( ground[table[index]-gnd], op1val*op1phase )
							) ;
					dumre[ig] += GSL_REAL( dum );
					dumim[ig] += GSL_IMAG( dum );
				}
			}
		}
	}
	for(int i=gnd; i<gnd+gndblock; i++) {
		int ig = i-gnd ; 
		GSL_REAL(opopval) += dumre[ig] ; 
		GSL_IMAG(opopval) += dumim[ig] ; 
	}
	delete []dumre ;
	delete []dumim ;
	return opopval ; 
}

gsl_complex calcOp( OpZ op1, gsl_complex *ground, typebasis **basis, int gnd, int gndblock , int *table) {
	gsl_complex opopval ;
	GSL_REAL(opopval) = 0;
	GSL_IMAG(opopval) = 0;
	double *dumre = new double[gndblock] ; 
	double *dumim = new double[gndblock] ; 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=1 ; 
		dumim[ig]=1 ; 
	}
	#pragma omp parallel for default(shared) 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=0 ; 
		dumim[ig]=0 ; 
		for( int bb=0 ; bb<op1.ndat ; bb++ )  {
			int mu1 = op1.cDat[bb].i ;
			int nu1 = op1.cDat[bb].j ;
			gsl_complex op1val = op1.cDat[bb].dat ;
			if( One(  basis[i][nu1%2], nu1/2) ) {
				typebasis d1bin[2] ;
				d1bin[0]	= basis[i][0];
				d1bin[1]	= basis[i][1];
				int d1phase	= permu( d1bin, nu1/2, nu1%2);  d1bin[nu1%2] = Turnoff( d1bin[nu1%2], nu1/2 );
				if( Zero( d1bin[mu1%2], mu1/2) ) {
					typebasis op1bin[2] ;
					op1bin[0]	= d1bin[0];
					op1bin[1]	= d1bin[1];
					int op1phase	= d1phase * permu( op1bin, mu1/2, mu1%2);  op1bin[mu1%2] = Turnon(  op1bin[mu1%2], mu1/2 );
					int index = (op1bin[0]<<Ns) + op1bin[1];
					gsl_complex dum = gsl_complex_mul(
							gsl_complex_conjugate(ground[ig]),
							gsl_complex_mul( ground[table[index]-gnd], gsl_complex_mul_real(op1val,op1phase) )
							) ;
					dumre[ig] += GSL_REAL( dum );
					dumim[ig] += GSL_IMAG( dum );
				}
			}
		}
	}
	for(int i=gnd; i<gnd+gndblock; i++) {
		int ig = i-gnd ; 
		GSL_REAL(opopval) += dumre[ig] ; 
		GSL_IMAG(opopval) += dumim[ig] ; 
	}
	delete []dumre ;
	delete []dumim ;
	return opopval ; 
}

double calcDiagOpOp( Op op2, Op op1, gsl_complex *ground, typebasis **basis, int gnd, int gndblock ) {
	//printf( "calcDiagOpOp:" ) ;
	double opopval=0 ;
	double *dumre = new double[gndblock] ; 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=0 ; 
	}
	#pragma omp parallel for default(shared) 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=0 ; 
		typebasis bin[2] ;
		bin[0] = basis[i][0];
		bin[1] = basis[i][1];
		for( int aa=0 ; aa<op1.ndat ; aa++ )  {
			int mu1 = op1.cDat[aa].i ;
			//int nu1 = op1.cDat[aa].j ;
			if( One( basis[i][mu1%2], mu1/2) ) {
				double op1val = op1.cDat[aa].dat ;
				for( int bb=0 ; bb<op2.ndat ; bb++ )  {
					int mu2 = op2.cDat[bb].i ;
					//int nu2 = op2.cDat[bb].j ;
					if( One( basis[i][mu2%2], mu2/2) ) {
						double op2val = op2.cDat[bb].dat ;
						double dum =  op1val*op2val*gsl_complex_abs2(ground[ig]) ;
						dumre[ig] += dum ;
					}
				}
			}
		}
	}
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		opopval += dumre[ig] ;
	}
	delete []dumre ;
	return opopval ;
}

gsl_complex calcOffdiagOpOp( Op op2, Op op1, gsl_complex *ground, typebasis **basis, int gnd, int gndblock , int *table) {
	gsl_complex opopval ;
	GSL_REAL(opopval) = 0;
	GSL_IMAG(opopval) = 0;
	double *dumre = new double[gndblock] ; 
	double *dumim = new double[gndblock] ; 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=1 ; 
		dumim[ig]=1 ; 
	}
	#pragma omp parallel for default(shared) 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=0 ; 
		dumim[ig]=0 ; 
		for( int bb=0 ; bb<op1.ndat ; bb++ )  {
			int mu1 = op1.cDat[bb].i ;
			int nu1 = op1.cDat[bb].j ;
			double op1val = op1.cDat[bb].dat ;
			if( Zero( basis[i][mu1%2], mu1/2) && One( basis[i][nu1%2], nu1/2) ) {
				typebasis bin[2] ;
				bin[0] = basis[i][0];
				bin[1] = basis[i][1];
				int phase=0 ;
				phase  = permu( bin, nu1/2, nu1%2);  bin[nu1%2] = Turnoff( bin[nu1%2], nu1/2 );
				phase *= permu( bin, mu1/2, mu1%2);  bin[mu1%2] = Turnon(  bin[mu1%2], mu1/2 );
				for( int aa=0 ; aa<op2.ndat ; aa++ )  {
					int mu2 = op2.cDat[aa].i ;
					int nu2 = op2.cDat[aa].j ;
					double op2val = op2.cDat[aa].dat ;
					typebasis op1bin[2] ;
					op1bin[0] = bin[0];
					op1bin[1] = bin[1];
					int op1phase = phase ;
					if( Zero( op1bin[mu2%2], mu2/2) && One( op1bin[nu2%2], nu2/2) ) {
						op1phase *= permu( op1bin, nu2/2, nu2%2);  op1bin[nu2%2] = Turnoff( op1bin[nu2%2], nu2/2 );
						op1phase *= permu( op1bin, mu2/2, mu2%2);  op1bin[mu2%2] = Turnon(  op1bin[mu2%2], mu2/2 );
						int index = (op1bin[0]<<Ns) + op1bin[1];
						gsl_complex dum = gsl_complex_mul(
								gsl_complex_conjugate(ground[ig]),
								gsl_complex_mul_real( ground[table[index]-gnd], op1val*op2val*op1phase )
								) ;
						dumre[ig] += GSL_REAL( dum );
						dumim[ig] += GSL_IMAG( dum );
					}
				}
			}
		}
	}
	for(int i=gnd; i<gnd+gndblock; i++) {
		int ig = i-gnd ; 
		GSL_REAL(opopval) += dumre[ig] ; 
		GSL_IMAG(opopval) += dumim[ig] ; 
	}
	delete []dumre ;
	delete []dumim ;
	return opopval ; 
}

gsl_complex calcOpOp( Op op2, Op op1, gsl_complex *ground, typebasis **basis, int gnd, int gndblock , int *table) {
	gsl_complex opopval ;
	GSL_REAL(opopval) = 0;
	GSL_IMAG(opopval) = 0;
	double *dumre = new double[gndblock] ; 
	double *dumim = new double[gndblock] ; 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=1 ; 
		dumim[ig]=1 ; 
	}
	#pragma omp parallel for default(shared) 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=0 ; 
		dumim[ig]=0 ; 
		for( int bb=0 ; bb<op1.ndat ; bb++ )  {
			int mu1 = op1.cDat[bb].i ;
			int nu1 = op1.cDat[bb].j ;
			double op1val = op1.cDat[bb].dat ;
			if( One(  basis[i][nu1%2], nu1/2) ) {
				typebasis d1bin[2] ;
				d1bin[0]	= basis[i][0];
				d1bin[1]	= basis[i][1];
				int d1phase	= permu( d1bin, nu1/2, nu1%2);  d1bin[nu1%2] = Turnoff( d1bin[nu1%2], nu1/2 );
				if( Zero( d1bin[mu1%2], mu1/2) ) {
					typebasis op1bin[2] ;
					op1bin[0]	= d1bin[0];
					op1bin[1]	= d1bin[1];
					int op1phase	= d1phase * permu( op1bin, mu1/2, mu1%2);  op1bin[mu1%2] = Turnon(  op1bin[mu1%2], mu1/2 );
					for( int aa=0 ; aa<op2.ndat ; aa++ )  {
						int mu2 = op2.cDat[aa].i ;
						int nu2 = op2.cDat[aa].j ;
						double op2val = op2.cDat[aa].dat ;
						if( One(  op1bin[nu2%2], nu2/2) ) {
							typebasis d2op1bin[2] ;
							d2op1bin[0]	= op1bin[0];
							d2op1bin[1]	= op1bin[1];
							int d2op1phase	= op1phase * permu( d2op1bin, nu2/2, nu2%2);  d2op1bin[nu2%2] = Turnoff( d2op1bin[nu2%2], nu2/2 );
							if( Zero( d2op1bin[mu2%2], mu2/2) ) {
								typebasis op2op1bin[2] ;
								op2op1bin[0]	= d2op1bin[0];
								op2op1bin[1]	= d2op1bin[1];
								int op2op1phase = d2op1phase * permu( op2op1bin, mu2/2, mu2%2);  op2op1bin[mu2%2] = Turnon(  op2op1bin[mu2%2], mu2/2 );
								int index = (op2op1bin[0]<<Ns) + op2op1bin[1];
								gsl_complex dum = gsl_complex_mul(
										gsl_complex_conjugate(ground[ig]),
										gsl_complex_mul_real( ground[table[index]-gnd], op1val*op2val*op2op1phase )
										) ;
								dumre[ig] += GSL_REAL( dum );
								dumim[ig] += GSL_IMAG( dum );
							}
						}
					}
				}
			}
		}
	}
	for(int i=gnd; i<gnd+gndblock; i++) {
		int ig = i-gnd ; 
		GSL_REAL(opopval) += dumre[ig] ; 
		GSL_IMAG(opopval) += dumim[ig] ; 
	}
	delete []dumre ;
	delete []dumim ;
	return opopval ; 
}

int write_ABC( char pararesult[] , int ic , char datname[] , char writemode[] , gsl_complex A, gsl_complex B, gsl_complex C, char prefix[]  ){
	FILE *fw ; 
	char fname[1025] ;
	sprintf( fname, "%s/%s%s.dat", pararesult , prefix, datname ) ;
	fw = fopen( fname , writemode );
	fprintf( fw , "%d\t", ic ) ;
	fprintf( fw , "%19.16lf\t%19.16lf\t", GSL_REAL(A) , GSL_IMAG(A)	);
	fprintf( fw , "%19.16lf\t%19.16lf\t", GSL_REAL(B) , GSL_IMAG(B)	);
	fprintf( fw , "%19.16lf\t%19.16lf\t", GSL_REAL(C) , GSL_IMAG(C)	);
	fprintf( fw , "\n") ;
	fclose( fw ) ;
	return 0 ;
}
int write_ABC( char pararesult[] , int ic , char datname[] , char writemode[] , gsl_complex A, gsl_complex B, gsl_complex C  ){
    return write_ABC( pararesult , ic , datname , writemode , A, B, C , "" );
}

int write_rdmorb( char pararesult[] , int ic , gsl_complex **target , char datname[]  , char writemode[], char prefix[] ){
	FILE *fw ; 
	char fname[1025] ;
       	sprintf( fname, "%s/%s%s.dat", pararesult , prefix, datname ) ;
	fw = fopen( fname , writemode );
	fprintf( fw , "%d\t", ic ) ;
	int NN = 4 ;
	for(int mu=0; mu<NN; mu++) for(int nu=0; nu<NN; nu++) fprintf( fw , "%19.16lf\t%19.16lf\t\t", GSL_REAL(target[mu][nu]), GSL_IMAG(target[mu][nu])   ) ;
	fprintf( fw , "\n") ;

	fclose(fw) ;
	return 0 ;
}
int write_rdmorb( char pararesult[] , int ic , gsl_complex **target , char datname[]  , char writemode[] ){
    return write_rdmorb( pararesult , ic , target , datname  , writemode, "" );
}

int write_rdmorbdiag( char pararesult[] , int ic , gsl_complex **target , char datname[]  , char writemode[], char prefix[] ){
	FILE *fw ; 
	char fname[1025] ;
       	sprintf( fname, "%s/%s%s.dat", pararesult , prefix, datname ) ;
	fw = fopen( fname , writemode );
	fprintf( fw , "%d\t", ic ) ;
	int NN = 4 ;
	for(int mu=0; mu<NN; mu++) fprintf( fw , "%19.16lf\t%19.16lf\t\t", GSL_REAL(target[mu][mu]), GSL_IMAG(target[mu][mu])   ) ;
	fprintf( fw , "\n") ;

	fclose(fw) ;
	return 0 ;
}
int write_rdmorbdiag( char pararesult[] , int ic , gsl_complex **target , char datname[]  , char writemode[] ){
    return write_rdmorbdiag( pararesult , ic , target , datname  , writemode , "" );
}
