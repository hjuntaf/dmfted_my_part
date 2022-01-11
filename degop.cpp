extern int beta_real;

int gscmatrix_mulrealconst( gsl_complex ** mat, double C , int dim ) {
	for( int d1=0 ; d1<dim ; d1++ ){
		for( int d2=0 ; d2<dim ; d2++ ){
			mat[d1][d2] = gsl_complex_mul_real( mat[d1][d2] , C )  ;
		}
	}
	return 0 ;
}
int gscmatrix_mulimagconst( gsl_complex ** mat, double C , int dim ) {
	for( int d1=0 ; d1<dim ; d1++ ){
		for( int d2=0 ; d2<dim ; d2++ ){
			mat[d1][d2] = gsl_complex_mul_imag( mat[d1][d2] , C )  ;
		}
	}
	return 0 ;
}

gsl_complex calcOp( Op op1, gsl_complex *ground1, int gnd1, int gndblock1, gsl_complex *ground2, int gnd2, int gndblock2, typebasis **basis, int *table) {
	gsl_complex opval ;
	GSL_REAL(opval) = 0;
	GSL_IMAG(opval) = 0;
	if( gnd1 != gnd2 )		{ }//printf("Warnning :: 'gnd' from ground states are different.\n") ;	}
	if( gndblock1 != gndblock2 )	{ }//printf("Warnning :: 'gndblock' from ground states are different.\n") ;}
	int gnd		= gnd2 ;
	int gndblock	= gndblock2 ;
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
					int tindexgnd1 = table[ (op1bin[0]<<Ns) + op1bin[1] ] - gnd1;
					gsl_complex dum ;
					if( (-1<tindexgnd1) && (tindexgnd1<gndblock1)  )
						dum = gsl_complex_mul(
								gsl_complex_conjugate(ground1[tindexgnd1]),
								gsl_complex_mul_real( ground2[ig], op1val*op1phase )
								) ;
					else { GSL_REAL(dum)=0.; GSL_IMAG(dum)=0.; }
					dumre[ig] += GSL_REAL( dum );
					dumim[ig] += GSL_IMAG( dum );
				}
			}
		}
	}
	for(int i=gnd; i<gnd+gndblock; i++) {
		int ig = i-gnd ; 
		GSL_REAL(opval) += dumre[ig] ; 
		GSL_IMAG(opval) += dumim[ig] ; 
	}
	delete []dumre ;
	delete []dumim ;
	return opval ; 
}

gsl_complex calcOp( OpZ op1, gsl_complex *ground1, int gnd1, int gndblock1, gsl_complex *ground2, int gnd2, int gndblock2, typebasis **basis, int *table) {
	gsl_complex opval ;
	GSL_REAL(opval) = 0;
	GSL_IMAG(opval) = 0;
	if( gnd1 != gnd2 )		{ }//printf("Warnning :: 'gnd' from ground states are different.\n") ;		}
	if( gndblock1 != gndblock2 )	{ }//printf("Warnning :: 'gndblock' from ground states are different.\n") ;	}
	int gnd		= gnd2 ;
	int gndblock	= gndblock2 ;
	double *dumre = new double[gndblock] ; 
	double *dumim = new double[gndblock] ; 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=0. ; 
		dumim[ig]=0. ; 
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
					int tindexgnd1 = table[ (op1bin[0]<<Ns) + op1bin[1] ]-gnd1 ;
					gsl_complex dum ;
					if( (-1<tindexgnd1) && (tindexgnd1<gndblock1)  )
						dum = gsl_complex_mul(
								gsl_complex_conjugate(ground1[tindexgnd1]),
								gsl_complex_mul( ground2[ig], gsl_complex_mul_real(op1val,op1phase) )
								) ;
					else { GSL_REAL(dum)=0.; GSL_IMAG(dum)=0.; }
					dumre[ig] += GSL_REAL( dum );
					dumim[ig] += GSL_IMAG( dum );
				}
			}
		}
	}
	for(int i=gnd; i<gnd+gndblock; i++) {
		int ig = i-gnd ; 
		GSL_REAL(opval) += dumre[ig] ; 
		GSL_IMAG(opval) += dumim[ig] ; 
	}
	delete []dumre ;
	delete []dumim ;
	return opval ; 
}

gsl_complex calcOpOp( Op op2, Op op1, gsl_complex *ground1, int gnd1, int gndblock1 , gsl_complex *ground2, int gnd2, int gndblock2 , typebasis **basis, int *table) {
	gsl_complex opopval ;
	GSL_REAL(opopval) = 0;
	GSL_IMAG(opopval) = 0;
	//Calculating  <g1|op1op2|g2>
	if( gnd1 != gnd2 )		{ }//printf("Warnning :: 'gnd' from ground states are different.\n") ;	}
	if( gndblock1 != gndblock2 )	{ }//printf("Warnning :: 'gndblock' from ground states are different.\n") ;	}
	int gnd		= gnd2 ;
	int gndblock	= gndblock2 ;
	double *dumre = new double[gndblock] ; 
	double *dumim = new double[gndblock] ; 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=0 ; 
		dumim[ig]=0 ; 
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
								int tindexgnd1 = table[ (op2op1bin[0]<<Ns) + op2op1bin[1] ] - gnd1;
								gsl_complex dum ;
								if( (-1<tindexgnd1) && (tindexgnd1<gndblock1)  )
									dum = gsl_complex_mul(
											gsl_complex_conjugate(ground1[tindexgnd1]),
											gsl_complex_mul_real( ground2[ig], op1val*op2val*op2op1phase )
											) ;
								else { GSL_REAL(dum)=0.; GSL_IMAG(dum)=0.; }
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

gsl_complex calcOpOp( OpZ op2, OpZ op1, gsl_complex *ground1, int gnd1, int gndblock1 , gsl_complex *ground2, int gnd2, int gndblock2 , typebasis **basis, int *table) {
	gsl_complex opopval ;
	GSL_REAL(opopval) = 0;
	GSL_IMAG(opopval) = 0;
	//Calculating  <g1|op1op2|g2>
	if( gnd1 != gnd2 )		{ }//printf("Warnning :: 'gnd' from ground states are different.\n") ;	}
	if( gndblock1 != gndblock2 )	{ }//printf("Warnning :: 'gndblock' from ground states are different.\n") ;	}
	int gnd		= gnd2 ;
	int gndblock	= gndblock2 ;
	double *dumre = new double[gndblock] ; 
	double *dumim = new double[gndblock] ; 
	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 
		dumre[ig]=0 ; 
		dumim[ig]=0 ; 
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
					for( int aa=0 ; aa<op2.ndat ; aa++ )  {
						int mu2 = op2.cDat[aa].i ;
						int nu2 = op2.cDat[aa].j ;
						gsl_complex op2val = op2.cDat[aa].dat ;
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
								int tindexgnd1 = table[ (op2op1bin[0]<<Ns) + op2op1bin[1] ] - gnd1;
								gsl_complex dum ;
								if( (-1<tindexgnd1) && (tindexgnd1<gndblock1)  )
									dum = gsl_complex_mul(
											gsl_complex_conjugate(ground1[tindexgnd1]),
											gsl_complex_mul( ground2[ig], gsl_complex_mul_real(gsl_complex_mul(op1val,op2val),op2op1phase) )
											) ;
								else { GSL_REAL(dum)=0.; GSL_IMAG(dum)=0.; }
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

gsl_complex **calcOpOpMat( Op op2, Op op1, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opopmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ) {
		gsl_complex *ground1	= groundarr[d1] ;
		int gnd1		= gndarr[d1] ;
		int gndblock1		= gndblockarr[d1] ;
		opopmat[d1][d1] = calcOpOp( op2, op1, ground1, basis, gnd1, gndblock1, table) ;
		for( int d2=d1+1 ; d2<degeneracy ; d2++ ) {
			gsl_complex *ground2	= groundarr[d2] ;
			int gnd2		= gndarr[d2] ;
			int gndblock2		= gndblockarr[d2] ;
			opopmat[d1][d2] = calcOpOp( op2, op1, ground1, gnd1, gndblock1, ground2, gnd2, gndblock2, basis, table) ;
		}
		for( int d2=0 ; d2<d1 ; d2++ ) {
			gsl_complex *ground2	= groundarr[d2] ;
			int gnd2		= gndarr[d2] ;
			int gndblock2		= gndblockarr[d2] ;
			opopmat[d1][d2] = calcOpOp( op2, op1, ground1, gnd1, gndblock1, ground2, gnd2, gndblock2, basis, table) ;
		}
	}
	return opopmat ;
}

gsl_complex **calcOpOpMat( OpZ op2, OpZ op1, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opopmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ) {
		gsl_complex *ground1	= groundarr[d1] ;
		int gnd1		= gndarr[d1] ;
		int gndblock1		= gndblockarr[d1] ;
		//opopmat[d1][d1] = calcOpOp( op2, op1, ground1, basis, gnd1, gndblock1, table) ;
        opopmat[d1][d1] = calcOpOp( op2, op1, ground1, gnd1, gndblock1, ground1, gnd1, gndblock1, basis, table) ;
		for( int d2=d1+1 ; d2<degeneracy ; d2++ ) {
			gsl_complex *ground2	= groundarr[d2] ;
			int gnd2		= gndarr[d2] ;
			int gndblock2		= gndblockarr[d2] ;
			opopmat[d1][d2] = calcOpOp( op2, op1, ground1, gnd1, gndblock1, ground2, gnd2, gndblock2, basis, table) ;
		}
		for( int d2=0 ; d2<d1 ; d2++ ) {
			gsl_complex *ground2	= groundarr[d2] ;
			int gnd2		= gndarr[d2] ;
			int gndblock2		= gndblockarr[d2] ;
			opopmat[d1][d2] = calcOpOp( op2, op1, ground1, gnd1, gndblock1, ground2, gnd2, gndblock2, basis, table) ;
		}
	}
	return opopmat ;
}

gsl_complex **calcOpMat( Op op1, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ) {
		gsl_complex *ground1	= groundarr[d1] ;
		int gnd1		= gndarr[d1] ;
		int gndblock1		= gndblockarr[d1] ;
		opmat[d1][d1] = calcOp( op1, ground1, basis, gnd1, gndblock1, table) ;
		for( int d2=d1+1 ; d2<degeneracy ; d2++ ) {
			gsl_complex *ground2	= groundarr[d2] ;
			int gnd2		= gndarr[d2] ;
			int gndblock2		= gndblockarr[d2] ;
			opmat[d1][d2] = calcOp( op1, ground1, gnd1, gndblock1, ground2, gnd2, gndblock2, basis, table) ;
			opmat[d2][d1] = calcOp( op1, ground2, gnd2, gndblock2, ground1, gnd1, gndblock1, basis, table) ;
		}
	}
	return opmat ;
}

gsl_complex **calcOpMat( OpZ op1, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ) {
		gsl_complex *ground1	= groundarr[d1] ;
		int gnd1		= gndarr[d1] ;
		int gndblock1		= gndblockarr[d1] ;
		opmat[d1][d1] = calcOp( op1, ground1, basis, gnd1, gndblock1, table) ;
		for( int d2=d1+1 ; d2<degeneracy ; d2++ ) {
			gsl_complex *ground2	= groundarr[d2] ;
			int gnd2		= gndarr[d2] ;
			int gndblock2		= gndblockarr[d2] ;
			opmat[d1][d2] = calcOp( op1, ground1, gnd1, gndblock1, ground2, gnd2, gndblock2, basis, table) ;
			opmat[d2][d1] = calcOp( op1, ground2, gnd2, gndblock2, ground1, gnd1, gndblock1, basis, table) ;
		}
	}
	return opmat ;
}

int write_gscmat( char pararesult[] , int ic , char datname[] , char writemode[] , gsl_complex **A, int dim  ){
	FILE *fw ; 
	char fname[1025] ;
	sprintf( fname, "%s/%s.dat", pararesult , datname ) ;
	fw = fopen( fname , writemode );
	fprintf( fw , "%d\t", ic ) ;
	for( int d1=0 ; d1<dim ; d1++ ) for( int d2=0 ; d2<dim ; d2++ ) fprintf( fw , "%19.16lf\t%19.16lf\t", GSL_REAL(A[d1][d2]) , GSL_IMAG(A[d1][d2])	);
	fprintf( fw , "\n") ;
	fclose( fw ) ;
	return 0 ;
}

gsl_complex **calcVopSquareMat( Vecop vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **jzjzmat	= calcOpOpMat(  vop.zj, vop.zj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **jpjmmat	= calcOpOpMat(  vop.pj, vop.mj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **jmjpmat	= calcOpOpMat(  vop.mj, vop.pj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = GSL_REAL(jzjzmat[d1][d2]) + 0.5*( GSL_REAL(jpjmmat[d1][d2]) + GSL_REAL(jmjpmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = GSL_IMAG(jzjzmat[d1][d2]) + 0.5*( GSL_IMAG(jpjmmat[d1][d2]) + GSL_IMAG(jmjpmat[d1][d2]) ) ;
		}
	}
	return opmat ;
}

gsl_complex **calcVopSquareTransverseMat( Vecop vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **jpjmmat	= calcOpOpMat(  vop.pj, vop.mj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **jmjpmat	= calcOpOpMat(  vop.mj, vop.pj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = 0.5*( GSL_REAL(jpjmmat[d1][d2]) + GSL_REAL(jmjpmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = 0.5*( GSL_IMAG(jpjmmat[d1][d2]) + GSL_IMAG(jmjpmat[d1][d2]) ) ;
		}
	}
	return opmat ;
}

gsl_complex **calcVopSquareMat( VecopZ vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **jzjzmat	= calcOpOpMat(  vop.zj, vop.zj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **jpjmmat	= calcOpOpMat(  vop.pj, vop.mj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **jmjpmat	= calcOpOpMat(  vop.mj, vop.pj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = GSL_REAL(jzjzmat[d1][d2]) + 0.5*( GSL_REAL(jpjmmat[d1][d2]) + GSL_REAL(jmjpmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = GSL_IMAG(jzjzmat[d1][d2]) + 0.5*( GSL_IMAG(jpjmmat[d1][d2]) + GSL_IMAG(jmjpmat[d1][d2]) ) ;
		}
	}
	return opmat ;
}

gsl_complex **calcVopSquareTransverseMat( VecopZ vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **jpjmmat	= calcOpOpMat(  vop.pj, vop.mj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **jmjpmat	= calcOpOpMat(  vop.mj, vop.pj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = 0.5*( GSL_REAL(jpjmmat[d1][d2]) + GSL_REAL(jmjpmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = 0.5*( GSL_IMAG(jpjmmat[d1][d2]) + GSL_IMAG(jmjpmat[d1][d2]) ) ;
		}
	}
	return opmat ;
}

gsl_complex **calcVopXMat( Vecop vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **jpmat	= calcOpMat(  vop.pj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **jmmat	= calcOpMat(  vop.mj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = 0.5*( GSL_REAL(jpmat[d1][d2]) + GSL_REAL(jmmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = 0.5*( GSL_IMAG(jpmat[d1][d2]) + GSL_IMAG(jmmat[d1][d2]) ) ;
		}
	}
	return opmat ;
}

gsl_complex **calcVopXMat( VecopZ vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **jpmat	= calcOpMat(  vop.pj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **jmmat	= calcOpMat(  vop.mj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = 0.5*( GSL_REAL(jpmat[d1][d2]) + GSL_REAL(jmmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = 0.5*( GSL_IMAG(jpmat[d1][d2]) + GSL_IMAG(jmmat[d1][d2]) ) ;
		}
	}
	return opmat ;
}

gsl_complex **calcVopYMat( Vecop vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **jpmat	= calcOpMat(  vop.pj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **jmmat	= calcOpMat(  vop.mj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) =  0.5*( GSL_IMAG(jpmat[d1][d2]) - GSL_IMAG(jmmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = -0.5*( GSL_REAL(jpmat[d1][d2]) - GSL_REAL(jmmat[d1][d2]) ) ;
		}
	}
	return opmat ;
}

gsl_complex **calcVopYMat( VecopZ vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **jpmat	= calcOpMat(  vop.pj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **jmmat	= calcOpMat(  vop.mj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) =  0.5*( GSL_IMAG(jpmat[d1][d2]) - GSL_IMAG(jmmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = -0.5*( GSL_REAL(jpmat[d1][d2]) - GSL_REAL(jmmat[d1][d2]) ) ;
		}
	}
	return opmat ;
}

gsl_complex **calcOpOpCommuteMat( Op op1, Op op2, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **op1op2mat	= calcOpOpMat(  op1, op2, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **op2op1mat	= calcOpOpMat(  op2, op1, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = ( GSL_REAL(op1op2mat[d1][d2]) + GSL_REAL(op2op1mat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = ( GSL_IMAG(op1op2mat[d1][d2]) + GSL_IMAG(op2op1mat[d1][d2]) ) ;
		}
	}
	freegscmatrixd( op1op2mat, degeneracy )  ;
	freegscmatrixd( op2op1mat, degeneracy )  ;
	return opmat ;
}

gsl_complex **calcOpOpAnticommuteMat( Op op1, Op op2, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **op1op2mat	= calcOpOpMat(  op1, op2, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **op2op1mat	= calcOpOpMat(  op2, op1, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = ( GSL_REAL(op1op2mat[d1][d2]) - GSL_REAL(op2op1mat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = ( GSL_IMAG(op1op2mat[d1][d2]) - GSL_IMAG(op2op1mat[d1][d2]) ) ;
		}
	}
	freegscmatrixd( op1op2mat, degeneracy )  ;
	freegscmatrixd( op2op1mat, degeneracy )  ;
	return opmat ;
}

gsl_complex **calcVopQx2y2Mat( Vecop vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **VpVpmat	= calcOpOpMat(  vop.pj, vop.pj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **VmVmmat	= calcOpOpMat(  vop.mj, vop.mj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = ( GSL_REAL(VpVpmat[d1][d2]) + GSL_REAL(VmVmmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = ( GSL_IMAG(VpVpmat[d1][d2]) + GSL_IMAG(VmVmmat[d1][d2]) ) ;
		}
	}
	freegscmatrixd(VpVpmat, degeneracy) ;
	freegscmatrixd(VmVmmat, degeneracy) ;

	gscmatrix_mulrealconst( opmat, 0.5, degeneracy ) ;
	return opmat ;
}

gsl_complex **calcVopQz2r2Mat( Vecop vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **VzVzmat	= calcOpOpMat(		vop.zj, vop.zj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ; gscmatrix_mulrealconst( VzVzmat,   2. /sqrt(3), degeneracy ) ;
	gsl_complex **VpVmCmat	= calcOpOpCommuteMat(	vop.pj, vop.mj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ; gscmatrix_mulrealconst( VpVmCmat, -0.5/sqrt(3), degeneracy ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = ( GSL_REAL(VzVzmat[d1][d2]) + GSL_REAL(VpVmCmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = ( GSL_IMAG(VzVzmat[d1][d2]) + GSL_IMAG(VpVmCmat[d1][d2]) ) ;
		}
	}
	freegscmatrixd(VzVzmat, degeneracy) ;
	freegscmatrixd(VpVmCmat, degeneracy) ;

	return opmat ;
}

gsl_complex **calcVopQxyMat( Vecop vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **VpVpmat	= calcOpOpMat(  vop.pj, vop.pj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **VmVmmat	= calcOpOpMat(  vop.mj, vop.mj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = ( GSL_REAL(VpVpmat[d1][d2]) - GSL_REAL(VmVmmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = ( GSL_IMAG(VpVpmat[d1][d2]) - GSL_IMAG(VmVmmat[d1][d2]) ) ;
		}
	}
	freegscmatrixd(VpVpmat, degeneracy) ;
	freegscmatrixd(VmVmmat, degeneracy) ;

	gscmatrix_mulimagconst( opmat, -0.5, degeneracy ) ;
	return opmat ;
}

gsl_complex **calcVopQyzMat( Vecop vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **VpVzCmat	= calcOpOpCommuteMat(	vop.pj, vop.zj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **VmVzCmat	= calcOpOpCommuteMat(	vop.mj, vop.zj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = ( GSL_REAL(VpVzCmat[d1][d2]) - GSL_REAL(VmVzCmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = ( GSL_IMAG(VpVzCmat[d1][d2]) - GSL_IMAG(VmVzCmat[d1][d2]) ) ;
		}
	}
	freegscmatrixd(VpVzCmat, degeneracy) ;
	freegscmatrixd(VmVzCmat, degeneracy) ;

	gscmatrix_mulimagconst( opmat, -0.5, degeneracy ) ;
	return opmat ;
}

gsl_complex **calcVopQzxMat( Vecop vop, gsl_complex **groundarr, int *gndarr, int *gndblockarr , int degeneracy, typebasis **basis, int *table) {
	gsl_complex **opmat = mkgscmatrixd( degeneracy, degeneracy ) ;
	gsl_complex **VpVzCmat	= calcOpOpCommuteMat(	vop.pj, vop.zj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	gsl_complex **VmVzCmat	= calcOpOpCommuteMat(	vop.mj, vop.zj, groundarr, gndarr, gndblockarr, degeneracy, basis, table ) ;
	for( int d1=0 ; d1<degeneracy ; d1++ ){
		for( int d2=0 ; d2<degeneracy ; d2++ ){
			GSL_REAL(opmat[d1][d2]) = ( GSL_REAL(VpVzCmat[d1][d2]) + GSL_REAL(VmVzCmat[d1][d2]) ) ;
			GSL_IMAG(opmat[d1][d2]) = ( GSL_IMAG(VpVzCmat[d1][d2]) + GSL_IMAG(VmVzCmat[d1][d2]) ) ;
		}
	}
	freegscmatrixd(VpVzCmat, degeneracy) ;
	freegscmatrixd(VmVzCmat, degeneracy) ;

	gscmatrix_mulrealconst( opmat, 0.5, degeneracy ) ;
	return opmat ;
}

int write_gscmateig( char pararesult[] , int ic , char datname[] , char writemode[] , gsl_complex **A, int dim  ){
	FILE *fw ; 
	char fname[1025] ;
	sprintf( fname, "%s/%s.dat", pararesult , datname ) ;
	fw = fopen( fname , writemode );

	gsl_complex **opevec	= mkgscmatrixd(dim,dim) ;
	double *opeval		= mkvectord(dim) ;
	eigen_lapack( A, opeval, opevec, dim ,0) ;

	fprintf( fw , "%d\t", ic ) ;
	for( int d1=0 ; d1<dim ; d1++ ) fprintf( fw , "%19.16lf\t", opeval[d1] ) ;
	fprintf( fw , "\n") ;
	fclose( fw ) ;

	sprintf( fname, "%s/%sevec.dat", pararesult , datname ) ;
	fw = fopen( fname , writemode );
	fprintf( fw , "%d\t", ic ) ;
	for( int d1=0 ; d1<dim ; d1++ ) for( int d2=0 ; d2<dim ; d2++ ) fprintf( fw , "%19.16lf\t%19.16lf\t", GSL_REAL(opevec[d1][d2]), GSL_IMAG(opevec[d1][d2]) ) ;
	fprintf( fw , "\n") ;
	fclose( fw ) ;

	return 0 ;
}
int write_doublearr( char pararesult[] , int ic , char datname[] , char writemode[] , double *A, int dim  ){
	FILE *fw ; 
	char fname[1025] ;
	sprintf( fname, "%s/%s.dat", pararesult , datname ) ;
	fw = fopen( fname , writemode );

	fprintf( fw , "%d\t", ic ) ;
	for( int d1=0 ; d1<dim ; d1++ ) fprintf( fw , "%19.16lf\t", A[d1] ) ;
	fprintf( fw , "\n") ;
	fclose( fw ) ;

	return 0;
}

int write_gscmateigTemp( char pararesult[] , int ic , char datname[] , char writemode[] , gsl_complex **A, int dim  , double *weight ){
	FILE *fw ; 
	char fname[1025] ;
	sprintf( fname, "%s/%s.dat", pararesult , datname ) ;
	fw = fopen( fname , writemode );

	gsl_complex **opevec	= mkgscmatrixd(dim,dim) ;
	double *opeval		= mkvectord(dim) ;
	eigen_lapack( A, opeval, opevec, dim ,0) ;

	gsl_complex **right	= mkgscmatrixd(dim,dim) ;
	gsl_complex **result	= mkgscmatrixd(dim,dim) ;

        for( int mu=0; mu<dim; mu++) for( int nu=0; nu<dim; nu++) { right[mu][nu] = zero;  result[mu][nu] = zero; }
        for( int mu=0; mu<dim; mu++) for( int k=0; k<dim; k++) for( int nu=0; nu<dim; nu++) right[mu][nu] = gsl_complex_add( right[mu][nu], gsl_complex_mul( A[mu][k], opevec[nu][k] )); 
        for( int mu=0; mu<dim; mu++) for( int nu=0; nu<dim; nu++) opevec[mu][nu] = gsl_complex_mul_real( gsl_complex_conjugate(opevec[mu][nu]) , weight[nu] );
        for( int mu=0; mu<dim; mu++) for( int k=0; k<dim; k++) for( int nu=0; nu<dim; nu++) result[mu][nu] = gsl_complex_add( result[mu][nu], gsl_complex_mul( opevec[mu][k], right[k][nu] )); 

	fprintf( fw , "%d\t", ic ) ;
	for( int d1=0 ; d1<dim ; d1++ ) fprintf( fw , "%19.16lf\t", GSL_REAL(result[d1][d1]) ) ;
	fprintf( fw , "\n") ;
	fclose( fw ) ;

	free(opeval);
	freegscmatrixd(opevec,dim);
	freegscmatrixd(right,dim);
	freegscmatrixd(result,dim);

	return 0 ;
}

int write_longtimecorrFromGscmatTemp( char pararesult[] , int ic , char datname[] , char writemode[] , gsl_complex **A, int dim  , double *weight, double *energy ){
	FILE *fw ; 
	char fname[1025] ;
	sprintf( fname, "%s/%s.dat", pararesult , datname ) ;
	fw = fopen( fname , writemode );

	double longtimecorrbeta = 0.;
        for( int mu=0; mu<dim; mu++) for( int nu=0; nu<dim; nu++) {
		double transitweight = exp( -beta_real *(energy[nu]-energy[mu])/2. );
		longtimecorrbeta += gsl_complex_abs2( A[mu][nu] ) * weight[mu] * transitweight;
	}

	fprintf( fw , "%d\t", ic ) ;
	fprintf( fw , "%19.16lf\t", longtimecorrbeta );
	fprintf( fw , "\n") ;
	fclose( fw ) ;

	return 0 ;
}

int write_doubleMat( char pararesult[] , int ic , char datname[] , char writemode[] , double **A, int dim1, int dim2  ) {
	FILE *fw ; 
	char fname[1025] ;
	sprintf( fname, "%s/%s.dat", pararesult , datname ) ;
	fw = fopen( fname , writemode );

	fprintf( fw , "%d\t", ic ) ;
	for( int d1=0 ; d1<dim1 ; d1++ ) for( int d2=0 ; d2<dim2 ; d2++ ) fprintf( fw , "%19.16lf\t", A[d1][d2] ) ;
	fprintf( fw , "\n") ;
	fclose( fw ) ;
	return 0 ;
}

void gscmatrixTrReCartesian( double *Vec, gsl_complex **Ax, gsl_complex **Ay, gsl_complex **Az, int dim  , double *weight ){
        for( int mu=0; mu<3; mu++) Vec[mu] = 0.;
        for( int mu=0; mu<dim; mu++) for( int nu=0; nu<dim; nu++) {
		Vec[0] += GSL_REAL(Ax[mu][mu]) * weight[nu];
		Vec[1] += GSL_REAL(Ay[mu][mu]) * weight[nu];
		Vec[2] += GSL_REAL(Az[mu][mu]) * weight[nu];
	}
}


int read_eigvec(char *save_directory,  int degeneracy, gsl_complex **eigvec , char *LocationTeigen){
	char para[1024];
	sprintf(para, "%s/%s", save_directory, LocationTeigen );
	FILE *fd;
	fd = fopen(para, "r");
	if( fd == NULL ){
		printf("%s does not exist!\n", para);
		return 1;
	}
	printf("Eigvec \"%s\" :\n", para );
	for( int mu=0; mu<degeneracy; mu++ ) {
		for( int nu=0; nu<degeneracy; nu++ ) {
			double dumre, dumim ;
			fscanf(fd, "%lf", &dumre );
			fscanf(fd, "%lf", &dumim );
			eigvec[mu][nu] = gsl_complex_rect( dumre, dumim );
			printf("%19.16f %19.16f\t", GSL_REAL(eigvec[mu][nu]), GSL_IMAG(eigvec[mu][nu]) );
		}
		printf("\n");
	}
	fclose(fd);
	fflush(stdout);
	return 0;
}

int print_eigvec(int degeneracy, gsl_complex **eigvec ){
	printf("Eigvec : \n") ;
	for( int mu=0; mu<degeneracy; mu++ ) {
		for( int nu=0; nu<degeneracy; nu++ ) {
			printf("%19.16f %19.16f\t", GSL_REAL(eigvec[mu][nu]), GSL_IMAG(eigvec[mu][nu]) );
		}
		printf("\n");
	}
	fflush(stdout);
	return 0;
}


int obtainGroundTarr( gsl_complex **eigvec, int degeneracy, gsl_complex **groundarr, gsl_complex **groundTarr , int *gndarr, int *gndblockarr ){	// You must hold the block structure in 'eigvec'.
	for(int dd=0; dd<degeneracy; dd++ ) {
		int gnd = gndarr[dd];
		int gndblock = gndblockarr[dd];
		for(int i=0; i<gndblock; i++) groundTarr[dd][i] = zero ;
		for(int degen=0; degen<degeneracy; degen++ ) {
			if( abs(gnd - gndarr[degen])>0 ){}
			else { 
#pragma omp parallel for default(shared) 
				for(int i=0; i<gndblock; i++){
					groundTarr[dd][i] = gsl_complex_add(
							groundTarr[dd][i] ,
							gsl_complex_mul( eigvec[dd][degen] , groundarr[degen][i] )
							);
				}
			}
		}
	}
	return 0;
}

int obtainGroundTarr( int degeneracy, gsl_complex **groundarr, gsl_complex **groundTarr , int *gndarr, int *gndblockarr, typebasis **basis, int *table , char opname[] ){
	Op op1(opname) ;  //Trev
	//op1.mkOp(6,6);	// !! Complex conjugate the coefficients AND the followings (its order should be kept. It is very important.)
	//if( fermionsgn ) {
	//	op1.cDat[0].setDat( 3, 0,  1 ); // 1 
	//	op1.cDat[1].setDat( 2, 1, -1 ); // 1 
	//	op1.cDat[2].setDat( 1, 2,  1 ); // 1 
	//	op1.cDat[3].setDat( 0, 3, -1 ); // 1 
	//	op1.cDat[4].setDat( 5, 4,  1 ); // 1 
	//	op1.cDat[5].setDat( 4, 5, -1 ); // 1 
	//}
	//else { 
	//	op1.cDat[0].setDat( 3, 0,  1 ); // 1 
	//	op1.cDat[1].setDat( 2, 1,  1 ); //-1 
	//	op1.cDat[2].setDat( 1, 2,  1 ); // 1 
	//	op1.cDat[3].setDat( 0, 3,  1 ); //-1 
	//	op1.cDat[4].setDat( 5, 4,  1 ); // 1 
	//	op1.cDat[5].setDat( 4, 5,  1 ); //-1 
	//}

	for(int dd=0; dd<degeneracy; dd++ ) {
		int gndblock	= gndblockarr[dd];
		for(int i=0; i<gndblock; i++) groundTarr[dd][i] = zero ;
	}

	for(int dd=0; dd<degeneracy; dd++ ) {
		int gnd		= gndarr[dd];
		int gndblock	= gndblockarr[dd];
		double Peven=0. ;
		double Podd =0. ;
		double PevenTrev=0. ;
		double PoddTrev =0. ;
		for(int i=gnd; i<gnd+gndblock; i++){
			int ig = i-gnd ; 

			typebasis bin[2] ;
			bin[0]	= basis[i][0];
			bin[1]	= basis[i][1];

			typebasis Trevbin[2] ;
			Trevbin[0]	= basis[i][0];
			Trevbin[1]	= basis[i][1];
			for( int mu=0 ; mu<tNC; mu++ ) 
				if( One(Trevbin[mu%2], mu/2) ) Trevbin[mu%2] = Turnoff( Trevbin[mu%2], mu/2 );
			//printf("coeff,bin,Trevbin : %12d %19.16f+i%19.16f ",i, GSL_REAL(groundarr[dd][ig]), GSL_IMAG(groundarr[dd][ig]) );
			//for( int mu=0 ; mu<tNC; mu++ ){ printf("%d", One( bin[mu%2], mu/2 ) );} for( int mu=0 ; mu<tNC; mu++ ){ printf("%d", One( Trevbin[mu%2], mu/2 ) );}
			int phaseAnn = 1;
			int phaseCre = 1;
			double op1valprod = 1.;
			int ncorr=0;
			for( int bb=0 ; bb<op1.ndat ; bb++ ) {
				int mu1 = op1.cDat[bb].i ;
				int nu1 = op1.cDat[bb].j ;
				if( One( basis[i][nu1%2], nu1/2) ) {
					op1valprod *= op1.cDat[bb].dat ;
					phaseAnn *= permu(bin, nu1/2, nu1%2); bin[nu1%2] = Turnoff(bin[nu1%2], nu1/2 );
					if( Zero(Trevbin[mu1%2], mu1/2) ) {
						phaseCre *= permu(Trevbin, mu1/2, mu1%2); Trevbin[mu1%2] = Turnon(Trevbin[mu1%2], mu1/2 );
					}
					else if( One(Trevbin[mu1%2], mu1/2) ) {
						printf("ERROR :: In 'obtainGroundTrevarr()', the transform does not have the one-to-one correspondence.\n");
						//char aa[1024];
						//sprintf(aa,"%d-[%d,%d](%d):", i,mu1,nu1, omp_get_thread_num() );
						//printf("%s bO : %d\tbO : %d\n", aa, One(bin[mu1%2], mu1/2), One(bin[nu1%2], nu1/2) );
						//printf("%s TO : %d\tTO : %d\n", aa, One(Trevbin[mu1%2], mu1/2), One(Trevbin[nu1%2], nu1/2) );
						//bin[nu1%2] = Turnon(bin[nu1%2], nu1/2 );
						//printf("%s bO': %d\tbO': %d\n", aa, One(bin[mu1%2], mu1/2), One(bin[nu1%2], nu1/2) );
						//printf("%s TO': %d\tTO': %d\n", aa, One(Trevbin[mu1%2], mu1/2), One(Trevbin[nu1%2], nu1/2) );
						exit(1);
					}
					else {
						printf("ERROR :: In 'obtainGroundTrevarr()', the transform is doing nothing (something wrong).\n");
						//char aa[1024];
						//sprintf(aa,"%d-[%d,%d]:\n", i,mu1,nu1);
						//printf("%s bO : %d\tbO : %d\n", aa, One(bin[mu1%2], mu1/2), One(bin[nu1%2], nu1/2) );
						//printf("%s TO : %d\tTO : %d\n", aa, One(Trevbin[mu1%2], mu1/2), One(Trevbin[nu1%2], nu1/2) );
						exit(1);
					}
					ncorr++;
				}
				//printf("After munu=%d,%d:\n", mu1,nu1 );
				//printf("bbin : "); for( int mu=0 ; mu<tNC; mu++ ) printf("%d\t", One( Trevbin[mu%2], mu/2 ) );
				//printf("Trev : "); for( int mu=0 ; mu<tNC; mu++ ) printf("%d\t", One( bin[mu%2], mu/2 ) ); printf("\n");
			}
			if( phaseAnn < 0  ) {
				printf("ERROR :: In 'obtainGroundTrevarr()', the order of the annihilation is incorrect in transforming.\n");
				exit(1);
			}
			//phaseCre *= op1valprod ;
			double fact = op1valprod*phaseCre ;
			int index = (Trevbin[0]<<Ns) + Trevbin[1];
			groundTarr[dd][table[index]-gnd] = gsl_complex_rect( GSL_REAL(groundarr[dd][ig])*fact, -GSL_IMAG(groundarr[dd][ig])*fact );
			//printf("| %2d*%4.1f*%19.16f+i%19.16f  ",phaseCre,op1valprod,GSL_REAL(groundTarr[dd][table[index]-gnd]), GSL_IMAG(groundTarr[dd][table[index]-gnd]) );
			//for( int mu=0 ; mu<tNC; mu++ ){ printf("%d", One( bin[mu%2], mu/2 ) );} for( int mu=0 ; mu<tNC; mu++ ){ printf("%d", One( Trevbin[mu%2], mu/2 ) );}
			//printf("\n");
			//if ( ig > 1000 ) exit(1);
			double prob = gsl_complex_abs2( groundarr[dd][ig] ),
			       probT= gsl_complex_abs2( groundTarr[dd][table[index]-gnd] );
			if( ncorr%2 )	{ Podd += prob; PoddTrev += probT; }
			else		{ Peven+= prob; PevenTrev+= probT; }
		}
		printf("P(    phi[%d]){even,odd} : %19.16f %19.16f\n", dd, Peven,	  Podd     );
		printf("P(Trevphi[%d]){even,odd} : %19.16f %19.16f\n", dd, PevenTrev, PoddTrev );
	}
	op1.freeOp() ;
	return 0;
}

int obtainGroundTarr( int degeneracy, gsl_complex **groundarr, gsl_complex **groundTarr , int *gndarr, int *gndblockarr, typebasis **basis, int *table , OpZ op1 ){
	for(int dd=0; dd<degeneracy; dd++ ) {
		int gndblock	= gndblockarr[dd];
		for(int i=0; i<gndblock; i++) groundTarr[dd][i] = zero ;
	}

	for(int dd=0; dd<degeneracy; dd++ ) {
		int gnd		= gndarr[dd];
		int gndblock	= gndblockarr[dd];
		double Peven=0. ;
		double Podd =0. ;
		double PevenTrev=0. ;
		double PoddTrev =0. ;
		for(int i=gnd; i<gnd+gndblock; i++){
			int ig = i-gnd ; 

			typebasis bin[2] ;
			bin[0]	= basis[i][0];
			bin[1]	= basis[i][1];

			typebasis Trevbin[2] ;
			Trevbin[0]	= basis[i][0];
			Trevbin[1]	= basis[i][1];
			for( int mu=0 ; mu<tNC; mu++ ) 
				if( One(Trevbin[mu%2], mu/2) ) Trevbin[mu%2] = Turnoff( Trevbin[mu%2], mu/2 );
			//printf("coeff,bin,Trevbin : %12d %19.16f+i%19.16f ",i, GSL_REAL(groundarr[dd][ig]), GSL_IMAG(groundarr[dd][ig]) );
			//for( int mu=0 ; mu<tNC; mu++ ){ printf("%d", One( bin[mu%2], mu/2 ) );} for( int mu=0 ; mu<tNC; mu++ ){ printf("%d", One( Trevbin[mu%2], mu/2 ) );}
			int phaseAnn = 1;
			int phaseCre = 1;
			gsl_complex op1valprod = gsl_complex_rect(1.,0);
			int ncorr=0;
			for( int bb=0 ; bb<op1.ndat ; bb++ ) {
				int mu1 = op1.cDat[bb].i ;
				int nu1 = op1.cDat[bb].j ;
				if( One( basis[i][nu1%2], nu1/2) ) {
					op1valprod = gsl_complex_mul( op1valprod, op1.cDat[bb].dat );
					phaseAnn *= permu(bin, nu1/2, nu1%2); bin[nu1%2] = Turnoff(bin[nu1%2], nu1/2 );
					if( Zero(Trevbin[mu1%2], mu1/2) ) {
						phaseCre *= permu(Trevbin, mu1/2, mu1%2); Trevbin[mu1%2] = Turnon(Trevbin[mu1%2], mu1/2 );
					}
					else if( One(Trevbin[mu1%2], mu1/2) ) {
						printf("ERROR :: In 'obtainGroundTrevarr()', the transform does not have the one-to-one correspondence.\n");
						//char aa[1024];
						//sprintf(aa,"%d-[%d,%d](%d):", i,mu1,nu1, omp_get_thread_num() );
						//printf("%s bO : %d\tbO : %d\n", aa, One(bin[mu1%2], mu1/2), One(bin[nu1%2], nu1/2) );
						//printf("%s TO : %d\tTO : %d\n", aa, One(Trevbin[mu1%2], mu1/2), One(Trevbin[nu1%2], nu1/2) );
						//bin[nu1%2] = Turnon(bin[nu1%2], nu1/2 );
						//printf("%s bO': %d\tbO': %d\n", aa, One(bin[mu1%2], mu1/2), One(bin[nu1%2], nu1/2) );
						//printf("%s TO': %d\tTO': %d\n", aa, One(Trevbin[mu1%2], mu1/2), One(Trevbin[nu1%2], nu1/2) );
						exit(1);
					}
					else {
						printf("ERROR :: In 'obtainGroundTrevarr()', the transform is doing nothing (something wrong).\n");
						//char aa[1024];
						//sprintf(aa,"%d-[%d,%d]:\n", i,mu1,nu1);
						//printf("%s bO : %d\tbO : %d\n", aa, One(bin[mu1%2], mu1/2), One(bin[nu1%2], nu1/2) );
						//printf("%s TO : %d\tTO : %d\n", aa, One(Trevbin[mu1%2], mu1/2), One(Trevbin[nu1%2], nu1/2) );
						exit(1);
					}
					ncorr++;
				}
				//printf("After munu=%d,%d:\n", mu1,nu1 );
				//printf("bbin : "); for( int mu=0 ; mu<tNC; mu++ ) printf("%d\t", One( Trevbin[mu%2], mu/2 ) );
				//printf("Trev : "); for( int mu=0 ; mu<tNC; mu++ ) printf("%d\t", One( bin[mu%2], mu/2 ) ); printf("\n");
			}
			if( phaseAnn < 0  ) {
				printf("ERROR :: In 'obtainGroundTrevarr()', the order of the annihilation is incorrect in transforming.\n");
				exit(1);
			}
			//phaseCre *= op1valprod ;
			//double fact = op1valprod*phaseCre ;
			int index = (Trevbin[0]<<Ns) + Trevbin[1];
			groundTarr[dd][table[index]-gnd] = gsl_complex_mul( gsl_complex_conjugate(groundarr[dd][ig]) , gsl_complex_mul_real(op1valprod,phaseCre) );
			//printf("| %2d*%4.1f*%19.16f+i%19.16f  ",phaseCre,op1valprod,GSL_REAL(groundTarr[dd][table[index]-gnd]), GSL_IMAG(groundTarr[dd][table[index]-gnd]) );
			//for( int mu=0 ; mu<tNC; mu++ ){ printf("%d", One( bin[mu%2], mu/2 ) );} for( int mu=0 ; mu<tNC; mu++ ){ printf("%d", One( Trevbin[mu%2], mu/2 ) );}
			//printf("\n");
			//if ( ig > 1000 ) exit(1);
			double prob = gsl_complex_abs2( groundarr[dd][ig] ),
			       probT= gsl_complex_abs2( groundTarr[dd][table[index]-gnd] );
			if( ncorr%2 )	{ Podd += prob; PoddTrev += probT; }
			else		{ Peven+= prob; PevenTrev+= probT; }
		}
		printf("P(    phi[%d]){even,odd} : %19.16f %19.16f\n", dd, Peven,	  Podd     );
		printf("P(Trevphi[%d]){even,odd} : %19.16f %19.16f\n", dd, PevenTrev, PoddTrev );
	}
	op1.freeOp() ;
	return 0;
}

int obtainGroundT( gsl_complex *ground, gsl_complex *groundT , int gnd, int gndblock, typebasis **basis, int *table ) {
#pragma omp parallel for default(shared)
	for(int i=0; i<gndblock; i++) groundT[i] = zero ;
	gsl_complex **T = mkgscmatrixd(NU, NU);
	initgscmatrixd(T, NU, NU);
	init_transform(T);

	for(int i=gnd; i<gnd+gndblock; i++){
		int ig = i-gnd ; 

		int nuon[tNC], nlatt=0, lattiter=1;
		for( int nu1=0; nu1<tNC ; nu1++ ) {
			if( One( basis[i][nu1%2], nu1/2) )  {
				nuon[nu1] = 1;
				nlatt	 += 1;
				lattiter *= tNC;
			}
			else	nuon[nu1] = 0;
		}
		gsl_complex coeffprod = gsl_complex_rect(1.,0.);
		int nuarr[nlatt], ii=0;
		for( int nu1=0; nu1<tNC ; nu1++ ) 
			if( nuon[nu1]>0 ) {
				nuarr[ii]=nu1;
				ii++;
			}

		double coeffground = gsl_complex_abs2( ground[ig] );
		if( coeffground>0 ){
			for( int imu=0; imu<lattiter ; imu++ ) {
				int muarr[nlatt];
				for( int nuind=0; nuind<nlatt ; nuind++ ) {
					muarr[nuind] = imu;
					for( int ndeno=0; ndeno<nuind ; ndeno++ )
						muarr[nuind] /=tNC;
					muarr[nuind] %=tNC;
				}
				typebasis Tbin[2] ;
				Tbin[0]	= basis[i][0];
				Tbin[1]	= basis[i][1];
				GSL_REAL(coeffprod) = GSL_REAL( ground[ig] );
				GSL_IMAG(coeffprod) = GSL_IMAG( ground[ig] );
				int phaseCre = 1;
				for( int nuind=0; nuind<nlatt ; nuind++ ) {
					int mu=nuarr[nuind];
					Tbin[mu%2] = Turnoff( Tbin[mu%2], mu/2 );
				}

				for( int nuind=0 ; nuind<nlatt; nuind++ ) {
					int mu = muarr[nuind],
					    nu = nuarr[nuind];
					coeffprod = gsl_complex_mul( coeffprod, T[mu][nu] ); 
				}
				double coeffreabs = fabs(GSL_REAL(coeffprod)),
				       coeffimabs = fabs(GSL_IMAG(coeffprod));
				bool diffmu = 1 ;
				for( int nuind=0 ; nuind<nlatt; nuind++ ) for( int nuind2=nuind+1 ; nuind2<nlatt; nuind2++ ) {
					int dmuabs = abs(muarr[nuind]-muarr[nuind2]);
					if( dmuabs!=0 ) {}
					else		diffmu=0;
				}
				if( diffmu!=0 && coeffreabs>0 && coeffimabs>0 ) {
					for( int nuind=0 ; nuind<nlatt; nuind++ ) {
						int mu = muarr[nuind];
						if( Zero( Tbin[mu%2], mu/2) ) {
							phaseCre *= permu(Tbin, mu/2, mu%2);
							Tbin[mu%2] = Turnon( Tbin[mu%2], mu/2 );
						}
						else  {
							printf("ERROR :: In 'obtainGroundT()', something wrong with acting creation operator.\n");
							exit(1);
						}
					}
					int index = (Tbin[0]<<Ns) + Tbin[1];
					//int indexorig = (basis[i][0]<<Ns) + basis[i][1];
					//printf("m");
					//for( int nuind=0; nuind<nlatt ; nuind++ ) printf("%d", muarr[nuind] );
					//printf("n");
					//for( int nuind=0; nuind<nlatt ; nuind++ ) printf("%d", nuarr[nuind] );
					//printf("mn");
					//for( int nuind=0; nuind<nlatt ; nuind++ ) printf("%d%d", muarr[nuind], nuarr[nuind] );
					//printf(" \n");
					//for( int nuind=0 ; nuind<nlatt; nuind++ ) {
					//	int mu = muarr[nuind],
					//	    nu = nuarr[nuind];
					//	printf("T[%d][%d]=%10.6f +i %10.6f\n", mu,nu, GSL_REAL(T[mu][nu]), GSL_IMAG(T[mu][nu]) ); 
					//}
					//for( int mu=0 ; mu<tNC; mu++ )
					//	printf("%d",One( basis[i][mu%2], mu/2 ) );
					//printf("  ");
					//for( int mu=0 ; mu<tNC; mu++ )
					//	printf("%d",One( Tbin[mu%2], mu/2 ) );
					//printf("  ");
					//printf( "nlatt=%d i=%d ind=%d\tgnd=%d\tindorig=%d\tgndblock=%d\tindg=%d\tindoirgg=%d coeff=%19.16f\t+i%19.16f\tcoeff=%19.16f\t+i%19.16f\n", nlatt, i, table[index], gnd, table[indexorig] , gndblock , table[index]-gnd , table[indexorig]-gnd ,
					//		GSL_REAL(coeffprod), GSL_IMAG(coeffprod), coeffreabs, coeffimabs 
					//      );
					groundT[table[index]-gnd] = gsl_complex_add( groundT[table[index]-gnd], gsl_complex_mul_real(coeffprod,phaseCre) );
				}
				//printf("i=%10d imu=%10d lattiter=%10dnlatt=%10d\r",i,imu,lattiter,nlatt);fflush(stdout);
			}
			//printf("\rig=%10d nlatt=%10d",ig,nlatt);fflush(stdout);
		}
	}
	//printf("\n");
	double prob = 0.;
#pragma omp parallel for default(shared) reduction(+:prob)
	for(int i=0; i<gndblock; i++) prob += gsl_complex_abs2( groundT[i] );
	printf("Transformed :: Prob = %f\n", prob );
	fflush(stdout);
	
	freegscmatrixd(T, NU);
	return 0;
}
int obtainGroundT( gsl_complex *ground, gsl_complex *groundT , int gnd, int gndblock, typebasis **basis, int *table , OpZ op1 ){
    return obtainGroundT( ground, groundT , gnd, gndblock, basis, table ) ;
}

int write_entropy( char pararesult[] , int ic , double *weightarr , char writemode[] , int ndeg, char prefix[] ){
    FILE *fw ; 
    char fname[1025] ;
    sprintf( fname, "%s/%sentropy.dat", pararesult, prefix ) ;
    fw = fopen( fname , writemode );
    fprintf( fw , "%d\t", ic ) ;
    double entropy = 0 ; 
    for(int mu=0; mu<ndeg; mu++)  {
        double evalmu = weightarr[mu];
        if( fabs(evalmu)<1e-15 ) evalmu += 1e-15;
        entropy -= log( evalmu ) * evalmu ;
    }
    fprintf( fw , "%19.16lf\n", entropy ) ;
    printf("entropy : %19.16lf\n", entropy  ) ;

    fclose(fw) ;
    return 0 ;
}
int write_entropy( char pararesult[] , int ic , double *weightarr , char writemode[], int ndeg ){
    return write_entropy( pararesult, ic, weightarr, writemode, ndeg, "" );
}

int write_EEall( char pararesult[] , int ic , double *weightarr , double *EEarr, char writemode[] , int ndeg, char prefix[] ){
    FILE *fw ; 
    char fname[1025] ;
    sprintf( fname, "%s/%sEEall.dat", pararesult, prefix ) ;
    fw = fopen( fname , writemode );
    fprintf( fw , "%d\t", ic ) ;
    double EEall = 0 ; 
    for(int mu=0; mu<ndeg; mu++)  {
        double evalmu = weightarr[mu] * EEarr[mu];
        EEall += evalmu ;
    }
    fprintf( fw , "%19.16lf\n", EEall ) ;
    printf("EEall : %19.16lf\n", EEall  ) ;

    fclose(fw) ;
    return 0 ;
}
int write_EEall( char pararesult[] , int ic , double *weightarr , double *EEarr, char writemode[], int ndeg ){
    return write_EEall( pararesult, ic, weightarr, EEarr, writemode, ndeg, "" );
}
