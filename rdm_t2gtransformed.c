#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
#include "matrix.h"
#include "basic.h"
#include "lattice.h"
#include "ls_basis.c"
#include "ls.c"
#include "opclass.cpp"
#include "rdmclass.cpp"
#include "exciton_static.cpp"
#include "mps.cpp"
#include "degop.cpp"
#include "convert_ground.cpp"

#define DIMd	64						//for rdm_static.c
#define DIMd4	15						//for rdm_static.c
#define gcr(zre, zim)	gsl_complex_rect( zre, zim )		//for rdm_static.c
void compute_static(typebasis **basis, int count, char *save_directory);

//long seed = -134551;
double J, MU, DELTA, Uinter, JTD, SOC, Hz;
int beta, Nmax, nctot, nb, tnb, Spinmax, Blocks,Powns, Powns2, myrank=0, beta_real, Conversion_required, Nonzero_SOC ;
gsl_complex zero, ***Greennew, **Ecluster;
Hamiltonianclass Hc ;

int main(int argc, char *argv[]){
	if(argc != 11){
		printf("Usage :: %s <input.mak> <Hz> <JTD> <SOC> <J> <U_init> <U_final> <DELTA> <save_directory> <beta_real>\n", argv[0]);
		exit(1);
	}
	sectorinfo **ginfo;
	int count = 0;
	double U_init, U_final;

	FILE *fd; 
	char pathinform[1024], save_directory[1024];

	nctot = NC;
	nb = NB;
	tnb = 2*NB;
	Powns	= (int)pow(2,Ns);
	Powns2	= (int)pow(2,2*Ns);
	Blocks	= 2*Ns+1;
	Spinmax	= 2*Ns+1;
	zero = gsl_complex_rect(0,0);
	Hz	= atof( argv[2] );
	JTD	= atof( argv[3] );
	SOC	= atof( argv[4] );
	J	= atof( argv[5] );
	U_init	= atof( argv[6] );
	U_final	= atof( argv[7] );
	DELTA	= atof( argv[8] );
	beta_real	= atoi( argv[10] );

	Nonzero_SOC = 1;

	if( fabs(SOC) < 1e-6 )	Conversion_required = 1;
	else			Conversion_required = 0;

	sprintf(save_directory, "%s", argv[9]);

	fd = fopen(argv[1], "r");
	nofile(fd, argv[1]);
	fscanf(fd, "%d", &count);
	fscanf(fd, "%lf", &Uinter);
	fclose(fd);

	typebasis **basis;
	basis = mkmatrixb(Powns2, 2);
	ginfo = (sectorinfo **) malloc( Ni * sizeof(sectorinfo*));
	for(int iimp=0; iimp<Ni; iimp++){
		ginfo[iimp] = (sectorinfo *) malloc( Blocks*sizeof(sectorinfo) );
		build_bitbasis(Ns, basis, ginfo[iimp]);
		print_ginfo(Ns, ginfo[iimp]);
	}
	printf("U = %.2lf, MU = %.2lf (of UF%.2lf UT%.2lf) from '%s'\n", Uinter, MU, U_init, U_final, argv[1]);

	Ecluster = mkgscmatrixd( NU, NU);
	char paraEcluster[1024];
	sprintf(paraEcluster, "%s/Ecluster.dat", save_directory);
	load_Ecluster(paraEcluster, Ecluster);
	update_chem(Ecluster);

	Hc.allocgsl() ;
	Hc.eigsystem( Ecluster ) ;
	Hc.dotransform_Ecluster( Hc.evecEcluster, Ecluster ) ;
	freegscmatrixd( Ecluster , NU ) ;

	sprintf(pathinform, "%s/result/%c%.3lf", save_directory, control_char, control_para);
	mkdirectory(pathinform);

	compute_static(basis, count, save_directory);

	for(int iimp=0; iimp<Powns2; iimp++)	free(basis[iimp]);		free(basis);
	int i;
	for(int iimp=0; iimp<Ni; iimp++){
		for(i=0; i<Blocks; i++)
			freegscmatrixd(ginfo[iimp][i].ground, Nvec);
		free( ginfo[iimp] );
	}
	free(ginfo);
	Hc.freegsl() ;

	return 0;
}

int read_mak(PNSpara *egv, int *degeneracy, double ***ginformation, char *save_directory, int count){
	int tol, i, j, k;
	char para[1024];
	double re, im;
	FILE *fd;
	sprintf(para, "%s/%c%.2lf_%dth.mak", save_directory, control_char, control_para, count);
	printf("Reading .mak : %s\n", para);

	fd = fopen(para, "r");
	if( fd == NULL ){
		printf("%s does not exist!\n", para);
		return 1;
	}
	fscanf(fd, "%d", &count);
	fscanf(fd, "%lf", &Uinter);
	fscanf(fd, "%lf", &MU);
	fscanf(fd, "%d", &tol);
	fscanf(fd, "%d", &beta);
	fscanf(fd, "%d", &Nmax);

	for(int iimp=0; iimp<Ni; iimp++){
		fscanf(fd, "%d", &degeneracy[iimp]);

		for(i=0; i<degeneracy[iimp]; i++) for(j=0; j<3; j++)	fscanf(fd, "%lf", &ginformation[iimp][i][j]);

		for(k=0; k<tnb; k++){
			fscanf(fd, "%lf", &re);
			egv[iimp].egbath[k] = re;
		}
		for(j=0; j<tNC; j++) for(k=0; k<tnb; k++){
			fscanf(fd, "%lf %lf", &re, &im);
			egv[iimp].hybrid[j][k] = gsl_complex_rect( re, im );
		}
	}
	fclose(fd);

	return 0;
}

void compute_static(typebasis **basis, int count, char *save_directory){
	PNSpara egv[Ni];
	typebasis bin[2];
	int degen, gndblock, ic, degeneracy[Ni], mu, nu, i, gnd, phase, *table, index, final_nonzero=0;
	//char paramag[1024], paradouble[1024], parafilling[1024] ;
	char paragnd[1024], paraext[1024], 
	     parardmeval[1024] , parasusCC[1024] , parancrdm[1024] , parancrdmevec[1024] , parafillingdeg[1024] ,
	     parardminfo[1024] , parardmdiag[1024] ,
	     pararesult[1024] ;
	//FILE *fmag, *ffil, *fdoc ;
	FILE *fground, *fext,
	     *frdmeval , *fsusCC, *fncrdm, *fncrdmevec , *ffilldeg , *frdminfo ,
	     *frdmdiag ; 
	gsl_complex *ground;
	gsl_complex ***excitonic  = mkgsctritensord(Ni, tNC, tNC);
	gsl_complex ***excitonicT = mkgsctritensord(Ni, tNC, tNC);
	gsl_complex ***excitonict2g = mkgsctritensord(Ni, tNC, tNC);
	double doccu[Ni][NC], filling[Ni][NC][2], mag[Ni][NC], sum[3], magfactor[tNC]={1.5,0.5,-0.5,-1.5,0.5,-0.5};
	double jsqfactor[tNC]={3.75,3.75,3.75,3.75,0.75,0.75,};
	gsl_complex *jsq = mkgscvectord( Ni ) ;
	double ***ginformation;
	ginformation = mktritensord(Ni, 50*SDmax, 3);
	int Powtnb = int( pow(2,2*NB) ) ; 
	int Powtnc = int( pow(2, tNC) ) ; 
	gsl_complex ***rdm = mkgsctritensord(Ni, Powtnc, Powtnc);
	gsl_complex ***mps = mkgsctritensord(Ni, Powtnb, Powtnc);
	gsl_complex ****mps3orb = mkgsctetratensord(4,4,4, Powtnb);
	gsl_complex **susCC= mkgscmatrixd(tNC, tNC);
	gsl_complex **transform = mkgscmatrixd(NU, NU); init_transform(transform) ;

	//sprintf(paramag, "%s/result/%c%.3lf/mag.dat", save_directory, control_char, control_para);
	//sprintf(paradouble, "%s/result/%c%.3lf/double.dat", save_directory, control_char, control_para);
	//sprintf(parafilling, "%s/result/%c%.3lf/filling.dat", save_directory, control_char, control_para);
	sprintf(paraext, "%s/result/%c%.3lf/excitonic.dat", save_directory, control_char, control_para);
	sprintf(parafillingdeg, "%s/result/%c%.3lf/fillingdeg.dat", save_directory, control_char, control_para);
	sprintf(parardminfo, "%s/result/%c%.3lf/rdminfo.dat", save_directory, control_char, control_para);
	sprintf(pararesult, "%s/result/%c%.3lf/", save_directory, control_char, control_para);
	char parardmjjJ[1024] ; sprintf(parardmjjJ, "%s/rdmjjJd4diag.dat",pararesult) ; 
	char parardmLSJ[1024] ; sprintf(parardmLSJ, "%s/rdmLSJd4diag.dat",pararesult) ; 
	char paramagdeg[1024] ; sprintf( paramagdeg, "%s/magdeg.dat", pararesult ) ;
	//fmag = fopen(paramag, "w");
	//fdoc = fopen(paradouble, "w");
	//ffil = fopen(parafilling, "w");
	fext = fopen(paraext, "w");
	ffilldeg = fopen(parafillingdeg, "w");
	frdminfo = fopen(parardminfo, "w");
	Vecop Jop("Jop") ; 
	FILE *frdmjjJdiag = fopen( parardmjjJ , "w" ) ;
	FILE *frdmLSJdiag = fopen( parardmLSJ , "w" ) ;
	FILE *fmagdeg     = fopen( paramagdeg , "w" ) ;
	fclose(frdmjjJdiag) ;
	fclose(frdmLSJdiag) ;
	fclose(fmagdeg    ) ;


	table = mkvectori( Powns2 );
	for(i=0; i<Powns2; i++){
		index = (basis[i][0]<<Ns) + basis[i][1];
		table[index] = i;
	}

	for(int iimp=0; iimp<Ni; iimp++){
		egv[iimp].egbath = mkvectord(tnb);	
		egv[iimp].hybrid = mkgscmatrixd(tNC, tnb);
	}

	for(mu=0; mu<Powtnc; mu++) {
		unsigned int  cup = (mu>>NC) ,
			      cdn = mu - (cup<<NC) ;
		int ncdum = 0 ;
		for(nu=0; nu<NC; nu++){
			int ncupmu = (((cup)>>nu)&0x01) ; 
			int ncdnmu = (((cdn)>>nu)&0x01) ;
			fprintf(frdminfo, "%d%d ", ncupmu, ncdnmu ) ;
			ncdum += ncupmu + ncdnmu ;
		}
		fprintf(frdminfo, "%d\n", ncdum ) ;
	}

    typebasis **basis_t2g ; basis_t2g       = mkmatrixb(Powns2, 2);
	sectorinfo **ginfo_t2g ;
	int Blocks_t2g = (Ns+1)*(Ns+1);

	ginfo_t2g = (sectorinfo **) malloc( Ni * sizeof(sectorinfo*));
	for(int iimp=0; iimp<Ni; iimp++){
		ginfo_t2g[iimp] = (sectorinfo *) malloc( Blocks_t2g*sizeof(sectorinfo) );
		build_bitbasis(Ns, basis_t2g, ginfo_t2g[iimp], 0 );
		print_ginfo(Ns, ginfo_t2g[iimp]);
	}
	for(int iimp=0; iimp<Ni; iimp++){
		for(int iblock=0; iblock<Blocks_t2g; iblock++)
			freegscmatrixd(ginfo_t2g[iimp][iblock].ground, Nvec);
		free( ginfo_t2g[iimp] );
	}
	free(ginfo_t2g);



	char writemode[1024]		; sprintf(writemode,"w") ;
	char writemodedegen[1024]	; sprintf(writemodedegen,"w") ;
	for(ic=count; ic<count+1; ic++){
		if( read_mak(egv, degeneracy, ginformation, save_directory, ic) ){
			if( ic == count ){
				printf("no mak file exists!\n");
				exit(1);
			}
			continue;
		}

		//fprintf(fmag, "%d\t", ic);
		//fprintf(fdoc, "%d\t", ic);
		//fprintf(ffil, "%d\t", ic);
		fprintf(fext, "%d\t", ic);
		//fprintf(fsusCC, "%d\t", ic);
		sum[0] = 0;
		sum[1] = 0;
		sum[2] = 0;
		for(int iimp=0; iimp<Ni; iimp++){
			char prefixiimp[1024]; sprintf( prefixiimp, "imp%d_", iimp );
            sprintf(parardmeval, "%s/result/%c%.3lf/%srdmeval.dat", save_directory, control_char, control_para, prefixiimp);
            sprintf(parardmdiag, "%s/result/%c%.3lf/%srdmdiag.dat", save_directory, control_char, control_para, prefixiimp);
            sprintf(parasusCC, "%s/result/%c%.3lf/%ssusCC.dat", save_directory, control_char, control_para, prefixiimp);
            sprintf(parancrdm, "%s/result/%c%.3lf/%sncrdm.dat", save_directory, control_char, control_para, prefixiimp);
            sprintf(parancrdmevec, "%s/result/%c%.3lf/%sncrdmevec.dat", save_directory, control_char, control_para, prefixiimp);
            char parardmt2gdiag[1024]; sprintf(parardmt2gdiag, "%s/result/%c%.3lf/%srdmt2gdiag.dat", save_directory, control_char, control_para, prefixiimp);
            frdmeval = fopen(parardmeval, writemode);
            frdmdiag = fopen(parardmdiag, writemode);
            //fsusCC = fopen(parasusCC, writemode);
            fncrdm = fopen(parancrdm, writemode);
            fncrdmevec = fopen(parancrdmevec, writemode);
            FILE *frdmt2gdiag = fopen(parardmt2gdiag, writemode);

			sprintf(paragnd, "%s/u%.2lf_%dth_n%d.gnd", save_directory, Uinter, ic, iimp);
			fground = fopen(paragnd, "rb");
			nofile(fground, paragnd);

			for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	{
				excitonic[iimp][mu][nu] = zero;
				excitonicT[iimp][mu][nu] = zero;
				excitonict2g[iimp][mu][nu] = zero;
			}
			for(mu=0; mu<nctot; mu++){
				doccu[iimp][mu] = 0;
				filling[iimp][mu][0] = 0;
				filling[iimp][mu][1] = 0;
				mag[iimp][mu] = 0;
			}

            double *myenergyarr     = mkvectord( degeneracy[iimp] ) ;
            double *myweightarr     = mkvectord( degeneracy[iimp] ) ;
            double *myEEarr         = mkvectord( degeneracy[iimp] ) ;
            printf("Ground state(s) initilized :\n") ;
            for(degen=0; degen<degeneracy[iimp]; degen++){
                myenergyarr[degen]  =      ginformation[iimp][degen][0];
                printf("ic=%d, n=%d, degen=%3d/%3d\t, energy=%19.16f\t;\n",
                        ic, iimp, degen, degeneracy[iimp], myenergyarr[degen]
                      );
            }
            double gm = myenergyarr[0],
                   partition = 0.,
                   weightsum = 0.;
            printf("parition weight :\n");
            for(degen=0; degen<degeneracy[iimp]; degen++){
                if( beta_real ) myweightarr[degen] = exp( - ( myenergyarr[degen] - gm )*beta_real );
                else        myweightarr[degen] = 1.;
                partition += myweightarr[degen];
            }
            for(degen=0; degen<degeneracy[iimp]; degen++){
                printf("degen%3d : %19.16f\t%19.16f\n", degen, myweightarr[degen], myweightarr[degen]/partition );
                myweightarr[degen]  /= partition;
                weightsum += myweightarr[degen];
            }
            printf("sum : %19.16f\t%19.16f\n", partition, weightsum );
            write_doublearr( pararesult, ic, "boltzweight", writemode, myweightarr, degeneracy[iimp] );
            write_entropy( pararesult, ic, myweightarr, writemode, degeneracy[iimp], prefixiimp );


			for(degen=0; degen<degeneracy[iimp]; degen++){
				fprintf(frdmeval, "%d\t", ic);
				fprintf(frdmdiag, "%d\t", ic);
				fprintf(fncrdm, "%d\t", ic);
				fprintf(fncrdmevec, "%d\t", ic);
				fprintf(ffilldeg, "%d\t", ic);
				fprintf(frdmt2gdiag, "%d\t", ic);
				jsq[iimp] = zero ;
				for(mu=0; mu<nctot; mu++){
					for(nu=0; nu<nctot; nu++){
						susCC[mu][nu] = zero ; 
					}
				}
				for(mu=0; mu<Powtnc; mu++) {
					for(nu=0; nu<Powtnb; nu++)  mps[iimp][nu][mu] = zero;
					for(nu=0; nu<Powtnc; nu++)  rdm[iimp][mu][nu] = zero;
				}
				for(int mu0=0; mu0<4; mu0++) for(int mu1=0; mu1<4; mu1++) for(int mu2=0; mu2<4; mu2++) for(nu=0; nu<Powtnb; nu++) 
					mps3orb[mu0][mu1][mu2][nu] = zero;
				double fillingdeg[tNC] ;
				for(mu=0; mu<tNC; mu++) {
					//printf("%d/%d\n", mu,tNC) ; fflush(stdout) ;
					fillingdeg[mu]=0 ;
				}

				gnd = (int)ginformation[iimp][degen][1];
				gndblock = (int)ginformation[iimp][degen][2];

				ground = mkgscvectord( gndblock );
				fread(ground, sizeof(gsl_complex), gndblock, fground);
				printf("ic=%d, n=%d, gnd=%d, gndblock=%d, degen=%d\t;\n", ic, iimp, gnd, gndblock, degen);

                double prob = 0.;
#pragma omp parallel for default(shared) reduction(+:prob)
                for(int i=0; i<gndblock; i++) prob += gsl_complex_abs2( ground[i] );
                printf("Original :: Prob = %f\n", prob );

				double mpssum_t2g=0., rdmsum_t2g=0.;
				printf( "startMPS:" ) ;
				mpssum_t2g = obtainMPS(	 mps[iimp], ground, basis_t2g, gnd, gndblock )	; printf( "Done.\n" ) ; printf( "mpssum_t2g=%19.16f\n", mpssum_t2g ) ;
				obtainRDM( rdm[iimp], mps[iimp], Powtnb, Powtnc ) ; 
				for(mu=0; mu<Powtnc; mu++)  rdmsum_t2g += GSL_REAL(rdm[iimp][mu][mu]) ; 
				for(mu=0; mu<Powtnc; mu++)  { fprintf( frdmt2gdiag , "%19.16lf\t", GSL_REAL(rdm[iimp][mu][mu]) ) ; } fprintf( frdmt2gdiag , "\n") ;
				printf("rdmdiagsum_t2g   = %19.16lf\n", rdmsum_t2g ) ;
				for(mu=0; mu<Powtnc; mu++) {
					for(nu=0; nu<Powtnb; nu++)  mps[iimp][nu][mu] = zero;
					for(nu=0; nu<Powtnc; nu++)  rdm[iimp][mu][nu] = zero;
				}


                int gindex ;
                if( Conversion_required )
                    ground = convert_ground(ginformation[iimp][degen], ground, &gindex, &gnd, &gndblock);

				for(mu=0; mu<nctot; mu++){
					for(i = gnd; i<gnd+gndblock; i++){	
						doccu[iimp][mu] += gsl_complex_abs2(ground[i-gnd])* One(basis[i][0], mu) * One(basis[i][1], mu);
						mag[iimp][mu] += gsl_complex_abs2(ground[i-gnd])* ( (magfactor[2*mu]* (int)One(basis[i][0], mu)) + (magfactor[2*mu+1]* (    int)One(basis[i][1], mu)) );
						filling[iimp][mu][0] += gsl_complex_abs2(ground[i-gnd])* ((int)One(basis[i][0], mu) );
						filling[iimp][mu][1] += gsl_complex_abs2(ground[i-gnd])* ((int)One(basis[i][1], mu) );
					}
				}
				for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
					for(i = gnd; i<gnd+gndblock; i++){	
						if( Zero( basis[i][mu%2], mu/2) && One( basis[i][nu%2], nu/2 ) ){
							bin[0] = basis[i][0];
							bin[1] = basis[i][1];
							phase  = permu( bin, nu/2, nu%2);	bin[nu%2] = Turnoff( bin[nu%2], nu/2 );
							phase *= permu( bin, mu/2, mu%2);	bin[mu%2] = Turnon( bin[mu%2], mu/2 );
							index = (bin[0]<<Ns) + bin[1];
							//excitonic[iimp][mu][nu]
							//	= gsl_complex_add(
							//		excitonic[iimp][mu][nu],
							//		gsl_complex_mul( gsl_complex_conjugate(ground[table[index]-gnd]), gsl_complex_mul_real( ground[i-gnd], phase ) )
							//	);
						}
					}
				}

				//excitonclass excc ;
				//excc.allocgsl() ;
				//excc.obtain_exc( gnd, gndblock, basis, table, ground )   ;
				//excc.write_excdeg( pararesult, ic , writemodedegen ) ;
				//excc.dotransforminv_excorder( transform ) ;
				//excc.write_excTdeg_datname(	pararesult, ic , "excitonict2gdeg.dat" , writemodedegen ) ;
				//excc.write_excTdegdiag_datname(	pararesult, ic , "fillingt2gdeg.dat" , writemodedegen ) ;
				//for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) excitonict2g[iimp][mu][nu] = gsl_complex_add( excitonict2g[iimp][mu][nu], gsl_complex_rect( GSL_REAL(excc.excorderT[mu][nu]), GSL_IMAG(excc.excorderT[mu][nu]) ) ) ;
				//excc.dotransform_excorder( Hc.evecEcluster ) ;
				//excc.write_excTdeg( pararesult, ic , writemodedegen ) ;
				//for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
				//	excitonicT[iimp][mu][nu] = gsl_complex_add( excitonicT[iimp][mu][nu], gsl_complex_rect( GSL_REAL(excc.excorderT[mu][nu]), GSL_IMAG(excc.excorderT[mu][nu]) ) ) ;
				//}
				//excc.freegsl() ;

				double magdeg = 0 ; 
				for(mu=0; mu<tNC; mu++) {
					for(i = gnd; i<gnd+gndblock; i++)	fillingdeg[mu] += gsl_complex_abs2(ground[i-gnd])* ((int)One(basis[i][mu%2], mu/2) );
					fprintf(ffilldeg, "%19.16lf\t", fillingdeg[mu] ) ;
				}
				fprintf(ffilldeg, "\n") ;

				for(mu=0; mu<NC; mu++) 
					for(i = gnd; i<gnd+gndblock; i++)	magdeg += gsl_complex_abs2(ground[i-gnd])* ( (magfactor[2*mu]* (int)One(basis[i][0], mu)) + (magfactor[2*mu+1]* (    int)One(basis[i][1], mu)) );
				fmagdeg = fopen( paramagdeg, writemodedegen ) ;
				fprintf( fmagdeg,  "%d\t%19.16lf\n", ic , magdeg ) ;
				fclose(fmagdeg ) ;


				printf( "db%d,dc%d:", Powtnb,Powtnc ) ;
				//printf( "startJsq:\n" ) ;
				//printf( "JJ JaJa : \n" ) ;
				//double jzjzval = calcDiagOpOp( Jop.zj, Jop.zj, ground, basis, gnd, gndblock ) ;
				//gsl_complex jzjzcalc	= calcOpOp( 		Jop.zj, Jop.zj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex jpjmcalc 	= calcOffdiagOpOp(	Jop.pj, Jop.mj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex jmjpcalc 	= calcOffdiagOpOp(	Jop.mj, Jop.pj, ground, basis, gnd, gndblock, table ) ;
				////printf( "jzjzval\t: %19.16f\n", jzjzval ) ;
				////printf( "jpjm\t: %19.16f +i %19.16f\n", GSL_REAL(jpjmcalc) , GSL_IMAG(jpjmcalc) ) ;
				////printf( "jmjp\t: %19.16f +i %19.16f\n", GSL_REAL(jmjpcalc) , GSL_IMAG(jmjpcalc) ) ;
				//printf( "jzjz\t: %19.16f +i %19.16f\n", GSL_REAL(jzjzcalc) , GSL_IMAG(jzjzcalc) ) ;
				//jsq[iimp] = gsl_complex_add(	jpjmcalc, jmjpcalc	) ;
				//jsq[iimp] = gsl_complex_mul_real(	jsq[iimp], 0.5		) ;		gsl_complex jxyjxycalc = gsl_complex_rect( GSL_REAL(jsq[iimp]), GSL_IMAG(jsq[iimp]) ) ;
				//jsq[iimp] = gsl_complex_add(	jsq[iimp], jzjzcalc	) ;
				//printf("jxyjxy\t: %19.16f +i %19.16f\n",	GSL_REAL(jxyjxycalc) ,	GSL_IMAG(jxyjxycalc) 	);
				//printf(    "jj\t: %19.16f +i %19.16f\n",	GSL_REAL(jsq[iimp]) , 	GSL_IMAG(jsq[iimp]) 	);

				Vecop Jaop("Jaop") ; 
				//gsl_complex jazjazcalc = calcOpOp( Jaop.zj, Jaop.zj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex japjamcalc = calcOpOp( Jaop.pj, Jaop.mj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex jamjapcalc = calcOpOp( Jaop.mj, Jaop.pj, ground, basis, gnd, gndblock, table ) ;
				////printf( "japjam\t: %19.16f +i %19.16f\n", GSL_REAL(japjamcalc) , GSL_IMAG(japjamcalc) ) ;
				////printf( "jamjap\t: %19.16f +i %19.16f\n", GSL_REAL(jamjapcalc) , GSL_IMAG(jamjapcalc) ) ;
				//printf( "jazjaz\t: %19.16f +i %19.16f\n", GSL_REAL(jazjazcalc) , GSL_IMAG(jazjazcalc) ) ;
				//gsl_complex jajacalc ;
				//jajacalc = gsl_complex_add( 		japjamcalc, jamjapcalc	) ;
				//jajacalc = gsl_complex_mul_real(	jajacalc, 0.5		) ;		gsl_complex jaxyjaxycalc = gsl_complex_rect( GSL_REAL(jajacalc), GSL_IMAG(jajacalc) ) ;
				//jajacalc = gsl_complex_add(		jajacalc, jazjazcalc	) ;
				//printf( "jaxyjaxy\t: %19.16f +i %19.16f\n", GSL_REAL(jaxyjaxycalc) , 	GSL_IMAG(jaxyjaxycalc)	) ;
				//printf(     "jaja\t: %19.16f +i %19.16f\n", GSL_REAL(jajacalc) , 	GSL_IMAG(jajacalc) 	) ;
				//write_ABC(	pararesult, ic ,   "JJdeg" , writemodedegen , jsq[iimp] ,   jzjzcalc ,   jxyjxycalc		) ;
				//write_ABC(	pararesult, ic , "JaJadeg" , writemodedegen , jajacalc , jazjazcalc , jaxyjaxycalc		) ;

				//printf( "J Ja : \n" ) ;
				//gsl_complex jzcalc  = calcOp(  Jop.zj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex jazcalc = calcOp( Jaop.zj, ground, basis, gnd, gndblock, table ) ;
				//printf(  "jzcalc\t: %19.16f +i %19.16f\n", GSL_REAL(jzcalc) ,  GSL_IMAG(jzcalc) ) ;
				//printf( "jazcalc\t: %19.16f +i %19.16f\n", GSL_REAL(jazcalc) , GSL_IMAG(jazcalc) ) ;
				//gsl_complex jpcalc  = calcOp(  Jop.pj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex jmcalc  = calcOp(  Jop.mj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex jxcalc  = gsl_complex_add( jpcalc, jmcalc	) ;	jxcalc = gsl_complex_mul_real( jxcalc, 0.5	) ;
				//gsl_complex jycalc  = gsl_complex_sub( jpcalc, jmcalc	) ;	jycalc = gsl_complex_mul_imag( jycalc,-0.5	) ;
				//gsl_complex japcalc  = calcOp(  Jaop.pj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex jamcalc  = calcOp(  Jaop.mj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex jaxcalc  = gsl_complex_add( japcalc, jamcalc	) ;	jaxcalc = gsl_complex_mul_real( jaxcalc, 0.5	) ;
				//gsl_complex jaycalc  = gsl_complex_sub( japcalc, jamcalc	) ;	jaycalc = gsl_complex_mul_imag( jaycalc,-0.5	) ;
				//printf(  "jxcalc\t: %19.16f +i %19.16f\n", GSL_REAL(jxcalc) ,  GSL_IMAG(jxcalc)		) ;
				//printf( "jaxcalc\t: %19.16f +i %19.16f\n", GSL_REAL(jaxcalc) , GSL_IMAG(jaxcalc)	) ;
				//printf(  "jycalc\t: %19.16f +i %19.16f\n", GSL_REAL(jycalc) ,  GSL_IMAG(jycalc) 	) ;
				//printf( "jaycalc\t: %19.16f +i %19.16f\n", GSL_REAL(jaycalc) , GSL_IMAG(jaycalc) 	) ;
				//write_ABC(	pararesult, ic ,   "Jzxydeg" , writemodedegen , jzcalc ,	jxcalc ,	jycalc	) ;
				//write_ABC(	pararesult, ic ,  "Jazxydeg" , writemodedegen , jazcalc ,	jaxcalc ,	jaycalc	) ;

				Vecop Lop("Lop") ;
				Vecop Sop("Sop") ;
				//printf( "L S : \n" ) ;
				//gsl_complex lzcalc  = calcOp( Lop.zj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex szcalc  = calcOp( Sop.zj, ground, basis, gnd, gndblock, table ) ;
				//printf( "lzcalc\t: %19.16f +i %19.16f\n", GSL_REAL(lzcalc) , GSL_IMAG(lzcalc) ) ;
				//printf( "szcalc\t: %19.16f +i %19.16f\n", GSL_REAL(szcalc) , GSL_IMAG(szcalc) ) ;
				//gsl_complex lpcalc  = calcOp(  Lop.pj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex lmcalc  = calcOp(  Lop.mj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex lxcalc  = gsl_complex_add( lpcalc, lmcalc	) ;	lxcalc = gsl_complex_mul_real( lxcalc, 0.5	) ;
				//gsl_complex lycalc  = gsl_complex_sub( lpcalc, lmcalc	) ;	lycalc = gsl_complex_mul_imag( lycalc,-0.5	) ;
				//gsl_complex spcalc  = calcOp(  Sop.pj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex smcalc  = calcOp(  Sop.mj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex sxcalc  = gsl_complex_add( spcalc, smcalc	) ;	sxcalc = gsl_complex_mul_real( sxcalc, 0.5	) ;
				//gsl_complex sycalc  = gsl_complex_sub( spcalc, smcalc	) ;	sycalc = gsl_complex_mul_imag( sycalc,-0.5	) ;
				//printf(  "lxcalc\t: %19.16f +i %19.16f\n", GSL_REAL(lxcalc) , 	GSL_IMAG(lxcalc) ) ;
				//printf(  "sxcalc\t: %19.16f +i %19.16f\n", GSL_REAL(sxcalc) ,	GSL_IMAG(sxcalc) ) ;
				//printf(  "lycalc\t: %19.16f +i %19.16f\n", GSL_REAL(lycalc) ,	GSL_IMAG(lycalc) ) ;
				//printf(  "sycalc\t: %19.16f +i %19.16f\n", GSL_REAL(sycalc) ,	GSL_IMAG(sycalc) ) ;
				//write_ABC(	pararesult, ic ,   "Lzxydeg" , writemodedegen , lzcalc ,	lxcalc ,	lycalc	, prefixiimp) ;
				//write_ABC(	pararesult, ic ,   "Szxydeg" , writemodedegen , szcalc ,	sxcalc ,	sycalc	, prefixiimp) ;

				//printf( "LL SS : \n" ) ;
				//gsl_complex lzlzcalc = calcOpOp( Lop.zj, Lop.zj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex lplmcalc = calcOpOp( Lop.pj, Lop.mj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex lmlpcalc = calcOpOp( Lop.mj, Lop.pj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex llcalc ;
				//llcalc = gsl_complex_add( 	lplmcalc, lmlpcalc	) ;
				//llcalc = gsl_complex_mul_real(	llcalc, 0.5		) ;		gsl_complex lxylxycalc = gsl_complex_rect( GSL_REAL(llcalc), GSL_IMAG(llcalc) ) ;
				//llcalc = gsl_complex_add(	llcalc, lzlzcalc 	) ;

				//gsl_complex szszcalc = calcOpOp( Sop.zj, Sop.zj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex spsmcalc = calcOpOp( Sop.pj, Sop.mj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex smspcalc = calcOpOp( Sop.mj, Sop.pj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex sscalc ;
				//sscalc = gsl_complex_add( 	spsmcalc, smspcalc 	) ;
				//sscalc = gsl_complex_mul_real(	sscalc, 0.5		) ;		gsl_complex sxysxycalc = gsl_complex_rect( GSL_REAL(sscalc), GSL_IMAG(sscalc) ) ;
				//sscalc = gsl_complex_add(	sscalc, szszcalc 	) ;
				//printf(   "lzlzcalc\t: %19.16f +i %19.16f\n", GSL_REAL(lzlzcalc) ,	 GSL_IMAG(lzlzcalc) 	) ;
				//printf(   "szszcalc\t: %19.16f +i %19.16f\n", GSL_REAL(szszcalc) ,	 GSL_IMAG(szszcalc) 	) ;
				//printf( "lxylxycalc\t: %19.16f +i %19.16f\n", GSL_REAL(lxylxycalc) ,	 GSL_IMAG(lxylxycalc) 	) ;
				//printf( "sxysxycalc\t: %19.16f +i %19.16f\n", GSL_REAL(sxysxycalc) ,	 GSL_IMAG(sxysxycalc) 	) ;
				//printf(     "llcalc\t: %19.16f +i %19.16f\n", GSL_REAL(llcalc) ,	 GSL_IMAG(llcalc) 	) ;
				//printf(     "sscalc\t: %19.16f +i %19.16f\n", GSL_REAL(sscalc) ,	 GSL_IMAG(sscalc) 	) ;
				//write_ABC(	pararesult, ic ,   "LLdeg" , writemodedegen , llcalc ,	lzlzcalc ,	lxylxycalc	) ;
				//write_ABC(	pararesult, ic ,   "SSdeg" , writemodedegen , sscalc ,	szszcalc ,	sxysxycalc	) ;

				Vecop Mop("Mop") ;
				//printf( "M : \n" ) ;
				//gsl_complex mzcalc  = calcOp( Mop.zj, ground, basis, gnd, gndblock, table ) ;
				//printf( "mzcalc\t: %19.16f +i %19.16f\n", GSL_REAL(mzcalc) , GSL_IMAG(mzcalc) ) ;
				//gsl_complex mpcalc  = calcOp(  Mop.pj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex mmcalc  = calcOp(  Mop.mj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex mxcalc  = gsl_complex_add( mpcalc, mmcalc	) ;	mxcalc = gsl_complex_mul_real( mxcalc, 0.5	) ;
				//gsl_complex mycalc  = gsl_complex_sub( mpcalc, mmcalc	) ;	mycalc = gsl_complex_mul_imag( mycalc,-0.5	) ;
				//printf(  "mxcalc\t: %19.16f +i %19.16f\n", GSL_REAL(mxcalc) ,	GSL_IMAG(mxcalc) ) ;
				//printf(  "mycalc\t: %19.16f +i %19.16f\n", GSL_REAL(mycalc) ,	GSL_IMAG(mycalc) ) ;
				//write_ABC(	pararesult, ic ,   "Mzxydeg" , writemodedegen , mzcalc ,	mxcalc ,	mycalc	) ;

				//printf( "MM : \n" ) ;
				//gsl_complex mzmzcalc = calcOpOp( Mop.zj, Mop.zj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex mpmmcalc = calcOpOp( Mop.pj, Mop.mj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex mmmpcalc = calcOpOp( Mop.mj, Mop.pj, ground, basis, gnd, gndblock, table ) ;
				//gsl_complex MMcalc ;
				//MMcalc = gsl_complex_add( 	mpmmcalc, mmmpcalc 	) ;
				//MMcalc = gsl_complex_mul_real(	MMcalc, 0.5		) ;		gsl_complex mxymxycalc = gsl_complex_rect( GSL_REAL(MMcalc), GSL_IMAG(MMcalc) ) ;
				//MMcalc = gsl_complex_add(	MMcalc, mzmzcalc 	) ;
				//printf(   "mzmzcalc\t: %19.16f +i %19.16f\n", GSL_REAL(mzmzcalc) ,	 GSL_IMAG(mzmzcalc) 	) ;
				//printf( "mxymxycalc\t: %19.16f +i %19.16f\n", GSL_REAL(mxymxycalc) ,	 GSL_IMAG(mxymxycalc) 	) ;
				//printf(     "MMcalc\t: %19.16f +i %19.16f\n", GSL_REAL(MMcalc) ,	 GSL_IMAG(MMcalc) 	) ;
				//write_ABC(	pararesult, ic ,   "MMdeg" , writemodedegen , MMcalc ,	mzmzcalc ,	mxymxycalc	) ;

				Jaop.freeVecop() ;
				Lop.freeVecop() ;
				Sop.freeVecop() ;
				Mop.freeVecop() ;

				double mpssum;
				printf( "startMPS:" ) ;
				mpssum = obtainMPS(	 mps[iimp], ground, basis, gnd, gndblock )	; printf( "Done.\n" ) ; printf( "mpssum=%19.16f\n", mpssum ) ;
				/*
				printf( "startMPSorb:" ) ;
				mpssum = obtainMPS3orb(	mps3orb, ground, basis, gnd, gndblock )	; printf( "Done.\n" ) ; printf( "mpssum=%19.16f\n", mpssum ) ;

				for( int mu=0 ; mu<NC ; mu++ ){
					gsl_complex **rdmorb= mkgscmatrixd(4,4);
					for( int mu0=0 ; mu0<4 ; mu0++ ) for( int mu1=0 ; mu1<4 ; mu1++ )
						rdmorb[mu0][mu1] = zero;
					obtainRDMorb( rdmorb, mps3orb, Powtnb, mu );
					char prdm[1024]	;
					sprintf(prdm,"rdmorb%d",mu)	; write_rdmorb(	    pararesult , ic , rdmorb, prdm	, writemodedegen ) ;
					sprintf(prdm,"rdmorb%ddiag",mu)	; write_rdmorbdiag( pararesult , ic , rdmorb, prdm	, writemodedegen ) ;
				}
				*/

				printf("startRDM:\n" )  ;
				gsl_complex	**rdmevec	= mkgscmatrixd(Powtnc,Powtnc) ;
				double		*rdmeval	= mkvectord(Powtnc) ;
				double		jzrdm[Powtnc] ;
				int   		ncrdm[Powtnc] ;
				obtainRDM( rdm[iimp], mps[iimp], Powtnb, Powtnc ) ; 
				calcNcRDMJzRDM( ncrdm , jzrdm , Powtnc ) ;

				for(mu=0; mu<Powtnc; mu++)  { fprintf(fncrdm, "%4d\t", ncrdm[mu] ) ; } fprintf(fncrdm, "\n") ; //fflush(fncrdm) ;

				double	rdmsum  = 0 ;
				double	rdmsumd4= 0 ;
				double	jzsum   = 0 ;
				double	jzsumd4 = 0 ;
				printf("d4_mu_index : " )  ;
				for(mu=0; mu<Powtnc; mu++)  rdmsum += GSL_REAL(rdm[iimp][mu][mu]) ; 
				for(mu=0; mu<Powtnc; mu++)  jzsum  += GSL_REAL(rdm[iimp][mu][mu]) * jzrdm[mu]  ; 
				for(mu=0; mu<Powtnc; mu++)  if( ncrdm[mu] == 4 ) { printf("%d ",mu) ; jzsumd4  += GSL_REAL(rdm[iimp][mu][mu]) * jzrdm[mu]  ; }  printf("\n") ;
				for(mu=0; mu<Powtnc; mu++)  if( ncrdm[mu] == 4 ) rdmsumd4  += GSL_REAL(rdm[iimp][mu][mu]) ;
				for(mu=0; mu<Powtnc; mu++)  { fprintf( frdmdiag , "%19.16lf\t", GSL_REAL(rdm[iimp][mu][mu]) ) ; } fprintf( frdmdiag , "\n") ;
				printf("rdmdiagsum   = %19.16lf\n", rdmsum ) ;
				printf("rdmdiagsumd4 = %19.16lf\n", rdmsumd4 ) ;
				printf("jzsum        = %19.16lf\n", jzsum  ) ;
				printf("jzsumd4      = %19.16lf\n", jzsumd4) ;

				rdmclass rc ; 
				rc.allocgsl() ;
				rc.init_transformMatd4jjJ() ;
				rc.init_transformMatd4LSJ() ;
				rc.init_transformMatd4LSJcubic() ;
				rc.init_transformMatd4t2g() ;
				rc.init_transformMatd4Localdiag() ;
				rc.dotransform_jjJd4( rdm[iimp] ) ;
				gsl_complex **rdmd4jjJevec	= mkgscmatrixd(DIMd4,DIMd4);
				double *rdmd4jjJeval		= mkvectord(DIMd4) ;
				eigen_lapack( rc.rdmjjJ, rdmd4jjJeval, rdmd4jjJevec, DIMd4 ,0) ; printf("eigen_lapack rdmd4jjJ completed.\n" ) ; fflush(stdout) ; // Row-wise vectors, i.e. rdmevec[mu][] is the mu-th eigenvectors.
				rc.dotransform_LSJd4( rdm[iimp] ) ;
				rc.dotransform_LSJcubicd4( rdm[iimp] ) ;
				rc.dotransform_t2gd4( rdm[iimp] ) ;
				rc.dotransform_Localdiagd4( rdm[iimp] ) ;
				rc.checkrdmHerimitian( rdm[iimp] ) ;
				rc.write_jjJd4( pararesult , ic , writemodedegen, prefixiimp ) ;
				rc.write_LSJd4( pararesult , ic , writemodedegen, prefixiimp ) ;
				rc.write_LSJcubicd4( pararesult , ic , writemodedegen, prefixiimp ) ;
				rc.write_t2gd4( pararesult , ic , writemodedegen, prefixiimp ) ;
				rc.write_Localdiagd4( pararesult , ic , writemodedegen, prefixiimp ) ;
				//rc.obtain_jzmat_LSJd4() ;
				rc.write_target(  pararesult , ic , rdmd4jjJevec, "rdmd4jjJevec" , writemodedegen, prefixiimp ) ;
				rc.write_targetd( pararesult , ic , rdmd4jjJeval, "rdmd4jjJeval" , writemodedegen, prefixiimp ) ;
				freegscmatrixd(rdmd4jjJevec,DIMd4) ;
				free(rdmd4jjJeval) ;

				eigen_lapack( rdm[iimp], rdmeval, rdmevec, Powtnc ,0) ; printf("eigen_lapack completed.\n" ) ; //fflush(stdout) ; // Row-wise vectors, i.e. rdmevec[mu][] is the mu-th eigenvectors.
				myEEarr[degen] = rc.write_EE( pararesult , ic , rdmeval , writemodedegen, prefixiimp ) ;

				gsl_complex **rdmd4	= mkgscmatrixd(DIMd4,DIMd4);
				gsl_complex **rdmevecd4	= mkgscmatrixd(DIMd4,DIMd4);
				double *rdmevald4	= mkvectord(DIMd4) ;
				int rdmvecIndd4[DIMd4] =  { 51, 53, 54, 60, 29, 45, 43, 39, 30, 57, 27, 23, 15, 46, 58 };
				for(mu=0; mu<DIMd4; mu++) for(nu=0; nu<DIMd4; nu++) rdmd4[mu][nu] = rdm[iimp][ rdmvecIndd4[mu] ][ rdmvecIndd4[nu] ] ;printf("rdm is copied to rdmd4.\n" ) ; fflush(stdout) ; 
				for(mu=0; mu<DIMd4; mu++) rdmevald4[mu] = GSL_REAL(rdmd4[mu][mu]) ;
				rc.write_targetd( pararesult , ic , rdmevald4, "rdmdiagd4" , writemodedegen, prefixiimp ) ; for(mu=0; mu<DIMd4; mu++) rdmevald4[mu] = 0 ;
				eigen_lapack( rdmd4, rdmevald4, rdmevecd4, DIMd4 ,0) ; printf("eigen_lapack rdmd4 completed.\n" ) ; fflush(stdout) ; // Row-wise vectors, i.e. rdmevec[mu][] is the mu-th eigenvectors.
				rc.write_targetd( pararesult , ic , rdmevald4, "rdmevald4" , writemodedegen, prefixiimp ) ;
				rc.write_target(  pararesult , ic , rdmevecd4, "rdmevecd4" , writemodedegen, prefixiimp ) ;
				printf("TESTrdmevecd4" ) ;
				rc.dotransform_testd4unit( rdmevecd4 ) ;

				gsl_complex **rdmevecd4jjJ = mkgscmatrixd(DIMd4,DIMd4) ;
				printf("TESTtransformMatd4jjJ" ) ;
				rc.dotransform_testd4unit( rc.transformMatd4jjJ ) ;
				rc.dotransformevecd4( rdmevecd4, rdmevecd4jjJ, rc.transformMatd4jjJ ) ;
				printf("TESTrdmevecd4jjJ" ) ;
				rc.dotransform_testd4unit( rdmevecd4jjJ ) ;
				rc.write_target( pararesult , ic , rdmevecd4jjJ, "rdmevecd4jjJ" , writemodedegen, prefixiimp ) ;
				freegscmatrixd(rdmevecd4jjJ,DIMd4) ;

				printf("TESTtransformMatd4LSJ" ) ;
				rc.dotransform_testd4unit( rc.transformMatd4LSJ ) ;
				gsl_complex **rdmevecd4LSJ = mkgscmatrixd(DIMd4,DIMd4) ;
				rc.dotransformevecd4( rdmevecd4, rdmevecd4LSJ, rc.transformMatd4LSJ ) ;
				printf("TESTrdmevecd4LSJ" ) ;
				rc.dotransform_testd4unit( rdmevecd4LSJ ) ;
				rc.write_target( pararesult , ic , rdmevecd4LSJ, "rdmevecd4LSJ" , writemodedegen, prefixiimp ) ;
				freegscmatrixd(rdmevecd4LSJ,DIMd4) ;

				printf("TESTtransformMatd4t2g" ) ;
				rc.dotransform_testd4unit( rc.transformMatd4t2g ) ;
				gsl_complex **rdmevecd4t2g = mkgscmatrixd(DIMd4,DIMd4) ;
				rc.dotransformevecd4( rdmevecd4, rdmevecd4t2g, rc.transformMatd4t2g ) ;
				printf("TESTrdmevecd4t2g" ) ;
				rc.dotransform_testd4unit( rdmevecd4t2g ) ;
				rc.write_target( pararesult , ic , rdmevecd4t2g, "rdmevecd4t2g" , writemodedegen, prefixiimp ) ;
				freegscmatrixd(rdmevecd4t2g,DIMd4) ;

				freegscmatrixd(rdmd4,DIMd4) ;
				freegscmatrixd(rdmevecd4,DIMd4) ;
				free(rdmevald4) ;
				rc.freegsl() ;


				rdmsum =0 ;
				printf("rdmeval:\n") ;
				for(mu=0; mu<Powtnc; mu++) {
					printf("%19.16lf ", rdmeval[mu] ) ;
					fprintf( frdmeval, "%19.16lf\t", rdmeval[mu] ) ;
					rdmsum += rdmeval[mu] ;
				}
				printf("\nrdmevalsum = %19.16lf\n", rdmsum ) ;
				fprintf( frdmeval, "%19.16lf\n", rdmsum ) ;
				double ncrdmevec[Powtnc] ;
				for(mu=0; mu<Powtnc; mu++) {
					double probnu=0 ; 
					for(nu=0; nu<Powtnc; nu++) {
						probnu += gsl_complex_abs2( rdmevec[mu][nu] ) * ncrdm[nu] ;
					}
					fprintf( fncrdmevec, "%19.16lf\t", probnu ) ;
					ncrdmevec[mu] = probnu ;
				}
				double ncimprdm = 0 ;
				for(mu=0; mu<Powtnc; mu++) {
					ncimprdm += rdmeval[mu]*ncrdmevec[mu] ;
					//printf( "%19.16lf ; %19.16lf ; %19.16lf \n", rdmeval[mu], ncrdmevec[mu], ncimprdm ) ;
				}
				printf( "ncimprdm : %19.16lf\n", ncimprdm ) ;
				fprintf( fncrdmevec, "\n" ) ; fflush(fncrdmevec) ;

				printf( "jsqrdm : \n" ) ;
				gsl_complex jsqrdm=zero ; 
				for(mu=0; mu<Powtnc; mu++) {
					printf( "\rMU : %d/%d ", mu, Powtnc ) ; //fflush(stdout) ;
					unsigned int  cup = (mu>>NC) ,
						      cdn = mu - (cup<<NC) ;
					for( int aa=0 ; aa<Jop.zj.ndat ; aa++ )  {
						int mu1 = Jop.zj.cDat[aa].i ;
						int nu1 = Jop.zj.cDat[aa].j ;
						double jz1 = Jop.zj.cDat[aa].dat ;
						for( int bb=0 ; bb<Jop.zj.ndat ; bb++ )  {
							int mu2 = Jop.zj.cDat[bb].i ;
							int nu2 = Jop.zj.cDat[bb].j ;
							double jz2 = Jop.zj.cDat[bb].dat ;
							bin[0] = cup ; 
							bin[1] = cdn ; 
							if( One( bin[nu2%2], nu2/2 ) ) {
								phase  = permu( bin, nu2/2, nu2%2);  bin[nu2%2] = Turnoff( bin[nu2%2], nu2/2 );
								if( Zero( bin[mu2%2], mu2/2 ) ) {  
									phase *= permu( bin, mu2/2, mu2%2);  bin[mu2%2] = Turnon(  bin[mu2%2], mu2/2 );
									if( One( bin[nu1%2], nu1/2 ) ) {
										phase *= permu( bin, nu1/2, nu1%2);  bin[nu1%2] = Turnoff( bin[nu1%2], nu1/2 );
										if( Zero( bin[mu1%2], mu1/2 ) ) {  
											phase *= permu( bin, mu1/2, mu1%2);  bin[mu1%2] = Turnon(  bin[mu1%2], mu1/2 );
											nu = (bin[0]<<NC) + (bin[1]) ;
											jsqrdm = gsl_complex_add( jsqrdm ,
													gsl_complex_mul_real( rdm[iimp][nu][mu] , phase*jz1*jz2 ) 
													) ;
											//printf("jz(%d,%d,%d,%d,%.2g)",mu1,nu1,mu2,nu2,GSL_IMAG(rdm[iimp][nu][mu]) );
										}
									}
								}
							}
						}
					}
					for( int aa=0 ; aa<Jop.pj.ndat ; aa++ )  {
						int mu1 = Jop.pj.cDat[aa].i ;
						int nu1 = Jop.pj.cDat[aa].j ;
						double jp1 = Jop.pj.cDat[aa].dat ;
						for( int bb=0 ; bb<Jop.mj.ndat ; bb++ )  {
							int mu2 = Jop.mj.cDat[bb].i ;
							int nu2 = Jop.mj.cDat[bb].j ;
							double jm2 = Jop.mj.cDat[bb].dat ;
							bin[0] = cup ; 
							bin[1] = cdn ; 
							if( One( bin[nu2%2], nu2/2 ) ) {
								phase  = permu( bin, nu2/2, nu2%2);  bin[nu2%2] = Turnoff( bin[nu2%2], nu2/2 );
								if( Zero( bin[mu2%2], mu2/2 ) ) {  
									phase *= permu( bin, mu2/2, mu2%2);  bin[mu2%2] = Turnon(  bin[mu2%2], mu2/2 );
									if( One( bin[nu1%2], nu1/2 ) ) {
										phase *= permu( bin, nu1/2, nu1%2);  bin[nu1%2] = Turnoff( bin[nu1%2], nu1/2 );
										if( Zero( bin[mu1%2], mu1/2 ) ) {  
											phase *= permu( bin, mu1/2, mu1%2);  bin[mu1%2] = Turnon(  bin[mu1%2], mu1/2 );
											nu = (bin[0]<<NC) + (bin[1]) ;
											jsqrdm = gsl_complex_add_real( jsqrdm ,
													GSL_REAL(rdm[iimp][nu][mu]) * phase*jp1*jm2  // 0.5*2Re=1Re
													) ;
											//printf("jpm(%d,%d,%.2g)",mu,nu,rdm[iimp][nu][mu] );
										}
									}
								}
							}
						}
					}
				}
				printf("\r%19.16lf +i %19.16lf\n", GSL_REAL(jsqrdm) , GSL_IMAG(jsqrdm) ) ;

				printf( "jsqdiagrdmevec : \n" ) ;
				double jsqdiagrdmevec[Powtnc] ;
				double jzsqdiagrdmevec[Powtnc] ;
				double jzdiagrdmevec[Powtnc] ;
				for(mu=0; mu<Powtnc; mu++) {
					jsqdiagrdmevec[mu] = 0 ;
					jzsqdiagrdmevec[mu] = 0 ;
					jzdiagrdmevec[mu] = 0 ;
					for(nu=0; nu<Powtnc; nu++) {
						unsigned int  cup = (nu>>NC) ,
							      cdn = nu - (cup<<NC) ;
						for( int aa=0 ; aa<Jop.zj.ndat ; aa++ )  {
							int mu1 = Jop.zj.cDat[aa].i ;
							int nu1 = Jop.zj.cDat[aa].j ;
							double jz1 = Jop.zj.cDat[aa].dat ;
							bin[0] = cup ; 
							bin[1] = cdn ; 
							if( One( bin[nu1%2], nu1/2 ) ) {
								phase = permu( bin, nu1/2, nu1%2);  bin[nu1%2] = Turnoff( bin[nu1%2], nu1/2 );
								if( Zero( bin[mu1%2], mu1/2 ) ) {  
									phase *= permu( bin, mu1/2, mu1%2);  bin[mu1%2] = Turnon(  bin[mu1%2], mu1/2 );
									int ou = (bin[0]<<NC) + (bin[1]) ;
									gsl_complex dd = gsl_complex_mul_real( 
											gsl_complex_mul( rdmevec[mu][nu] ,
												gsl_complex_conjugate( rdmevec[mu][ou] ) 
												) ,
											phase*jz1
											) ;
									jzdiagrdmevec[mu] += GSL_REAL( dd ) ;
								}
							}
							for( int bb=0 ; bb<Jop.zj.ndat ; bb++ )  {
								int mu2 = Jop.zj.cDat[bb].i ;
								int nu2 = Jop.zj.cDat[bb].j ;
								double jz2 = Jop.zj.cDat[bb].dat ;
								bin[0] = cup ; 
								bin[1] = cdn ; 
								if( One( bin[nu2%2], nu2/2 ) ) {
									phase  = permu( bin, nu2/2, nu2%2);  bin[nu2%2] = Turnoff( bin[nu2%2], nu2/2 );
									if( Zero( bin[mu2%2], mu2/2 ) ) {  
										phase *= permu( bin, mu2/2, mu2%2);  bin[mu2%2] = Turnon(  bin[mu2%2], mu2/2 );
										if( One( bin[nu1%2], nu1/2 ) ) {
											phase *= permu( bin, nu1/2, nu1%2);  bin[nu1%2] = Turnoff( bin[nu1%2], nu1/2 );
											if( Zero( bin[mu1%2], mu1/2 ) ) {  
												phase *= permu( bin, mu1/2, mu1%2);  bin[mu1%2] = Turnon(  bin[mu1%2], mu1/2 );
												int ou = (bin[0]<<NC) + (bin[1]) ;
												gsl_complex dd = gsl_complex_mul_real( 
														gsl_complex_mul( rdmevec[mu][nu] ,
															gsl_complex_conjugate( rdmevec[mu][ou] ) 
															) ,
														phase*jz1*jz2
														) ;
												jsqdiagrdmevec[mu] += GSL_REAL( dd ) ;
											}
										}
									}
								}
							}
						}
						jzsqdiagrdmevec[mu] = jsqdiagrdmevec[mu] ;
						for( int aa=0 ; aa<Jop.pj.ndat ; aa++ )  {
							int mu1 = Jop.pj.cDat[aa].i ;
							int nu1 = Jop.pj.cDat[aa].j ;
							double jp1 = Jop.pj.cDat[aa].dat ;
							for( int bb=0 ; bb<Jop.mj.ndat ; bb++ )  {
								int mu2 = Jop.mj.cDat[bb].i ;
								int nu2 = Jop.mj.cDat[bb].j ;
								double jm2 = Jop.mj.cDat[bb].dat ;
								bin[0] = cup ; 
								bin[1] = cdn ; 
								if( One( bin[nu2%2], nu2/2 ) ) {
									phase  = permu( bin, nu2/2, nu2%2);  bin[nu2%2] = Turnoff( bin[nu2%2], nu2/2 );
									if( Zero( bin[mu2%2], mu2/2 ) ) {  
										phase *= permu( bin, mu2/2, mu2%2);  bin[mu2%2] = Turnon(  bin[mu2%2], mu2/2 );
										if( One( bin[nu1%2], nu1/2 ) ) {
											phase *= permu( bin, nu1/2, nu1%2);  bin[nu1%2] = Turnoff( bin[nu1%2], nu1/2 );
											if( Zero( bin[mu1%2], mu1/2 ) ) {  
												phase *= permu( bin, mu1/2, mu1%2);  bin[mu1%2] = Turnon(  bin[mu1%2], mu1/2 );
												int ou = (bin[0]<<NC) + (bin[1]) ;
												gsl_complex dd = gsl_complex_mul_real( 
														gsl_complex_mul( rdmevec[mu][nu] ,
															gsl_complex_conjugate( rdmevec[mu][ou] ) 
															) ,
														phase*jp1*jm2
														) ;
												jsqdiagrdmevec[mu] += GSL_REAL( dd ) ;
											}
										}
									}
								}
							}
						}
					}
					printf( "%19.16lf(%19.16lf,%.2g)\t", jsqdiagrdmevec[mu], jzsqdiagrdmevec[mu] , ncrdmevec[mu] ) ;
				}
				printf( "\n" ) ;
				printf( "jzdiagrdmevec : \n" ) ;
				jzsum = 0 ;
				for(mu=0; mu<Powtnc; mu++) { 
					printf( "%19.16lf\t", jzdiagrdmevec[mu] ) ;
					jzsum += jzdiagrdmevec[mu] ;
				}
				printf( "\n" ) ;
				printf( "jzsumrdm : %19.16lf\n" , jzsum ) ;
				printf( "\n" ) ;

				printf( "mpssum : %19.16lf\n", mpssum ) ;
				for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
					//printf( "MUNU : %d %d\r", mu,nu ) ; //fflush(stdout) ;
					for(i = gnd; i<gnd+gndblock; i++){
						bin[0] = basis[i][0];
						bin[1] = basis[i][1];
						unsigned int  tstate = (bin[0]<<Ns) + bin[1]  ;
						unsigned int  bup = (bin[0]>>NC) ,
							      bdn = (bin[1]>>NC) ;
						unsigned int  bstate = (bup<<NB) + bdn ;
						unsigned int  cup = bin[0] - (bup<<NC) ,
							      cdn = bin[1] - (bdn<<NC) ;
						unsigned int  cstate = (cup<<NC) + cdn ;
						unsigned int  CCtstate ;
						unsigned int  CCbstate ;
						unsigned int  CCcstate ;
						if( One( basis[i][mu%2], mu/2) && One( basis[i][nu%2], nu/2 ) && (mu!=nu) ) {
							bin[0] = basis[i][0];
							bin[1] = basis[i][1];
							phase  = permu( bin, nu/2, nu%2);	bin[nu%2] = Turnoff( bin[nu%2], nu/2 );
							phase *= permu( bin, mu/2, mu%2);	bin[mu%2] = Turnoff( bin[mu%2], mu/2 );
							CCtstate = (bin[0]<<Ns) + bin[1]  ;
							unsigned int  CCbup = (bin[0]>>NC) ,
								      CCbdn = (bin[1]>>NC) ;
							CCbstate = (CCbup<<NB) + CCbdn ;
							if( bstate!=CCbstate ) printf("bath isn't conserved.\n") ;
							unsigned int  CCcup = bin[0] - (bup<<NC) ,
								      CCcdn = bin[1] - (bdn<<NC) ;
							CCcstate = (CCcup<<NC) + CCcdn ;
							susCC[mu][nu] = gsl_complex_add_real(
									susCC[mu][nu] , 
									gsl_complex_abs2( mps[iimp][bstate][CCcstate] )
									//gsl_complex_mul(
									//	mps[iimp][bstate][cstate] ,
									//	gsl_complex_mul_real( 
									//		gsl_complex_conjugate(mps[iimp][bstate][CCcstate]) ,
									//		phase )
									//	)
									) ;
						}
					}
				}
				for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
					double im = GSL_IMAG(susCC[mu][nu]) ;
					//fprintf(fsusCC, "%lf\t", GSL_REAL(susCC[mu][nu]) );
					if( fabs(im) > 1e-5 )  printf( "Non-Hermitian susCC obtained.:im=%19.16lf\n",im ) ;
				}
				//fprintf(fsusCC, "\n") ;

				freegscmatrixd( rdmevec ,Powtnc ) ;
				free( rdmeval ) ;
				free(ground);
				sprintf(writemodedegen, "a") ;
			}
            write_EEall( pararesult, ic, myweightarr, myEEarr, writemode, degeneracy[iimp], prefixiimp );
            free(myenergyarr);
            free(myweightarr);
            free(myEEarr    );

			fclose(fground);
			
			for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	excitonic[iimp][mu][nu]	= gsl_complex_div_real( excitonic[iimp][mu][nu],		degeneracy[iimp] );
			for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	excitonicT[iimp][mu][nu]	= gsl_complex_div_real( excitonicT[iimp][mu][nu],		degeneracy[iimp] );
			for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	excitonict2g[iimp][mu][nu]	= gsl_complex_div_real( excitonict2g[iimp][mu][nu],	degeneracy[iimp] );
			for(mu=0; mu<nctot; mu++){
				mag[iimp][mu] /= degeneracy[iimp];
				doccu[iimp][mu] /= degeneracy[iimp];
				filling[iimp][mu][0] /= degeneracy[iimp];
				filling[iimp][mu][1] /= degeneracy[iimp];
				sum[0] += mag[iimp][mu];
				sum[1] += doccu[iimp][mu];
				sum[2] += filling[iimp][mu][0] + filling[iimp][mu][1];
			}

			//fprintf(fmag, "%lf\t\t", sum[0]);
			//fprintf(fdoc, "%lf\t\t", sum[1]);
			//fprintf(ffil, "%lf\t\t", sum[2]);
			//for(mu=0; mu<nctot; mu++) fprintf(fmag, "%lf\t", mag[iimp][mu]);
			//for(mu=0; mu<nctot; mu++) fprintf(fdoc, "%lf\t", doccu[iimp][mu]);
			//for(mu=0; mu<nctot; mu++) fprintf(ffil, "%lf\t%lf\t\t", filling[iimp][mu][0], filling[iimp][mu][1]);
			for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
				fprintf(fext, "%lf\t%lf\t\t", GSL_REAL(excitonic[iimp][mu][nu]), GSL_IMAG(excitonic[iimp][mu][nu]));
				if( ic==count && gsl_complex_abs2( excitonic[iimp][mu][nu] ) > 1e-6 )	final_nonzero = 1;
			}
			write_excT(	pararesult, ic , excitonicT[iimp] , writemode ) ;
			write_excname(	pararesult, ic , excitonict2g[iimp] , "excitonict2g" , writemode ) ;
			sprintf(writemode, "a") ;

            fclose(frdmeval);
            fclose(frdmdiag);
            //fclose(fsusCC);
            fclose(fncrdm);
            fclose(fncrdmevec);
            fclose(frdmt2gdiag);
		}
		//fprintf(fmag, "\n");
		//fprintf(fdoc, "\n");
		//fprintf(ffil, "\n");
		fprintf(fext, "\n");
	}
	if(final_nonzero) for(int iimp=0; iimp<Ni; iimp++)	print_gscmatrixd("excitonic",  excitonic[iimp], tNC, tNC);
	for(int iimp=0; iimp<Ni; iimp++)			print_gscmatrixd("excitonicT", excitonicT[iimp], tNC, tNC);
	freetritensord(ginformation, Ni, 50*SDmax);
	//fclose(fmag);
	//fclose(fdoc);
	//fclose(ffil);
	fclose(frdminfo);
	fclose(fext);
	freegsctritensord(excitonic,	Ni, tNC);
	freegsctritensord(excitonicT,	Ni, tNC);
	freegsctritensord(excitonict2g,	Ni, tNC);
	freegsctritensord(rdm , Ni, Powtnc );
	freegsctritensord(mps , Ni, Powtnb );
	freegsctetratensord(mps3orb , 4,4,4 );
	freegscmatrixd(susCC , tNC ) ;
	freegscmatrixd(transform, NU);
	free(jsq) ;
	for(int iimp=0; iimp<Ni; iimp++){
		free( egv[iimp].egbath );
		freegscmatrixd( egv[iimp].hybrid, tNC );
	}
	free(table);
    Jop.freeVecop() ;
	for(int iimp=0; iimp<Powns2; iimp++)	free(basis_t2g[iimp]);		free(basis_t2g);
}

