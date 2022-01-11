#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <unistd.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sort.h>

#include "f2c.h"
#include "ran2.c"

#include "basic.h"
#include "matrix.h"
#include "lattice.h"
#include "inverse.c"
#include "hop.c"
#include "hamil.h"
#include "arp.h"
#include "weighting.cpp"
#include "hybblock.cpp"


double find_gm(sectorinfo ginfo[], int argument[], int forced_gindex, int forced_ndeg);
void build_ginfo(int id_imp, char *save_directory, int count, typebasis **basis, PNSpara egv, sectorinfo *ginfo, int *laniter, int target_sector, int first);
void build_argument(int argument[], sectorinfo ginfo[], double gm, int *degenerac, int forced_gindex, int forced_ndegy);

double find_ground(int id_imp, char *save_directory, int count, gsl_complex **groundt, typebasis **basis, sectorinfo ginfo[], int gindex, PNSpara egv, int *lancount, int tol, int CALL_ARPACK);
int init_table(Glist *table);
int print_table(Glist *table);
int free_table(Glist *table);
void compute_self(typebasis **basis, int argument[], sectorinfo ginfo[], PNSpara egv, int degeneracy, int count, char *save_directory, gsl_complex ***Self, Glist *table, int n);
int compute_gnew(gsl_complex *hopmatrix, gsl_complex ****Self, int count, char *save_directory, Quad *quadrature);
int build_bath(int n, gsl_complex ***Green_bath, PNSpara egv);
int continued_fraction(gsl_complex ***Green, int mu, int nu, int direction, double *ai, double *bi, double gm, double p0inner, int dimension);

double func_zeroSOC(double p[], int n);
double my_f_zeroSOC(const gsl_vector *v, void *params);
void my_df_zeroSOC(const gsl_vector *v, void *params, gsl_vector *df);
void my_fdf_zeroSOC(const gsl_vector *x, void *params, double *f, gsl_vector *df);

double func(double p[], int n);
double my_f(const gsl_vector *v, void *params);
void my_df(const gsl_vector *v, void *params, gsl_vector *df);
void my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df);

int corcd(typebasis **basis, sectorinfo ginfo[], int gindex, int direction, int mu, int spin, gsl_complex *p0, gsl_complex coeff, gsl_complex *ground);
int build_p0_gen(typebasis **basis, sectorinfo ginfo[], int gindex, gsl_complex *ground, gsl_complex *p0, double *p0inner, Glist *table);
int tosector(typebasis **basis, sectorinfo ginfo[], int gindex, int spin, int direction);


void minimize(int count, PNSpara egv, double *minimal, int *iter, double converge, int n);
void minimize_zeroSOC(int count, PNSpara egv, double *minimal, int *iter, double converge, int n, int sp);

int save_aibi(Glist *table, char *save_directory, int count, int degeneracy, int n);
void save_parameters(sectorinfo **ginfo, int **argument, PNSpara *egv, int count, int *degeneracy, int tol, FILE *fsave);

int print_gnew(FILE *file, gsl_complex ***Green);
int print_green(FILE *file, gsl_complex ***Green);
int print_self(FILE *file, gsl_complex ***Green, int start, int end);
int print_lattice(FILE *file, gsl_complex ***Green, int start, int end);
double normalize_gsc(gsl_complex *vec, int block);
void print_complex(gsl_complex cmnumber);

long ranseed = -1;
//long seed = -134551;
double SOC, J, DELTA, Uinter, MU, Interval, Hz, JTD;
gsl_complex **Ecluster, cone;
double *Matsu;
double filling[Ni], doccu[Ni];
int beta, Nmax, Lanmax, Blocks, Spinmax, nset, nctot, Powns, Powns2, myrank, size, nb, tnb, Njj, Ndd, tol_small, Nonzero_SOC, ROTATE, beta_real = 0, Vmax = SDmax;
gsl_complex ****Self, ****Greennew, zero, *hopmatrix;
MPI_Datatype MPI_GSL_COMPLEX;
PNSpara fegv;
gsl_complex ***Tbath;

int main(int argc, char *argv[]){
	int i, j, k, n, count_init, count, titer, iter[Ni] = {0}, cont=0;
	int degeneracy[Ni]={0}, tol, Repeat;
	int **laniter;
	int forced_gindex=0, forced_ndeg=0;

	typebasis **basis;
	//sectorinfo ginfo[Ni][2*Ns+1];
	sectorinfo **ginfo;
	PNSpara egv[Ni];

	double U_init, U_final, difference[Ni] = {0}, dummy, re, im, Converge, gm, ratio, Efermi;
	//Converge = tolerance for the function 'minimize'
	double *p = mkvectord(Ni*Np);

	FILE *fd, *fp=NULL, *fmain=NULL, *fsave, *fresult; 
	char save_directory[1024]={0}, prmain[1024], parameters[1024], para[1024], saveEV[1024], paratemp[1024];


	gsl_complex ****Greenold;

	MPI_Datatype itype[2] = {MPI_DOUBLE, MPI_DOUBLE};
	MPI_Aint idisp[2] = {0, sizeof(double)};
	int iblock[2] = {1, 1};//real and imaginary part, respectively

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Type_create_struct(2, iblock, idisp, itype, &MPI_GSL_COMPLEX);
	MPI_Type_commit( &MPI_GSL_COMPLEX );
	if(argc != 16 && argc != 18){
		printf("Usage :: %s <inputfilename> <Hz> <JTD> <SOC> <J> <U_init> <U_final> <DELTA> <Lanmax> <beta> <Nmax> <ratio new/old> <Repeat> <Interval> <tol_small> <optional: forced_gindex> <optional: forced_ndeg>\n", argv[0]);
		printf("Usage :: %s <inputfilename> <Hz> <JTD> <SOC> <J> <U_init> <U_final> <DELTA> <Lanmax> <beta> <Nmax> <ratio new/old> <Repeat> <Interval> <tol_small>  <optional: beta_real> <optional: number_max_states>\n", argv[0]);
		exit(1);
	}
	zero = gsl_complex_rect(0,0);
	cone = gsl_complex_rect(1,0);

	printf("OMPNUMTHREADS initilized : %d\n", OMPNUMTHREADS ) ;

	//parameters
	nb=NB;
	tnb=2*NB;
	nctot=NC;
	Spinmax	= 2*Ns+1;
	Powns	= (int)pow(2,Ns);
	Powns2	= (int)pow(2,2*Ns);
	tol = 8;
	Converge = pow( 10, -tol );				//tolerance for the function 'minimize'
	//nset	= nb/4;
	if(Ns != nctot+nb ){
		printf("'Ns' should be equal to 'nctot+nb'!\n");
		exit(1);
	}

	//read parameters
	Hz	= atof( argv[2] );
	JTD	= atof( argv[3] );
	SOC	= atof( argv[4] );
	J	= atof( argv[5] );
	U_init	= atof( argv[6] );
	U_final	= atof( argv[7] );
	DELTA	= atof( argv[8] );
	Lanmax	= atoi( argv[9] );
	beta	= atoi( argv[10] );
	Nmax	= atoi( argv[11] );
	ratio	= atof( argv[12] );
	Repeat	= atoi( argv[13] );
	Interval = atof( argv[14] );
	tol_small = atof( argv[15] );
	if( argc == 18 ){
		i = atoi( argv[17] );
		if( i > 10 ){
			beta_real = atoi( argv[16] );
			Vmax = i;
		}
		else{
			forced_gindex = atoi( argv[16] );
			forced_ndeg = atoi( argv[17] );
		}
	}
	int **argument = mkmatrixi(Ni, Vmax);
	for(n=0; n<Ni; n++) for(i=0; i<Vmax; i++)	argument[n][i] = 0;
	//end of read parameters

	char Shead;
	if( fabs(SOC) < 1e-6 ){
		ROTATE = 0;
		Nonzero_SOC = 0;
		Shead = 'o';
		Blocks = (Ns+1)*(Ns+1);
		if( PARAMAGNETIC )	Shead = 'n';
	}
	else{
		ROTATE = 1;
		Nonzero_SOC = 1;
		Shead = 'i';
		Blocks	= 2*Ns+1;
		if( PARAMAGNETIC )	Shead = 'j';
	}

	//initialize
	fd = fopen(argv[1], "r");
	nofile(fd, argv[1]);
	fscanf(fd, "%d", &count_init);
	fscanf(fd, "%lf", &Uinter);
	fscanf(fd, "%lf", &MU);
	fscanf(fd, "%d", &i);
	fscanf(fd, "%d", &n);
	fscanf(fd, "%d", &n);

	for(n=0; n<Ni; n++){
		egv[n].egbath = mkvectord(tnb);	
		egv[n].hybrid = mkgscmatrixd(tNC, tnb);
		fscanf(fd, "%d", &degeneracy[n]);
		for(i=0; i<degeneracy[n]; i++) for(j=0; j<3; j++)	fscanf(fd, "%lf", &dummy);
		for(k=0; k<tnb; k++){
			fscanf(fd, "%lf", &dummy);
			egv[n].egbath[k] = dummy;
		}
		for(j=0; j<tNC; j++) for(k=0; k<tnb; k++){
			fscanf(fd, "%lf %lf", &re, &im);
			egv[n].hybrid[j][k] = gsl_complex_rect( re, im );
		}
	}
	fclose(fd);
	//end of initialize


	Tbath = mkgsctritensord(Nmax, tNC, tNC);
	fegv.egbath = mkvectord(tnb);
	fegv.hybrid = mkgscmatrixd(tNC, tnb);

	Ecluster = mkgscmatrixd( NU, NU);
	basis = mkmatrixb(Powns2, 2);
	laniter = mkmatrixi(Ni, Blocks);
	Greennew = mkgsctetratensord(Ni, Nmax, tNC, tNC);
	Greenold = mkgsctetratensord(Ni, Nmax, tNC, tNC);
	Self	 = mkgsctetratensord(Ni, Nmax, tNC, tNC);

	Matsu = mkvectord(Nmax);
	for(i=0; i<Nmax; i++)	wn = (double)( (2*i+1) * PI/beta );
	initWwn( Matsu, Nmax );

	if( USEHYBBLOCK ) {
		infoHybblock iHyb( "inputs/HYBBLOCK" );	if(myrank<1) iHyb.show();
		int nblock = iHyb.nblock;
		int nSpinFit = -Nonzero_SOC+2;
		for(n=0; n<Ni; n++){
			PNSpara egvtmp;
			egvtmp.egbath = mkvectord(tnb);
			egvtmp.hybrid = mkgscmatrixd(tNC, tnb);
			for(k=0; k<tnb; k++){
				egvtmp.egbath[k] = 0.;
				for(j=0; j<tNC; j++) {
					GSL_REAL( egvtmp.hybrid[j][k] )  = 0.;
					GSL_IMAG( egvtmp.hybrid[j][k] )  = 0.;
				}
			}
			for(int sp=0; sp<nSpinFit; sp++) {
				for(int iblock=0; iblock<nblock; iblock++){
					int nInd	= iHyb.blockIndArr[iblock].size();
					int indArr[tNC];
					for(i=0; i<nInd; i++ ) indArr[i] = iHyb.blockIndArr[iblock][i] + sp ;

					int nIndBath	= iHyb.nbathBlockArr[iblock];
					int indBathArr[Np];
					for(i=0; i<nIndBath; i++ ) indBathArr[i] = iHyb.bathBlockIndArr[iblock][i]*(iHyb.ispin+1) + sp ;
					convert_egv_to_p_blockSym( egv[n], p, sp, iblock,
						       	nInd, indArr, nIndBath, indBathArr
							);
					convert_p_to_egv_blockSym( egvtmp, p, sp, iblock,
						       	nInd, indArr, nIndBath, indBathArr
							);
				}
			}
			for(k=0; k<tnb; k++){
				egv[n].egbath[k] = egvtmp.egbath[k] ;
				for(j=0; j<tNC; j++) {
					GSL_REAL( egv[n].hybrid[j][k] )  = GSL_REAL( egvtmp.hybrid[j][k] );
					GSL_IMAG( egv[n].hybrid[j][k] )  = GSL_IMAG( egvtmp.hybrid[j][k] );
				}
			}
			free(egvtmp.egbath);
			freegscmatrixd(egvtmp.hybrid, tNC);
			if( myrank == 0 ) print_egv("Convert (USEHYBBLOCK=1)", egv[n]);
		}
	}
	else { 
		//apply constraint
		for(n=0; n<Ni; n++){
			convert_egv_to_p(egv[n], p);
			convert_p_to_egv(egv[n], p);
			if( myrank == 0 ) print_egv("Convert", egv[n]);
		}

		//for(n=0; n<Ni; n++){
		//	convert_egv_to_p(egv[n], p);
		//	convert_p_to_egv(egv[n], p);
		//	if( myrank == 0 ) print_egv("Convert", egv[n]);
		//}
		//end of constraint
	}
	char tails[1024];
	int mytail;

	if( forced_gindex ){
		sprintf(tails,
				"%cDir_OPT%dVec%d_UF%.2lf_UT%.2lf_ns%d_nc%dnb%d_beta%d_Nmax%d_Hz%.3lf_JT%.3lf_S%.3lf_J%.3lf_D%.3lf_mix%.2lf_tol%d_Nint%d_n%d_R%d_ts%d_fg%d_fd%d",
				Shead, Tolgreen, Nvec, U_init, U_final, Ns, nctot, nb, beta, Nmax, Hz, JTD, SOC, J, DELTA, ratio, tol, Nintx, size, ROTATE, tol_small, forced_gindex, forced_ndeg
		       );
		sprintf(save_directory, "%s/%s", Path_save, tails);
	}
	else{
		sprintf(tails,
				"%cDir_OPT%dVec%d_UF%.2lf_UT%.2lf_ns%d_nc%dnb%d_beta%d_Nmax%d_Hz%.3lf_JT%.3lf_S%.3lf_J%.3lf_D%.3lf_mix%.2lf_tol%d_Nint%d_n%d_R%d_ts%d_IT%d",
				Shead, Tolgreen, Nvec, U_init, U_final, Ns, nctot, nb, beta, Nmax, Hz, JTD, SOC, J, DELTA, ratio, tol, Nintx, size, ROTATE, tol_small, beta_real
		       );
		sprintf(save_directory, "%s/%s", Path_save, tails);
	}

	if( myrank == 0 ){
		mytail = mk_savedir(save_directory);
		if( mytail )	sprintf(tails, "%s_%d", tails, mytail);
		if( fabs(Uinter-U_init) < 0.0001 && mytail )	cont = 1;
		else						cont = 0;
	}
	MPI_Bcast(save_directory, 1024, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(&cont, 1, MPI_INT, 0, MPI_COMM_WORLD);

	Quad quadrature;
	hopmatrix = mkgscvectord(Nintx*Ninty*Nintz*NU*NU);
	Efermi = compute_hopmatrix_all(hopmatrix, &quadrature);
	MU = 0;

	char paraEcluster[1024];
	gsl_complex **transform = mkgscmatrixd(NU, NU);
	initgscmatrixd(transform, NU, NU);
	for(i=0; i<NU; i++)	transform[i][i] = gsl_complex_rect(1,1);
	if( ROTATE ){
		init_transform(transform);
		//init_Ecluster( Ecluster, hopmatrix, &quadrature );
		//sprintf(paraEcluster, "%s/Ecluster_before_rotation.dat", save_directory);
		//save_Ecluster(paraEcluster, Ecluster);
		//iunitary(Ecluster, Ecluster, transform, NU);
		//sprintf(paraEcluster, "%s/Ecluster_independent_rotation.dat", save_directory);
		//save_Ecluster(paraEcluster, Ecluster);

		rotate_hopmatrix(hopmatrix, transform);
	}
	init_Ecluster( Ecluster, hopmatrix, &quadrature );
	//symmetrize bath parameters
	for(n=0; n<Ni; n++){
		if( !Nonzero_SOC )	set_inter_spin_zero(&egv[n]);
		//symmetrize_egv(&egv[n], transform);
	}
	freegscmatrixd(transform, NU);

	sprintf(paraEcluster, "%s/Ecluster.dat", save_directory);
	save_Ecluster(paraEcluster, Ecluster);
	update_chem(Ecluster);

	MPI_Barrier( MPI_COMM_WORLD);
	if( myrank == 0 ){
		if( ratio > 1 || ratio < 0 ){
			printf("ratio must be in the range [0:1]!\n");
			exit(1);
		}
		sprintf(parameters, "%s/datafile", save_directory);
		mkdirectory(parameters);
		sprintf(parameters, "%s/ongoing", save_directory);
		mkdirectory(parameters);
		sprintf(parameters, "%s/ongoing/weightwn.dat", save_directory);
		print_doubleArr( parameters , Matsu, Wwn, 0, Nmax );
		if( Debug ){
			printf("Debug output ON\n");
			sprintf(parameters, "%s/debug", save_directory);
			mkdirectory(parameters);
		}
	}
	ginfo = (sectorinfo **) malloc( Ni * sizeof(sectorinfo*));
	for(n=0; n<Ni; n++){
		ginfo[n] = (sectorinfo *) malloc( Blocks*sizeof(sectorinfo) );
		build_bitbasis(Ns, basis, ginfo[n], Nonzero_SOC);
		construct_mapnumber(Ecluster, basis, ginfo[n], 0, Blocks, &Ndd, &Njj);
		if(myrank<1) print_ginfo(Ns, ginfo[n], Blocks);
	}

	int iu, iumax, sign, sp;
	iumax = (U_final - U_init)/Interval;
	sign = iumax > 0 ? 1 : -1;
	iumax = abs(iumax)+1;
	//control_para = U_init;

	Glist table[Ni][8*nctot*nctot];
	for(n=0; n<Ni; n++){
		init_table(table[n]);
		//if( myrank == 0 )	print_table(table[n]);
	}
	for(iu=0; iu<iumax; iu++){
		Uinter = U_init + sign*Interval*iu;
		//MU = Efermi + Uinter/2. + DELTA;
		MU = Efermi + 2*Uinter - 4*J + DELTA;
		//MU = Uinter/2 + DELTA;
		init_Ecluster( Ecluster, hopmatrix, &quadrature );
		update_chem(Ecluster);
		if(myrank<1) print_gscmatrixd("Ecluster: after update_chem", Ecluster, NU, NU);

		if( cont != 1 )	count = 0;
		else {
			count = count_init;
			Repeat += count;
		}
		sprintf(parameters,	"%s/parameters_u%.2lf.dat",	save_directory, Uinter);
		sprintf(prmain,		"%s/output_u%.2lf.txt",		save_directory, Uinter);

		if( myrank == 0 ){
			fmain = fopen(prmain, "w"); nofile(fmain, "fmain");
			fp = fopen(parameters, "w"); nofile(fp, "fp");

			fprintf(fmain, "#Initial data from %s\n", argv[1]);
			fprintf(fp, "#count\t");
			for(n=0; n<Ni; n++)
				fprintf(fp, "chi2_%d\t\titer_%d\tdegen_%d\tlanc_%d\tgm_%d\t\tgptl_%d\tnpos_%d\tdouble_%d\ttfilling_%d\t", n, n, n, n, n, n, n, n, n);
			fprintf(fp, "\n");
			fflush(fp);
		}

		do{
			count++;

			int already_done, gnode;
			for(n=0; n<Ni; n++){
				build_ginfo(n, save_directory, count, basis, egv[n], ginfo[n], laniter[n], forced_gindex, count == count_init+1);
				//gm = ginfo[Ns].energy;
				gm = find_gm(ginfo[n], argument[n], forced_gindex, forced_ndeg);							//Find lowest eigenvalue.
				build_argument(argument[n], ginfo[n], gm, &degeneracy[n], forced_gindex, forced_ndeg);
				for(i=0; i<degeneracy[n]; i++){
					already_done = 0;
					for(j=0; j<i; j++) if( argument[n][i] == argument[n][j] ) already_done = 1;
					gnode	= argument[n][i]%size;
					if( !already_done && myrank == gnode ){
						find_ground(n, save_directory, count, ginfo[n][ argument[n][i] ].ground, basis, ginfo[n], argument[n][i], egv[n], &j, Tolgreen, 1);	//machin epsilon ~ 1.11e-16
					}
				}
				for(i=0; i<degeneracy[n]; i++){
					gnode	= argument[n][i]%size;
					MPI_Bcast(&ginfo[n][argument[n][i]].energy, Dmax, MPI_DOUBLE, gnode, MPI_COMM_WORLD);
					MPI_Bcast(&ginfo[n][argument[n][i]].ndeg, 1, MPI_INT, gnode, MPI_COMM_WORLD);
				}
				gm = ginfo[n][argument[n][0]].energy[0];
				if(myrank == 0)	printf("after arpack:\n");
				build_argument(argument[n], ginfo[n], gm, &degeneracy[n], forced_gindex, forced_ndeg);

				omp_set_num_threads( OMPNUMTHREADS );
				compute_self(basis, argument[n], ginfo[n], egv[n], degeneracy[n], count, save_directory, Self[n], table[n], n);
			}
			omp_set_num_threads( OMPNUMTHREADS );
			compute_gnew(hopmatrix, Self, count, save_directory, &quadrature);

			sprintf(saveEV, "u%.2lf_%dth.mak", Uinter, count);
			if( myrank == 0 ){
				if( count == 1 ) for(n=0; n<Ni; n++) for(i=0; i<Nmax; i++) for(j=0; j<tNC; j++) for(k=0; k<tNC; k++)	Greenold[n][i][j][k] = Greennew[n][i][j][k];

				for(n=0; n<Ni; n++) for(i=0; i<Nmax; i++) for(j=0; j<tNC; j++) for(k=0; k<tNC; k++){
					Greennew[n][i][j][k] = gsl_complex_add( gsl_complex_mul_real( Greenold[n][i][j][k], 1.-ratio ), gsl_complex_mul_real( Greennew[n][i][j][k], ratio ) );
					Greenold[n][i][j][k] = Greennew[n][i][j][k];
				}

				sprintf(para, "%s/%s", save_directory, saveEV);
				fsave = fopen(para, "w");
				save_parameters(ginfo, argument, egv, count, degeneracy, tol, fsave);
				fclose(fsave);

				printf("U : %.2lf, %dth, dw = %lf, omegamax = %lf\n", Uinter, count, delta_w, (2*Nmax)*delta_w);
				for(n=0; n<Ni; n++){
					convert_egv_to_p(egv[n], p);
					difference[n] = func(p, n);
					printf("chisquare of %dth site: = %lg, iter = %d\n", n, difference[n], iter[n]);
					sprintf(paratemp, "egv[%d]", n);
					print_egv(paratemp, egv[n]);
				}
				//fprint_egv(fmain, "main", egv);
			}
			titer = 1;
			if( USEHYBBLOCK ){
				infoHybblock iHyb( "inputs/HYBBLOCK" );	if(myrank<1) iHyb.show();
				int nblock = iHyb.nblock;
				int nSpinFit = -Nonzero_SOC+2;
				for(n=0; n<Ni; n++){
					for(sp=0; sp<nSpinFit; sp++) {
						for(int iblock=0; iblock<nblock; iblock++){
							omp_set_num_threads( OMPNUMTHREADS );
							minimize_block(count, egv[n], &difference[n], &iter[n], Converge, n, sp, iblock);
							titer *= iter[n];
						}
					}
				}
			}
			else if( Nonzero_SOC ){
				for(n=0; n<Ni; n++){
					omp_set_num_threads( OMPNUMTHREADS );
					minimize(count, egv[n], &difference[n], &iter[n], Converge, n);
					titer *= iter[n];
				}
			}
			else{
				for(n=0; n<Ni; n++){
					omp_set_num_threads( OMPNUMTHREADS );
					for(sp=0; sp<2; sp++)	minimize_zeroSOC(count, egv[n], &difference[n], &iter[n], Converge, n, sp);
					titer *= iter[n];
				}
			}
			if( myrank == 0 ){
				fprintf(fmain, "\n---- < U : %.2lf, %dth >---- titer = %d, dw = %lf, omegamax = %lf\n", Uinter, count, titer, delta_w, (2*Nmax)*delta_w);
				sprintf(para, "%s/next.mak", save_directory);

				fprintf(fp, "%d\t", count);
				for(n=0; n<Ni; n++){
					i=0;
					for(k=0; k<tnb; k++)	if( egv[n].egbath[k] > 0 )	i++;
					k = argument[n][0];
					fprintf(fp, "%.10lg\t%5d\t%d\t%5d\t%.10lf\t%d\t%d\t%.10lf\t%.10lf\t", difference[n], iter[n], degeneracy[n], laniter[n][k], ginfo[n][k].energy[0], ginfo[n][k].ptl, i, doccu[n], filling[n]);
				}
				fprintf(fp, "\n");
				fflush(fp);

				fsave = fopen(para, "w");
				save_parameters(ginfo, argument, egv, count+1, degeneracy, tol, fsave);
				fclose(fsave);
				fflush(fmain);
				fflush(fp);
			}

			MPI_Bcast(iter, Ni, MPI_INT, 0, MPI_COMM_WORLD);
			for(n=0; n<Ni; n++)
				convert_egv_to_p(egv[ind(n)], &p[n*Np]);
			MPI_Bcast(p, Ni*Np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			for(n=0; n<Ni; n++)
				convert_p_to_egv(egv[n], &p[n*Np]);

			if( titer == 1 || count == Repeat){
				free(hopmatrix);
				if( myrank == 0 ){
					sprintf(para, "%s/mkresult.bat", save_directory);

					fresult = fopen(para, "a");
					fprintf(fresult, "./$1 \t%s/%s \t%.3lf\t \t%.3lf \t %.6lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%s \t%d\n", tails, saveEV, Hz, JTD, SOC, J, U_init, U_final, DELTA, tails, beta_real);
					fclose( fresult );
					fclose(fmain);
					fclose(fp);
					if( Post_process ){
						printf("static: Started\n");	fflush(stdout);
						sprintf(para, "./static %s/%s \t%.3lf\t %.3lf %.6lf %.3lf %.3lf %.3lf %.3lf %s %d", tails, saveEV, Hz, JTD, SOC, J, U_init, U_final, DELTA, tails, beta_real);
						system(para);
						printf("static: Done\n");	fflush(stdout);

						printf("vert_new: Started\n");	fflush(stdout);
						sprintf(para, "./vert_new %s/%s \t%.3lf\t %.3lf %.6lf %.3lf %.3lf %.3lf %.3lf %s %d", tails, saveEV, Hz, JTD, SOC, J, U_init, U_final, DELTA, tails, beta_real);
						system(para);
						printf("vert_new: Done\n");	fflush(stdout);

						//printf("local: Started\n");	fflush(stdout);
						//sprintf(para, "./local %s/%s \t%.3lf\t %.3lf %.6lf %.3lf %.3lf %.3lf %.3lf %s %d", tails, saveEV, Hz, JTD, SOC, J, U_init, U_final, DELTA, tails, beta_real);
						//system(para);
						//printf("local: Done\n");	fflush(stdout);

						//sprintf(para, "./opt %s/%s \t%.3lf\t %.3lf %.6lf %.3lf %.3lf %.3lf %.3lf %s 64 256", tails, saveEV, Hz, JTD, SOC, J, U_init, U_final, DELTA, tails);
						//system(para);
					}
				}
				break;
			}
			else{
				if( myrank == 0 ){
					sprintf(para, "%s/temp.bat", save_directory);

					fresult = fopen(para, "w");
					fprintf(fresult, "./$1 \t%s/%s \t%.3lf \t%.3lf %.6lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%s \t%d\n", tails, saveEV, Hz, JTD, SOC, J, U_init, U_final, DELTA, tails, beta_real);
					fclose( fresult );
				}
			}
			MPI_Barrier( MPI_COMM_WORLD );
		}while( count < Repeat );
	}
	for(i=0; i<Powns2; i++)	free(basis[i]);		free(basis);
	freegscmatrixd(Ecluster, NU);
	for(n=0; n<Ni; n++){
		for(i=0; i<Blocks; i++)
			freegscmatrixd(ginfo[n][i].ground, Nvec);
		free( ginfo[n] );
		free_table(table[n]);
		free( egv[n].egbath );
		freegscmatrixd( egv[n].hybrid, tNC );
	}
	free(p);
	free(ginfo);
	freematrixi(laniter, Ni);

	freegsctetratensord( Greennew, Ni, Nmax, tNC);
	freegsctetratensord( Greenold, Ni, Nmax, tNC);
	freegsctetratensord( Self, Ni, Nmax, tNC);
	freegsctritensord(Tbath, Nmax, tNC);
	free(fegv.egbath);
	freegscmatrixd(fegv.hybrid, tNC);
	rmdir(para);
	free( Wwn );

	MPI_Finalize();
	return 0;
}

/*
void H_vector_map(int *istart, double *Hdia, gsl_complex *Hamil, int *column, double *vectorin, double *vectorout, int block){
	int i, j;

#pragma omp parallel for default(none)	\
	private(i) shared(block, vectorout, Hdia, vectorin)
	for(i=0; i<block; i++){
		vectorout[i] = Hdia[i] * vectorin[i];
	}
#pragma omp parallel for default(none)	\
	private(i,j) shared(istart, column, Hamil, vectorin, vectorout, block)
	for(j=0; j<block; j++){
		for(i=istart[j]; i<istart[j+1]; i++){
			vectorout[j] += Hamil[i] * vectorin[column[i]];
		}
	}
}//end of H_vector_map
*/

void build_ginfo(int id_imp, char *save_directory, int count, typebasis **basis, PNSpara egv, sectorinfo *ginfo, int *laniter, int target_sector, int first){
	int gindex=0, refer=0, block=0, bn, rank, tol=Tolgreen, ED_iter;
	double gm;

	if( first || !target_sector ){
		for(bn=myrank; bn<Blocks; bn+=size){
			gindex = bn;

			gm = find_ground(id_imp, save_directory, count, ginfo[gindex].ground, basis, ginfo, gindex, egv, &ED_iter, tol, 0);
			laniter[gindex] = ED_iter;

			block = ginfo[gindex].block;
			refer = ginfo[gindex].refer;

			printf("rank%2d: %5dth\t%15d\t%10d\t%lf\t%d\n", myrank, gindex, refer, block, gm, laniter[gindex]);
			fflush(stdout);
		}
	}
	else{
		gindex = target_sector;

		gm = find_ground(id_imp, save_directory, count, ginfo[gindex].ground, basis, ginfo, gindex, egv, &ED_iter, tol, 0);
		laniter[gindex] = ED_iter;

		block = ginfo[gindex].block;
		refer = ginfo[gindex].refer;

		printf("rank%2d: %5dth\t%15d\t%10d\t%lf\t%d\n", myrank, gindex, refer, block, gm, laniter[gindex]);
		fflush(stdout);
	}
	for(bn=0; bn<Blocks; bn++){
		rank = bn%size;
		MPI_Bcast(&ginfo[bn].energy, Dmax, MPI_DOUBLE, rank, MPI_COMM_WORLD);
		MPI_Bcast(&ginfo[bn].ndeg, 1, MPI_INT, rank, MPI_COMM_WORLD);
		MPI_Bcast(&laniter[bn], 1, MPI_INT, rank, MPI_COMM_WORLD);
		//printf("Bcast from rank %d completed!\n", rank);
	}
}//end of build_ginfo

double find_gm(sectorinfo ginfo[], int argument[], int forced_gindex, int forced_ndeg){
	int i=0; 
	double gm = ginfo[i].energy[0];
	for(i=1; i<Blocks; i++)	if( ginfo[i].energy[0] < gm ){
		gm = ginfo[i].energy[0];
		argument[0] = i;
	}
	if( forced_gindex ){
		i = forced_gindex;
		gm = ginfo[i].energy[0];
		argument[0] = i;
	}

	return gm;
}


void build_argument(int argument[], sectorinfo ginfo[], double gm, int *degeneracy, int forced_gindex, int forced_ndeg){
	double Tolerance, minvalue;
	int i, j, itemp, num_visit[Blocks], arg_firsts[Blocks], num_blocks;
	if( beta_real )	Tolerance = -log(0.001)/beta_real;
	else		Tolerance = pow(10, -tol_small);
	*degeneracy = 0;

	num_blocks = 0;
	for(i=0; i<Blocks; i++){
		//if( fabs(Hz) < 1e-6 && fabs(SOC) > 1e-6 && fabs(SOC) < 0.01 &&fabs(JTD) < 1e-6 && ginfo[i].ndeg > 2 )	ginfo[i].ndeg = 2;
		num_visit[i] = 0;
		if( fabs( gm-ginfo[i].energy[0] ) < Tolerance ){	//find the valid sectors
			arg_firsts[num_blocks] = i;
			num_blocks++;
		}
	}

	for(i=0; i<num_blocks; i++){					//sorting
		minvalue = ginfo[arg_firsts[i]].energy[0];
		for(j=i; j<num_blocks; j++){
			if( ginfo[arg_firsts[j]].energy[0] < minvalue ){
				minvalue = ginfo[arg_firsts[j]].energy[0];
				itemp = arg_firsts[i];
				arg_firsts[i] = arg_firsts[j];
				arg_firsts[j] = itemp;
			}
		}
	}
	int nstates = 0;
	for(i=0; i<num_blocks; i++){
		for(j=0; j<Dmax; j++){
			if( fabs(ginfo[arg_firsts[i]].energy[j] - gm) < Tolerance ){
				nstates++;
			}
		}
	}

	if( nstates >= Vmax ){
		printf("build_argument: Too many states in %d blocks, nstates %d > %d", num_blocks, nstates, Vmax);	fflush(stdout);
		exit(1);
	}

	for(i=0; i<num_blocks; i++){
		for(j=0; j<Dmax; j++){
			if( fabs(ginfo[arg_firsts[i]].energy[j] - gm) < Tolerance ){
				argument[*degeneracy] = arg_firsts[i];
				*degeneracy = *degeneracy + 1;
			}
		}
	}
	if( forced_gindex ){
		*degeneracy = forced_ndeg;
		for(i=0; i<forced_ndeg; i++)	argument[i] = forced_gindex;
		ginfo[forced_gindex].ndeg = forced_ndeg;
	}

	if( myrank == 0 ){
		printf("Tolerance in build_argument: %lf\n", Tolerance);
		for(i=0; i<*degeneracy; i++){
			printf("argument[%d] = %d\t%.14lf\n", i, argument[i], ginfo[argument[i]].energy[num_visit[argument[i]]]);
			num_visit[argument[i]] += 1;
		}
		if( *degeneracy < 1 ){
			puts("build_argument: at least a state is required");	fflush(stdout);
			exit(1);
		}
	}
}//end of build_argument



double find_ground(int id_imp, char *save_directory, int count, gsl_complex **ground, typebasis **basis, sectorinfo ginfo[], int gindex, PNSpara egv, int *lancount, int tol, int CALL_ARPACK){
	omp_set_num_threads( OMPNUMTHREADS );
	int i, block, *istart;
	int mapnumber, *column;
	double *Hdia;
	gsl_complex *Hamil;
	double Eg=0, Egt=0, overlap, Tolerance = pow(10, -tol);	

	block = ginfo[gindex].block;

	if( block == 1 ){
		compute_Hdia(id_imp, &Eg, ginfo, gindex, basis, egv);
		*lancount = 0;
		ground[0][0] = gsl_complex_rect(1,0);
		return Eg;
	}

	istart	= mkvectori(block+1);
	mapnumber = ginfo[gindex].mapnumber;
	Hdia	= mkvectord(block);
	Hamil	= mkgscvectord(mapnumber);
	column	= mkvectori(mapnumber);
	read_map(id_imp, basis, istart, Hamil, column, ginfo, gindex, egv);
	if( Debug )
		print_map(save_directory, count, id_imp, basis, istart, ginfo, gindex, egv);
	compute_Hdia(id_imp, Hdia, ginfo, gindex, basis, egv);

	int k, ndeg, iter_jd = 0;
	double res[Nvec]={10}, split=100, inn_r, inn_i, values[SDmax]={0};
	gsl_complex **gtemp, inner;
	gtemp = mkgscmatrixd( Nvec, block );


	FILE *ftemp;
	char paratemp[1024];
	if( Save_davidson ){
		sprintf(paratemp, "%s/iter/n%d_da_gindex%d.txt", save_directory, id_imp, gindex);
		ftemp = fopen(paratemp, "a");
	}
	if( CALL_ARPACK ){
		int nconv;
		int *sort_perm = mkvectori(SDmax+1);

		double *value_temp = mkvectord(SDmax+1);
		Eg = arpack_areig_Hc(istart, Hdia, Hamil, column, block, ground, &overlap, lancount, value_temp, res, Tolerance, &nconv);

		indexxd(nconv, value_temp-1, sort_perm-1);
		for(i=0; i<nconv; i++){
			values[i] = value_temp[sort_perm[i]];
			for(k=0; k<block; k++)
				gtemp[i][k] = ground[sort_perm[i]][k];
		}
		if( Save_davidson ){
			fprintf(ftemp, "%d\t%d\t", count, Lanmax);
			for(k=0; k<nconv; k++)
				fprintf(ftemp, "%.18lf\t%.18lf\t\t", values[k], res[k]);
			fprintf(ftemp, "\n\n");
			fclose(ftemp);
		}
		free(sort_perm);
		free(value_temp);
	}
	else{
		for(k=0; k<Nvec; k++){
			for(i=0; i<block; i++){
				gtemp[k][i] = gsl_complex_rect( ran2(&ranseed)-0.5, ran2(&ranseed)-0.5  );
			}
			normalize_gsc(gtemp[k], block);
		}
		*lancount = 0;
		Egt = Davidson(istart, Hdia, Hamil, column, block, gtemp, ground, &overlap, lancount, values, res, Tolerance);
		for(k=0; k<Nvec; k++) for (i=0; i<block; i++)	gtemp[k][i] = ground[k][i];
		if( Save_davidson )	fprintf(ftemp, "%d\t%d\t%.18lf\t%.18lf\t\t%.18lf\n", count, *lancount, Egt, overlap, res[0]);

		if(block >= 100 ){
			do{
				Eg = Egt;
				Egt = MLanczos_Nvector(istart, Hdia, Hamil, column, block, gtemp[0], ground[0], &overlap, lancount, values, res);
				for (i=0; i<block; i++)	gtemp[0][i] = ground[0][i];
				//for(k=0; k<Nvec; k++) for (i=0; i<block; i++)	gtemp[k][i] = ground[k][i];
				if( Save_davidson )	fprintf(ftemp, "%d\t%d\t%.18lf\t%.18lf\t\t%.18lf\t%.18lf\t%.18lf\t%d\t%d\t%d\t%lf\n", count, *lancount, Egt, overlap, fabs((Egt-Eg)/Egt), fabs((fabs(overlap)-1)), res[0], (fabs((Eg-Egt)/Egt) > Tolerance), (fabs((fabs(overlap)-1)) > Tolerance), res[0] > Tolerance, values[1]);
				//printf("Eg %.18lf, Egt %.18lf: cri %.18lf, %.18lf\n", Eg, Egt, fabs((Eg-Egt)/Eg), (fabs(overlap)-1));
			}while( ( (fabs((Eg-Egt)/Egt) > sqrt(Tolerance)) | (fabs((fabs(overlap)-1)) > sqrt(Tolerance))) && *lancount < N_pre && res[0] > sqrt(Tolerance) );

			do{
				Eg = Egt;
				Egt = Davidson(istart, Hdia, Hamil, column, block, gtemp, ground, &overlap, lancount, values, res, Tolerance);
				for(k=0; k<Nvec; k++) for (i=0; i<block; i++)	gtemp[k][i] = ground[k][i];
				if( Save_davidson )	fprintf(ftemp, "%d\t%d\t%.18lf\t%.18lf\t\t%.18lf\t%.18lf\t%.18lf\t%d\t%d\t%d\t%lf\n", count, *lancount, Egt, overlap, fabs((Egt-Eg)/Egt), fabs((fabs(overlap)-1)), res[0], (fabs((Eg-Egt)/Egt) > Tolerance), (fabs((fabs(overlap)-1)) > Tolerance), res[0] > Tolerance, values[1]);
				//printf("Eg %.18lf, Egt %.18lf: cri %.18lf, %.18lf\n", Eg, Egt, fabs((Eg-Egt)/Eg), (fabs(overlap)-1));
			}while( ( (fabs((Eg-Egt)/Egt) > Tolerance) | (fabs((fabs(overlap)-1)) > Tolerance)) && *lancount < Lanmax && res[0] > Tolerance );
		}

		if( Save_davidson ){
			fclose(ftemp);

			sprintf(paratemp, "%s/iter/n%d_last_ptl%d_gindex%d.txt", save_directory, id_imp, ginfo[gindex].ptl, gindex);
			ftemp = fopen(paratemp, "a");
			fprintf(ftemp, "%d\t%d\t%.18lf\t%.18lf\t%.18lf\t%.18lf\t", count, *lancount, Egt, overlap, Egt-Eg, res[0]);
			for(k=1; k<Nvec; k++)	fprintf(ftemp, "%.18lf\t%.18lf\t\t", values[k], res[k]);
			fprintf(ftemp, "\n");
			fclose(ftemp);
		}
	}

	ndeg = 1;
	printf("%d:\t%.18lf, ", gindex, values[0]);
	for(k=1; k<Nvec; k++){
		split = values[k] - values[0];
		inn_r = 0; inn_i = 0;

#pragma omp parallel for default(none) private(i) shared(k, block, gtemp) reduction(+:inn_r, inn_i)
			for(i=0; i<block; i++){
				inn_r += gsl_complex_in_r( gtemp[k][i], gtemp[0][i] );
				inn_i += gsl_complex_in_i( gtemp[k][i], gtemp[0][i] );
			}
			inner = gsl_complex_rect( inn_r, inn_i );
		//if( fabs((values[k]-values[0])/values[0]) < 1e-4 && gsl_complex_abs2(inner) < 1e-6 ) 	ndeg++;
		if( fabs((values[k]-values[0])/values[0]) < pow(10, -tol_small)  && gsl_complex_abs2(inner) < 1e-6 ) 	ndeg++;
		else if( fabs((values[k]-values[0])/values[0]) < 1e-6 && (gsl_complex_abs2(inner) > 1e-6 && gsl_complex_abs2(inner) < 0.1) ){
			Gram_Schmidt_gsl(gtemp, k+1, block, 1);
			ndeg++;
		}
		printf("%.18lf (%.18lf) inner %.10lf, ", values[k], split, gsl_complex_abs2(inner));
	}
	if( CALL_ARPACK)	printf("AR");
	printf("Final ndeg= %d (after jd %d)\n", ndeg, iter_jd); fflush(stdout);
	if( fabs(values[0])<1e-5 ) 
		printf("Warning :: ARPACK gives something wrong.\n");
	else { 
		ginfo[gindex].ndeg = ndeg;
		for(k=0; k<Nvec; k++)	ginfo[gindex].energy[k] = values[k];
		for(k=0; k<Nvec; k++) for(i=0; i<block; i++)
			ground[k][i] = gtemp[k][i];

        	sprintf(paratemp, "%s/iter/energyall.txt", save_directory );
        	ftemp = fopen(paratemp, "a");
        	char method;
        	if( CALL_ARPACK ) method='A';
        	else              method='D';
        	fprintf(ftemp,"%3d %3d %3c\t", count, gindex, method );
        	for(k=0; k<SDmax; k++) fprintf(ftemp, "%19.16lf\t", values[k]);
        	fprintf(ftemp,"\n");
        	fclose(ftemp);
	}

	/*
	for(k=0; k<SDmax; k++){
		for(j=0; j<SDmax; j++){
			inn_r = 0; inn_i = 0;
#pragma omp parallel for default(none) private(i) shared(j, k, block, gtemp) reduction(+:inn_r, inn_i)
			for(i=0; i<block; i++){
				inn_r += gsl_complex_in_r( gtemp[k][i], gtemp[j][i] );
				inn_i += gsl_complex_in_i( gtemp[k][i], gtemp[j][i] );
			}
			printf("(%d,%d) (%.18lf, %.18lf)\t", k, j, inn_r, inn_i);
		}
		printf("\n");
	}
	*/

	freegscmatrixd(gtemp, Nvec);
	free(istart); free(Hdia); free(Hamil); free(column);
	return Egt;
}//end of find_ground

int continued_fraction(gsl_complex ***Green, int mu, int nu, int direction, double *ai, double *bi, double gm, double p0inner, int dimension){
	int i, k;
	gsl_complex tempg;
	for(i=0; i<Nmax; i++) {
		tempg = zero;
		for(k=dimension-1; k>0; k--){
			tempg = gsl_complex_div(
					gsl_complex_rect(bi[k], 0),
					gsl_complex_sub(
						gsl_complex_rect(-1*ai[k] * direction + gm * direction, wn ),
						tempg
					)
				);
		}
		tempg = gsl_complex_sub(
				gsl_complex_rect(-1*ai[0] * direction + gm * direction, wn ),
				tempg
			);
		Green[mu][nu][i]
			= gsl_complex_add(
					Green[mu][nu][i],
					gsl_complex_mul_real(
						gsl_complex_inverse(tempg),
						p0inner 
					)
			);
	}
	return 0;
}//end of continued_fraction

int build_green(typebasis **basis, sectorinfo ginfo[], int gindex, gsl_complex *ground, PNSpara egv, Glist *table, int degen, int n){
	int p0index, grblock, mapnumber, *column, dimension, *istart, k, direction=table->direction, toindex1, toindex2;
	double *ai, *bi, *Hdia, p0inner, *aibi = (table->aibi) + degen*2*Upper;
	gsl_complex *Hamil, *p0_complex;

	p0index = tosector(basis, ginfo, gindex, table->sp[0], direction);

	if( table->sp[0] > -1 && table->sp[1] > -1 ){
		toindex1 = tosector(basis, ginfo, gindex, table->sp[0], direction);
		toindex2 = tosector(basis, ginfo, gindex, table->sp[1], direction);
		
		if( toindex1 != toindex2 )	p0index = Blocks;
	}
	if( p0index >= Blocks || p0index < 0 ){
		grblock = 0;
		for(k=0; k<2*Upper; k++)	aibi[k] = 0;
		table->dimension = 0;
		return 0;
	}
	else{
		grblock = ginfo[p0index].block;
	}

	istart	= mkvectori(grblock+1);
	mapnumber = ginfo[p0index].mapnumber;
	Hdia	= mkvectord(grblock);
	Hamil	= mkgscvectord(mapnumber);
	column	= mkvectori(mapnumber);
	read_map(n, basis, istart, Hamil, column, ginfo, p0index, egv);
	compute_Hdia(n, Hdia, ginfo, p0index, basis, egv);

	ai = mkvectord(grblock+1);	bi = mkvectord(grblock+1);

	p0_complex = mkgscvectord( grblock );
	dimension = ( grblock < Upper ? grblock : Upper );

	build_p0_gen(basis, ginfo, gindex, ground, p0_complex, &p0inner, table);
	lanczos_vector_complex(istart, Hdia, Hamil, column, &dimension, p0_complex, ai, bi, grblock, dimension, 0);
	for(k=0; k<dimension; k++){
		aibi[k] = ai[k];
		aibi[k+Upper] = bi[k];
	}
	aibi[Upper] = p0inner;
	table->dimension = dimension;

	//if( table->type == 'g' && table->mu[0] == table->mu[1] && direction==-1 )	filling[n] += p0inner;
	free(istart);	free(Hdia);	free(Hamil);	free(column);
	free(ai);	free(bi);
	free(p0_complex);

	return 0;
}//end of build_green

int init_table(Glist *table){
	int i, direc, spin, mu, nu, k;

	i=0;
	for(direc=-1; direc<2; direc+=2){
		for(spin=0; spin<Spmax; spin++){
			for(mu=0; mu<nctot; mu++) for(nu=mu+1; nu<nctot; nu++){
				table[i].direction = direc;
				table[i].mu[0] = mu;
				table[i].mu[1] = nu;
				table[i].sp[0] = spin;
				table[i].sp[1] = spin;
				table[i].type = 'i';
				i++;
			}
			for(mu=0; mu<nctot; mu++) for(nu=mu; nu<nctot; nu++){
				table[i].direction = direc;
				table[i].mu[0] = mu;
				table[i].sp[0] = spin;
				table[i].mu[1] = nu;
				if( mu == nu ){
					table[i].sp[1] = -1;
				}
				else{
					table[i].sp[1] = spin;
				}
				table[i].type = 'g';
				i++;
			}
		}
		for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++){
			table[i].direction = direc;
			table[i].mu[0] = mu;
			table[i].mu[1] = nu;
			table[i].sp[0] = 0;
			table[i].sp[1] = 1;
			table[i].type = 't';
			i++;
		}

		for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++){
			table[i].direction = direc;
			table[i].mu[0] = mu;
			table[i].mu[1] = nu;
			table[i].sp[0] = 0;
			table[i].sp[1] = 1;
			table[i].type = 'k';
			i++;
		}
	}
	for(i=0; i<8*nctot*nctot; i++){
		table[i].aibi = mkvectord(Vmax*2*Upper);
		for(k=0; k<Vmax*2*Upper; k++)	table[i].aibi[k] = 0;
		table[i].co[0] = cone;
		if( table[i].type == 'i' || table[i].type == 't' )
			table[i].co[1] = gsl_complex_rect(0, -table[i].direction);
		else	table[i].co[1] = cone;
	}
	return 0;
}
int print_table(Glist *table){
	int i;
	for(i=0; i<8*nctot*nctot; i++){
		printf("i%d: %d, %d, %d, %c\n", i, table[i].direction, table[i].mu[0], table[i].mu[1], table[i].type);
	}
	return 0;
}
int free_table(Glist *table){
	int i;
	for(i=0; i<8*nctot*nctot; i++)
		free(table[i].aibi);
	return 0;
}

void compute_self(typebasis **basis, int argument[], sectorinfo ginfo[], PNSpara egv, int degeneracy, int count, char *save_directory, gsl_complex ***Self, Glist *table, int n){
	//compute G(w) : on-site Green function, compute G0^Ns(iw) from G(iw) 
	int mu, nu, degen, i, spin;
	int gnd, gndblock, gnode;
	double gm, p0inner, myenergy, myweight, partition=0;
	gsl_complex *ground;

	FILE *fground;
	char paragnd[1024];
	gsl_complex ***kai, ***tildekai, ***Ginverse, ****Green_on, ***Green_total, ***Green_bath;
	gsl_complex ****Giplus, **Grplus, onepi=gsl_complex_rect(1,1), onemi=gsl_complex_rect(1,-1);


	Ginverse	= mkgsctritensord(Nmax, tNC, tNC);
	Green_total	= mkgsctritensord(Nmax, tNC, tNC);
	Green_bath	= mkgsctritensord(Nmax, tNC, tNC);

	Green_on	= mkgsctetratensord(Spmax, nctot, nctot, Nmax);
	Giplus		= mkgsctetratensord(Spmax, nctot, nctot, Nmax);
	tildekai	= mkgsctritensord(nctot, nctot, Nmax);
	kai		= mkgsctritensord(nctot, nctot, Nmax);

	Grplus		= mkgscmatrixd(Spmax, Nmax);

	for(i=0; i<Nmax; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)
		Green_total[i][mu][nu] = zero;

	int myorder;
	doccu[n]=0; filling[n]=0;
	gm = ginfo[argument[0]].energy[0];
	for(degen = 0; degen < degeneracy; degen++){
		gnd	= ginfo[argument[degen]].refer;
		gndblock= ginfo[argument[degen]].block;

		myorder = 0;
		for(i=0; i<degen; i++) if( argument[i] == argument[degen] )	myorder++;
		ground 		= ginfo[argument[degen]].ground[myorder];
		myenergy	= ginfo[argument[degen]].energy[myorder];
		if( beta_real )	myweight	= exp( - ( myenergy - gm )*beta_real );
		else		myweight	= 1;
		partition +=	myweight;

		gnode	= argument[degen]%size;
		if( myrank == gnode ) printf("gnode = %d\n", gnode);
		MPI_Bcast(ground, gndblock, MPI_GSL_COMPLEX, gnode, MPI_COMM_WORLD);

		if( myrank == 0 ){
			sprintf(paragnd, "%s/u%.2lf_%dth_n%d.gnd", save_directory, Uinter, count, n);
			if( degen == 0 )	fground = fopen(paragnd, "wb");
			else			fground = fopen(paragnd, "ab");
			fwrite(ground, sizeof(gsl_complex), gndblock, fground);
			fclose(fground);
			for(mu=0; mu<nctot; mu++) {					//compute elements of on-site Green function.
				for(i = gnd; i<gnd+gndblock; i++){	
					doccu[n] += myweight * gsl_complex_abs2(ground[i-gnd])* One(basis[i][0], mu) * One(basis[i][1], mu);
					filling[n] += myweight * gsl_complex_abs2(ground[i-gnd])* ((int)One(basis[i][0], mu) + (int)One(basis[i][1], mu));
				}
			}
		}

		for(spin=0; spin<Spmax; spin++) {
			for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nmax; i++){
				Green_on[spin][mu][nu][i] = zero;
				Giplus[spin][mu][nu][i] = zero;
			}
		}
		for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nmax; i++){
			kai[mu][nu][i] = zero;
			tildekai[mu][nu][i] = zero;
		}

		for(i=myrank; i<8*nctot*nctot; i+=size)
			build_green( basis, ginfo, argument[degen], ground, egv, &table[i], degen, n);

		for(i=0; i<8*nctot*nctot; i++){
			MPI_Bcast( table[i].aibi+degen*2*Upper, 2*Upper, MPI_DOUBLE, i%size, MPI_COMM_WORLD);
			MPI_Bcast( &(table[i].dimension), 1, MPI_INT, i%size, MPI_COMM_WORLD);
			p0inner = table[i].aibi[Upper+degen*2*Upper];
			if( table[i].type == 'i' && fabs(p0inner) > 1e-10 )
				continued_fraction(Giplus[table[i].sp[0]], table[i].mu[0], table[i].mu[1], table[i].direction, table[i].aibi+degen*2*Upper, &table[i].aibi[Upper+degen*2*Upper], myenergy, p0inner, Upper);
			if( table[i].type == 'g' && fabs(p0inner) > 1e-10 )
				continued_fraction(Green_on[table[i].sp[0]], table[i].mu[0], table[i].mu[1], table[i].direction, table[i].aibi+degen*2*Upper, &table[i].aibi[Upper+degen*2*Upper], myenergy, p0inner, Upper);
			if( Nonzero_SOC ){
				if( table[i].type == 't' && fabs(p0inner) > 1e-10 )
					continued_fraction(tildekai, table[i].mu[0], table[i].mu[1], table[i].direction, table[i].aibi+degen*2*Upper, &table[i].aibi[Upper+degen*2*Upper], myenergy, p0inner, Upper);
				if( table[i].type == 'k' && fabs(p0inner) > 1e-10 )
					continued_fraction(kai, table[i].mu[0], table[i].mu[1], table[i].direction, table[i].aibi+degen*2*Upper, &table[i].aibi[Upper+degen*2*Upper], myenergy, p0inner, Upper);
			}
		}

		for(spin=0; spin<Spmax; spin++) {
			for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nmax; i++){
				Giplus[spin][mu][nu][i] = gsl_complex_mul( Giplus[spin][mu][nu][i], gsl_complex_rect(0, -1 ) );
			}
			for(mu=0; mu<nctot; mu++) for(nu=mu+1; nu<nctot; nu++) for(i=0; i<Nmax; i++){
				Grplus[spin][i] = Green_on[spin][mu][nu][i];
				Green_on[spin][mu][nu][i]
					= gsl_complex_sub(
							Green_on[spin][mu][nu][i],
							gsl_complex_add(
								gsl_complex_mul( Green_on[spin][mu][mu][i], onemi ),
								gsl_complex_mul( Green_on[spin][nu][nu][i], onemi )
								)
							);
				Green_on[spin][mu][nu][i] = gsl_complex_add( Green_on[spin][mu][nu][i], Giplus[spin][mu][nu][i] );
				Green_on[spin][mu][nu][i] = gsl_complex_div_real( Green_on[spin][mu][nu][i], 2 );

				Green_on[spin][nu][mu][i]
					= gsl_complex_sub(
							Grplus[spin][i],
							gsl_complex_add(
								gsl_complex_add( Green_on[spin][mu][mu][i], Green_on[spin][nu][nu][i] ),
								Green_on[spin][mu][nu][i]
								)
							);
			}
			for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nmax; i++){
				Green_total[i][2*mu+spin][2*nu+spin] = gsl_complex_add( Green_total[i][2*mu+spin][2*nu+spin], gsl_complex_mul_real(Green_on[spin][mu][nu][i], myweight) );
			}
		}
		if( Nonzero_SOC ) for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++) for(i=0; i<Nmax; i++){
			Green_total[i][2*mu][2*nu+1]
				= gsl_complex_add(
					Green_total[i][2*mu][2*nu+1],
					gsl_complex_mul_real(
						gsl_complex_sub(
							gsl_complex_add(
								kai[mu][nu][i],
								gsl_complex_mul(
									tildekai[mu][nu][i],
									gsl_complex_rect(0,-1)
								)
							),
							gsl_complex_mul(
								gsl_complex_add(
									Green_on[0][mu][mu][i],
									Green_on[1][nu][nu][i]
								),
								onemi
							)
						),
						myweight/2.
					)
				);
			Green_total[i][2*nu+1][2*mu]
				= gsl_complex_add(
					Green_total[i][2*nu+1][2*mu],
					gsl_complex_mul_real(
						gsl_complex_sub(
							gsl_complex_add(
								kai[mu][nu][i],
								gsl_complex_mul(
									tildekai[mu][nu][i],
									gsl_complex_rect(0,1)
								)
							),
							gsl_complex_mul(
								gsl_complex_add(
									Green_on[0][mu][mu][i],
									Green_on[1][nu][nu][i]
								),
								onepi
							)
						),
						myweight/2.
					)
				);
		}
	}
	for(i=0; i<Nmax; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) {
		Green_total[i][mu][nu] = gsl_complex_div_real(Green_total[i][mu][nu], partition);	//up to here, Green_total = on-site Green function(not inverse).
	}
	inverse_complex_matrix_all_omp( Green_total, Ginverse, 0, Nmax, tNC);

	doccu[n] /= nctot*partition;
	double red_filling;
	MPI_Reduce(&filling[n], &red_filling, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	filling[n] = red_filling/partition;

	build_bath(n, Green_bath, egv);

	for(i=0; i<Nmax; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) {
		Self[i][mu][nu] = gsl_complex_sub( Green_bath[i][mu][nu], Ginverse[i][mu][nu] );		//Self[spin][mu][nu] = self energy of cluster(Sigma_munu^c)
	}
	if( PARAMAGNETIC ){
		int pairtnc[tNC] =  { 1, 0, 3, 2, 5, 4 }; 
		int angletnc[tNC]=  {-1, 1,-1, 1,-1, 1 }; 
		if( Nonzero_SOC ){
			int zpairtnc[tNC] =  { 3, 2, 1, 0, 5, 4 }; 
			int zangletnc[tNC]=  { 1,-1, 1,-1, 1,-1 };       // MTrev_{mu,pairtnc[mu}
			for( int aa=0; aa<tNC; aa++ ){
				pairtnc[aa]	= zpairtnc[aa];
				angletnc[aa]	= zangletnc[aa];
			}
		}
#pragma omp parallel for default(shared) private(mu,nu) shared(i,Self)
		for(i=0; i<Nmax; i++) {
			gsl_complex **selftrev;
			selftrev = mkgscmatrixd( tNC, tNC );
			for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)  selftrev[mu][nu] = zero;
			for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)  {
				int ggamma = pairtnc[nu];
				int ddelta = pairtnc[mu];
				int ThetaNuGammaT=-angletnc[nu];        // Theta^transpose = -Theta
				int ThetaMuDelta = angletnc[mu];
				GSL_REAL(selftrev[mu][nu]) = -GSL_REAL(Self[i][ggamma][ddelta]) * ThetaNuGammaT * ThetaMuDelta;
				GSL_IMAG(selftrev[mu][nu]) = -GSL_IMAG(Self[i][ggamma][ddelta]) * ThetaNuGammaT * ThetaMuDelta;
			}    
			for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)  {
				GSL_REAL(Self[i][mu][nu]) = ( GSL_REAL(Self[i][mu][nu]) + GSL_REAL(selftrev[mu][nu]) )/2.;
				GSL_IMAG(Self[i][mu][nu]) = ( GSL_IMAG(Self[i][mu][nu]) + GSL_IMAG(selftrev[mu][nu]) )/2.;
			}    
			freegscmatrixd( selftrev, tNC );
		}   
	}
	///////////////////////
	FILE *ftemp;
	char order[1024];

	if( myrank == 0 ){
		save_aibi(table, save_directory, count, degeneracy, n);

		sprintf(order, "%s/ongoing/self%d_u%.2lf_%dth.txt", save_directory, n, Uinter, count);
		ftemp = fopen(order, "w");
		print_self(ftemp, Self, 0, Nmax);
		fclose(ftemp);
		sprintf(order, "%s/ongoing/ginverse%d_u%.2lf_%dth.txt", save_directory, n, Uinter, count);
		ftemp = fopen(order, "w");
		print_self(ftemp, Ginverse, 0, Nmax);
		fclose(ftemp);

		sprintf(order, "%s/ongoing/bathinverse%d_u%.2lf_%dth.txt", save_directory, n, Uinter, count);
		ftemp = fopen(order, "w");
		print_self(ftemp, Green_bath, 0, Nmax);
		fclose(ftemp);

		inverse_complex_matrix_all_omp( Green_bath, Ginverse, 0, Nmax, tNC);
		sprintf(order, "%s/ongoing/green%d_u%.2lf_%dth.txt", save_directory, n, Uinter, count);
		ftemp = fopen(order, "w");
		print_self(ftemp, Green_total, 0, Nmax);
		fclose(ftemp);

		sprintf(order, "%s/ongoing/bath%d_u%.2lf_%dth.txt", save_directory, n, Uinter, count);
		ftemp = fopen(order, "w");
		print_self(ftemp, Ginverse, 0, Nmax);
		fclose(ftemp);

	}
	MPI_Barrier( MPI_COMM_WORLD);
	///////////////////////

	freegsctritensord( Ginverse,		Nmax, tNC);
	freegsctritensord( Green_total,		Nmax, tNC);
	freegsctritensord( Green_bath,		Nmax, tNC);
	freegsctritensord( tildekai,		nctot, nctot);
	freegsctritensord( kai,			nctot, nctot);

	freegsctetratensord( Green_on, Spmax, nctot, nctot);
	freegsctetratensord( Giplus, Spmax, nctot, nctot);

	freegscmatrixd( Grplus, Spmax);
}

int compute_gnew(gsl_complex *hopmatrix, gsl_complex ****Self, int count, char *save_directory, Quad *quadrature){
	int mu, nu, i, kx, ky, kz, same, n, Omp_num = omp_get_max_threads();
	int mystart, myend, myblock;
	gsl_complex ***Glocalinverse, ***Glocal, *oneline, ***myG;
	
	mystart	=  myrank	* ((Nmax)/size);
	myend	= (myrank+1)	* ((Nmax)/size);
	if( myrank+1 == size )	myend = Nmax;
	myblock = myend - mystart;

	//for(n=0; n<Ni; n++) for(i=0; i<Nmax; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	Self[n][i][mu][nu] = zero;
	Glocalinverse	= mkgsctritensord(Nmax, NU, NU);
	myG		= mkgsctritensord(Omp_num, NU, NU);
	Glocal		= mkgsctritensord(Nmax, NU, NU);
	oneline		= mkgscvectord(Nmax*NU*NU);

	double diag[NU];
	update_diag(diag);

	for(i=mystart; i<myend; i++){
		for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)	Glocal[i][mu][nu] = zero;

#pragma omp parallel private(mu,nu,kx,ky,kz,n,same) shared(i,oneline,Matsu,hopmatrix,Ecluster,Glocal,quadrature,zero,diag)
		{
			int mythread = omp_get_thread_num();
			gsl_complex **Glocalt;
			Glocalt = mkgscmatrixd(NU, NU);
			for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)	myG[mythread][mu][nu] = zero;
		#pragma omp for
			for(kx=0; kx<Nintx; kx++) for(ky=0; ky<Ninty; ky++) for(kz=0; kz<Nintz; kz++){
				for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++) {
					same = (int)(mu==nu);
					Glocalt[mu][nu]
						= gsl_complex_rect(
								//- GSL_REAL( hopmatrix(kx,ky,kz,mu,nu) ) - GSL_REAL(Ecluster[mu][nu])*same ,
								- GSL_REAL( hopmatrix(kx,ky,kz,mu,nu) ) - diag[mu]*same ,
								- GSL_IMAG( hopmatrix(kx,ky,kz,mu,nu) ) + wn*same  
								);
				}
				for(n=0; n<Ni; n++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) 
					Glocalt[tNC*n+mu][tNC*n+nu]
						= gsl_complex_rect(
								 GSL_REAL( Glocalt[tNC*n+mu][tNC*n+nu] ) - GSL_REAL(Self[n][i][mu][nu]) ,
								 GSL_IMAG( Glocalt[tNC*n+mu][tNC*n+nu] ) - GSL_IMAG(Self[n][i][mu][nu]) 
								);
				inverse_complex_matrix_lapack( Glocalt, Glocalt, NU);

				for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++){
					myG[mythread][mu][nu]
						= gsl_complex_rect(
								GSL_REAL(myG[mythread][mu][nu]) + GSL_REAL(Glocalt[mu][nu]) *(quadrature->weight[0])[kx]*(quadrature->weight[1])[ky]*(quadrature->weight[2])[kz]/Area,
								GSL_IMAG(myG[mythread][mu][nu]) + GSL_IMAG(Glocalt[mu][nu]) *(quadrature->weight[0])[kx]*(quadrature->weight[1])[ky]*(quadrature->weight[2])[kz]/Area
								);
				}
			}
			freegscmatrixd(Glocalt, NU);
		}
		for(kx=0; kx<Omp_num; kx++) for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)
			Glocal[i][mu][nu]
				= gsl_complex_rect(
					GSL_REAL(Glocal[i][mu][nu]) + GSL_REAL(myG[kx][mu][nu]),
					GSL_IMAG(Glocal[i][mu][nu]) + GSL_IMAG(myG[kx][mu][nu])
				);
		for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)
			oneline[i*NU*NU+mu*NU+nu] = Glocal[i][mu][nu];
	}
	int *offset, *number;
	offset = mkvectori( size+1 );
	number = mkvectori( size );
	for(i=0; i<size; i++)
		offset[i] = NU*NU*i * (Nmax/size);
	offset[i] = NU*NU*Nmax;
	for(i=0; i<size; i++)
		number[i] = offset[i+1] - offset[i];

	MPI_Allgatherv(MPI_IN_PLACE, myblock*NU*NU, MPI_GSL_COMPLEX, oneline, number, offset, MPI_GSL_COMPLEX, MPI_COMM_WORLD);
#pragma omp parallel for default(none) private(i,mu,nu) shared(Nmax,Glocal,oneline)
	for(i=0; i<Nmax; i++) for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)	Glocal[i][mu][nu] = oneline[i*NU*NU+mu*NU+nu];

	gsl_complex ***sub;
	sub	= mkgsctritensord(Nmax, tNC, tNC);

	for(n=0; n<Ni; n++){
		for(i=0; i<Nmax; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	sub[i][mu][nu] = Glocal[i][tNC*n+mu][tNC*n+nu];
		inverse_complex_matrix_all_omp( sub, sub, 0, Nmax, tNC);
		
		for(i=0; i<Nmax; i++) for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	Greennew[n][i][mu][nu] = gsl_complex_add( Self[n][i][mu][nu], sub[i][mu][nu] );
	}
	freegsctritensord(sub, Nmax, tNC);
	free(offset);
	free(number);

	if( myrank == 0 ){
		char order[1024];
		FILE *ftemp;
		for(n=0; n<Ni; n++){
			sprintf(order, "%s/ongoing/greennew%d_u%.2lf_%dth.txt", save_directory, n, Uinter, count);
			ftemp = fopen(order, "w");
			print_self(ftemp, Greennew[n], 0, Nmax);
			fclose(ftemp);
		}

		for(i=0; i<Nmax; i++){
			inverse_complex_matrix_lapack( Glocal[i], Glocalinverse[i], NU );
		}
		sprintf(order, "%s/ongoing/Glocalinverse_u%.2lf_%dth.txt", save_directory, Uinter, count);
		ftemp = fopen(order, "w");
		print_lattice(ftemp, Glocalinverse, 0, Nmax);
		fclose(ftemp);

		sprintf(order, "%s/ongoing/Glocal_u%.2lf_%dth.txt", save_directory, Uinter, count);
		ftemp = fopen(order, "w");
		print_lattice(ftemp, Glocal, 0, Nmax);
		fclose(ftemp);
	}

	freegsctritensord( Glocalinverse,	Nmax, NU);
	freegsctritensord( myG,			Omp_num, NU);
	freegsctritensord( Glocal,		Nmax, NU);

	return 0;
}//end of compute_self

int build_bath(int n, gsl_complex ***Green_bath, PNSpara egv){
	gsl_complex tempGreen, inner, deno;
	int mu, nu, k, i;
#pragma omp parallel for default(none)	\
	private(tempGreen,mu,nu,k,inner,deno) shared(Green_bath,egv,Ecluster,zero,Matsu,Nmax,tnb, n)
	for(i=0; i<Nmax; i++) {
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++) {
			tempGreen = zero;
			for(k=0; k<tnb; k++){
				//tempGreen = gsl_complex_sub( tempGreen, gsl_complex_div( gsl_complex_mul( gsl_complex_conjugate( egv.hybrid[nu][k] ), egv.hybrid[mu][k] ), gsl_complex_rect( - egv.egbath[k], wn )));
					//(-) sign is included.
				inner = gsl_complex_rect(
						gsl_complex_in_r(egv.hybrid[nu][k], egv.hybrid[mu][k]),
						gsl_complex_in_i(egv.hybrid[nu][k], egv.hybrid[mu][k])
					);
				deno = gsl_complex_rect(
						egv.egbath[k]/(egv.egbath[k]*egv.egbath[k]+wn*wn),
						wn/(egv.egbath[k]*egv.egbath[k]+wn*wn)
					);
				tempGreen = gsl_complex_rect(
						GSL_REAL(tempGreen) + gsl_complex_mul_r(&inner, &deno),
						GSL_IMAG(tempGreen) + gsl_complex_mul_i(&inner, &deno)
					);
			}
			//Green_bath[i][mu][nu] = gsl_complex_add( gsl_complex_sub( gsl_complex_rect( 0, wn*(mu==nu) ), Ecluster[mu][nu]), tempGreen);
			Green_bath[i][mu][nu]
				= gsl_complex_rect(
					GSL_REAL(tempGreen) - GSL_REAL(Ecluster[tNC*n+mu][tNC*n+nu]),
					GSL_IMAG(tempGreen) - GSL_IMAG(Ecluster[tNC*n+mu][tNC*n+nu]) + wn*(mu==nu)
				);
		}
	}
	return 0;
}//end of build_bath

int save_aibi(Glist *table, char *save_directory, int count, int degeneracy, int n){
	FILE *faibi;
	char paradata[1024];
	int i, k, degen=0;

	if(nctot>1){
		sprintf(paradata, "%s/datafile/giplus_u%.2lf_%dth.ab", save_directory, Uinter, count);
		faibi = fopen(paradata, "a");
		for(degen=0; degen<degeneracy; degen++) for(i=0; i<8*nctot*nctot; i++) if( table[i].type == 'i' ){
			fprintf(faibi, "%d\n%d\t%d\t%d\t%d\t%d\t1\t\t", n, degen, table[i].sp[0], table[i].mu[0], table[i].mu[1], table[i].direction);
			fprintf(faibi, "%d\t%.18lf\n", table[i].dimension, table[i].aibi[degen*2*Upper+Upper]);
			for(k=0; k<table[i].dimension; k++)	fprintf(faibi, "%.18lf\t", table[i].aibi[degen*2*Upper+k]); fprintf(faibi, "\n");
			for(k=0; k<table[i].dimension; k++)	fprintf(faibi, "%.18lf\t", table[i].aibi[degen*2*Upper+Upper+k]); fprintf(faibi, "\n");
		}
		fclose(faibi);
	}

	sprintf(paradata, "%s/datafile/greenon_u%.2lf_%dth.ab", save_directory, Uinter, count);
	faibi = fopen(paradata, "a");
	for(degen=0; degen<degeneracy; degen++) for(i=0; i<8*nctot*nctot; i++) if( table[i].type == 'g' ){
		fprintf(faibi, "%d\n%d\t%d\t%d\t%d\t%d\t1\t\t", n, degen, table[i].sp[0], table[i].mu[0], table[i].mu[1], table[i].direction);
		fprintf(faibi, "%d\t%.18lf\n", table[i].dimension, table[i].aibi[degen*2*Upper+Upper]);
		for(k=0; k<table[i].dimension; k++)	fprintf(faibi, "%.18lf\t", table[i].aibi[degen*2*Upper+k]); fprintf(faibi, "\n");
		for(k=0; k<table[i].dimension; k++)	fprintf(faibi, "%.18lf\t", table[i].aibi[degen*2*Upper+Upper+k]); fprintf(faibi, "\n");
	}
	fclose(faibi);

	if( Nonzero_SOC ){
		sprintf(paradata, "%s/datafile/kai_u%.2lf_%dth.ab", save_directory, Uinter, count);
		faibi = fopen(paradata, "a");
		for(degen=0; degen<degeneracy; degen++) for(i=0; i<8*nctot*nctot; i++) if( table[i].type == 'k' ){
			fprintf(faibi, "%d\n%d\t%d\t%d\t%d\t1\t\t", n, degen, table[i].mu[0], table[i].mu[1], table[i].direction);
			fprintf(faibi, "%d\t%.18lf\n", table[i].dimension, table[i].aibi[degen*2*Upper+Upper]);
			for(k=0; k<table[i].dimension; k++)	fprintf(faibi, "%.18lf\t", table[i].aibi[degen*2*Upper+k]); fprintf(faibi, "\n");
			for(k=0; k<table[i].dimension; k++)	fprintf(faibi, "%.18lf\t", table[i].aibi[degen*2*Upper+Upper+k]); fprintf(faibi, "\n");
		}
		fclose(faibi);

		sprintf(paradata, "%s/datafile/tildekai_u%.2lf_%dth.ab", save_directory, Uinter, count);
		faibi = fopen(paradata, "a");
		for(degen=0; degen<degeneracy; degen++) for(i=0; i<8*nctot*nctot; i++) if( table[i].type == 't' ){
			fprintf(faibi, "%d\n%d\t%d\t%d\t%d\t1\t\t", n, degen, table[i].mu[0], table[i].mu[1], table[i].direction);
			fprintf(faibi, "%d\t%.18lf\n", table[i].dimension, table[i].aibi[degen*2*Upper+Upper]);
			for(k=0; k<table[i].dimension; k++)	fprintf(faibi, "%.18lf\t", table[i].aibi[degen*2*Upper+k]); fprintf(faibi, "\n");
			for(k=0; k<table[i].dimension; k++)	fprintf(faibi, "%.18lf\t", table[i].aibi[degen*2*Upper+Upper+k]); fprintf(faibi, "\n");
		}
		fclose(faibi);
	}

	return 0;
}
void save_parameters(sectorinfo **ginfo, int **argument, PNSpara *egv, int count, int *degeneracy, int tol, FILE *fsave){
	int n, k, mu, myorder;

	fprintf(fsave, "%d\t%lf\t%lf\t%d\n", count, Uinter, MU, tol);
	fprintf(fsave, "%d\t%d\n", beta, Nmax);
	for(n=0; n<Ni; n++){
		fprintf(fsave, "%d\n", degeneracy[n]);
		for(k=0; k<degeneracy[n]; k++){
			myorder = 0;
			for(int i=0; i<k; i++) if( argument[n][i] == argument[n][k] )	myorder++;
			fprintf(fsave, "%.17lf\t%d\t%d\n", ginfo[n][argument[n][k]].energy[myorder], ginfo[n][argument[n][k]].refer, ginfo[n][argument[n][k]].block);
		}

		for(k=0; k<tnb; k++)	fprintf(fsave, "%.17lf\t\t\t\t", egv[n].egbath[k]);
		fprintf(fsave, "\n\n");
		for(mu=0; mu<tNC; mu++){
			for(k=0; k<tnb; k++)	fprintf(fsave, "%20.17lf  %20.17lf\t", GSL_REAL(egv[n].hybrid[mu][k]), GSL_IMAG(egv[n].hybrid[mu][k]));
			fprintf(fsave, "\n");
		}
		fprintf(fsave, "\n\n");
	}
}// end of save_parameters

void minimize_zeroSOC(int count, PNSpara egv, double *minimal, int *iter, double Converge, int n, int sp){
	double p[Np_zeroSOC];

	/*
	gsl_complex **transform = mkgscmatrixd(NU, NU);
	if( ROTATE ){
		init_transform(transform);
		symmetrize_egv(&egv[n], transform);
	}
	freegscmatrixd(transform, NU);
	*/
	convert_egv_to_p_zeroSOC(egv, p, sp);

	int i, status, **par;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	gsl_vector *x;
	gsl_multimin_function_fdf my_func;

	int *offset, *number, bound[2], resi=Np_zeroSOC%size;
	offset = mkvectori( size+1 );
	number = mkvectori( size );
	if( size > Np_zeroSOC ){
		for(i=0; i<Np_zeroSOC; i++)
			offset[i] = i;
		for(i=Np_zeroSOC; i<size+1; i++)
			offset[i] = Np_zeroSOC;
		bound[0] = myrank < Np_zeroSOC ? myrank : Np_zeroSOC;
		bound[1] = myrank < Np_zeroSOC ? myrank+1 : Np_zeroSOC;
		//printf("at rank%02d myblock = %d-%d=%d\n", myrank, bound[1], bound[0], bound[1]-bound[0]);
	}
	else{
		for(i=0; i<size; i++)
			offset[i] = i * (Np_zeroSOC/size);
		for(i=0; i<resi+1; i++)
			offset[i] += i;
		for(i=resi+1; i<size; i++)
			offset[i] += resi;
		offset[i] = Np_zeroSOC;
		bound[0] = offset[myrank];
		if( myrank+1 == size ){
			bound[1] = Np_zeroSOC;
		}
		else	bound[1] = offset[myrank+1];
	}
	for(i=0; i<size; i++)
		number[i] = offset[i+1] - offset[i];
	//for(i=0; i<size; i++)	printf("rank%d: number[%d] = %d\n", myrank, i, number[i]);

	par = (int **) malloc( 5*sizeof(int*) );
	par[0] = offset;
	par[1] = number;
	par[2] = bound;
	par[3] = &n;

	par[4] = &sp;

	my_func.n = Np_zeroSOC;
	my_func.f = my_f_zeroSOC;
	my_func.df = my_df_zeroSOC;
	my_func.fdf = my_fdf_zeroSOC;
	my_func.params = par;

	x = gsl_vector_alloc(Np_zeroSOC);
	for(i=0; i<Np_zeroSOC; i++)	gsl_vector_set(x, i, p[i]);

	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc(T, Np_zeroSOC);

	int ex = ( count > 7 ) ? 10 : count+2; 
	double step = pow(10,-ex);
	gsl_multimin_fdfminimizer_set(s, &my_func, x, step, 2e-4);
	*iter = 0;
	int find = 1, allfind;

	do{
		*iter = *iter + 1;
		status = gsl_multimin_fdfminimizer_iterate (s);
		status = gsl_multimin_test_gradient (s->gradient, Converge);
		if( status == GSL_CONTINUE )	find = 0;

		for(i=0; i<Np_zeroSOC; i++)	p[i] = gsl_vector_get(s->x, i);
		MPI_Bcast(p, Np_zeroSOC, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for(i=0; i<Np_zeroSOC; i++)	gsl_vector_set(s->x, i, p[i]);
		MPI_Allreduce(&find, &allfind, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if( allfind == size )	break;
		else if( allfind > 0 ){
			printf("iter%d rank%d find%d allfind%d\n", *iter, myrank, find, allfind);
			printf("Inconsistency in minimization!\n");	fflush(stdout);
			exit(1);
		}
		else find = 1;
	}while(status == GSL_CONTINUE && *iter < Maxmin);

	for(i=0; i<Np_zeroSOC; i++)	p[i] = gsl_vector_get(s->x, i);
	*minimal = s->f;

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);

	convert_p_to_egv_zeroSOC(egv, p, sp);
	sort_bath_zeroSOC(&egv, sp);
	char paratemp[1024];
	if( 1 || myrank == 0 ){
		sprintf(paratemp, "rank%d minimize_zeroSOC %d sp %d (iter %d)", myrank, n, sp, *iter);
		print_egv(paratemp, egv);
	}

	free(par);
	free(offset);
	free(number);

}//end of minimize_zeroSOC



void minimize(int count, PNSpara egv, double *minimal, int *iter, double Converge, int n){
	double p[Np];

	/*
	gsl_complex **transform = mkgscmatrixd(NU, NU);
	if( ROTATE ){
		init_transform(transform);
		symmetrize_egv(&egv[n], transform);
	}
	freegscmatrixd(transform, NU);
	*/
	convert_egv_to_p(egv, p);

	int i, status, **par;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	gsl_vector *x;
	gsl_multimin_function_fdf my_func;

	int *offset, *number, bound[2], resi=Np%size;
	offset = mkvectori( size+1 );
	number = mkvectori( size );
	if( size > Np ){
		for(i=0; i<Np; i++)
			offset[i] = i;
		for(i=Np; i<size+1; i++)
			offset[i] = Np;
		bound[0] = myrank < Np ? myrank : Np;
		bound[1] = myrank < Np ? myrank+1 : Np;
		//printf("at rank%02d myblock = %d-%d=%d\n", myrank, bound[1], bound[0], bound[1]-bound[0]);
	}
	else{
		for(i=0; i<size; i++)
			offset[i] = i * (Np/size);
		for(i=0; i<resi+1; i++)
			offset[i] += i;
		for(i=resi+1; i<size; i++)
			offset[i] += resi;
		offset[i] = Np;
		bound[0] = offset[myrank];
		if( myrank+1 == size ){
			bound[1] = Np;
		}
		else	bound[1] = offset[myrank+1];
	}
	for(i=0; i<size; i++)
		number[i] = offset[i+1] - offset[i];
	//for(i=0; i<size; i++)	printf("rank%d: number[%d] = %d\n", myrank, i, number[i]);

	par = (int **) malloc( 4*sizeof(int*) );
	par[0] = offset;
	par[1] = number;
	par[2] = bound;
	par[3] = &n;

	my_func.n = Np;
	my_func.f = my_f;
	my_func.df = my_df;
	my_func.fdf = my_fdf;
	my_func.params = par;

	x = gsl_vector_alloc(Np);
	for(i=0; i<Np; i++)	gsl_vector_set(x, i, p[i]);

	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc(T, Np);

	int ex = ( count > 7 ) ? 10 : count+2; 
	double step = pow(10,-ex);
	gsl_multimin_fdfminimizer_set(s, &my_func, x, step, 2e-4);
	*iter = 0;
	int find = 1, allfind;

	do{
		*iter = *iter + 1;
		status = gsl_multimin_fdfminimizer_iterate (s);
		status = gsl_multimin_test_gradient (s->gradient, Converge);
		if( status == GSL_CONTINUE )	find = 0;

		for(i=0; i<Np; i++)	p[i] = gsl_vector_get(s->x, i);
		MPI_Bcast(p, Np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for(i=0; i<Np; i++)	gsl_vector_set(s->x, i, p[i]);
		MPI_Allreduce(&find, &allfind, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if( allfind == size )	break;
		else if( allfind > 0 ){
			printf("iter%d rank%d find%d allfind%d\n", *iter, myrank, find, allfind);
			printf("Inconsistency in minimization!\n");	fflush(stdout);
			exit(1);
		}
		else find = 1;
	}while(status == GSL_CONTINUE && *iter < Maxmin);

	for(i=0; i<Np; i++)	p[i] = gsl_vector_get(s->x, i);
	*minimal = s->f;

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);

	convert_p_to_egv(egv, p);
	sort_bath(&egv);
	char paratemp[1024];
	if( 1 || myrank == 0 ){
		sprintf(paratemp, "rank%d minimize %d (iter %d)", myrank, n, *iter);
		print_egv(paratemp, egv);
	}

	free(par);
	free(offset);
	free(number);

}//end of minimize

void print_complex(gsl_complex cmnumber){
	printf("%lf + %lf i\n", GSL_REAL(cmnumber), GSL_IMAG(cmnumber) );
}

int corcd(typebasis **basis, sectorinfo ginfo[], int gindex, int direction, int mu, int spin, gsl_complex *p0, gsl_complex coeff, gsl_complex *ground){
	int nbasis, phase, opern, occupied, toindex, torefer, fromrefer, fromblock, index, *table;
	typebasis grstate[2];

	fromrefer = ginfo[gindex].refer;
	fromblock = ginfo[gindex].block;

	toindex = tosector(basis, ginfo, gindex, spin, direction);

	torefer = ginfo[toindex].refer;
	
	table = mkvectori( Powns2 );
	for(nbasis=0; nbasis<Powns2; nbasis++){
		index = (basis[nbasis][0]<<Ns) + basis[nbasis][1];
		table[index] = nbasis;
	}

	occupied = (1-direction)/2;

	if( occupied == 0 ){
		for(nbasis=fromrefer; nbasis<fromrefer+fromblock; nbasis++){				//compute elements of |p0>
			if( Zero(basis[nbasis][spin], mu ) ){
				grstate[0] = basis[nbasis][0];		grstate[1] = basis[nbasis][1];
				phase = permu(grstate, mu, spin);	grstate[spin] = Turnon(grstate[spin], mu);
				index = (grstate[0]<<Ns) + grstate[1];
				opern = (int)table[index] - torefer;
				p0[opern] = gsl_complex_add(
						p0[opern],
						gsl_complex_mul(
							ground[nbasis-fromrefer],
							gsl_complex_mul_real(coeff, phase )
						)
					);

			}
		}
	}
	else if( occupied == 1 ){
		for(nbasis=fromrefer; nbasis<fromrefer+fromblock; nbasis++){				//compute elements of |p0>
			if( One(basis[nbasis][spin], mu ) ){
				grstate[0] = basis[nbasis][0];		grstate[1] = basis[nbasis][1];
				phase = permu(grstate, mu, spin);	grstate[spin] = Turnoff(grstate[spin], mu);
				index = (grstate[0]<<Ns) + grstate[1];
				opern = (int)table[index] - torefer;
				p0[opern] = gsl_complex_add(
						p0[opern],
						gsl_complex_mul(
							ground[nbasis-fromrefer],
							gsl_complex_mul_real(coeff, phase )
						)
					);
			}
		}
	}
	else{
		printf("Error in corcd!\n");
		exit(1);
	}
	free(table);
	return 0;
}//end of corcd

int build_p0_gen(typebasis **basis, sectorinfo ginfo[], int gindex, gsl_complex *ground, gsl_complex *p0, double *p0inner, Glist *table){
	int imax=2, i, k, toindex, grblock, direction = table->direction, toindex1, toindex2;


	toindex = tosector(basis, ginfo, gindex, table->sp[0], direction);
	grblock = ginfo[toindex].block;

	initgscvectord(p0, grblock);

	if( table->sp[0] > -1 && table->sp[1] > -1 ){
		toindex1 = tosector(basis, ginfo, gindex, table->sp[0], direction);
		toindex2 = tosector(basis, ginfo, gindex, table->sp[1], direction);
		
		if( toindex1 != toindex2 )	imax = 0;
	}

	for(i=0; i<imax; i++){
		if( table->sp[i] < 0 )	continue;
		corcd(basis, ginfo, gindex, direction, table->mu[i], table->sp[i], p0, table->co[i], ground);
	}

	*p0inner = 0.;
	for(k=0; k<grblock; k++)	*p0inner += gsl_complex_abs2( p0[k] );
	if( *p0inner > 1e-10 ) for(k=0; k<grblock; k++)	p0[k] = gsl_complex_div_real( p0[k], sqrt(*p0inner) );

	if(direction == -1)	printf("%c: mu = %d, nu = %d\tp0inner\t= %.10lf\n", table->type, table->mu[0], table->mu[1], *p0inner);
	return 0;
}//end of build_p0_gen

int tosector(typebasis **basis, sectorinfo ginfo[], int gindex, int spin, int direction){
	int fromrefer,ptl,netspin, ptl_found, netspin_found;
	fromrefer = ginfo[gindex].refer;
	ptl	= count_bit(Ns, basis[fromrefer][1]) + count_bit(Ns, basis[fromrefer][0]);
	netspin	= count_bit(Ns, basis[fromrefer][1]) - count_bit(Ns, basis[fromrefer][0]);

	int theindex = find_index(Ns, ptl+direction, netspin+direction*(2*spin-1), Blocks);
	if( Nonzero_SOC ) return ptl+direction;

	if( ptl+direction > 2*Ns || abs(netspin+direction*(2*spin-1)) > Ns - abs(ptl+direction-Ns) )	return Blocks+1;
	ptl_found	= count_bit(Ns, basis[ginfo[theindex].refer][1]) + count_bit(Ns, basis[ginfo[theindex].refer][0]);
	netspin_found	= count_bit(Ns, basis[ginfo[theindex].refer][1]) - count_bit(Ns, basis[ginfo[theindex].refer][0]);
	if( ptl+direction != ptl_found || netspin+direction*(2*spin-1) != netspin_found ){
		printf("ptl+direction %d != ptl_found %d, or netspin+direction*(2*spin-1) %d != netspin_found %d\n", ptl+direction, ptl_found, netspin+direction*(2*spin-1), netspin_found);
		exit(1);
	}

	return	theindex;
}

int print_green(FILE *file, gsl_complex ***Green){
	int mu, nu, i;
	for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++)	fprintf(file, "%d,%d\t\t\t\t\t", mu, nu );
	fprintf(file, "\n");
	for(i=0; i<Nmax; i++){
		for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++)	fprintf(file, "%.10lf\t%.10lf\t\t", GSL_REAL( Green[mu][nu][i] ), GSL_IMAG( Green[mu][nu][i] ) );
		fprintf(file, "\n");
	}
	return 0;
}

int print_gnew(FILE *file, gsl_complex ***Green){
	int mu, nu, i;
	fprintf(file, "#wn\t\t");
	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	fprintf(file, "%d,%d\t\t\t\t\t", mu, nu );
	fprintf(file, "\n");
	for(i=0; i<Nmax; i++){
		fprintf(file, "%.10lf\t", wn);
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	fprintf(file, "%.10lf\t%.10lf\t\t", GSL_REAL( Green[mu][nu][i] ), GSL_IMAG( Green[mu][nu][i] ) );
		fprintf(file, "\n");
	}
	return 0;
}
int print_lattice(FILE *file, gsl_complex ***Green, int start, int end){
	int mu, nu, i;
	fprintf(file, "#wn\t\t");
	for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)	fprintf(file, "%d,%d\t\t\t\t\t", mu, nu );
	fprintf(file, "\n");
	for(i=start; i<end; i++){
		fprintf(file, "%.10lf\t", wn);
		for(mu=0; mu<NU; mu++) for(nu=0; nu<NU; nu++)	fprintf(file, "%.10lf\t%.10lf\t\t", GSL_REAL( Green[i][mu][nu] ), GSL_IMAG( Green[i][mu][nu] ) );
		fprintf(file, "\n");
	}
	return 0;
}
int print_self(FILE *file, gsl_complex ***Green, int start, int end){
	int mu, nu, i;
	fprintf(file, "#wn\t\t");
	for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	fprintf(file, "%d,%d\t\t\t\t\t", mu, nu );
	fprintf(file, "\n");
	for(i=start; i<end; i++){
		fprintf(file, "%.10lf\t", wn);
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++)	fprintf(file, "%.10lf\t%.10lf\t\t", GSL_REAL( Green[i][mu][nu] ), GSL_IMAG( Green[i][mu][nu] ) );
		fprintf(file, "\n");
	}
	return 0;
}

double func(double p[], int n){
	int i, mu, nu;
	double chisquare = 0, ctemp;

	convert_p_to_egv(fegv, p);

	build_bath(n, Tbath, fegv);

#pragma omp parallel for ordered private(ctemp, i,mu,nu) shared(Greennew, Tbath, Nmax, n, Wwn) default(none) reduction(+: chisquare)
	for(i=0; i<Nmax; i++){
		ctemp = 0;
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
			ctemp +=
				(GSL_REAL(Greennew[n][i][mu][nu])-GSL_REAL(Tbath[i][mu][nu]))*(GSL_REAL(Greennew[n][i][mu][nu])-GSL_REAL(Tbath[i][mu][nu])) +
				(GSL_IMAG(Greennew[n][i][mu][nu])-GSL_IMAG(Tbath[i][mu][nu]))*(GSL_IMAG(Greennew[n][i][mu][nu])-GSL_IMAG(Tbath[i][mu][nu]));
		}
#pragma omp ordered
		chisquare += ctemp * WWN;
	}

	return chisquare / (Nmax) / tNC/tNC;
}//end of func

double func_zeroSOC(double p[], int n, int sp){
	int i, mu, nu;
	double chisquare = 0, ctemp;

	convert_p_to_egv_zeroSOC(fegv, p, sp);

	build_bath(n, Tbath, fegv);

#pragma omp parallel for ordered private(i,mu,nu,ctemp) shared(Greennew, Tbath, Nmax, n, nctot, sp, Wwn) default(none) reduction(+: chisquare)
	for(i=0; i<Nmax; i++){
		ctemp = 0;
		for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++){
			ctemp += (GSL_REAL(Greennew[n][i][2*mu+sp][2*nu+sp])-GSL_REAL(Tbath[i][2*mu+sp][2*nu+sp]))*(GSL_REAL(Greennew[n][i][2*mu+sp][2*nu+sp])-GSL_REAL(Tbath[i][2*mu+sp][2*nu+sp])) +
				 (GSL_IMAG(Greennew[n][i][2*mu+sp][2*nu+sp])-GSL_IMAG(Tbath[i][2*mu+sp][2*nu+sp]))*(GSL_IMAG(Greennew[n][i][2*mu+sp][2*nu+sp])-GSL_IMAG(Tbath[i][2*mu+sp][2*nu+sp]));
		}
#pragma omp ordered
		chisquare += ctemp *WWN;
	}

	return chisquare / (Nmax) / NC/NC;
}//end of func_zeroSOC


double my_f_zeroSOC(const gsl_vector *v, void *params){
	int **par = (int **) params;
	int i, mu, nu, n = par[3][0], sp = par[4][0];
	double chisquare = 0, ctemp=0;
	double p[Np_zeroSOC];
	for(i=0; i<Np_zeroSOC; i++)	p[i] = gsl_vector_get(v, i);

	convert_p_to_egv_zeroSOC(fegv, p, sp);

	build_bath(n, Tbath, fegv);

#pragma omp parallel for ordered private(i,mu,nu,ctemp) shared(Greennew, Tbath, Nmax, n, nctot, sp, Wwn) default(none) reduction(+: chisquare)
	for(i=0; i<Nmax; i++){
		ctemp = 0;
		for(mu=0; mu<nctot; mu++) for(nu=0; nu<nctot; nu++){
			ctemp += (GSL_REAL(Greennew[n][i][2*mu+sp][2*nu+sp])-GSL_REAL(Tbath[i][2*mu+sp][2*nu+sp]))*(GSL_REAL(Greennew[n][i][2*mu+sp][2*nu+sp])-GSL_REAL(Tbath[i][2*mu+sp][2*nu+sp])) +
				 (GSL_IMAG(Greennew[n][i][2*mu+sp][2*nu+sp])-GSL_IMAG(Tbath[i][2*mu+sp][2*nu+sp]))*(GSL_IMAG(Greennew[n][i][2*mu+sp][2*nu+sp])-GSL_IMAG(Tbath[i][2*mu+sp][2*nu+sp]));
		}
#pragma omp ordered
		chisquare += ctemp * WWN;
	}

	return chisquare / Nmax / NC/NC;
}//end of my_f_zeroSOC

void my_df_zeroSOC(const gsl_vector *v, void *params, gsl_vector *df){
	int i,j;
	double p2plus[Np_zeroSOC], pplus[Np_zeroSOC], pminus[Np_zeroSOC], p2minus[Np_zeroSOC], h=diff/2.0L;
	double p[Np_zeroSOC], dp[Np_zeroSOC];
	for(i=0; i<Np_zeroSOC; i++)	p[i] = gsl_vector_get(v, i);

	int **par = (int **)params;
	int *offset = par[0];
	int *number = par[1];
	int mystart = par[2][0], myend = par[2][1];
	int n = par[3][0];
	int sp = par[4][0];

	for(i=mystart; i<myend; i++){
		for(j=0; j<Np_zeroSOC; j++){
			p2plus[j] = p[j];
			pplus[j] = p[j];
			pminus[j] = p[j];
			p2minus[j] = p[j];
		}
		p2plus[i] += 2*h;
		pplus[i] += h;
		pminus[i] -= h;
		p2minus[i] -= 2*h;
		dp[i] = (-func_zeroSOC(p2plus, n, sp) + 8*func_zeroSOC(pplus, n, sp) - 8*func_zeroSOC(pminus, n, sp) + func_zeroSOC(p2minus, n, sp) ) / (12*h);
	}
	MPI_Allgatherv(MPI_IN_PLACE, myend-mystart, MPI_DOUBLE, dp, number, offset, MPI_DOUBLE, MPI_COMM_WORLD);

	for(i=0; i<Np_zeroSOC; i++){
		gsl_vector_set( df,  i, dp[i] );
	}
	//printf("myrank = %d, [%d:%d] of %d\n", myrank, mystart, myend, Np_zeroSOC);
}//end of my_df_zeroSOC

void my_fdf_zeroSOC(const gsl_vector *x, void *params, double *f, gsl_vector *df){

	*f = my_f_zeroSOC(x, params);
	my_df_zeroSOC(x, params, df);
}

double my_f(const gsl_vector *v, void *params){
	int **par = (int **) params;
	int i, mu, nu, n = par[3][0];
	double chisquare = 0, ctemp=0;
	double p[Np];
	for(i=0; i<Np; i++)	p[i] = gsl_vector_get(v, i);

	convert_p_to_egv(fegv, p);

	build_bath(n, Tbath, fegv);

#pragma omp parallel for ordered private(i,mu,nu,ctemp) shared(Greennew, Tbath, Nmax, n, Wwn) default(none) reduction(+: chisquare)
	for(i=0; i<Nmax; i++){
		ctemp = 0;
		for(mu=0; mu<tNC; mu++) for(nu=0; nu<tNC; nu++){
			ctemp += (GSL_REAL(Greennew[n][i][mu][nu])-GSL_REAL(Tbath[i][mu][nu]))*(GSL_REAL(Greennew[n][i][mu][nu])-GSL_REAL(Tbath[i][mu][nu])) +
				(GSL_IMAG(Greennew[n][i][mu][nu])-GSL_IMAG(Tbath[i][mu][nu]))*(GSL_IMAG(Greennew[n][i][mu][nu])-GSL_IMAG(Tbath[i][mu][nu]));
		}
#pragma omp ordered
		chisquare += ctemp * WWN ;
	}

	return chisquare / Nmax / tNC/tNC;
}//end of my_f

void my_df(const gsl_vector *v, void *params, gsl_vector *df){
	int i,j;
	double p2plus[Np], pplus[Np], pminus[Np], p2minus[Np], h=diff/2.0L;
	double p[Np], dp[Np];
	for(i=0; i<Np; i++)	p[i] = gsl_vector_get(v, i);

	int **par = (int **)params;
	int *offset = par[0];
	int *number = par[1];
	int mystart = par[2][0], myend = par[2][1];
	int n = par[3][0];

	for(i=mystart; i<myend; i++){
		for(j=0; j<Np; j++){
			p2plus[j] = p[j];
			pplus[j] = p[j];
			pminus[j] = p[j];
			p2minus[j] = p[j];
		}
		p2plus[i] += 2*h;
		pplus[i] += h;
		pminus[i] -= h;
		p2minus[i] -= 2*h;
		dp[i] = (-func(p2plus, n) + 8*func(pplus, n) - 8*func(pminus, n) + func(p2minus, n) ) / (12*h);
	}
	MPI_Allgatherv(MPI_IN_PLACE, myend-mystart, MPI_DOUBLE, dp, number, offset, MPI_DOUBLE, MPI_COMM_WORLD);

	for(i=0; i<Np; i++){
		gsl_vector_set( df,  i, dp[i] );
	}
	//printf("myrank = %d, [%d:%d] of %d\n", myrank, mystart, myend, Np);
}//end of my_df

void my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df){

	*f = my_f(x, params);
	my_df(x, params, df);
}


