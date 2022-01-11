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

#define DIMd	64						//for rdm_static.c
#define DIMd4	15						//for rdm_static.c
#define gcr(zre, zim)	gsl_complex_rect( zre, zim )		//for rdm_static.c

#include "opclass.cpp"
#include "rdmclass.cpp"
#include "mps.cpp"
#include "degop.cpp"

void compute_static(typebasis **basis, int count, char *save_directory);

//long seed = -134551;
double J, MU, DELTA, Uinter, JTD, SOC, Hz;
int beta, Nmax, nctot, nb, tnb, Spinmax, Blocks,Powns, Powns2, myrank=0, beta_real;
gsl_complex zero, ***Greennew, **Ecluster;

int main(int argc, char *argv[]){
	if(argc != 11){
		printf("Usage :: %s <input.mak> <Hz> <JTD> <SOC> <J> <U_init> <U_final> <DELTA> <save_directory> <beta_real>\n", argv[0]);
		exit(1);
	}
	sectorinfo **ginfo;
	int n, count = 0;
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
	beta_real = atoi( argv[10] );

	sprintf(save_directory, "%s", argv[9]);

	fd = fopen(argv[1], "r");
	nofile(fd, argv[1]);
	fscanf(fd, "%d", &count);
	fscanf(fd, "%lf", &Uinter);
	fclose(fd);

	typebasis **basis;
	basis = mkmatrixb(Powns2, 2);
	ginfo = (sectorinfo **) malloc( Ni * sizeof(sectorinfo*));
	for(n=0; n<Ni; n++){
		ginfo[n] = (sectorinfo *) malloc( Blocks*sizeof(sectorinfo) );
		build_bitbasis(Ns, basis, ginfo[n]);
		print_ginfo(Ns, ginfo[n]);
	}
	printf("U = %.2lf, MU = %.2lf (of UF%.2lf UT%.2lf) from '%s'\n", Uinter, MU, U_init, U_final, argv[1]);

	Ecluster = mkgscmatrixd( NU, NU);
	char paraEcluster[1024];
	sprintf(paraEcluster, "%s/Ecluster.dat", save_directory);
	load_Ecluster(paraEcluster, Ecluster);
	update_chem(Ecluster);

	freegscmatrixd( Ecluster , NU ) ;

	sprintf(pathinform, "%s/result/%c%.3lf", save_directory, control_char, control_para);
	mkdirectory(pathinform);

	compute_static(basis, count, save_directory);

	for(n=0; n<Powns2; n++)	free(basis[n]);		free(basis);
	int i;
	for(n=0; n<Ni; n++){
		for(i=0; i<Blocks; i++)
			freegscmatrixd(ginfo[n][i].ground, Nvec);
		free( ginfo[n] );
	}
	free(ginfo);

	return 0;
}

int read_mak(PNSpara *egv, int *degeneracy, double ***ginformation, char *save_directory, int count){
	int tol, n, i, j, k;
	char para[1024];
	double re, im;
	FILE *fd;
	sprintf(para, "%s/%c%.2lf_%dth.mak", save_directory, control_char, control_para, count);
	//printf("trying to read %s\n", para);

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

	for(n=0; n<Ni; n++){
		fscanf(fd, "%d", &degeneracy[n]);

		for(i=0; i<degeneracy[n]; i++) for(j=0; j<3; j++)	fscanf(fd, "%lf", &ginformation[n][i][j]);

		for(k=0; k<tnb; k++){
			fscanf(fd, "%lf", &re);
			egv[n].egbath[k] = re;
		}
		for(j=0; j<tNC; j++) for(k=0; k<tnb; k++){
			fscanf(fd, "%lf %lf", &re, &im);
			egv[n].hybrid[j][k] = gsl_complex_rect( re, im );
		}
	}
	fclose(fd);

	return 0;
}

void compute_static(typebasis **basis, int count, char *save_directory){
	PNSpara egv[Ni];
	int n, degen, ic, degeneracy[Ni], mu, nu, i, *table, index;
	//char paramag[1024], paradouble[1024], parafilling[1024] ;
	char paragnd[1024], 
	     pararesult[1024] ;
	//FILE *fmag, *ffil, *fdoc ;
	FILE *fground ;
	double ***ginformation;
	ginformation = mktritensord(Ni, 100, 3);
	//int Powtnb = int( pow(2,2*NB) ) ; 
	int Powtnc = int( pow(2, tNC) ) ; 
	gsl_complex **transform = mkgscmatrixd(NU, NU); init_transform(transform) ;

	sprintf(pararesult, "%s/result/%c%.3lf/", save_directory, control_char, control_para);

	table = mkvectori( Powns2 );
	for(i=0; i<Powns2; i++){
		index = (basis[i][0]<<Ns) + basis[i][1];
		table[index] = i;
	}

	for(n=0; n<Ni; n++){
		egv[n].egbath = mkvectord(tnb);	
		egv[n].hybrid = mkgscmatrixd(tNC, tnb);
	}

	for(mu=0; mu<Powtnc; mu++) {
		unsigned int  cup = (mu>>NC) ,
			      cdn = mu - (cup<<NC) ;
		int ncdum = 0 ;
		for(nu=0; nu<NC; nu++){
			int ncupmu = (((cup)>>nu)&0x01) ; 
			int ncdnmu = (((cdn)>>nu)&0x01) ;
			ncdum += ncupmu + ncdnmu ;
		}
	}

	char writemode[1024]		; sprintf(writemode,"w") ;
	for(ic=count; ic<count+1; ic++){
		if( read_mak(egv, degeneracy, ginformation, save_directory, ic) ){
			if( ic == count ){
				printf("no mak file exists!\n");
				exit(1);
			}
			continue;
		}

		double **Jvecimp = mkmatrixd( Ni, 3 ),
		       **Lvecimp = mkmatrixd( Ni, 3 ),
		       **Svecimp = mkmatrixd( Ni, 3 ),
		       **Mvecimp = mkmatrixd( Ni, 3 );
		for(n=0; n<Ni; n++) for(mu=0; mu<3; mu++) {
			Jvecimp[n][mu] = 0.;
			Lvecimp[n][mu] = 0.;
			Svecimp[n][mu] = 0.;
			Mvecimp[n][mu] = 0.;
		}

		for(n=0; n<Ni; n++){
			sprintf(paragnd, "%s/u%.2lf_%dth_n%d.gnd", save_directory, Uinter, ic, n);
			fground = fopen(paragnd, "rb");
			nofile(fground, paragnd);

			int *gndarr			= mkvectori( degeneracy[n] ) ;
			int *gndblockarr		= mkvectori( degeneracy[n] ) ;
			double *myenergyarr		= mkvectord( degeneracy[n] ) ;
			double *myweightarr		= mkvectord( degeneracy[n] ) ;
			gsl_complex	**groundarr	= (gsl_complex **)malloc( sizeof(gsl_complex*) * degeneracy[n] ) ;

			printf("Ground state(s) initilized :\n") ;
			for(degen=0; degen<degeneracy[n]; degen++){
				myenergyarr[degen]	=      ginformation[n][degen][0];
				gndarr[degen]		= (int)ginformation[n][degen][1];
				gndblockarr[degen]	= (int)ginformation[n][degen][2];
				groundarr[degen]	= mkgscvectord( gndblockarr[degen] );
				fread(groundarr[degen], sizeof(gsl_complex), gndblockarr[degen], fground);
				printf("ic=%d, n=%d, gnd=%d, gndblock=%d, degen=%3d/%3d\t, energy=%19.16f\t;\n",
						ic, n, gndarr[degen], gndblockarr[degen], degen, degeneracy[n], myenergyarr[degen]
						);
			}
			double gm = myenergyarr[0],
			       partition = 0.,
			       weightsum = 0.;
			printf("parition weight :\n");
			for(degen=0; degen<degeneracy[n]; degen++){
				if( beta_real ) myweightarr[degen] = exp( - ( myenergyarr[degen] - gm )*beta_real );
				else		myweightarr[degen] = 1.;
				partition += myweightarr[degen];
			}
			for(degen=0; degen<degeneracy[n]; degen++){
				printf("degen%3d : %19.16f\t%19.16f\n", degen, myweightarr[degen], myweightarr[degen]/partition );
				myweightarr[degen]  /= partition;
				weightsum += myweightarr[degen];
			}
			printf("sum : %19.16f\t%19.16f\n", partition, weightsum );
			write_doublearr( pararesult, ic, "boltzweight", writemode, myweightarr, degeneracy[n] );
			printf("Start calculating : ") ;

			Vecop Jop("Jop") ; 
			Vecop Jaop("Jaop") ; 
			Vecop Lop("Lop") ; 
			Vecop Sop("Sop") ; 
			Vecop Mop("Mop") ; 

			gsl_complex	**jzmat		= calcOpMat(  Jop.zj, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;
			gsl_complex	**lzmat		= calcOpMat(  Lop.zj, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;
			gsl_complex	**szmat		= calcOpMat(  Sop.zj, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;
			gsl_complex	**mzmat		= calcOpMat(  Mop.zj, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;

			gsl_complex	**jxmat		= calcVopXMat(   Jop, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;
			gsl_complex	**jymat		= calcVopYMat(   Jop, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;
			gsl_complex	**lxmat		= calcVopXMat(   Lop, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;
			gsl_complex	**lymat		= calcVopYMat(   Lop, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;
			gsl_complex	**sxmat		= calcVopXMat(   Sop, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;
			gsl_complex	**symat		= calcVopYMat(   Sop, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;
			gsl_complex	**mxmat		= calcVopXMat(   Mop, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;
			gsl_complex	**mymat		= calcVopYMat(   Mop, groundarr, gndarr, gndblockarr, degeneracy[n], basis, table ) ;

			gscmatrixTrReCartesian(	Jvecimp[n], jxmat, jymat, jzmat,  degeneracy[n] , myweightarr );
			gscmatrixTrReCartesian(	Svecimp[n], sxmat, symat, szmat,  degeneracy[n] , myweightarr );
			gscmatrixTrReCartesian(	Lvecimp[n], lxmat, lymat, lzmat,  degeneracy[n] , myweightarr );
			gscmatrixTrReCartesian(	Mvecimp[n], mxmat, mymat, mzmat,  degeneracy[n] , myweightarr );


			for(degen=0; degen<degeneracy[n]; degen++) free(groundarr[degen]);
			free( gndarr ) ; 
			free( gndblockarr ) ; 
			free( groundarr ) ; 
			free( myenergyarr ) ; 
			free( myweightarr ) ; 
			fclose(fground);
			Jop.freeVecop() ;
			Jaop.freeVecop() ;
			Lop.freeVecop() ;
			Sop.freeVecop() ;
			Mop.freeVecop() ;

			printf("Done.\n") ;
		}
		write_doubleMat(	pararesult, ic ,      "Jvec" , writemode ,    Jvecimp , Ni, 3 ) ;
		write_doubleMat(	pararesult, ic ,      "Lvec" , writemode ,    Lvecimp , Ni, 3 ) ;
		write_doubleMat(	pararesult, ic ,      "Svec" , writemode ,    Svecimp , Ni, 3 ) ;
		write_doubleMat(	pararesult, ic ,      "Mvec" , writemode ,    Mvecimp , Ni, 3 ) ;
		sprintf(writemode, "a") ;

		freematrixd( Jvecimp , Ni );
		freematrixd( Lvecimp , Ni );
		freematrixd( Svecimp , Ni );
		freematrixd( Mvecimp , Ni );
	}
	freetritensord(ginformation, Ni, 100);
	freegscmatrixd(transform, NU);
	for(n=0; n<Ni; n++){
		free( egv[n].egbath );
		freegscmatrixd( egv[n].hybrid, tNC );
	}
	free(table);
}

