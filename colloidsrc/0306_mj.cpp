//20181107大きさを加えた複数粒子
//20181108粒子のサイズごとにファイル出力,3次元化
//20181124相図作成プログラム
//20181211パラメータ改定
//20181214powの廃止
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
//#include <random>
#include "ran.h"
#include "dSFMT.h"
#include "gasdev.h"
#include <time.h>
//#include <cmath>
//#include <limits.h>



	/* ファイルの入力関数 */
FILE *fopen_wb( char filename[] ){

	FILE *fp;
	if(  (fp=fopen(filename,"wb")) == NULL  ){
		printf( "cannot open output file\n" );
		exit(1);
	}
	return fp;
}


void change_to_coordinates( double *px, double *py, double *pz, double *X, 
		int p_num, double *rad, int ixmax, int iymax, double enlx, double enly)
{
	int i, j, num;
	//配列初期化
	for( i = 0; i < ixmax; i++ ){
		for( j = 0; j < iymax; j++ ){
			X[i*(iymax)+j] = 0.0;
		}
	}

	for(num=0; num<p_num; num++){
		for(i=(int)ceil((px[num]-rad[num])*enlx); i<(int)((px[num]+rad[num])*enlx)+1; i++){
			for(j=(int)ceil((py[num]-rad[num])*enly); j<(int)((py[num]+rad[num])*enly)+1; j++){
				if(((i-px[num]*enlx)*(i-px[num]*enlx) + (j-py[num]*enly)*(j-py[num]*enly)) <= rad[num]*enlx*rad[num]*enly){
					if(i>=0 && i<ixmax && j>=0 && j<iymax){
						X[i*(iymax)+j] += pz[num];
					}
				}
			}
		}
	}
}
void output_raw(double *u, char fname[], int ixmax, int iymax){
	FILE *fp;
	fp = fopen_wb(fname);
	fwrite(u, sizeof(double), (ixmax)*(iymax)+1, fp);
	fclose (fp);
}

#define PI M_PI

int main(int argc, char* argv[]){
	printf("program start\n");
	const double dt=0.001 ;
	const double tmax = 10000;//10000;
	const double xmax =120.;
	const double ymax =120.;
	const double enl=5.;//raw file enlarge
	int ixmax = (int)(xmax*enl);
	int iymax = (int)(ymax*enl);

	if(argc!=10){
		std::cout<<"error!!\n";
		std::cout<<"this file require 5 arguments (n2:256) (n3:0) (alpha:1) (beta:10-1000) (zeta:0.1) (kaappa:5) (g=1.) (kz=0.01) (kz2=0.01_but3to1)\n";
		return 1;
	}
	int N = atoi(argv[1])+atoi(argv[2]);
	int N1 = atoi(argv[1]);
	double alpha = atof(argv[3]); //5;
	double beta = atof(argv[4]);//2000;
	double zeta = atof(argv[5]);//0.1;
	double kappa = atof(argv[6]); //5.;
	double g = atof(argv[7]);
	double kz= atof(argv[8]);
	double kz2= atof(argv[9]);
	clock_t start,end;
	double stime;
	time_t timer;

	const double tau = 1;

	FILE *file1;
	std::ofstream os;
	int i,j,n;

	double *x,*y,*z,*xnew, *ynew,*znew, *l, *rawXY;
	x = new double[N];
	y = new double[N];
	z = new double[N];
	xnew = new double[N];
	ynew = new double[N];
	znew = new double[N];
	rawXY = new double[ixmax*iymax];
	l = new double[N];

	double t,r;//,g=0.01;
	double ex,ey,ez;

	double kx=1.0;//,kz=0.01,kz2=0.01;//EHDの強さの比をいじる
	double kf=1000.;
	double kw=1000.;
	double F_int_x,F_x;
	double F_int_y,F_y;
	double F_int_z,F_z,F_gravity,F_floor;

	char filename[301];

	int nmax=(int)(tmax/dt);
	int mabiki = 10000;//100000;

	double F_wall_x0,F_wall_x10,F_wall_y0,F_wall_y10;//壁の設定

	srand((unsigned int)time(NULL));

	std::ostringstream command,name1;//, name3;

	command.str("");
	command<<"mkdir -p dat";
	system((command.str()).c_str());

	//kz=0.01;
	//kz2=0.01;
	name1.str("");
	name1<<N1<<"_"<<(N-N1)<<"_"<<alpha<<"_"<<beta<<"_"<<zeta<<"_"<<kappa<<"_"<<g<<"_"<<kz<<"_"<<kz2;

	command.str("");
	command<<"mkdir -p text_"<<name1.str();
	system((command.str()).c_str());

	command.str("");
	command<<"mkdir -p raw_"<<ixmax<<"_"<<iymax<<"_"<<name1.str();
	system((command.str()).c_str());

	sprintf(filename,"text_%s/data_0000.txt",(name1.str()).c_str());
	file1=fopen(filename,"w");

	for(i=0;i<N;i++){
		x[i]=rand()/(RAND_MAX+1.0)*(xmax-4.)+2.;
		y[i]=rand()/(RAND_MAX+1.0)*(ymax-4.)+2.;
		z[i]=3.;//+0.1*rnz[i];
		//大きさの代入
		if(i<N1){
			l[i]=0.005;
		}
		else{
			l[i]=1.0;
		}
		fprintf(file1, "%lf\t%lf\t%lf\t%lf\t%d\n",x[i],y[i],z[i],l[i],i);
	}
	fclose(file1);
	start=clock();
	for(n=1;n<=nmax;n++){ 
		for(i=0;i<N;i++){
			F_x=0;
			F_y=0;
			F_z=0;
			F_int_x=0;
			F_int_y=0;
			F_int_z=0;
			F_floor=0;
			F_wall_x0=0;
			F_wall_x10=0;
			F_wall_y0=0;
			F_wall_y10=0;

			//wall_x
			if(x[i]<l[i]){
				F_wall_x0=kw*(l[i]-x[i]);
			}
			if(x[i]>xmax-l[i]){
				F_wall_x10=kw*(xmax-l[i]-x[i]);
			}

			//wall_y
			if(y[i]<l[i]){
				F_wall_y0=kw*(l[i]-y[i]);
			}
			if(y[i]>ymax-l[i]){
				F_wall_y10=kw*(ymax-l[i]-y[i]);
			}

			//F_gravity
			F_gravity=-g*l[i]*l[i];//);

			//F_floor
			if(z[i]<l[i]){
				F_floor=kf*(l[i]-z[i]);
			}
			if(n<(int)(10./dt)){
				for(j=0;j<N;j++){
					if(i!=j){
						//r=1.0;
						r=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]));
						ex=(x[i]-x[j])/r;
						ey=(y[i]-y[j])/r;
						ez=(z[i]-z[j])/r;


						//F_int
						if(r<(l[i]+l[j])){
							F_int_x+=kf*((l[i]+l[j])-r)*ex;
							F_int_y+=kf*((l[i]+l[j])-r)*ey;
							F_int_z+=kf*((l[i]+l[j])-r)*ez;
						}
					}//if終了
				}//j終了
				xnew[i]=x[i]+(F_int_x*n*dt*0.1+F_wall_x0+F_wall_x10)*dt+(zeta/pow(l[i],1.5)*gasdev())*sqrt(dt);
				ynew[i]=y[i]+(F_int_y*n*dt*0.1+F_wall_y0+F_wall_y10)*dt+(zeta/pow(l[i],1.5)*gasdev())*sqrt(dt);
				znew[i]=z[i]+(F_int_z*n*dt*0.1+F_floor)*dt+(zeta/pow(l[i],1.5)*gasdev())*sqrt(dt);
				//znew[i]=10;
			}
			else{
				for(j=0;j<N;j++){
					if(i!=j){
						//r=1.0;
						r=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]));
						ex=(x[i]-x[j])/r;
						ey=(y[i]-y[j])/r;
						ez=(z[i]-z[j])/r;

						//F_dipole
						F_x+=(alpha/pow(r,5)*(x[i]-x[j])*(1.-3.*((z[i]-z[j])/r)*((z[i]-z[j])/r))*pow(l[i],2)*pow(l[j],3)
								+beta/kappa/pow((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+z[j]*z[j],2.5)*(-kx*(x[i]-x[j])*z[j])*exp(-kappa*(z[i]-l[i]))*pow(l[j],3)/l[i]);
						F_y+=(alpha/pow(r,5)*(y[i]-y[j])*(1.-3.*((z[i]-z[j])/r)*((z[i]-z[j])/r))*pow(l[i],2)*pow(l[j],3)
								+beta/kappa/pow((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+z[j]*z[j],2.5)*(-kx*(y[i]-y[j])*z[j])*exp(-kappa*(z[i]-l[i]))*pow(l[j],3)/l[i]);
						F_z+=(alpha/pow(r,5)*(z[i]-z[j])*(3.-5.*((z[i]-z[j])/r)*((z[i]-z[j])/r))*pow(l[i],2)*pow(l[j],3)
								+beta/kappa/pow((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])
								+z[j]*z[j],2.5)*kz*(2-5*((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))/((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+z[j]*z[j]))*exp(-kappa*(z[i]-l[i]))*pow(l[j],3)/l[i]);

						if(r<(l[i]+l[j])){
							F_int_x+=kf*((l[i]+l[j])-r)*ex;
							F_int_y+=kf*((l[i]+l[j])-r)*ey;
							F_int_z+=kf*((l[i]+l[j])-r)*ez;
						}
					}
				}//j終了
				F_z+=2.*beta/kappa/(z[i]*z[i]*z[i]*z[i])*kz2*exp(-kappa*(z[i]-l[i]))*l[i]*l[i];
				xnew[i]=x[i]+(F_x+F_int_x+F_wall_x0+F_wall_x10)*dt+(zeta/pow(l[i],1.5)*gasdev())*sqrt(dt);
				ynew[i]=y[i]+(F_y+F_int_y+F_wall_y0+F_wall_y10)*dt+(zeta/pow(l[i],1.5)*gasdev())*sqrt(dt);
				znew[i]=z[i]+(F_z+F_int_z+F_gravity+F_floor)*dt+(zeta/pow(l[i],1.5)*gasdev())*sqrt(dt); 
			}
		}


		//transfer data from new one to original

		for(i=0;i<N;i++){
			x[i]=xnew[i];
			y[i]=ynew[i];
			z[i]=znew[i];
		}//iend
		if(n%mabiki==0){
			sprintf(filename,"text_%s/data_%04d.txt",(name1.str()).c_str(),n/mabiki);
			file1=fopen(filename,"w");
			for(i=0;i<N;i++){
				fprintf(file1,"%lf\t%lf\t%lf\t%lf\t%d \n",x[i],y[i],z[i],l[i],i);
			}
			fclose(file1);

			sprintf(filename,"raw_%d_%d_%s/%s_%04d.raw",ixmax, iymax, (name1.str()).c_str(),(name1.str()).c_str(), n/mabiki);
			change_to_coordinates(x, y, z, rawXY, N, l, ixmax, iymax, enl, enl);
			output_raw(rawXY, filename, ixmax, iymax);

		}
	}//n終了
	end = clock();
	stime = (double)(end-start)/(double)CLOCKS_PER_SEC;
	time(&timer);
	command.str("");
	command<<"log.txt";
	os.open((command.str()).c_str(), std::ios::out | std::ios::app);
	os<<name1.str()<<"\ttime: "<<stime<<"\tlocal time: "<<ctime(&timer);
	os.close();

	command.str("");
	command<<"dat_"<<name1.str()<<".dat";
	os.open((command.str()).c_str(), std::ios::binary);
	os.write((char*)x, sizeof(double)* N);
	os.write((char*)y, sizeof(double)* N);
	os.write((char*)z, sizeof(double)* N);
	os.close();

	delete[] l,x,y,z,xnew,ynew,znew, rawXY;

	return (0);
}
