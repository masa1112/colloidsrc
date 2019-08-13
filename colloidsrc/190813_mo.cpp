
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "ran.h"
#include "dSFMT.h"
#include "gasdev.h"
#include <time.h>

//global

const double kx=1.0;
const double kf=1000.;
const double kw=1000.;
double F_int_x,F_x;
double F_int_y,F_y;
double F_int_z,F_z,F_gravity,F_floor;
double F_wall_x0,F_wall_x10,F_wall_y0,F_wall_y10;
double forx,fory,forz;
double x_ij,y_ij,z_ij,l_ij;
double F_calc1,F_calc2;

const double dt=0.01;
const double tmax = 10000;//10000;
const double xmax =240.;
const double ymax =240.;
const double enl=5.;//raw file enlarge
int nmax=(int)(tmax/dt);
const int mabiki = 1000;//100000;
double sqdt = sqrt(dt);
int ixmax = (int)(xmax*enl);
int iymax = (int)(ymax*enl);

typedef struct {
    double x, y, z;
    } Force;
Force force;
/*
Force getforce(double forx,double fory,double forz){
    return (Force){forx,fory,forz};
}
*/
typedef struct {
    int c1,c2;
    double c3,c4,c5,c6,c7,c8,c9,c10;
} Cnst;

Cnst const_value;

Cnst getcnst(int c1,int c2,double c3,double c4,double c5,double c6,double c7,double c8,double c9,double c10){
    return (Cnst){c1,c2,c3,c4,c5,c6,c7,c8,c9,c10};
}


//open file for writing binary style
FILE *fopen_wb( char filename[] ){

	FILE *fp;
	if(  (fp=fopen(filename,"wb")) == NULL  ){
		printf( "cannot open output file\n" );
		exit(1);
	}
	return fp;
}

//draw colloid motion in two dimensional manner height is super imposed
void change_to_coordinates( double *px, double *py, double *pz, double *X,
		int p_num, double *rad, int ixmax, int iymax, double enlx, double enly) {
	int i, j, num;

	for( i = 0; i < ixmax; i++ ){
		for( j = 0; j < iymax; j++ ){
			X[i*(iymax)+j] = 0.0;
		}
	}

	for(num=0; num<p_num; num++){
		for(i=(int)ceil((px[num]-rad[num])*enlx); i<(int)((px[num]+rad[num])*enlx)+1; i++){
			for(j=(int)ceil((py[num]-rad[num])*enly); j<(int)((py[num]+rad[num])*enly)+1; j++){
				//the rasius caluculation is done only when i and j is within a square around the particle center
				//this will reduce the distance calculation
				if(((i-px[num]*enlx)*(i-px[num]*enlx) + (j-py[num]*enly)*(j-py[num]*enly)) <= rad[num]*enlx*rad[num]*enly){
					if(i>=0 && i<ixmax && j>=0 && j<iymax){
						
						if(X[i*(iymax)+j] < pz[num]){
							X[i*(iymax)+j] = pz[num];
						//drawing region must be in a reserved array
						}
					}
				}
			}
		}
	}
}

inline double pow_i(double x, int imax){
	double output = 1.;
	for (int i=0; i<imax; i++){
		output *=x;
	}
	return (output);
}

//output raw file
void output_raw(double *u, char fname[], int ixmax, int iymax){
	FILE *fp;
	fp = fopen_wb(fname);
	fwrite(u, sizeof(double), (ixmax)*(iymax)+1, fp);
	fclose (fp);
}

Force Force_all(double *x,double *y,double *z,double *l,bool Q,int i,int n,Cnst all){
	//typedef struct { double x, y, z; } Force;
	/*Force getforce(double forx,double fory,double forz){
        return (Force){forx,fory,forz};
    };*/

	double forx = 0;
	double fory = 0;
	double forz = 0;
    

	double F_x=0;
	double F_y=0;
	double F_z=0;
	double F_int_x=0;
	double F_int_y=0;
	double F_int_z=0;
	double F_floor=0;
	double F_wall_x0=0;
	double F_wall_x10=0;
	double F_wall_y0=0;
	double F_wall_y10=0;
	
	double x_ij = 0;
	double y_ij = 0;
	double z_ij = 0;
	double l_ij = 0;

	int j = 0;
	double r = 0;
	double rho2 = 0;
	double r2 = 0;
	double br2 = 0;
	double br = 0;
    /*
    int N = all.c1;
    int N1 = all.c2;
    double alpha = all.c3; //5;
    double beta = all.c4;//2000;
    double sigma = all.c5;//0.1;
    double kappa = all.c6; //5.;
    double lambda = all.c7;
    double g = all.c8;
    double kz= all.c9;
    double kz2= all.c10;
	*/
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
	F_gravity=- all.c8 * l[i] * l[i];//);

	//F_floor
	if(z[i]<l[i]){
		F_floor=kf*(l[i]-z[i]);
	}

	if(Q == 0){
		for(j=0;j<all.c1;j++){
			if(i!=j){
				if(r<(l[i]+l[j])){

					x_ij=(x[i]-x[j]);
					y_ij=(y[i]-y[j]);
					z_ij=(z[i]-z[j]);
					r = sqrt((x_ij * x_ij)+(y_ij * y_ij)+(z_ij * z_ij));
					l_ij = (l[i]+l[j]) / r;

					F_int_x += kf * ((l_ij) - 1) * x_ij;
					F_int_y += kf * ((l_ij) - 1) * y_ij;
					F_int_z += kf * ((l_ij) - 1) * z_ij;
				}
			}//if終了
		}//j終了
		forx = F_int_x * n * dt * 0.01 + F_wall_x0 + F_wall_x10;
		fory = F_int_y * n * dt * 0.01 + F_wall_y0 + F_wall_y10;
        forz = F_int_z * n * dt * 0.01 + F_floor;

        return (Force){forx,fory,forz};
	}
	else{
		for(j=0;j<all.c1;j++){
			if(i!=j){

				x_ij=(x[i]-x[j]);
				y_ij=(y[i]-y[j]);
				z_ij=(z[i]-z[j]);
				
				rho2 = (x_ij * x_ij) + (y_ij * y_ij);
				r2 = rho2 + (z_ij * z_ij);
				br2 = rho2 + pow_i(z[j],2);
				r=sqrt(r2);
				l_ij = (l[i]+l[j]) / r;
				br=sqrt(br2);

				F_calc1 = all.c3 * (l[i] * l[i]) * pow_i(l[j],3) / pow_i(r,5) * exp(- all.c6 * r);
				F_calc2 = F_calc1 * (1. - 3. * (z_ij * z_ij) / r2) - all.c4 * pow_i(l[j],3) * l[i] / (z[i] - l[i] + all.c7) * z[j] / pow_i(br,5) * exp(- all.c6 *br);

				F_x += F_calc2 * x_ij;
				F_y += F_calc2 * y_ij;
				F_z += F_calc1 * (3. - 5. * (z_ij * z_ij) / r2) * z_ij + all.c4 * pow_i(l[j],3) * l[i] * z[j] * (2. - 5. * rho2 / br2 - all.c6 * rho2 / br) / pow_i(br,5) * exp(- all.c6 *br) * all.c9 / (z[i] - l[i] + all.c7);

				if(r<(l[i]+l[j])){
					
					F_int_x += kf * ((l_ij)-1) * x_ij;
					F_int_y += kf * ((l_ij)-1) * y_ij;
					F_int_z += kf * ((l_ij)-1) * z_ij;
				}
			}
		}//j終了
		F_z += all.c4 * pow_i(l[i] / z[i],4) * exp(- all.c6 * z[i]) * all.c10 / (z[i]-l[i] + all.c7 ); 	

		forx = F_x + F_int_x * 0.01 + F_wall_x0 + F_wall_x10;
		fory = F_y + F_int_y * 0.01 + F_wall_y0 + F_wall_y10;
		forz = F_z + F_int_z * 0.01 + F_gravity + F_floor;

        return (Force){forx,fory,forz};
	}
}

int main(int argc, char* argv[]){
	printf("program start\n");
	
	if(argc!=10){
		std::cout<<"error!!\n";
		std::cout<<"this file require 5 arguments (n2:256) (n3:0) (alpha:1) (beta:10-1000)"
			<<"(sigma:0.1) (kaappa:5) (g=1.) (kz=0.01) (kz2=0.01_but3to1)"<<"\n";
		return 1;
	}
	int N = atoi(argv[1])+atoi(argv[2]);
	int N1 = atoi(argv[1]);
	double alpha = atof(argv[3]); //5;
	double beta = atof(argv[4]);//2000;
	double sigma = atof(argv[5]);//0.1;
	double kappa = atof(argv[6]); //5.;
	double lambda = 1./kappa;
	double g = atof(argv[7]);
	double kz= atof(argv[8]);
	double kz2= atof(argv[9]);
    
    Cnst const_value = getcnst(N,N1,alpha,beta,sigma,kappa,lambda,g,kz,kz2);
    
	clock_t start,end;
	double stime;
	time_t timer;

	FILE *file1;
	std::ofstream os;
	int i,j,n;

	double *x,*y,*z,*v,*xnew, *ynew,*znew, *l, *rawXY;
	x = new double[N];
	y = new double[N];
	z = new double[N];
	v = new double[N];
	xnew = new double[N];
	ynew = new double[N];
	znew = new double[N];
	rawXY = new double[ixmax*iymax];
	l = new double[N];

	double t,r, r2, br, br2, rho2;//,g=0.01;
	
	char filename[301];

	bool Q;
    double *x_new1,*y_new1,*z_new1;
    x_new1 = new double[N];
    y_new1 = new double[N];
    z_new1 = new double[N];

    int syoki = 10./dt;

	srand((unsigned int)time(NULL));

	std::ostringstream command,name1;//, name3;

	command.str("");
	command<<"mkdir -p dat";
	system((command.str()).c_str());

	name1.str("");
	name1<<N1<<"_"<<(N-N1)<<"_"<<alpha<<"_"<<beta<<"_"<<sigma
		<<"_"<<kappa<<"_"<<g<<"_"<<kz<<"_"<<kz2;

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
			l[i]=1.0;
		}
		else{
			l[i]=1.5;
		}
		fprintf(file1, "%lf\t%lf\t%lf\t%lf\t%d\n",x[i],y[i],z[i],l[i],i);
	}
	fclose(file1);
	start=clock();
	//typedef struct { double x, y, z; } Force;
	Force force1,force2;

	for(n=1;n<=nmax;n++){ 
		for(i=0;i<N;i++){
			if(n < syoki ){
				bool Q = 0;
			}
			else{
				Q = 1;
			}
			//force1 = Force_all(x,y,z,l,Q,i,n,const_value);

			//xnew[i] = x[i] + (force1.x ) * dt + (sigma/pow(l[i],0.5)*gasdev())*sqdt;
			//ynew[i] = y[i] + (force1.y ) * dt + (sigma/pow(l[i],0.5)*gasdev())*sqdt;
			//znew[i] = z[i] + (force1.z ) * dt + (sigma/pow(l[i],0.5)*gasdev())*sqdt;
			
			force1 = Force_all(x,y,z,l,Q,i,n,const_value);
			
			//x_new1[i] = x[i] + (force1.x)*dt + (sigma/pow(l[i],0.5)*gasdev())*sqdt;
			//y_new1[i] = y[i] + (force1.y)*dt + (sigma/pow(l[i],0.5)*gasdev())*sqdt;
			//z_new1[i] = z[i] + (force1.z)*dt + (sigma/pow(l[i],0.5)*gasdev())*sqdt;
			
			x_new1[i] = x[i] + (force1.x)*dt;
			y_new1[i] = y[i] + (force1.y)*dt;
			z_new1[i] = z[i] + (force1.z)*dt;
			
			force2 = Force_all(x_new1,y_new1,z_new1,l,Q,i,n+1,const_value);
			
			//xnew[i]=x[i]+(force1.x + force2.x)* dt/2;
			//ynew[i]=y[i]+(force1.y + force2.y)* dt/2;
			//znew[i]=z[i]+(force1.z + force2.z)* dt/2;
			
			xnew[i]=x[i]+(force1.x + force2.x)* dt/2 + (sigma/pow(l[i],0.5)*gasdev())*sqdt;
			ynew[i]=y[i]+(force1.y + force2.y)* dt/2 + (sigma/pow(l[i],0.5)*gasdev())*sqdt;
			znew[i]=z[i]+(force1.z + force2.z)* dt/2 + (sigma/pow(l[i],0.5)*gasdev())*sqdt;
				
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

	delete[] l,x,y,z,xnew,ynew,znew,x_new1,y_new1,z_new1,rawXY;

	return (0);
}
