#include <stdio.h>
#include <math.h>

/*
# Note
Last: 2/11/2019
motor cyclist under power-law air resistance
$ \vec{f} = - k v^{\beta} \vec{v} $

# Definitions
t        current time
T        final time
x,y      pos
vx,vy    vel
ax,ay    acc
*/

// physical params
#define beta    1./3  // power-law resistance
#define n       1000  // max time steps
#define g       9.80  // gravity
#define angle   42.5*M_PI/180 // shooting angle
#define speed   67    // launch speed
#define mass    250   // mass of cyclist
#define area    0.93  // area of cyclist
#define density 1.2
#define k       area*density/(2*mass)
#define dt      2*speed*sin(angle)/(g*n)
#define eps     1e-5  // small num

int    t,T;
double d,d2,d3,x[n],y[n],vx[n],vy[n],ax[n],ay[n];

void init(){
	t=0;

	/* init at t=0,1 */
	d=dt*dt/2,d2=2*dt,d3=dt/3,
	x[0]=0;
	y[0]=0;
	vx[0]=speed*cos(angle);
	vy[0]=speed*sin(angle);
	double v=sqrt(vx[0]*vx[0]+vy[0]*vy[0]);
	ax[0]=-k*pow(v,beta)*vx[0];
	ay[0]=-g-k*pow(v,beta)*vy[0];

	double p=vx[0]*ax[0]+vy[0]*ay[0];
	x[1]=x[0]+dt*vx[0]+d*ax[0];
	y[1]=y[0]+dt*vy[0]+d*ay[0];
	vx[1]=vx[0]+dt*ax[0]-d*k*(pow(v,beta)*ax[0]+beta*p*pow(v,beta-2)*vx[0]);
	vy[1]=vy[0]+dt*ay[0]-d*k*(pow(v,beta)*ay[0]+beta*p*pow(v,beta-2)*vy[0]);
	v=sqrt(vx[1]*vx[1]+vy[1]*vy[1]);
	ax[1]=-k*pow(v,beta)*vx[1];
	ay[1]=-g-k*pow(v,beta)*vy[1];
}

void update(){
	/* predictor */
	//x[t+2]=x[t]+d2*x[t+1];
	//y[t+2]=y[t]+d2*y[t+1];
	vx[t+2]=vx[t]+d2*vx[t+1];
	vy[t+2]=vy[t]+d2*vy[t+1];
	double v=sqrt(vx[t+2]*vx[t+2]+vy[t+2]*vy[t+2]);
	ax[t+2]=-k*pow(v,beta)*vx[t+2];
	ay[t+2]=-g-k*pow(v,beta)*vy[t+2];

	/* corrector: Simpson update */
	x[t+2]=x[t]+d3*(vx[t+2]+4*vx[t+1]+vx[t]);
	y[t+2]=y[t]+d3*(vy[t+2]+4*vy[t+1]+vy[t]);
	vx[t+2]=vx[t]+d3*(ax[t+2]+4*ax[t+1]+ax[t]);
	vy[t+2]=vy[t]+d3*(ay[t+2]+4*ay[t+1]+ay[t]);

	t++;
}

void iter(){
	/* iter till hitting ground or max steps */
	while(y[t]>-eps && t<n-2) update();
	T=t;
}

void data(){
	/* print data to file */
	FILE *f=fopen("data.csv","w");
	fprintf(f,"x,y,vx,vy,ax,ay\n");
	for(int t=0; t<T; t++)
		fprintf(f,"%lf,%lf,%lf,%lf,%lf,%lf\n",
		x[t],y[t],vx[t],vy[t],ax[t],ay[t]);
	fclose(f);
}

int main(){
	init();
	iter();
	data();
	return 0;
}
