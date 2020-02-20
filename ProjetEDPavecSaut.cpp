// EE/EI/CN nombre I nombre N
/*ce fichier génère un fichier nommé 'graph' que vous pouvez charger avec gnuplot avec la commande suivante :
plot 'graph' using 1:2 with lines, '' u 1:3 with lines
*/

#include </usr/include/eigen3/Eigen/LU>
#include <iostream>
#include <ctime>
#include <math.h>
// pour ptr_fun
#include <fstream>
#include <functional>
#include <algorithm>
using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;
float u0(float s,float K){
	return fmax(K-s,0.0);
}

float ul(float t,int K=100,float r=0.05,int Smin=80){
	return K*exp(-r*t)-Smin;
}

float ur(float t){
	return 0.0;
}

float u0k(float s){
	return u0(s,100);
}


int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

VectorXf q(const VectorXf a,const VectorXf b,int I,float t){
	VectorXf qq=VectorXf::Zero(I);
	qq(0)=(-a(0)+b(0))*ul(t);
	qq(I-1)=(-a(I-1)-b(I-1))*ur(t);
	return qq;

}

float stox(float s){
	return log(s/100);;
}


float Normal(float x){
 return	0.5*erf(x/sqrt(2))+0.5;
}


VectorXf BS(float t,float K,float sigma,float r,const VectorXf s){
	VectorXf PP;	
	float dm;
	float dp;
	if (t==0){
	PP=s.unaryExpr(&u0k);
	}
	else{
	PP=VectorXf::Ones(s.size())*K*exp(-r*t);
		float tau=pow(sigma,2)*t;
	for(int i=0;i<=s.size()-1;i++){
		if(s(i)>0){
			dm=(log(s(i)/K) + r*t - 0.5*tau)/ sqrt(tau);	
			dp=(log(s(i)/K) + r*t + 0.5*tau)/ sqrt(tau);
			PP(i)=K*exp(-r*t)*(Normal(-dm)) - s(i)*(Normal(-dp));
		}	
	
	
	}
	}
	return PP;
	}




float g(float x){
	return (exp((-pow(log(x)-0.90,2)/(2*pow(0.45,2))))/(sqrt(2*M_PI)*0.45*x))*x;
//	return (exp((-pow(log(x)-0.90,2)/(2*pow(0.45,2))))/(sqrt(2*M_PI)*0.45*x))*x;
}







int main(int argc, char* argv[])
{

float K, r, sigma, T, Smin, Smax,gamma,mu,lambda;
K=100;sigma=0.15;r=0.05;T=0.25;Smin=80;Smax=120;gamma=0.45;mu=-0.90;lambda=0.10;


int I=std::stoi(argv[2]);
int N=std::stoi(argv[3]);

const char* SCHEMA=argv[1];
const char* CENTRAGE="CENTRE";

float Xmin,Xmax,Ymin,Ymax,err_scale;
Xmin=Smin;Xmax=Smax;Ymin=-20;Ymax=K;
err_scale=0;
int deltan=N/10;
float dt=T/N;
float h=(Smax-Smin)/(I+1);
VectorXf s=VectorXf::LinSpaced(I, 80,120);
float cfl=(dt/pow(h,2))*pow(sigma*Smax,2);

MatrixXf A = MatrixXf::Zero(I,I);
VectorXf alpha=VectorXf::Zero(I);
VectorXf bet=VectorXf::Zero(I);

if (strcmp(CENTRAGE, "CENTRE") == 0) 
{

	alpha=((pow(sigma,2))/2)*(s.array().pow(2))/pow(h,2);
	bet=r*s.array()/(2*h);
	for (int i=0;i<=I-1;i++) A(i,i)=2*alpha(i)+(r-lambda*(exp(mu+pow(gamma,2)/2)-1));
	for (int i=1;i<=I-1;i++) A(i,i-1)=-alpha(i)+bet(i);
	for (int i=0;i<=I-2;i++) A(i,i+1)=-alpha(i)-bet(i);
	std::cout << "Here is the matrix A:\n" << A << std::endl;
} 
else if (strcmp(CENTRAGE, "DROIT") == 0)
{
  // do something else
}
/* more else if clauses */
else /* default: */
{
}

VectorXf P=s.unaryExpr(&u0k);


// Starting cputime counter
std::clock_t c_start = std::clock();

MatrixXf Id=MatrixXf::Identity(I,I); 
//approx intergral
int Pm=3*I/5;
VectorXf xx=VectorXf::LinSpaced((3*I)-2, 40,160);
VectorXf gvect0=xx.unaryExpr(&stox);
VectorXf gvect1=gvect0.array().exp();
VectorXf gvect2=gvect1.unaryExpr(&g);


std::cout << " matrice vect1 :\n" << gvect1 << std::endl;
std::cout << " matrice vect2 :\n" << gvect2 << std::endl;
MatrixXf G=MatrixXf::Zero(I,I);
MatrixXf D=MatrixXf::Zero(I,I);
VectorXf Vminus=VectorXf::Zero(I);
for (int i=0; i<=I-1;i++){
	for(int j=0;j<=I-1;j++){
		if(abs(i-j)<=Pm) G(i,j)=gvect2(I-1+j-i);
		if(j>=i && j<Pm) D(i,j)=gvect2(I-1-(Pm)+j-i);
	}
}
std::cout << " matrice G :\n" << G << std::endl;

std::cout << " matrice D :\n" << D << std::endl;



float t,t1,to;
VectorXf q0;
VectorXf q1;
VectorXf Pex;
VectorXf Pmerton;
float errLI;
VectorXf sn;
float sigman;
float rn;
int Nmerton=10;
MatrixXf Result(I,3);
for (int n=0 ; n<=N-1 ; n++)
{
	t=n*dt;

	for(int o=Pm;o>=0;o--) Vminus(o)=ul(t);

	if (strcmp(SCHEMA, "EE") == 0) 
	{

		// sans saut
		// 	P=(Id-dt*A)*P -dt*q(alpha,bet,I,t);
		//avec saut
		P=((1-dt*lambda)*Id-dt*((-h*lambda*G)+A))*P -dt*(q(alpha,bet,I,t) +dt*h*D*lambda*Vminus);
	}
	else if (strcmp(SCHEMA, "EI") == 0)
	{
		t1=t+dt; //Vminus en fonction de t1
		P=(Id+dt*A).lu().solve((1-dt*lambda)*P-dt*q(alpha,bet,I,t1)+dt*h*lambda*G*P+lambda*h*dt*D*Vminus);
	}

	/* more else if clauses */
	else	//Crank-nicholson
	{
		q0=q(alpha,bet,I,t);
		q1=q(alpha,bet,I,t+dt);
		P = (Id +dt/2*A).lu().solve((((1-lambda*dt)*Id) -dt/2*A) *P - dt*(q0+q1)/2 +dt/2*h*D*lambda*Vminus +dt/2*h*G*lambda*P);
	}

	if ((n+1)%deltan==0){
		//Pex=BS(t1,K,sigma,r,s);

		Pmerton=VectorXf::Zero(I);
		to=T-t;
		for(int z=0;z<=Nmerton;z++){
			rn=r-lambda*(exp((mu+pow(gamma,2))/2) -1) + z*(mu/t);
			sigman=sqrt(pow(sigma,2)+z*(pow(gamma,2)/t));
			sn=s*exp(z*pow(gamma,2)/2);

			Pmerton=Pmerton + exp(-lambda*to)*(pow((lambda*to),z)/factorial(z))*exp(rn*to)*BS(t,K,sigman,rn,sn);

		}
	Pex=exp(-r*to)*Pmerton;

	Result.col(0)<<s;
	Result.col(1)<<P;
	Result.col(2)<<Pex;
	std::cout << "S P Pex\n" <<Result<< std::endl;


	errLI=(P-Pex).lpNorm<Infinity>();	
	std::cout << "Here is the vector errLI:\n" << errLI << std::endl;
	}



}

Result.col(0)<<s;
Result.col(1)<<P;
Result.col(2)<<Pex;


std::cout << "Methode utilisé   " <<SCHEMA<< std::endl;
std::cout << "Cfl \n" <<cfl<< std::endl;
std::cout << "S P Pex\n" <<Result<< std::endl;


//creation graph file
std::ofstream outfile("graph");
outfile << "S P Pex\n" <<Result<< std::endl;
outfile.close();




std::clock_t c_end = std::clock();
long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";






}
