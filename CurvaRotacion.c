#include <math.h>
#include <random>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

#define SIZE_FILE 300
#define N 50000

using namespace std;

const double bb = 0.2497;
const double bd = 5.16;
const double ad = 0.3105;
const double ah = 64.3;
const double sigma = 10.0;

void load_data( double R[], double VR[] ){
	ifstream file("RadialVelocities.dat");

	string t;
	getline(file, t);

	int k=0, l=0;
	for(int i = 0; i < 2*SIZE_FILE; i++){
		if( i%2 == 0 ){
			file >> R[k];
			k += 1;
		}
		else{
			file >> VR[l];
			l += 1;
		}
	}
}

void my_model( double* R, double* VR, double Mb, double Md, double Mh){
	for(int i = 0; i < SIZE_FILE; i++ ){
		VR[i] = (sqrt(Mb)*R[i])/pow(R[i]*R[i] + bb*bb, 3./4.) + (sqrt(Md)*R[i])/pow( R[i]*R[i] + (bd+ad)*(bd+ad),3./4.)  + sqrt(Mh)/pow( R[i]*R[i] + ah*ah, 1./4.);
	}
}

double likelihood( double* y_model, double* y ){
	double ysum=0;
	for(int i = 0; i < SIZE_FILE; i++ ){
		ysum += (y_model[i] - y[i])*(y_model[i] - y[i]);
	}
	return 0.5*ysum;
}

int max_index( double *l ){
	vector<double> v(N);
	for(int i=0; i < N; i++)
		v[i] = l[i];

	vector<double>::const_iterator first = v.begin() + int(N/2);
	vector<double>::const_iterator last = v.end();
	vector<double> newVec(first, last);

	vector<double>::iterator max = max_element(newVec.begin(), newVec.end());
	return int(N/2) + std::distance( newVec.begin(), max);
}

int main(){
	double Mb_prime = 0.0, Md_prime = 0.0, Mh_prime = 0.0, l_prime = 0.0, l_init, alpha, beta;
	double R[SIZE_FILE], VR[SIZE_FILE], VR_init[SIZE_FILE], VR_prime[SIZE_FILE];
	double Mb_walk[N], Md_walk[N], Mh_walk[N], l_walk[N];
	default_random_engine generator;

	load_data( R, VR );

	Mb_walk[0] = ((double) rand()/(RAND_MAX));
	Md_walk[0] = ((double) rand()/(RAND_MAX));
	Mh_walk[0] = ((double) rand()/(RAND_MAX));

	my_model( R, VR_init, Mb_walk[0], Md_walk[0], Mh_walk[0] );
	l_walk[0] = likelihood( VR, VR_init );

	for( int i = 1; i < N; i++){
		normal_distribution<double> Ndistribution_Mb( Mb_walk[i-1], 10 );
		normal_distribution<double> Ndistribution_Md( Md_walk[i-1], 10 );
		normal_distribution<double> Ndistribution_Mh( Mh_walk[i-1], 10 );

		Mb_prime = Ndistribution_Mb(generator);
		Md_prime = Ndistribution_Md(generator);
		Mh_prime = Ndistribution_Mh(generator);

		my_model( R, VR_init, Mb_walk[i-1], Md_walk[i-1], Mh_walk[i-1] );
		my_model( R, VR_prime, Mb_prime, Md_prime, Mh_prime );

		l_prime = likelihood( VR, VR_prime);
		l_init = likelihood( VR, VR_init);

		alpha = exp(-(l_prime - l_init));
		if( alpha >= 1.0 ){
			Mb_walk[i] = Mb_prime;
			Md_walk[i] = Md_prime;
			Mh_walk[i] = Mh_prime;
			l_walk[i] = l_prime;
		} else {
			beta = ((double) rand()/(RAND_MAX));
			if( beta <= alpha ){
				Mb_walk[i] = Mb_prime;
				Md_walk[i] = Md_prime;
				Mh_walk[i] = Mh_prime;
				l_walk[i] = l_prime;
			} else {
				Mb_walk[i] = Mb_walk[i-1];
				Md_walk[i] = Md_walk[i-1];
				Mh_walk[i] = Mh_walk[i-1];
				l_walk[i] = l_init;
			}
		}
	}

	int max_idx = max_index( l_walk );
	cout << "Mb = " << Mb_walk[max_idx] << endl;
	cout << "Md = " << Md_walk[max_idx] << endl;
	cout << "Mh = " << Mh_walk[max_idx] << endl;

	ofstream OF(".parameters.dat");
	OF << Mb_walk[max_idx] << endl;
	OF << Md_walk[max_idx] << endl;
	OF << Mh_walk[max_idx] << endl;
	return 0;
}
