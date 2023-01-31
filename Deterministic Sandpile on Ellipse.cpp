#include <iostream>
#include<cmath>
#include<math.h>
#include<cstdlib>
#include <vector>
#include<fstream>
#include <ctime>

using namespace std;

void routine(int L, int max);
void ellipsematrixprintBOOLS(vector<vector<bool> > matrix, int R, int C);
void ellipsematrixprintINTS(vector<vector<int> > matrix, int R, int C);
void ellipsematrixprintDOUBLES(vector<vector<double> > matrix, int R, int C);
void Laplacian(vector<vector<int> > &L, vector<vector<int> > v, int n);
int x(int i, int half);
int y(int i, int half);
int xmax(vector<vector<int> > Corrected_O, int L, int max);
int ymax(vector<vector<int> > Corrected_O, int L, int max);
void Toppling(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<int> > Domain, int x_m, int y_m);
void BoundaryTopple(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<int> > Domain, vector<vector<int> > bf, int x_m, int y_m, int m);
bool RecurrencyTest (vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<int> > test, vector<vector<int> > Domain, vector<vector<int> > bf, int x_m, int y_m);
/*void ComponentwiseAddition(vector<vector<int> > &v3, vector<vector<int> > v1, vector<vector<int> > v2, int x_m, int y_m, int m);
void Toppling(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<bool> > Domain, int x_m, int y_m);
void BoundaryTopple(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<bool> > Domain, vector<vector<bool> > bf, int x_m, int y_m, int m);
bool RecurrencyTest (vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<int> > test, vector<vector<bool> > Domain, vector<vector<bool> > bf, int x_m, int y_m);
void Laplacian(vector<vector<int> > &L, vector<vector<int> > v, vector<vector<bool> > Domain, int x_m, int y_m);*/
void WriteMathematicaFile (vector<vector<int> > Function, int x_m, int y_m);
void WriteMathematicaFile2 (vector<vector<int> > Function, int x_m, int y_m);
void WriteMathematicaFile3 (vector<vector<int> > Function, int x_m, int y_m, int max);



int main () {

	//int currenttime;
	//cout << "Enter Current Time: " << endl;
	//cin >> currenttime;

	//cout << "ticks per sec: " << CLOCKS_PER_SEC << endl << endl;
	
	int i;

	vector<int> max (5);
	vector<int> L (5);

	max[0] = 5;
	max[1] = 4000;
	max[2] = 8000;
	max[3] = 16000;
	max[4] = 24000;

	L[0] = 19;
	L[1] = 399;
	L[2] = 599;
	L[3] = 799;
	L[4] = 999;

	for(i = 0; i < 1; i++) {
	 routine(L[i], max[i]);
	}	
	
return 0;
}

void routine (int L, int max) {
	clock_t t;
	t = clock();	

/*
	int L ;
	cout << "Input side length (must be odd number): " << endl;
	cin >> L;

	int max;
	cout << "Input level set value: " << endl;
	cin >> max;

	*/
	double a, b, c;
	int i, j;
	int half = (L - 1)/2;
	
	//n = 10;
	a = 1.25;
	b = 0.5;
	c = 1;
	
	//BUILD DOMAIN VECTOR

	vector<vector<int> > OtherDomain(L);

	
	for (i = 0 ; i < L ; i++ ) {
   		OtherDomain[i].resize(L);
	}

		
		

	int xx, yy;
	for (i = 0; i < L; i++) {
			xx = i - half;
		for (j = 0; j < L; j++) {
			yy = j - half;
		
			
			if ( a*xx*xx + 2*b*xx*yy + c*yy*yy < 2*max){
				OtherDomain[i][j] = 1;
	
			} else {
				OtherDomain[i][j] = 0;
			}
	
		}
	}

		cout << "get here?" << endl;
	


	/*

	cout << "Input matrix A" << endl
	<< "a:" << endl;
	cin >> a;

	cout << "b:" << endl;
	cin >> b;

	cout << "c:" << endl;
	cin >> c;*/

	

	

	vector<vector<int> > Uncorrected_O(L);
	vector<vector<int> > Corrected_O(L);

	for (i = 0 ; i < L ; i++) {
   		Uncorrected_O[i].resize(L);
	}

	for (i = 0 ; i < L ; i++) {
   		Corrected_O[i].resize(L);
	}


	/*Uncorrected_O[half][half] = 0;
	Uncorrected_O[half+1][half] = 0;	
	Uncorrected_O[half+2][half] = 0;
	Uncorrected_O[half+3][half] = 0;	
	Uncorrected_O[half+4][half] = 0;	
	Uncorrected_O[half+5][half] = 0;	
	Uncorrected_O[half+6][half] = 0;	
	Uncorrected_O[half+7][half] = 0;	
	Uncorrected_O[half+8][half] = 0;
	
	Uncorrected_O[half][half+1] = 0;
	Uncorrected_O[half+1][half+1] = 1;	
	Uncorrected_O[half+2][half+1] = 2;	
	Uncorrected_O[half+3][half+1] = 3;	
	Uncorrected_O[half+4][half+1] = 4;	
	Uncorrected_O[half+5][half+1] = 5;	
	Uncorrected_O[half+6][half+1] = 6;	
	Uncorrected_O[half+7][half+1] = 7;	
	Uncorrected_O[half+8][half+1] = 7;	

	Uncorrected_O[half][half+2] = 0;	
	Uncorrected_O[half+1][half+2] = 2;	
	Uncorrected_O[half+2][half+2] = 4;
	Uncorrected_O[half+3][half+2] = 6;
	Uncorrected_O[half+4][half+2] = 8;
	Uncorrected_O[half+5][half+2] = 10;
	Uncorrected_O[half+6][half+2] = 12;
	Uncorrected_O[half+7][half+2] = 13;
	Uncorrected_O[half+8][half+2] = 14;

	Uncorrected_O[half][half+3] = 0;
	Uncorrected_O[half+1][half+3] = 3;
	Uncorrected_O[half+2][half+3] = 6;
	Uncorrected_O[half+3][half+3] = 9;
	Uncorrected_O[half+4][half+3] = 12;
	Uncorrected_O[half+5][half+3] = 15;
	Uncorrected_O[half+6][half+3] = 17;
	Uncorrected_O[half+7][half+3] = 19;
	Uncorrected_O[half+8][half+3] = 21;	

	Uncorrected_O[half][half+4] = 0;
	Uncorrected_O[half+1][half+4] = 4;
	Uncorrected_O[half+2][half+4] = 8;
	Uncorrected_O[half+3][half+4] = 12;
	Uncorrected_O[half+4][half+4] = 16;
	Uncorrected_O[half+5][half+4] = 19;
	Uncorrected_O[half+6][half+4] = 22;
	Uncorrected_O[half+7][half+4] = 25;
	Uncorrected_O[half+8][half+4] = 28;

	Uncorrected_O[half][half+5] = 0;
	Uncorrected_O[half+1][half+5] = 5;
	Uncorrected_O[half+2][half+5] = 10;
	Uncorrected_O[half+3][half+5] = 15;
	Uncorrected_O[half+4][half+5] = 19;
	Uncorrected_O[half+5][half+5] = 23;
	Uncorrected_O[half+6][half+5] = 27;
	Uncorrected_O[half+7][half+5] = 31;
	Uncorrected_O[half+8][half+5] = 35;

	Uncorrected_O[half][half+6] = 0;
	Uncorrected_O[half+1][half+6] = 6;
	Uncorrected_O[half+2][half+6] = 12;
	Uncorrected_O[half+3][half+6] = 17;
	Uncorrected_O[half+4][half+6] = 22;
	Uncorrected_O[half+5][half+6] = 27;
	Uncorrected_O[half+6][half+6] = 32;
	Uncorrected_O[half+7][half+6] = 37;
	Uncorrected_O[half+8][half+6] = 42;

	Uncorrected_O[half][half+7] = 0;
	Uncorrected_O[half+1][half+7] = 7;
	Uncorrected_O[half+2][half+7] = 13;
	Uncorrected_O[half+3][half+7] = 19;
	Uncorrected_O[half+4][half+7] = 25;
	Uncorrected_O[half+5][half+7] = 31;
	Uncorrected_O[half+6][half+7] = 37;
	Uncorrected_O[half+7][half+7] = 43;
	Uncorrected_O[half+8][half+7] = 49;

	Uncorrected_O[half][half+8] = 0;
	Uncorrected_O[half+1][half+8] = 7;
	Uncorrected_O[half+2][half+8] = 14;
	Uncorrected_O[half+3][half+8] = 21;
	Uncorrected_O[half+4][half+8] = 28;
	Uncorrected_O[half+5][half+8] = 35;
	Uncorrected_O[half+6][half+8] = 42;
	Uncorrected_O[half+7][half+8] = 49;
	Uncorrected_O[half+8][half+8] = 56;*/

	Uncorrected_O[half][half] = 0;
	Uncorrected_O[half+1][half] = 0;	
	Uncorrected_O[half+2][half] = 0;
	Uncorrected_O[half+3][half] = 0;	
	Uncorrected_O[half+4][half] = 0;	
	Uncorrected_O[half+5][half] = 0;	
	Uncorrected_O[half+6][half] = 0;	
	Uncorrected_O[half+7][half] = 0;	
	Uncorrected_O[half+8][half] = 0;
	
	Uncorrected_O[half][half+1] = 0;
	Uncorrected_O[half+1][half+1] = 1;	
	Uncorrected_O[half+2][half+1] = 1;	
	Uncorrected_O[half+3][half+1] = 1;	
	Uncorrected_O[half+4][half+1] = 1;	
	Uncorrected_O[half+5][half+1] = 1;	
	Uncorrected_O[half+6][half+1] = 1;	
	Uncorrected_O[half+7][half+1] = 1;	
	Uncorrected_O[half+8][half+1] = 1;	

	Uncorrected_O[half][half+2] = 0;	
	Uncorrected_O[half+1][half+2] = 1;	
	Uncorrected_O[half+2][half+2] = 4;
	Uncorrected_O[half+3][half+2] = 6;
	Uncorrected_O[half+4][half+2] = 8;
	Uncorrected_O[half+5][half+2] = 10;
	Uncorrected_O[half+6][half+2] = 12;
	Uncorrected_O[half+7][half+2] = 2;
	Uncorrected_O[half+8][half+2] = 2;

	Uncorrected_O[half][half+3] = 0;
	Uncorrected_O[half+1][half+3] = 1;
	Uncorrected_O[half+2][half+3] = 6;
	Uncorrected_O[half+3][half+3] = 9;
	Uncorrected_O[half+4][half+3] = 12;
	Uncorrected_O[half+5][half+3] = 15;
	Uncorrected_O[half+6][half+3] = 17;
	Uncorrected_O[half+7][half+3] = 3;
	Uncorrected_O[half+8][half+3] = 3;	

	Uncorrected_O[half][half+4] = 0;
	Uncorrected_O[half+1][half+4] = 1;
	Uncorrected_O[half+2][half+4] = 8;
	Uncorrected_O[half+3][half+4] = 12;
	Uncorrected_O[half+4][half+4] = 16;
	Uncorrected_O[half+5][half+4] = 19;
	Uncorrected_O[half+6][half+4] = 22;
	Uncorrected_O[half+7][half+4] = 4;
	Uncorrected_O[half+8][half+4] = 4;

	Uncorrected_O[half][half+5] = 0;
	Uncorrected_O[half+1][half+5] = 1;
	Uncorrected_O[half+2][half+5] = 10;
	Uncorrected_O[half+3][half+5] = 15;
	Uncorrected_O[half+4][half+5] = 19;
	Uncorrected_O[half+5][half+5] = 23;
	Uncorrected_O[half+6][half+5] = 27;
	Uncorrected_O[half+7][half+5] = 5;
	Uncorrected_O[half+8][half+5] = 5;

	Uncorrected_O[half][half+6] = 0;
	Uncorrected_O[half+1][half+6] = 1;
	Uncorrected_O[half+2][half+6] = 12;
	Uncorrected_O[half+3][half+6] = 17;
	Uncorrected_O[half+4][half+6] = 22;
	Uncorrected_O[half+5][half+6] = 27;
	Uncorrected_O[half+6][half+6] = 32;
	Uncorrected_O[half+7][half+6] = 6;
	Uncorrected_O[half+8][half+6] = 6;

	Uncorrected_O[half][half+7] = 0;
	Uncorrected_O[half+1][half+7] = 1;
	Uncorrected_O[half+2][half+7] = 2;
	Uncorrected_O[half+3][half+7] = 3;
	Uncorrected_O[half+4][half+7] = 4;
	Uncorrected_O[half+5][half+7] = 5;
	Uncorrected_O[half+6][half+7] = 6;
	Uncorrected_O[half+7][half+7] = 7;
	Uncorrected_O[half+8][half+7] = 7;

	Uncorrected_O[half][half+8] = 0;
	Uncorrected_O[half+1][half+8] = 1;
	Uncorrected_O[half+2][half+8] = 2;
	Uncorrected_O[half+3][half+8] = 3;
	Uncorrected_O[half+4][half+8] = 4;
	Uncorrected_O[half+5][half+8] = 5;
	Uncorrected_O[half+6][half+8] = 6;
	Uncorrected_O[half+7][half+8] = 7;
	Uncorrected_O[half+8][half+8] = 8;
	


	for (i = half; i < L-8; i++) {
		Uncorrected_O[i+8][y(i, half)+1+half] = Uncorrected_O[i][y(i, half)+half] + x(i, half) + 7*y(i, half) + 7;
		Uncorrected_O[i+8][y(i, half)+2+half] = Uncorrected_O[i][y(i, half)+1+half] + x(i, half) + 7*(y(i, half)+1) + 7;
		Uncorrected_O[i+8][y(i, half)+3+half] = Uncorrected_O[i][y(i, half)+2+half] + x(i, half) + 7*(y(i, half)+2) + 7;
		Uncorrected_O[i+8][y(i, half)+4+half] = Uncorrected_O[i][y(i, half)+3+half] + x(i, half) + 7*(y(i, half)+3) + 7;
		Uncorrected_O[i+8][y(i, half)+5+half] = Uncorrected_O[i][y(i, half)+4+half] + x(i, half) + 7*(y(i, half)+4) + 7;
		Uncorrected_O[i+8][y(i, half)+6+half] = Uncorrected_O[i][y(i, half)+5+half] + x(i, half) + 7*(y(i, half)+5) + 7;
		Uncorrected_O[i+8][y(i, half)+7+half] = Uncorrected_O[i][y(i, half)+6+half] + x(i, half) + 7*(y(i, half)+6) + 7;
		Uncorrected_O[i+8][y(i, half)+8+half] = Uncorrected_O[i][y(i, half)+7+half] + x(i, half) + 7*(y(i, half)+7) + 7;
		Uncorrected_O[i+8][y(i, half)+9+half] = Uncorrected_O[i][y(i, half)+8+half] + x(i, half) + 7*(y(i, half)+8) + 7;
	}
	

	for (i = half+7; i > 7; i--){
		Uncorrected_O[i-8][y(i-8, half)+half] = Uncorrected_O[i][y(i-8, half)+half+1] - (x(i-8,half) + 7*y(i-8, half) + 7);
		Uncorrected_O[i-8][y(i-8, half)+1+half] = Uncorrected_O[i][y(i-8, half)+2+half] - (x(i-8,half) + 7*(y(i-8, half)+1)+7);	
		Uncorrected_O[i-8][y(i-8, half)+2+half] = Uncorrected_O[i][y(i-8, half)+half+3] - (x(i-8,half) + 7*(y(i-8, half)+2) + 7);
		Uncorrected_O[i-8][y(i-8, half)+3+half] = Uncorrected_O[i][y(i-8, half)+half+4] - (x(i-8,half) + 7*(y(i-8, half)+3) + 7);
		Uncorrected_O[i-8][y(i-8, half)+4+half] = Uncorrected_O[i][y(i-8, half)+half+5] - (x(i-8,half) + 7*(y(i-8, half)+4) + 7);
		Uncorrected_O[i-8][y(i-8, half)+5+half] = Uncorrected_O[i][y(i-8, half)+half+6] - (x(i-8,half) + 7*(y(i-8, half)+5) + 7);
		Uncorrected_O[i-8][y(i-8, half)+6+half] = Uncorrected_O[i][y(i-8, half)+half+7] - (x(i-8,half) + 7*(y(i-8, half)+6) + 7);
		Uncorrected_O[i-8][y(i-8, half)+7+half] = Uncorrected_O[i][y(i-8, half)+half+8] - (x(i-8,half) + 7*(y(i-8, half)+7) + 7);
		Uncorrected_O[i-8][y(i-8, half)+8+half] = Uncorrected_O[i][y(i-8, half)+half+9] - (x(i-8,half) + 7*(y(i-8, half)+8) + 7);
	}
		
	

	for (i = 0; i < L; i++) {
		j = y(i, half) + half;
		while (j + 8 < L) {
			Uncorrected_O[i][j+8] = Uncorrected_O[i][j] + 7*x(i, half);
			j++;
		}

		j = y(i, half) + 7 + half;
		while (j - 8 >= 0){
			Uncorrected_O[i][j-8] = Uncorrected_O[i][j] - 7*x(i, half);
			j--;
		}
	}



	vector<vector<int> > Lap_O(L);

	for (i = 0 ; i < L ; i++) {
   		Lap_O[i].resize(L);
	}



	
	

	for(i = 0; i < L; i++){
		for (j = 0; j< L; j++){
			Corrected_O[i][j] = Uncorrected_O[i][j] + (x(i,half) - (j-half) + 1)*(x(i,half) + (j-half))/2;
			//Corrected_O[i][j] = Uncorrected_O[i][j] + x(i,half)*(x(i,half) + 1);
			//Corrected_O[i][j] = Uncorrected_O[i][j] + (j-half)*((j-half)+1)/2+x(i,half)*(x(i,half) + 1)/2;
			
		}
	}	

	vector<vector<double> > LinTermsFinder(L);

	for (i = 0 ; i < L ; i++) {
   		LinTermsFinder[i].resize(L);
	}

	double f, g;

	for(i = 0; i < L; i++){
		f = i- half;
		for(j=0; j< L; j++){
			g = j - half;
			//LinTermsFinder[i][j] = (((65.0/64.0)*f*f+2*(7.0/8.0)*f*g + g*g)/2.0) + (0.4375*f + 0.5*g) + (Uncorrected_O[i][j] - (((1.0/64.0)*f*f+2*(7.0/8.0)*f*g)/2) + f/16.0);
			//LinTermsFinder[i][j] =  (Uncorrected_O[i][j] - (((1.0/3.0)*f*f+2*(7.0/8.0)*f*g + 1.0/4.0*g*g)/2) + f/16.0);
			LinTermsFinder[i][j] = Uncorrected_O[i][j] - ((1/64.0)*f*f+2*(7.0/8.0)*f*g)/2.0;
		}
	}

	Laplacian(Lap_O, Uncorrected_O, L);
	//WriteMathematicaFile (Lap_O, half, half);
	cout << "here's LapO" << endl;
	ellipsematrixprintINTS(Lap_O, L, L);
	//ellipsematrixprintINTS(Uncorrected_O, L, L);

	//cout << "here's p" << endl;
	//ellipsematrixprintDOUBLES(LinTermsFinder, L, L);



	vector<vector<int> > Domain(L);

	for (i = 0 ; i < L ; i++) {
   		Domain[i].resize(L);
	}

	for(i = 0; i < L; i++){
		for(j=0; j< L; j++){
			if (LinTermsFinder[i][j] < max) {
				Domain[i][j] = 1;
			} else {
				Domain[i][j] = 0;
			}
		}
	}


	bool DomainProblem = false;

	for (i = 0; i < L; i++) {
		if (Domain[i][0] == 1 || Domain[i][L-1] == 1 || Domain[0][i] == 1 || Domain[L-1][i] == 1) {
			DomainProblem = true;
			cout << "Domain Problem: Level set exceeds allotted arrayspace" << endl;
		}
	}

	vector<vector<int> > bf(L);
	for (i = 0 ; i < L ; i++ ) {
   		bf[i].resize(L);
	}

	for(i = 0; i < L; i++) {
		for (j = 0; j < L; j++){
			bf[j][i] = 0;
			if (Domain[j][i] == 0) {
				if (i < L-1 && Domain[j][i+1] == 1){
					bf[j][i] = 1;
				} else if (i > 0 && Domain[j][i-1] == 1){
					bf[j][i] = 1;
				} else if (j < L-1 && Domain[j+1][i] == 1){
					bf[j][i] = 1;
				} else if (j > 0 && Domain[j-1][i] == 1){
					bf[j][i] = 1;
				}
			}
		}
	}	

	

	/*double minboundaryvalue;
	bool minboundaryvaluefirsthit = false;
	
	for (i = 0; i < L; i++) {
		for (j = 0; j < L; j++) {
			if(minboundaryvaluefirsthit == true && minboundaryvalue > LinTermsFinder[i][j] && bf[i][j] == 1) {
				minboundaryvalue = LinTermsFinder[i][j];
			}
			if (minboundaryvaluefirsthit == false && bf[i][j] == 1) {
				minboundaryvalue = LinTermsFinder[i][j];
				minboundaryvaluefirsthit = true;
			}
		}
	}


	for (i = 0; i < L; i++) {
		for (j = 0; j < L; j++) {
			if (Domain[i][j] == 0 && bf[i][j] == 0) {
				LinTermsFinder[i][j] = 0;
			} else {
				LinTermsFinder[i][j] = LinTermsFinder[i][j] - max;
			}
		}
	}

	

	//cout << "Odometer is: " << endl;
	//ellipsematrixprintDOUBLES(LinTermsFinder, L, L);

//INITIALIZE ODOMETER
	
	vector<vector<int> > o(L);
	for (i = 0 ; i < L; i++ ) {
   		o[i].resize(L);
	}

	

	
//INITIALIZE SANDPILE
	
	vector<vector<int> > s(L);
	for (i = 0 ; i < L; i++ ) {
   		s[i].resize(L);
	}

int m =  16*L*L;
	//cout << "going into boundary topple" << endl;
	BoundaryTopple(s, o, Domain, bf, half, half, m);
	//cout << "done with bt, going into topple" << endl;
//TOPPLE INTERIOR UNTIL STABILITY
	Toppling(s, o, Domain, half, half);
	//cout << "done with topple" << endl;	
	//CHECK IF RECURRENT. IF NOT, TOPPLE THE BOUNDARY ONCE. CHECK AGAIN. CONTINUE UNTIL RECURRENT
	bool RecurrencyCheck = false;


	vector<vector<int> > test(L);
	
	for (i = 0 ; i < L ; i++ ) {
   		test[i].resize(L);
	}

	
	while(RecurrencyCheck == false) {
		RecurrencyCheck = RecurrencyTest(s, o, test, Domain, bf, half, half);
	} 
	
	
	int boundaryvalue;
	for (i = 0; i < L; i++) {
		for (j = 0; j < L; j++) {
			if (bf[i][j] == 1){
				boundaryvalue = o[i][j];
				o[i][j] = 0;
			} else if(Domain[i][j] == 1){
				o[i][j] = o[i][j] - boundaryvalue;
			}
		}
	}
		



	//Laplacian(Lap_O, Uncorrected_O, L);
	//cout << "v is: " << endl;
	//ellipsematrixprintINTS(o, L, L);
	//ellipsematrixprintINTS(Uncorrected_O, L, L);
	//ellipsematrixprintINTS(Corrected_O, L, L);
	//ellipsematrixprintINTS(Lap_O, L, L);
	//ellipsematrixprintINTS(s, L, L);
	//cout << "Domain is: " << endl;
	//ellipsematrixprintINTS(Domain, L, L);


	for (i = 0; i < L; i++) {
		for (j = 0; j < L; j++) {
			test[i][j] = Domain[i][j]-OtherDomain[i][j];
		}
	}

	//ellipsematrixprintINTS(bf, L, L);


	


	//cout << "the pair is: " << xmax(Corrected_O, L, max) << ", " << ymax(Corrected_O, L, max);	
	WriteMathematicaFile (Lap_O, half, half);
	WriteMathematicaFile2 (test, half, half); 
	WriteMathematicaFile3(bf, half, half, max);


	
	bool InequalityProblem = false;
	int maxdiff = 0;

	
	
		
	int xmaxdiff, ymaxdiff;
	
	for (i = 0; i < L; i++) {
		for (j = 0; j < L; j++) {
			if(LinTermsFinder[i][j] - o[i][j] > maxdiff && Domain[i][j] == 1) {
				maxdiff = LinTermsFinder[i][j] - o[i][j];
				xmaxdiff = i - half;
				ymaxdiff = j - half;
			}

			if (o[i][j] >= 1 + LinTermsFinder[i][j]) {
				InequalityProblem = true;
				cout << "InequalityProblem: BVP solution is greater than odometer. This is bad." << endl;
			}
		}
	}

	cout << "L is: " << L << endl;
	cout << "max is: " << max << endl;
	cout << "max difference is: " << maxdiff << endl;
	cout << "and it occured at: " << xmaxdiff << ", " << ymaxdiff << endl << endl;




	t = clock() - t;
	cout << "time to finish is " << t << endl << endl;
	*/
}



void Toppling(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<int> > Domain, int x_m, int y_m){
	int i, j;
	bool clear = false;
	while (clear == false) {
		for(i = 0; i < 2*x_m + 1; i++) {
			if (i == 0) {
				clear = true;
			}

			for (j = 0; j < 2*y_m + 1; j++){
				if (s[j][i] >= 4 && Domain[j][i] == 1){
					o[j][i]++;
					clear = false;
					s[j][i] = s[j][i] - 4;
					s[j+1][i]++;
					s[j-1][i]++;
					s[j][i+1]++;		
					s[j][i-1]++;
				}
			}
		}
	}
}		

void BoundaryTopple(vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<int> > Domain, vector<vector<int> > bf,  int x_m, int y_m, int m){
	
	int i,j;
	for (i = 0; i < 2*x_m + 1; i++) {
		for (j = 0; j < 2*y_m + 1; j++) {
			if (bf[j][i] == 1){
					
				o[j][i] = o[j][i] + m;				

				if (i < 2*x_m  && Domain[j][i+1] == 1) {
					s[j][i+1] = s[j][i+1] + m;
				}
				if (i > 0 && Domain[j][i-1] == 1) {
					s[j][i-1] = s[j][i-1] + m;
				}
				if (j < 2*y_m  && Domain[j+1][i] == 1) {
					s[j+1][i] = s[j+1][i] + m;
				}
				if (j > 0 && Domain[j-1][i] == 1) {
					s[j-1][i] = s[j-1][i] + m;
				}
			}
		}
	}
}



bool RecurrencyTest (vector<vector<int> > &s, vector<vector<int> > &o, vector<vector<int> > test, vector<vector<int> > Domain, vector<vector<int> > bf, int x_m, int y_m){
	
	//SAVE THE ORIGINAL SANDPILE
	
	int i, j;
	

	for (i = 0; i < 2*x_m + 1; i++) {
		for (j = 0; j < 2*y_m + 1; j++) {
			test[j][i] = s[j][i];
		}
	}
	
	
	//TOPPLE THE BOUNDARY, THEN EVERYTHING ELSE, TO GET A 'NEW' SANDPILE
	BoundaryTopple(s, o, Domain, bf, x_m, y_m, 1); 
	Toppling(s, o, Domain, x_m, y_m);

	//ellipsematrixprintINTS(s, 2*y_m +1, 2*x_m+1);
	//ellipsematrixprintINTS(test, 2*y_m +1, 2*x_m+1);

//cout << "hello" << endl;
	//NOW CHECK IF THIS NEW SANDPILE IS THE SAME AS THE ORIGINAL
	bool clear = true;
	for (i = 0; i < 2*x_m+1; i++) {
		for (j = 0; j < 2*y_m+1; j++) {
			if (Domain[j][i] == 1 && s[j][i] != test[j][i]){
				clear = false;
			}
		}
	}
	return clear; 

}	



int x(int i, int half) {
	int x = i - half;
	return x;
}

int y(int i, int half){
	int temp, y;
	temp = x(i, half);
	if (temp >= 0) {
		y = (temp/8);
	} else {
		y = ((temp + 1)/8) - 1;
	}
	return y;
}



void ellipsematrixprintINTS(vector<vector<int> > matrix, int R, int C){
		int matxpnt1;
		int matxpnt2;
		for(matxpnt1 = R-1; matxpnt1 >= 0; matxpnt1--){
			for(matxpnt2=0; matxpnt2 < C-1; matxpnt2++){
				cout << matrix[matxpnt2][matxpnt1];
				if (0 <= matrix[matxpnt2][matxpnt1] && matrix[matxpnt2][matxpnt1] <= 9) {
					cout << "     ";
				} else if ((10 <= matrix[matxpnt2][matxpnt1] && matrix[matxpnt2][matxpnt1] <= 99)||(0 >= matrix[matxpnt2][matxpnt1] && matrix[matxpnt2][matxpnt1] >= -9)) {
					cout << "    ";
				} else if (matrix[matxpnt2][matxpnt1] >= 100 || matrix[matxpnt2][matxpnt1]  <= -10 ){
					cout << "   ";
				}
			}
			cout << matrix[C-1][matxpnt1] <<endl;
			
		}
		cout << endl << endl;
}

void ellipsematrixprintBOOLS(vector<vector<bool> > matrix, int R, int C){
		int matxpnt1;
		int matxpnt2;
		for(matxpnt1 = R-1; matxpnt1 >= 0; matxpnt1--){
			for(matxpnt2=0; matxpnt2 < C-1; matxpnt2++){
				cout << matrix[matxpnt2][matxpnt1] << "                  ";
			}
			cout << matrix[C-1][matxpnt1] <<endl;
			
		}
		cout << endl << endl;
}

void ellipsematrixprintDOUBLES(vector<vector<double> > matrix, int R, int C){
		int matxpnt1;
		int matxpnt2;
		for(matxpnt1 = R-1; matxpnt1 >= 0; matxpnt1--){
			for(matxpnt2=0; matxpnt2 < C-1; matxpnt2++){
				cout << matrix[matxpnt2][matxpnt1];
				if (0 <= matrix[matxpnt2][matxpnt1] && matrix[matxpnt2][matxpnt1] <= 9) {
					cout << "     ";
				} else if ((10 <= matrix[matxpnt2][matxpnt1] && matrix[matxpnt2][matxpnt1] <= 99)||(0 >= matrix[matxpnt2][matxpnt1] && matrix[matxpnt2][matxpnt1] >= -9)) {
					cout << "    ";
				} else if (matrix[matxpnt2][matxpnt1] >= 100 || matrix[matxpnt2][matxpnt1]  <= -10 ){
					cout << "   ";
				}
			}
			cout << matrix[C-1][matxpnt1] <<endl;
			
		}
		cout << endl << endl;
}

void Laplacian(vector<vector<int> > &L, vector<vector<int> > v, int n) {
	
	int i,j;
	
	for (i = 1; i < n-1; i++) {
		for (j = 1; j < n-1; j++) {
				L[j][i] = -4*v[j][i] + v[j+1][i] + v[j-1][i] + v[j][i+1] + v[j][i-1];
		}
	}
}


int xmax(vector<vector<int> > Corrected_O, int L, int max) {
	int i,j;
	bool find = false;
	int xm;
	for (i = L-1; i >= 0; i--){
		for(j = 0; j < L; j++){
			if (Corrected_O[i][j] < max) {
				find = true;
				xm = x(i, (L-1)/2);
				break;
			}
		}
		if (find == true) {
			break;
		}
	}

	return xm;
}

int ymax(vector<vector<int> > Corrected_O, int L, int max){
	int i,j;
	bool find = false;
	int ym;
	for (j = L-1; j >= 0; j--){
		for(i = 0; i < L; i++){
			if (Corrected_O[i][j] < max) {
				find = true;
				ym = j - (L-1)/2;
				break;
			}
		}
		if (find == true) {
			break;
		}
	}

	return ym;
}
	


void WriteMathematicaFile (vector<vector<int> > Function, int x_m, int y_m){
	int i, j;
	
	ofstream outfile;
	outfile.open("Odometer.nb", ios::out);
	outfile << "CustomColorFunc[z_] := RGBColor[KroneckerDelta[z, 2-2] + KroneckerDelta[z, 3-2] + 0.75*KroneckerDelta[z, -1-2], KroneckerDelta[z, 1-2] + KroneckerDelta[z, 2-2] + 0.75*KroneckerDelta[z, -1-2], KroneckerDelta[z, 0-2] + KroneckerDelta[z, 1-2] + 0.75*KroneckerDelta[z, -1-2]]";

	outfile << endl;

	outfile << "Sandpile = ";

	for (j = 2*y_m; j >= 0; j--) {
		if (j == 2*y_m) {
			outfile << "{{";
		} else {
			outfile << "{";
		}

		for (i = 0; i < 2*x_m + 1; i++) {
			outfile << Function[i][j];
			if (i != 2*x_m){
				outfile << ",";
			} else if (j != 0){
				outfile << "},";
			} else if (j == 0){
				outfile << "}};";
			}
		}
	}


	outfile << endl;

	outfile << "ArrayPlot[Sandpile, ColorFunctionScaling -> False,  ColorFunction -> CustomColorFunc]";
	
	outfile.close();
}


void WriteMathematicaFile2 (vector<vector<int> > Function, int x_m, int y_m){
	int i, j;
	
	ofstream outfile;
	outfile.open("testDomainEllipseFinder2.nb", ios::out);
	outfile << "CustomColorFunc[z_] := RGBColor[KroneckerDelta[z, 2] + KroneckerDelta[z, 3] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 1] + KroneckerDelta[z, 2] + 0.75*KroneckerDelta[z, -1], KroneckerDelta[z, 0] + KroneckerDelta[z, 1] + 0.75*KroneckerDelta[z, -1]]";

	outfile << endl;

	outfile << "Sandpile = ";

	for (j = 0; j < 2*y_m+1; j++) {
		if (j == 0) {
			outfile << "{{";
		} else {
			outfile << "{";
		}

		for (i = 0; i < 2*x_m + 1; i++) {
			outfile << Function[j][i];
			if (i != 2*x_m){
				outfile << ",";
			} else if (j != 2*y_m){
				outfile << "},";
			} else if (j == 2*y_m){
				outfile << "}};";
			}
		}
	}

	outfile << endl;

	outfile << "ArrayPlot[Sandpile, ColorFunctionScaling -> False,  ColorFunction -> CustomColorFunc]";
	
	outfile.close();
}

void WriteMathematicaFile3 (vector<vector<int> > Function, int x_m, int y_m, int max){
	int i, j;
	
	ofstream outfile;
	outfile.open("BoundaryEllipseFinder2.nb", ios::out);
	
	outfile << "OuterBoundary = {";

	bool hit = false;

	for (i = 0; i < 2*x_m + 1; i++){
		for (j = 2*y_m ; j >= y_m; j--) {
			
			if (hit == true && Function[i][j] == 0) {
				break;
			}
						

			if ( Function[i][j] == 1 && !(j == y_m)) {
				outfile << "{" << i - x_m << ", " << j - y_m << "}, ";
				hit = true;
			}
			else if ( Function[i][j] == 1 && j == y_m){
				outfile << "{" << i - x_m << ", " << j - y_m << "}}";
				hit = true;
			}
			
			
		}
		hit = false;
	}


	outfile << endl;

	//outfile << "FindFit[OuterBoundary, {-b*x + Sqrt[(b^2 - a)*x^2 + 2*n^2]}, {a, b, n}, x]";
	outfile << "FindFit[OuterBoundary,  {-b/c*x + Sqrt[((b/c)^2 - a/c)*x^2 + 2*(" << max << ")/c]}, {a, b,c},x]";
	
	
	outfile.close();
}

