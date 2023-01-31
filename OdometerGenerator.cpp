#include <stack>
#include <vector>
#include <iostream>
#include<cmath>
#include<cstdlib>
#include <fstream>

using namespace std;

void matrixprint(vector<vector<int> > matrix, int R, int C);
int strongcounter (vector<bool> &istravatrightlevel, vector<vector<int> > strong, int level, int l, int travs);
void ChipCounter(vector<vector<int> > Function, vector<bool> istravatrightlevel, vector<vector<int> > endpoint, vector<vector<int> > &chips, int travs, int n, int l);
void WriteMathematicaFile (vector<vector<int> > &chips, int travs, int n);


int main () {

		
	//NODES, EDGES 
	int n, l, travs, level;
	/*n = 7;
	l = 17;
	travs = 2344;
	level = 0;

	n = 3;
	l = 3;
	travs = 3;
	level = 0;

	n = 4;
	l = 8;
	travs = 36;
	level = 0;*/

	n = 5;
	l = 12;
	travs = 96;
	level = 0;
	

	//vector<int> nodes(n);
	//vector<int> edges(l);

	int i,j,k,m,r;
	int q = 0;
	int u;

	/*for (i=0; i < n; i++){
		nodes[i] = i;
	}

	for (i=0; i < l; i++){
		edges[i] = i + n;
	}*/	

	//DEFINE ADJ MATRIX ("adj") AND ENDPOINTS ("endpoint")
	//vector<vector<int> > adj(n);
	vector<vector<int> > endpoint(l);
	/*for (i = 0 ; i < n; i++) {
   		adj[i].resize(n);
	}*/
	for (i = 0 ; i < l; i++) {
   		endpoint[i].resize(2);
	}
	
	
	/*endpoint[0][0] = 0;
	endpoint[0][1] = 1;
	endpoint[1][0] = 0;
	endpoint[1][1] = 2;
	endpoint[2][0] = 1;
	endpoint[2][1] = 2;*/


	endpoint[0][0] = 0;
	endpoint[0][1] = 1;
	endpoint[1][0] = 0;
	endpoint[1][1] = 1;

	endpoint[2][0] = 0;
	endpoint[2][1] = 2;
	endpoint[3][0] = 0;
	endpoint[3][1] = 2;

	endpoint[4][0] = 0;
	endpoint[4][1] = 3;
	endpoint[5][0] = 0;
	endpoint[5][1] = 3;

	endpoint[6][0] = 0;
	endpoint[6][1] = 4;
	endpoint[7][0] = 0;
	endpoint[7][1] = 4;

	endpoint[8][0] = 1;
	endpoint[8][1] = 2;
	endpoint[9][0] = 2;
	endpoint[9][1] = 4;
	endpoint[10][0] = 3;
	endpoint[10][1] = 4;
	endpoint[11][0] = 1;
	endpoint[11][1] = 3;

	/*endpoint[10][0] = 1;
	endpoint[10][1] = 2;
	endpoint[11][0] = 2;
	endpoint[11][1] = 4;
	endpoint[12][0] = 3;
	endpoint[12][1] = 4;
	endpoint[13][0] = 1;
	endpoint[13][1] = 3;
	endpoint[14][0] = 4;
	endpoint[14][1] = 6;
	endpoint[15][0] = 5;
	endpoint[15][1] = 6;
	endpoint[16][0] = 3;
	endpoint[16][1] = 5;

	endpoint[0][0] = 0;
	endpoint[0][1] = 2;
	endpoint[1][0] = 0;
	endpoint[1][1] = 2;
	endpoint[2][0] = 1;
	endpoint[2][1] = 3;

	endpoint[3][0] = 1;
	endpoint[3][1] = 3;
	endpoint[4][0] = 2;
	endpoint[4][1] = 3;
	endpoint[5][0] = 0;
	endpoint[5][1] = 1;

	endpoint[6][0] = 0;
	endpoint[6][1] = 3;
	endpoint[7][0] = 1;
	endpoint[7][1] = 2;
	endpoint[8][0] = 0;
	endpoint[8][1] = 3;
	
	endpoint[9][0] = 0;
	endpoint[9][1] = 4;
	endpoint[10][0] = 0;
	endpoint[10][1] = 4;
	endpoint[11][0] = 0;
	endpoint[11][1] = 4;

	endpoint[12][0] = 5;
	endpoint[12][1] = 1;
	endpoint[13][0] = 5;
	endpoint[13][1] = 2;
	endpoint[14][0] = 5;
	endpoint[14][1] = 3;
	endpoint[15][0] = 5;
	endpoint[15][1] = 4;

*/




	vector<int> degremaining(n,0);
	
	












	stack<int> element;
	stack<int> place;

	//WRITE THE TRAVERSALS
	vector<vector<int> > traversal(n+l);
	for (i = 0 ; i < n+l; i++) {
   		traversal[i].resize(l);
	}

	//DEFINE BOOL MATRIX
	vector<vector<int> > alreadyin(n+l);
	for (i = 0 ; i < n+l; i++) {
   		alreadyin[i].resize(l);
	}

	for (i = 0; i < n+l; i++){
		for (j = 0; j < l; j++){
			alreadyin[i][j] = 0;
		}
	}

	for (i = 0; i < l; i++){
		alreadyin[0][i] = 1;
	}

	vector<vector<int> > strong(l);
	for(i = 0; i < l; i++){
		strong[i].resize(l);
	}

	for (i = 0; i < l; i++){
		for (j = 0; j < l; j++){
			strong[i][j] = 0;
		}
	}

	bool highestedgeadded = false;
	



//START THE LOOP OVER ALL TRAVERSALS (i)
	for (i=0; i < travs; i++){
		cout << "i is: " << i << endl;
		traversal[0][i] = 0;
		j = 1;
		u = 0;
//INITIALIZE THE DEGREMAINING VECTOR BY DEGREE
		for (r = 0; r < n; r++) {
			for (k = 0; k < l; k++){
				if (endpoint[k][0] == r) {
					degremaining[r]++;
				}
				if (endpoint[k][1] == r) {
					degremaining[r]++;
				}
			}
		}

		

//START THE LOOP OVER ALL ELEMENTS (j)
		while (j < n+l) {

//ADD THE FIRST ELEMENTS ACCORDING TO WHAT'S WAITING ON THE QUEUE
			if (j == 1 && element.empty() == false) {
				
				for (k = 0; k < place.top(); k++){
					traversal[k][i] = traversal[k][i-1];
					alreadyin[traversal[k][i]][i] = 1;
					//cout << "hello" << endl;
					if (traversal[k][i] >= n){
						if (alreadyin[endpoint[traversal[k][i] - n][0]][i] && alreadyin[endpoint[traversal[k][i] - n][1]][i]){
							strong[u][i] = k;
							u++;
							//cout << k << endl;
							//cout << "hello" << endl;
						}
					}
					//cout << "hello" << endl;
				}
//cout << "hello" << endl;
				traversal[place.top()][i] = element.top();
				alreadyin[traversal[place.top()][i]][i] = 1;

				j = place.top() + 1;
//cout << "hello" << endl;
//AND TAKE OFF FROM THE DEGREMAININGS FOR THIS
				for (r = 0; r < n; r++) {
					for (k = 0; k < place.top() + 1; k++){

						if (traversal[k][i] >= n) {
							
							if (endpoint[traversal[k][i]-n][0] == r) {
								
								degremaining[r]--;
							}
							if (endpoint[traversal[k][i]-n][1] == r) {
								degremaining[r]--;
							}
			
						}	
					}
				}
				element.pop();
				place.pop();
			}

//IF, AFTER WE ADDED THE QUEUE AT THE BEGINNING, THE PREVIOUS ELEMENT WAS AN EDGE, PLACE NODE AT THE END OF THIS EDGE IN THE QUEUE
			

			if (traversal[j-1][i] >= n) {

				for (k = 0; k < n; k++) {

					if (alreadyin[k][i] == 0 && (endpoint[(traversal[j-1][i] - n)][0] == k || endpoint[(traversal[j-1][i] - n)][1] == k) ) {		
						element.push(k);
						place.push(j);
						break;
					}
				}
			}

//LOOP OVER EDGES FROM HIGHEST TO LOWEST
			for (k = n + l - 1; k >= n; k--) {
				
				for (m = 0; m < n; m++){

//ADD THE EDGE IF IT'S NOT ALREADY IN THERE AND IF IT'S CONNECTED TO WHAT WE ALREADY HAVE
					if (alreadyin[k][i] == 0 &&  alreadyin[m][i] == 1 && (endpoint[k-n][0] == m || endpoint[k-n][1] == m )) {
						if (alreadyin[endpoint[k-n][0]][i] && alreadyin[endpoint[k-n][1]][i]){
							strong[u][i] = k;
							u++;
						}

						traversal[j][i] = k;
						cout << "putting in an edge: " << traversal[j][i] << endl;
						alreadyin[k][i] = 1;
						for (r = 0; r < n; r++) {
							if (endpoint[k-n][0] == r) {
								degremaining[r]--;
							}
							if (endpoint[k-n][1] == r) {
								degremaining[r]--;
							}
						}
						highestedgeadded = true;
						break;
					}
				}
	
				if (highestedgeadded == true) {
					
					j++;
					highestedgeadded = false;
					break;
				}


//IF WE COULDN'T ADD AN EDGE, ADD THE TOP OF THE VERTEX LIST	
				if (k == n && highestedgeadded == false) {
					

					while (alreadyin[element.top()][i] == 1){
						element.pop();
						place.pop();
					}

					traversal[j][i] = element.top();
					cout << "putting in a vertex: " << traversal[j][i] << endl;
					alreadyin[element.top()][i] = 1;
					j++;
					element.pop();
					place.pop();
					
				}
			}	


		

//BEFORE RUNNING THE EDGE CHECK AGAIN, SEE IF ANY VERTICES NEED TO BE PUT IN BECAUSE ALL OF THEIR EDGES HAVE ALREADY BEEN ADDED
			for (k = 0; k < n; k++){
				if (degremaining[k] == 0 && alreadyin[k][i] == 0) {
					cout << "emergency putting in vertex: " << k << endl;
					traversal[j][i] = k;
					alreadyin[k][i] = 1;
					j++;	
				}
			}	

			if (j == n + l) {
				q++;
			}

		
			if ((i % l) == l - 1 && j == n + l && element.empty() == false) {
	
				for (r = 0; r < n + l; r++) {
					traversal[r].resize(l+q);
					alreadyin[r].resize(l+q);
				}
				
				for (r=0; r < l; r++){
					strong[r].resize(l+q);	
				}
				
				for (r = 0; r < n + l; r++) {
					for (k = q; k < q + l; k++) {
						alreadyin[r][k] = 0;
					}
				}
				for (r = 0; r < l; r++) {
					for (k = q; k < q + l; k++) {
						strong[r][k] = 0;
					}
				}
			}
			

		}
	

	cout << endl << endl << endl;
	}

	
	cout << endl;

	


		
	cout << "TRAVERSAL PRINT: " << endl;	
	matrixprint(traversal, n + l, travs);	
	cout << endl << endl;
	cout << "STRONG EDGES PRINT (SAME NUMBERING SYSTEM AS ABOVE)" << endl;
	matrixprint(strong, l, travs);
	cout << endl << endl;
	/*cout << "NUMBER OF TRAVERSALS OF GIVEN STRENGTH: " << endl;
	cout << strongcounter (istravstrong, strong, level, l, travs) << endl;*/


	int howmanyatcorrectlevel;
	vector<bool> istravatrightlevel(travs);
	howmanyatcorrectlevel = strongcounter(istravatrightlevel, strong, level, l, travs);
	cout << "Number of configs at level 0: " << howmanyatcorrectlevel << endl;
	vector<vector<int> > chips(n);
	for (i = 0; i < n; i++){
		chips[i].resize(howmanyatcorrectlevel);
	}

	ChipCounter(traversal, istravatrightlevel, endpoint, chips, travs, n, l);
	matrixprint(chips, n, howmanyatcorrectlevel);

	WriteMathematicaFile (chips, howmanyatcorrectlevel, n);


	

	return 0;


}


void matrixprint(vector<vector<int> > matrix, int R, int C){	
		int matxpnt1;
		int matxpnt2;
		for(matxpnt1 = 0; matxpnt1 < R; matxpnt1++){
			for(matxpnt2=0; matxpnt2 < (C-1); matxpnt2++){
				cout << matrix[matxpnt1][matxpnt2] << "      ";
			}
			cout << matrix[matxpnt1][C-1] <<endl;
		}


}
	






/*void PrintStack(const Stack<int>& list) {
	
	if (list.empty == true) {
		cout << "Stack empty" << endl;
	}
	
	while (list.empty() == false){
		cout << list.top() << endl;
		list.pop();
	}

	return 0;
}*/


int strongcounter (vector<bool> &istravatrightlevel, vector<vector<int> > strong, int level, int l, int travs) {

	int totalatlevel = 0;
	int strongedges;
	
	int i, j;

	for (i = 0; i < travs; i++) { 
		strongedges = 0;
		for (j = 0; j < l; j++) {
			if (strong[j][i] != 0) {
				strongedges++;
			}
		}
		if (strongedges == level){
			totalatlevel++;
			istravatrightlevel[i] = true;
		} else {
			istravatrightlevel[i] = false;
		}
	}

	return totalatlevel;
}			
	

void ChipCounter(vector<vector<int> > Function, vector<bool> istravatrightlevel, vector<vector<int> > endpoint, vector<vector<int> > &chips, int travs, int n, int l){

	int i, j, k, m, q;
	q = 0;

	int elements = n + l;

	for (i = 0; i < travs; i++){
		k = 0;
		if (istravatrightlevel[i] == true){
//cout << "hello" << endl;
			for (j = 0; j < elements; j++){
				if (Function[j][i] < n){
					chips[Function[j][i]][q] = 0;
					//see how many edges remain
					for(m = j; m < elements; m++){
//cout << "hello" << endl;
						
						if (Function[m][i] >= n) {
							
							if (endpoint[Function[m][i] - n][0] == Function[j][i]) {
								if (i != 0 && Function[j][i] == 1) {cout << "here and q is " << q << endl;}
								chips[Function[j][i]][q]++;
//cout << "hello" << endl;
							}
							if (endpoint[Function[m][i] - n][1] == Function[j][i]) {
								if (i != 0 && Function[j][i] == 1) {cout << "here and q is " << q << endl;}
								chips[Function[j][i]][q]++;
							}
//cout << "hello" << endl;
						}
					}
				}
			}
		q++;
		}
	}

	
}

void WriteMathematicaFile (vector<vector<int> > &chips, int howmanyatcorrectlevel, int n){
	int i, j;
	

	ofstream outfile;
	outfile.open("C3SolutionsTake2.nb", ios::out);
	

	outfile << "SolutionCandidates3 = {";
	

	for (i = 0; i < howmanyatcorrectlevel; i++) {
		
		if (i != 0){
			outfile << ",{";
		} else {
			outfile << "{";
		}
	/*
		outfile << "{3,3,3,";
		outfile << chips[3][i];	
		outfile	<< ",3,3,";
		outfile << chips[4][i] << "," << chips[5][i] << ","<< chips[2][i] << ",";
		

		outfile << "3," << chips[1][i] << ",3}";*/
	



		for (j = 1; j < n; j++){
			if ((j-1)%2 == 0) {
				outfile << "{3," << chips[j][i];
				
			} else if (j != (n-1)) {
				outfile << "," << chips[j][i] << "},";
				
			} else {
				outfile << "," << chips[j][i] << "},{3,3,3}}";
				
			}
		}

		
			


		if (i == howmanyatcorrectlevel - 1){
			outfile << "}";
		}
	}
	
			

					

		
	
	
	outfile.close();
}



/*
//PRINTSTACK

cout << "PRINT STACK" << endl;
r = 0;
vector<int> tempstack(n+l, -1);
if (element.empty() == true) {
	cout << "Stack empty" << endl;
} else {
	while( element.empty() == false) {
		tempstack[r] = element.top();
		r++;
		cout << element.top() << endl;
		element.pop();
	}

	for (r = n + l - 1; r >= 0; r--){
		if (tempstack[r] != -1) {
			element.push(tempstack[r]);
		}
	}	
}

cout << "Printstack done" << endl;
*/


		