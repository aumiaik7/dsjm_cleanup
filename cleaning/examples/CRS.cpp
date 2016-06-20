/* Step 5 */
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "algorithm2.h"
using namespace std;

void CRS(int *Ii, int *Jj, double *valc, double **spe, double **res, int nz, int N, int g, int *color, double *delta){
	int ia=0;
	int ai=-1;
	int l1=0;
	//for (int l=0;l<nz1;l++){
  		for (int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				if (spe[i][j]!=0){
				Jj[l1]=j;
				valc[l1]=spe[i][j];
					if(i>ai){
						Ii[ia]=l1;
						ia++;}
					ai=i;
					l1++;
				}
			}
		}
		Ii[ia]=nz;
ofstream out20;
		out20.open("IJ.txt"); 
out20<<"nz :"<<nz;
		out20<<"Ii: \n";
	for (int l = 0; l < ia+1; l++)
		 out20 << Ii[l] << ' ';
	out20<<"\nJi: \n"<<endl;
	for (int l = 0; l < nz; l++)
		 out20 << Jj[l] << ' ';
	out20<<"\nval: \n"<<endl;
	for (int l = 0; l < nz; l++)
		 out20 << valc[l] << ' ';
	out20<<endl;
	out20.close();

	double ee[N][g];
	double d[N][g];
		for (int m = 0; m < N; m++) 
			for (int n = 0; n < g; n++)
				ee[m][n]=0;
	cout<<"##############"<<endl;
		for (int pc = 0; pc < N; pc++) {
			//cout<<"pc: "<<pc<<"/ ";
			for (int y = 0; y < g; y++)
				if(color[pc]-1==y){
					//cout<<"y: "<<y<<"/ ";
					d[pc][y]=pc+1; 
					ee[pc][y]=delta[pc];}
				}
	cout<<"$$$$$$$$$$$$$$$"<<endl;
   		for (int pc = 0; pc < N; pc++) 
			for (int y = 0; y < g; y++)
      				for (int k = Ii[pc]; k < Ii[pc+1]; k++) 
         				res[pc][y] = res[pc][y] + valc[k]*ee[Jj[k]][y];

	//for (int pc = 0; pc < N; pc++){ 
			//for (int y = 0; y < g; y++)
				//cout<<d[pc][y]<<" ";
		//cout<<endl;}
}






