/* Step 10 */
#include <algorithm>
//#include "algorithm1p.h"
#include "algorithm2.h"
using namespace std;

/*Elimination of spurious zeros */
void procedure2(double **sz, double**St, double **rt, int N, int a, double ab, double *delta){
/*for (int i=0; i<N; i++)
		for (int j=0; j<N; j++){
			if (j>=i){
				int z = min(rt[i][j],rt[j][i]);
				if (z < 1)
					sz[i][j]=sz[j][i]=1;
			}
		}*/
for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			sz[i][j]=1;

	for (int i=0; i<N; i++)
		for (int j=0; j<a; j++)
			if(rt[i][j]!=0)	
				for (int k=0;k<N;k++)
					if(delta[k]==rt[i][j])
						sz[i][k]=0;
						
	
	/*Final S*/
	
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++){
				if (sz[i][j]==1){
					rt[i][j]=0;
					St[i][j]=0;}
				if (St[i][j]!=0)
					St[i][j]=1;			
				}
}
