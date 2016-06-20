#include "ProfitPartitionMatrix.hh"
#include "SimplePartitionedMatrix.hh"
#include "PartitionLoader.hh"
#include <iostream>
#include <cassert>
#include <fstream>


using namespace std;

/**
 * @Description : Constructor
 * @Creation Date: 17.08.2009
 */

ProfitPartitionMatrix::ProfitPartitionMatrix(int M, int N, int nz,bool value, PartitionLoader* partitionLoader)
    :SimplePartitionedMatrix(M,N,nz,value, partitionLoader)
{

}
/* ProfitPartitionMatrix() ENDS*/

/**
 * @Description : Destructor
 * @Creation Date: 17.08.2009
 */

ProfitPartitionMatrix::~ProfitPartitionMatrix()
{
    delete[] this->nnzTag;
}
/* ~ProfitPartitionMatrix() ENDS*/


/**
 * @Description : This function is a prototype for coloring in the
 * segmented graph.
 * @Creation Date: 05.08.2009
 * @Author: Mahmudul Hasan
 */
int
ProfitPartitionMatrix::cseg_profit(NNZTag* nnzTag)
{

    // Save this nnzTag for subsequent methods.
    this->nnzTag = nnzTag;
#ifdef DEBUG
    cout << "Matrix::cseg_profit() -- BEGIN" <<endl;
#endif
    // Our priority queue.

    int *head = new int[2*N];
    int *prev = new int[N+1];
    int *next = new int[N+1];
    int *list = new int[N+1];
    int *tag = new int[N+1]; // For Debugging Purpose.


    int *rowStatus = new int[M+1];

    const int already_conflicted = 1;
    const int free_to_use = 2;
    const int freshly_used = 3;

    const int NOT_COVERED = nnz + 1;

    for(int jcol =1; jcol <= N; jcol++)
    {
        for (int jp = jpntr[jcol] ; jp < jpntr[jcol+1] ; jp++)
        {
            int ir = indRow[jp];
            nnzTag->setValue(ir,jcol, NOT_COVERED);
        }

    }


    int maxdeg = 0;

    // Initialize the priority queue with number of non-zeroes in each
    // column.

    int groupNumber = 0;

    int c_coveredNNZ = 0;




    // Group Forming Loop.
    bool newGroupNeeded = true;
    bool newColumnNeeded = false;
    while(newGroupNeeded && (c_coveredNNZ != nnz) )
    {

        newGroupNeeded = false;
        groupNumber++;
        assert(groupNumber < NOT_COVERED);
        for(int i = 1; i <= M; i++)
        {
            rowStatus[i] = free_to_use;
        }
#ifdef DEBUG
        cout << "GroupNumber " << groupNumber << endl;
#endif
        newColumnNeeded = true;
        maxdeg = 0;
        int n_coveredNNZ = 0;
        int c_count = 0;
        for(int jcol = 1; jcol <= N; jcol++)
        {
            tag[jcol] = N;
            list[jcol] = 0;
            for(int jp= jpntr[jcol] ; jp < jpntr[jcol+1]; jp++)
            {
                int ir = indRow[jp];
                if(nnzTag->getValue(ir,jcol) == NOT_COVERED)
                {
                    list[jcol]++;
                    c_count++;
                    maxdeg = max(maxdeg,list[jcol]);
                }
                else
                {
                    n_coveredNNZ++;
                }
            }
        }
        assert(c_count+n_coveredNNZ == nnz);
        if(c_coveredNNZ != n_coveredNNZ)
        {
            cout << "c_coveredNNZ " << c_coveredNNZ << " n_coveredNNZ " << n_coveredNNZ << endl;
            cout << "GroupNumber " << groupNumber << endl;
        }
        assert(c_coveredNNZ == n_coveredNNZ);
        for(int deg = 0; deg < 2 * N; deg++)
        {
            head[deg] = 0;
        }

        for(int jcol = 1; jcol <= N; jcol++)
        {
            int numdeg = list[jcol];
            prev[jcol] = 0;
            next[jcol] = head[numdeg + N];
            if(head[numdeg+N] > 0)
            {
                prev[head[numdeg+N]] = jcol;
            }
            head[numdeg+N] = jcol;
        }

        assert(maxdeg > 0);
        while(newColumnNeeded)
        {
            // Step 1. get a jcol from the priority queue.
            // Step 2. Mark the conflicted rows for this for the group.
            // Step 3. Unmark all the tagged nnz in the conflicted
            // rows.

            int jcol;
            int numdeg;

            do
            {
                if(maxdeg < 0)
                {
                    newColumnNeeded = false;
                    break;
                }
                jcol = head[maxdeg+N];
                numdeg = maxdeg;
                if(jcol >0 )
                {
#ifdef DEBUG
                    cout << "Numdeg " << numdeg << endl ;
#endif
                    break;
                }
                maxdeg = maxdeg -1;
            }while(true);


            if(newColumnNeeded == false)
            {
                if(c_coveredNNZ != nnz)
                    newGroupNeeded = true;
#ifdef DEBUG
                cout << "c_coveredNNZ " << c_coveredNNZ << endl;
#endif
                break;
            }

            if(tag[jcol] != N )
            {
                cout << "Column " << jcol << " Tagged " << tag[jcol] << endl;
                cout << "Numdeg " << numdeg << endl;
            }
            assert(tag[jcol] == N);
            tag[jcol] = groupNumber;
            deleteColumn(head,next,prev,numdeg+N,jcol);

#ifdef DEBUG
            cout << "Column " << jcol << endl;
#endif
            // Traverse all the nonzeroes in the jcol
            for(int jp = jpntr[jcol] ; jp < jpntr[jcol+1] ; jp++)
            {
                int ir = indRow[jp];
                int nnzStatus = nnzTag->getValue(ir,jcol);
                switch(rowStatus[ir])
                {
                case free_to_use:
                    // Present
                    if(nnzStatus == NOT_COVERED)
                    {
                        rowStatus[ir] = freshly_used;
                        nnzTag->setValue(ir,jcol,groupNumber);
                        c_coveredNNZ++;

#ifdef DEBUG
                        cout << "c_coveredNNZ " << c_coveredNNZ << endl;
#endif

                        for(int ip = ipntr[ir] ; ip < ipntr[ir+1] ; ip++)
                        {
                            int ic = indCol[ip];
                            if(ic == jcol || tag[ic] == groupNumber)
                                continue;
                            int numdeg = list[ic];
                            list[ic]--;
                            if(nnzTag->getValue(ir,ic) == NOT_COVERED)
                            {
                                list[ic]--;
                            }
                            deleteColumn(head,next,prev,numdeg+N,ic);
                            addColumn(head,next,prev,list[ic]+N,ic);
                        }
                    }
                    else
                    {
                        rowStatus[ir] = already_conflicted;
                        for(int ip = ipntr[ir] ; ip < ipntr[ir+1]; ip++)
                        {
                            int ic = indCol[ip];
                            if(ic == jcol || tag[ic] == groupNumber)
                                continue;

                            if(nnzTag->getValue(ir,ic) == NOT_COVERED)
                            {
                                int numdeg = list[ic];
                                list[ic]--;
                                deleteColumn(head,next,prev,numdeg+N,ic);
                                addColumn(head,next,prev,list[ic]+N,ic);
                            }
                        }
                    }
                    // Future
                    break;
                case already_conflicted:
                    // Present
                    // Future
                    break;
                case freshly_used:
                    // Present
                    rowStatus[ir] = already_conflicted;
                    for (int ip = ipntr[ir] ; ip < ipntr[ir+1]; ip++)
                        // TODO: This loop can be avoided.
                    {
                        int ic = indCol[ip];
                        if(ic == jcol)
                            continue;
                        // cout << "nnzTag->getValue(ir,ic ) " << nnzTag->getValue(ir,ic) << endl;
                        if(nnzTag->getValue(ir,ic) == groupNumber)
                        {
                            nnzTag->setValue(ir,ic, NOT_COVERED);
                            c_coveredNNZ--;
#ifdef DEBUG
                            cout << "Decreasing: (groupNumber,  " << groupNumber << ") (ir,ic) = (" << ir << " , " << ic << " ) " << endl;
                            cout << "c_coveredNNZ " << c_coveredNNZ << endl;
#endif
                            //break; // We can only have one.
                        }
                        // Previously , we are doing it for all,
                        // With a simple else,
                        // now we are doing this such that the column
                        // is not in the current group.
                        else if ( tag[ic] != groupNumber)// if (nnzTag->getValue(ir,ic) == NOT_COVERED)
                        {
                            int numdeg = list[ic];
                            list[ic]++;
                            maxdeg = max(list[ic], maxdeg);
                            deleteColumn(head,next,prev,numdeg+N,ic);
                            addColumn(head,next,prev,list[ic]+N,ic);
                        }
                    }

                    // Future
                    break;
                default:
                    // error
                    break;
                }
            }
        }
    }


    delete[] head;
    delete[] next;
    delete[] prev;
    delete[] list;
    delete[] tag;

    generate_partition_from_nnz_color(nnzTag);

    return groupNumber;
}
/* cseg-profit() ENDS*/

/**
 * @Description : Genearte Paritioning Information from the NNZ
 * Color.
 * @Creation Date: 15.08.2009
 */
void
ProfitPartitionMatrix::generate_partition_from_nnz_color(NNZTag* nnzTag)
{
    int *row = new int[M+1];
    for (int ir = 1; ir <= M; ir++)
    {
        int ip = ipntr[ir];
        int ic = indCol[ip];
        int tval = nnzTag->getValue(ir,ic);
        row[ir] = tval;
#ifdef DEBUG
        cout << "Row[" << ir << "] = " << tval << endl;
#endif
    }

    delete[] row;
}
/* generate_partition_from_nnz_color() ENDS*/

/**
 * @Description : This functions
 * @Creation Date: 17.08.2009
 */
bool
ProfitPartitionMatrix::computedegree()
{

}
/* degr() ENDS*/

/**
 * @Description : This functions
 * @Creation Date: 17.08.2009
 */
bool
ProfitPartitionMatrix::slo(int* list)
{

}
/* slo() ENDS*/

/**
 * @Description : This functions
 * @Creation Date: 17.08.2009
 */
bool
ProfitPartitionMatrix::ido(int *order)
{

}
/* ido() ENDS*/

/**
 * @Description : This functions
 * @Creation Date: 17.08.2009
 */
bool
ProfitPartitionMatrix::lfo(int *order)
{

}
/* lfo() ENDS*/

/**
 * @Description : This functions
 * @Creation Date: 17.08.2009
 */
int
ProfitPartitionMatrix::rlf(int *ngrp)
{

}
/* rlf() ENDS*/

/**
 * @Description : This functions
 * @Creation Date: 17.08.2009
 */
int
ProfitPartitionMatrix::greedycolor(int *list, int *ngrp)
{

}
/* seq() ENDS*/


/**
 * @Description :
 * @Creation Date: 17.08.2009
 */
void
ProfitPartitionMatrix::writeColor(Result result, int *ngrp)
{
    ofstream out("out.color");
    out << result.totalColors << endl;
    for (int jcol = 1; jcol <= N; jcol++)
    {
        for (int jp = jpntr[jcol] ; jp < jpntr[jcol+1]; jp++)
        {
            int ir = indRow[jp];
            int segment = getSegment(ir,jcol);
            int tval = nnzTag->getValue(ir,jcol);
            out << segment << " " << tval << endl;
        }
    }
}
/* writeColor() ENDS*/














