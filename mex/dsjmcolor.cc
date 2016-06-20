#include "mex.h"
#include "matrix.h" // MATLAB Probably !!
#include "Matrix.hh"
#include "CLI.h"

void error_return(mxArray *plhs[])
{
    plhs[0] = mxCreateDoubleScalar(-1);
    double *err_val = mxGetPr(plhs[0]);
    *err_val = -1;
}


/**
 * @Description: Entry Point
 * @Creation Date: 30.07.2009
 */
void mexFunction(
                 int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mexPrintf("========================================\n");
    mexPrintf("dsjmcolor\n");
    int mrows, ncols;
    double *y;

    /* Check for Proper number of Arguments */
    if(nrhs != 2)
    {
        mexErrMsgTxt("Two input required");
    }
    else if (nlhs > 2)
    {
        mexErrMsgTxt("Too many output arguments");
    }

    /* Print Out the Dimension */
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);

    mexPrintf("MxN : %d x %d \n", mrows,ncols);


    // Print out the method name
    char *fcn;
    int   status;
    int   buflen=mxGetN(prhs[1])+1;

    fcn=(char *)mxCalloc(buflen,sizeof(char));
    status=mxGetString(prhs[1],fcn,buflen);

    // status=mexEvalString(fcn);
    printf("The Ordering Requested is %s\n",fcn);

    // Determine the oRdering method requested.
    CLI::ordering_method ordering_method = CLI::parseMethod(fcn);




    bool valid = true;

    if(mxIsDouble(prhs[0]))
    {
        mexPrintf("IsDouble = TRUE\n");
    }
    else
    {
        valid = false;
        mexPrintf("IsDouble == False\n");
        mexPrintf("Data has to be double\n");
    }

    if(mxIsComplex(prhs[0]))
    {
        valid = false;
        mexPrintf("IsComplex == True\n");
        mexPrintf("Complex data is not supported\n");
    }
    else
    {
        mexPrintf("IsComplex == False\n");
    }

    if ( mxIsSparse(prhs[0]) )
    {
        mexPrintf("IsSparse == True\n");
    }
    else
    {
        valid = false ;
        mexPrintf("IsSparse == False\n");
        mexPrintf("We operate only on Sparse Matrices\n");

        // Quick Fix:
        // TODO:
        // plhs[0] = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,
        // mxREAL);

        error_return(plhs);

        return ;
    }

    /**
     * Declare the Matrix Object
     */
    IMatrix *matrix;

    int nzmax = mxGetNzmax(prhs[0]);
    int *jc , *ir;
    ir = mxGetIr(prhs[0]);
    jc = mxGetJc(prhs[0]);
    mexPrintf("nzmax = %d\n",nzmax);
    mexPrintf("nnz = %d\n", jc[ncols]);

    matrix = new Matrix(mrows,ncols,jc[ncols],false);


    for (int j = 0, i = 1; j < ncols; j++)
    {
        for(int rowIndex = jc[j] ; rowIndex < jc[j+1] ; rowIndex++)
        {
            matrix->setIndRowEntry(i,ir[rowIndex] + 1);
            matrix->setIndColEntry(i, j + 1);
            // mexPrintf("Inputting (%d, %d)\n",ir[rowIndex] + 1, j+1);
            i++;
        }
    }
    mexPrintf("Input Given\n");

    matrix->computeCCS();

    int nnz = matrix->compress();

    matrix->computeCRS();

    matrix->computedegree();

    int maxgrp;

    int *color = new int[ncols+1];

    if (ordering_method == CLI::SLO )
    {
        int *order = new int[ncols+1];
        matrix->slo(order);


        maxgrp = matrix->greedycolor(order, color);
        delete[] order;
    }
    else if( ordering_method == CLI::IDO )
    {

        int *order = new int[ncols+1];
        matrix->ido(order);
        maxgrp = matrix->greedycolor(order,color);
        delete[] order;
    }
    else if ( ordering_method == CLI::LFO)
    {

        int *order = new int[ncols+1];
        matrix->lfo(order);
        maxgrp = matrix->greedycolor(order,color);
        delete[] order;
    }
    else if ( ordering_method == CLI::RLF)
    {
        maxgrp =  matrix->rlf(color);

    }
    else if ( ordering_method == CLI::SDO)
    {
        Matrix *nmatrix = dynamic_cast<Matrix*>(matrix);
        maxgrp = nmatrix->sdo(color);
    }
    else
    {
        mexPrintf("Unknown Method for Ordering\n");
        error_return(plhs);
        // TODO: Can't we implemente a Smart Pointer so
        // that matrix gets deleted automatically once
        // it is out of scope of this function. ?
        delete matrix;
        delete[] color;
        return;
    }

    //    maxgrp = matrix->greedycolor();

    mexPrintf("Colors = %d\n",maxgrp);



    plhs[0] = mxCreateDoubleMatrix(ncols,1,mxREAL);
    y = mxGetPr(plhs[0]);
    for (int i = 0; i < ncols; i++)
    {
        // mexPrintf("Coloring: %d ->
        // %d\n",i+1,matrix->getNgrpEntry(i+1));
        *(y+i) = color[i+1];
    }
    delete matrix;
}
/* mexFunction() ENDS */
