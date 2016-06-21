#include <iostream>
#include <cstdlib>
#include <cstring>
#include "mmio.h"
#include "Matrix.hh"


void display_usage(int argc, char *argv[])
{
    fprintf (stderr, "Usage: %s -i MatrixFile", argv[0]);
    exit(1);
}
int main(int argc, char *argv[])
{

    if(argc < 3)
    {
        fprintf (stderr, "Usage: %s -i MatrixFile", argv[0]);
        exit(1);
    }

    int opt = getopt(argc,argv, "i:h?");
    int n = -1;
    double d;
    char *inputfile = NULL;
    while(opt != -1)
    {
        switch(opt)
        {
        case 'i':
            if(optarg != NULL && strlen(optarg ) > 0)
            {
                if(inputfile != NULL)
                    delete[] inputfile ;
                inputfile = new char[strlen(optarg) + 1];
                strcpy(inputfile, optarg);
            }
            break;
        case 'd':
            d = atof(optarg);
            break;
        case 'h':
        case '?':
            display_usage(argc, argv);
        default:
            break;
        }
        opt = getopt(argc,argv,"n:x:d:h?");
    }

    MM_typecode matcode;
    int ret_code;

    FILE *f;

    f  = fopen(inputfile, "r");
    if (mm_read_banner (f, &matcode) != 0)
    {
        fprintf (stderr,"%s -> Could not process Matrix Market banner.\n",inputfile);
        exit (1);
    }
    if (mm_is_complex (matcode) && mm_is_matrix (matcode) &&
        mm_is_sparse (matcode))
    {
        printf ("Sorry, this application does not support ");
        printf ("Market Market type: [%s]\n", mm_typecode_to_str (matcode));
        exit (1);
    }

    int M,N, nz;

    if ((ret_code = mm_read_mtx_crd_size (f, &M, &N, &nz)) != 0)
        exit (1);

    int NZ = nz;

    bool is_symmetric = mm_is_symmetric(matcode);
    bool is_pattern = mm_is_pattern(matcode);
    if(is_symmetric)
        NZ = 2 * nz;

    Matrix *matrix;
    matrix = new Matrix(M,N,NZ);

    double value;

    for ( int i = 1,row,col ; i <= nz; i++ )
    {
        if (is_pattern)
        {
            fscanf (f, "%d %d\n", &row, &col);
        }
        else
        {
            fscanf (f, "%d %d %lg\n", &row, &col, &value);
        }

        matrix->setIndRowEntry(i,row);
        matrix->setIndColEntry(i,col);

        if(is_symmetric)
        {
            matrix->setIndRowEntry(i + nz,col);
            matrix->setIndColEntry(i + nz,row);
        }

    }

    if (f != stdin)
        fclose(f);

    matrix->srtdat();

    int nnz = matrix->compress();

    matrix->setr();
    matrix->degr();

    for(int ir = 1; ir <= M ; ir++)
    {
        int colind = -1;
        // cout << "Row : " << ir << endl ;
        for(int ip = matrix->getIpntrEntry(ir); ip < matrix->getIpntrEntry(ir+1); ip++)
        {
            int t = matrix->getIndColEntry(ip);
            // cout << t << endl;

            if (t >= colind )
            {
                colind = t;
            }
            else
            {
                assert( 0);
            }
        }

    }

    delete matrix;



    return 0;
}

