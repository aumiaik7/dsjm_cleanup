/* testsrtdat.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__9 = 9;
static integer c__3 = 3;
static integer c_b41 = 1302464;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1, i__2;
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *), s_wsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, n, jp, ir, iwa[1582459];
    static char rep[10];
    static integer nnz;
    static complex cval[1302464];
    static integer ival[1302464];
    static doublereal rval[1302464];
    static char symm[19], field[7];
    extern integer iargc_(void);
    static char ifile[32];
    static integer ncols, iunit, jpntr[263744], nrows;
    extern /* Subroutine */ int mmread_(integer *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, complex *, ftnlen, ftnlen, ftnlen);
    static integer indcol[1302464];
    extern /* Subroutine */ int getarg_(integer *, char *, ftnlen), mminfo_(
	    integer *, char *, char *, char *, integer *, integer *, integer *
	    , ftnlen, ftnlen, ftnlen), srtdat_(integer *, integer *, integer *
	    , integer *, integer *, integer *);
    static integer indrow[1302464];

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___28 = { 0, 6, 0, 0, 0 };
    static cilist io___29 = { 0, 6, 0, 0, 0 };



    if (iargc_() > 0) {
	getarg_(&c__1, ifile, (ftnlen)32);
	iunit = 8;
	o__1.oerr = 0;
	o__1.ounit = iunit;
	o__1.ofnmlen = 32;
	o__1.ofnm = ifile;
	o__1.orl = 0;
	o__1.osta = 0;
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    } else {
	iunit = 5;
    }
    s_wsle(&io___3);
    do_lio(&c__9, &c__1, "READING HEADER ONLY...", (ftnlen)22);
    e_wsle();
    mminfo_(&iunit, rep, field, symm, &nrows, &ncols, &nnz, (ftnlen)10, (
	    ftnlen)7, (ftnlen)19);
    s_wsle(&io___10);
    do_lio(&c__9, &c__1, "  MATRIX IS TYPE: ", (ftnlen)18);
    do_lio(&c__9, &c__1, rep, (ftnlen)10);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    do_lio(&c__9, &c__1, field, (ftnlen)7);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    do_lio(&c__9, &c__1, symm, (ftnlen)19);
    e_wsle();
    s_wsle(&io___11);
    do_lio(&c__9, &c__1, "  MATRIX SIZE: ", (ftnlen)15);
    do_lio(&c__3, &c__1, (char *)&nrows, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " BY ", (ftnlen)4);
    do_lio(&c__3, &c__1, (char *)&ncols, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " WITH ", (ftnlen)6);
    do_lio(&c__3, &c__1, (char *)&nnz, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " NONZEROS.", (ftnlen)10);
    e_wsle();

    s_wsle(&io___12);
    do_lio(&c__9, &c__1, "READING HEADER AND DATA...", (ftnlen)26);
    e_wsle();
    mmread_(&iunit, rep, field, symm, &nrows, &ncols, &nnz, &c_b41, indrow, 
	    indcol, ival, rval, cval, (ftnlen)10, (ftnlen)7, (ftnlen)19);
    s_wsle(&io___18);
    do_lio(&c__9, &c__1, "  MATRIX IS TYPE: ", (ftnlen)18);
    do_lio(&c__9, &c__1, rep, (ftnlen)10);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    do_lio(&c__9, &c__1, field, (ftnlen)7);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    do_lio(&c__9, &c__1, symm, (ftnlen)19);
    e_wsle();
    s_wsle(&io___19);
    do_lio(&c__9, &c__1, "  MATRIX SIZE: ", (ftnlen)15);
    do_lio(&c__3, &c__1, (char *)&nrows, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " BY ", (ftnlen)4);
    do_lio(&c__3, &c__1, (char *)&ncols, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " WITH ", (ftnlen)6);
    do_lio(&c__3, &c__1, (char *)&nnz, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " NONZEROS.", (ftnlen)10);
    e_wsle();
    srtdat_(&ncols, &nnz, indrow, indcol, jpntr, iwa);
    i__1 = nrows;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwa[i__ - 1] = 0;
/* L20: */
    }
    nnz = 1;
    i__1 = ncols;
    for (j = 1; j <= i__1; ++j) {
	k = nnz;
	i__2 = jpntr[j] - 1;
	for (jp = jpntr[j - 1]; jp <= i__2; ++jp) {
	    ir = indrow[jp - 1];
	    if (iwa[ir - 1] != j) {
		indrow[nnz - 1] = ir;
		++nnz;
		iwa[ir - 1] = j;
	    }
/* L30: */
	}
	jpntr[j - 1] = k;
/* L45: */
    }
    jpntr[n] = nnz;
/*      CALL SETR(NROWS,NCOLS,INDROW,JPNTR,INDCOL,IPNTR,IWA) */
/*      CALL DEGR(NCOLS,INDROW,JPNTR,INDCOL,IPNTR,IWA(5*NCOLS+1) */
/*     *     ,IWA(NCOLS+1)) */
    i__1 = nnz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsle(&io___28);
	do_lio(&c__3, &c__1, (char *)&indrow[i__ - 1], (ftnlen)sizeof(integer)
		);
	do_lio(&c__9, &c__1, ",", (ftnlen)1);
	do_lio(&c__3, &c__1, (char *)&indcol[i__ - 1], (ftnlen)sizeof(integer)
		);
	e_wsle();
/* L40: */
    }
    i__1 = ncols + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*         PRINT *,JPNTR(I),',',IPNTR(I),',',IWA(5*NCOLS+1+I) */
	s_wsle(&io___29);
	do_lio(&c__3, &c__1, (char *)&jpntr[i__ - 1], (ftnlen)sizeof(integer))
		;
	e_wsle();
/* L50: */
    }
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Main program alias */ int main_ () { MAIN__ (); return 0; }
