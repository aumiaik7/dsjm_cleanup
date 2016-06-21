/* testio.f -- translated by f2c (version 20050501).
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
static integer c__5 = 5;
static integer c__6 = 6;
static integer c_b197 = 263744;
static integer c__2 = 2;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1, i__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_wsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_wsle(void), s_cmp(char *, char *, ftnlen, 
	    ftnlen), f_clos(cllist *);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int ido_(integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *);
    static integer iwa[1582459];
    static char rep[10];
    static integer nnz;
    extern /* Subroutine */ int degr_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static complex cval[1302464];
    static integer ival[1302464];
    static doublereal rval[1302464];
    extern /* Subroutine */ int setr_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static char symm[19], field[7];
    extern integer iargc_(void);
    static char ifile[32], ofile[32];
    static integer ncols, iunit, ipntr[263744], jpntr[263744], ounit__, nrows;
    extern /* Subroutine */ int mmread_(integer *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, complex *, ftnlen, ftnlen, ftnlen);
    static integer indcol[1302464];
    extern /* Subroutine */ int getarg_(integer *, char *, ftnlen);
    static integer maxclq;
    extern /* Subroutine */ int mminfo_(integer *, char *, char *, char *, 
	    integer *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer mingrp;
    extern /* Subroutine */ int srtdat_(integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static integer indrow[1302464];
    extern /* Subroutine */ int mmwrite_(integer *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, complex *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, 0, 0 };
    static cilist io___27 = { 0, 6, 0, 0, 0 };
    static cilist io___28 = { 0, 6, 0, 0, 0 };
    static cilist io___29 = { 0, 6, 0, 0, 0 };
    static cilist io___30 = { 0, 6, 0, 0, 0 };
    static cilist io___31 = { 0, 6, 0, 0, 0 };
    static cilist io___32 = { 0, 6, 0, 0, 0 };
    static cilist io___33 = { 0, 6, 0, 0, 0 };
    static cilist io___34 = { 0, 6, 0, 0, 0 };
    static cilist io___35 = { 0, 6, 0, 0, 0 };
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___38 = { 0, 6, 0, 0, 0 };
    static cilist io___41 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___49 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };
    static cilist io___55 = { 0, 6, 0, 0, 0 };


/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */

/* USAGE:  TESTIO < DATAFILE.MTX */
/*         TESTIO DATAFILE.MTX */
/*         GUNZIP -C DATAFILE.MTX.Z | TESTIO */

/* NOTE:   MAKE 'TESTIO' EXECUTABLE WITH:   F77 -O TESTIO TESTIO.F MMIO.F */

/* THIS SAMPLE DRIVER TAKES A MATRIX MARKET FILE AS STANDARD INPUT, */
/* OR OPENS A THE SPECIFIED FILE IF A FILENAME ARGUMENT IS PRESENT, */
/* AND CALLS TWO ROUTINES: */

/*     CALL MMINFO(IUNIT,REP,FIELD,SYMM,NROWS,NCOLS,NNZ) */
/* AND */
/*     CALL MMREAD(IUNIT,REP,FIELD,SYMM,NROWS,NCOLS,NNZ,NNZMAX, */
/*    *             INDX,INDCOL,IVAL,RVAL,CVAL) */

/* SUBROUTINE MMINFO JUST PARSES HEADER INFORMATION, WHILE MMREAD PARSES */
/* THE HEADER AND READS THE NUMERICAL DATA INTO THE APPROPRIATE ARRAYS. */
/* EACH OF THESE REQUIRES A PRE-OPENED READ UNIT (IUNIT), AND */
/* REWINDS THE UNIT PRIOR TO RETURN. (SO MMINFON DOESN'T INTERFERE */
/* WITH MMREAD AND VICE VERSA). */

/* 18-OCT-96  KARIN A. REMINGTON, NIST ACMD (KARIN@CAM.NIST.GOV) */
/* 30-OCT-96  MINOR CHANGE TO CALLING SEQUENCES TO ACCOMMODATE NEW */
/*              INTEGER VALUE PARAMETER (IVAL) */
/*            + ASSOCIATED CHANGES TO CHECK FOR INTEGER VALUE DATA */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */

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

    if (s_cmp(rep, "ARRAY", (ftnlen)10, (ftnlen)5) == 0) {
	s_wsle(&io___20);
	do_lio(&c__9, &c__1, "  DENSE ARRAY.", (ftnlen)14);
	e_wsle();
	s_wsle(&io___21);
	do_lio(&c__9, &c__1, "  FIRST TWO ENTRIES:", (ftnlen)20);
	e_wsle();
	if (s_cmp(field, "INTEGER", (ftnlen)7, (ftnlen)7) == 0) {
	    s_wsle(&io___22);
	    for (i__ = 1; i__ <= 2; ++i__) {
		do_lio(&c__3, &c__1, (char *)&ival[i__ - 1], (ftnlen)sizeof(
			integer));
	    }
	    e_wsle();
	} else if (s_cmp(field, "REAL", (ftnlen)7, (ftnlen)4) == 0) {
	    s_wsle(&io___24);
	    for (i__ = 1; i__ <= 2; ++i__) {
		do_lio(&c__5, &c__1, (char *)&rval[i__ - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsle();
	} else if (s_cmp(field, "COMPLEX", (ftnlen)7, (ftnlen)7) == 0) {
	    s_wsle(&io___25);
	    for (i__ = 1; i__ <= 2; ++i__) {
		do_lio(&c__6, &c__1, (char *)&cval[i__ - 1], (ftnlen)sizeof(
			complex));
	    }
	    e_wsle();
	}
    } else {
	s_wsle(&io___26);
	do_lio(&c__9, &c__1, "  SPARSE (COORDINATE) ARRAY.", (ftnlen)28);
	e_wsle();
	s_wsle(&io___27);
	do_lio(&c__9, &c__1, "  FIRST TWO ENTRIES:", (ftnlen)20);
	e_wsle();
	if (s_cmp(field, "INTEGER", (ftnlen)7, (ftnlen)7) == 0) {
	    s_wsle(&io___28);
	    do_lio(&c__3, &c__1, (char *)&indrow[0], (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&indcol[0], (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&ival[0], (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___29);
	    do_lio(&c__3, &c__1, (char *)&indrow[1], (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&indcol[1], (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&ival[1], (ftnlen)sizeof(integer));
	    e_wsle();
	} else if (s_cmp(field, "REAL", (ftnlen)7, (ftnlen)4) == 0) {
	    s_wsle(&io___30);
	    do_lio(&c__3, &c__1, (char *)&indrow[0], (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&indcol[0], (ftnlen)sizeof(integer));
	    do_lio(&c__5, &c__1, (char *)&rval[0], (ftnlen)sizeof(doublereal))
		    ;
	    e_wsle();
	    s_wsle(&io___31);
	    do_lio(&c__3, &c__1, (char *)&indrow[1], (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&indcol[1], (ftnlen)sizeof(integer));
	    do_lio(&c__5, &c__1, (char *)&rval[1], (ftnlen)sizeof(doublereal))
		    ;
	    e_wsle();
	} else if (s_cmp(field, "COMPLEX", (ftnlen)7, (ftnlen)7) == 0) {
	    s_wsle(&io___32);
	    do_lio(&c__3, &c__1, (char *)&indrow[0], (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&indcol[0], (ftnlen)sizeof(integer));
	    do_lio(&c__6, &c__1, (char *)&cval[0], (ftnlen)sizeof(complex));
	    e_wsle();
	    s_wsle(&io___33);
	    do_lio(&c__3, &c__1, (char *)&indrow[1], (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&indcol[1], (ftnlen)sizeof(integer));
	    do_lio(&c__6, &c__1, (char *)&cval[1], (ftnlen)sizeof(complex));
	    e_wsle();
	} else if (s_cmp(field, "PATTERN", (ftnlen)7, (ftnlen)7) == 0) {
	    s_wsle(&io___34);
	    do_lio(&c__3, &c__1, (char *)&indrow[0], (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&indcol[0], (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___35);
	    do_lio(&c__3, &c__1, (char *)&indrow[1], (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&indcol[1], (ftnlen)sizeof(integer));
	    e_wsle();
	}
    }
    s_wsle(&io___36);
    do_lio(&c__9, &c__1, "BEFORE RUNNING", (ftnlen)14);
    e_wsle();
    i__1 = nnz;
    for (k = 1; k <= i__1; ++k) {
	s_wsle(&io___38);
	do_lio(&c__9, &c__1, "(", (ftnlen)1);
	do_lio(&c__3, &c__1, (char *)&indrow[k - 1], (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, ",", (ftnlen)1);
	do_lio(&c__3, &c__1, (char *)&indcol[k - 1], (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, ")", (ftnlen)1);
	e_wsle();
/* L10: */
    }
    srtdat_(&ncols, &nnz, indrow, indcol, jpntr, iwa);
    s_wsle(&io___41);
    do_lio(&c__9, &c__1, "WHERE IS THE SEGMENTATION ERROR", (ftnlen)31);
    e_wsle();
    setr_(&nrows, &ncols, indrow, jpntr, indcol, ipntr, iwa);
    s_wsle(&io___43);
    do_lio(&c__9, &c__1, "AFTER RUNNING", (ftnlen)13);
    e_wsle();
    i__1 = nnz;
    for (k = 1; k <= i__1; ++k) {
	s_wsle(&io___44);
	do_lio(&c__9, &c__1, "(", (ftnlen)1);
	do_lio(&c__3, &c__1, (char *)&indrow[k - 1], (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, ",", (ftnlen)1);
	do_lio(&c__3, &c__1, (char *)&indcol[k - 1], (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, ")", (ftnlen)1);
	e_wsle();
/* L110: */
    }
    mingrp = 0;
    s_wsle(&io___46);
    do_lio(&c__9, &c__1, "MINGRP ", (ftnlen)7);
    do_lio(&c__3, &c__1, (char *)&mingrp, (ftnlen)sizeof(integer));
    e_wsle();
    for (i__ = 1; i__ <= 263743; ++i__) {
/* Computing MAX */
	i__1 = mingrp, i__2 = ipntr[i__] - ipntr[i__ - 1];
	mingrp = max(i__1,i__2);
/* L20: */
    }
    s_wsle(&io___47);
    do_lio(&c__9, &c__1, "MINGRP ", (ftnlen)7);
    do_lio(&c__3, &c__1, (char *)&mingrp, (ftnlen)sizeof(integer));
    e_wsle();
    s_wsle(&io___48);
    do_lio(&c__3, &c_b197, (char *)&ipntr[0], (ftnlen)sizeof(integer));
    e_wsle();
    degr_(&ncols, indrow, jpntr, indcol, ipntr, &iwa[1318715], &iwa[263743]);
    for (i__ = 1; i__ <= 263743; ++i__) {
	s_wsle(&io___49);
	do_lio(&c__9, &c__1, "DEG(", (ftnlen)4);
	do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, ") : ", (ftnlen)4);
	do_lio(&c__3, &c__1, (char *)&iwa[i__ + 1318714], (ftnlen)sizeof(
		integer));
	e_wsle();
/* L30: */
    }
/*      CALL SLO(NCOLS,INDROW,JPNTR,INDCOL,IPNTR,IWA(5*N+1),IWA(4*N+1), */
/*     *     MAXCLQ,IWA(1),IWA(N+1),IWA(2*N+1), IWA(3*N+1)) */
/*     CALL IDO IN THE ALTERNATIVE PROGRAM. */
    ido_(&nrows, &ncols, indrow, jpntr, indcol, ipntr, &iwa[1318715], &iwa[
	    1054972], &maxclq, iwa, &iwa[263743], &iwa[527486], &iwa[791229]);
    for (i__ = 1; i__ <= 263743; ++i__) {
	s_wsle(&io___51);
	do_lio(&c__9, &c__1, "LIST;", (ftnlen)5);
	do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, ";", (ftnlen)1);
	do_lio(&c__3, &c__1, (char *)&iwa[i__ + 1054971], (ftnlen)sizeof(
		integer));
	e_wsle();
/* L40: */
    }
    if (iargc_() > 1) {
	getarg_(&c__2, ofile, (ftnlen)32);
	ounit__ = 9;
	o__1.oerr = 0;
	o__1.ounit = ounit__;
	o__1.ofnmlen = 32;
	o__1.ofnm = ofile;
	o__1.orl = 0;
	o__1.osta = 0;
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	s_wsle(&io___54);
	do_lio(&c__9, &c__1, "WRITING HEADER AND DATA...", (ftnlen)26);
	e_wsle();
	mmwrite_(&ounit__, rep, field, symm, &nrows, &ncols, &nnz, indrow, 
		indcol, ival, rval, cval, (ftnlen)10, (ftnlen)7, (ftnlen)19);
	cl__1.cerr = 0;
	cl__1.cunit = 9;
	cl__1.csta = 0;
	f_clos(&cl__1);
	s_wsle(&io___55);
	do_lio(&c__9, &c__1, "DONE WRITING TO FILE: ", (ftnlen)22);
	do_lio(&c__9, &c__1, ofile, (ftnlen)32);
	e_wsle();
    }
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Main program alias */ int main_ () { MAIN__ (); return 0; }
