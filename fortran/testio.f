      PROGRAM MAIN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C USAGE:  TESTIO < DATAFILE.MTX
C         TESTIO DATAFILE.MTX
C         GUNZIP -C DATAFILE.MTX.Z | TESTIO
C
C NOTE:   MAKE 'TESTIO' EXECUTABLE WITH:   F77 -O TESTIO TESTIO.F MMIO.F
C
C THIS SAMPLE DRIVER TAKES A MATRIX MARKET FILE AS STANDARD INPUT,
C OR OPENS A THE SPECIFIED FILE IF A FILENAME ARGUMENT IS PRESENT,
C AND CALLS TWO ROUTINES:
C
C     CALL MMINFO(IUNIT,REP,FIELD,SYMM,NROWS,NCOLS,NNZ)
C AND 
C     CALL MMREAD(IUNIT,REP,FIELD,SYMM,NROWS,NCOLS,NNZ,NNZMAX,
C    *             INDX,INDCOL,IVAL,RVAL,CVAL)
C 
C SUBROUTINE MMINFO JUST PARSES HEADER INFORMATION, WHILE MMREAD PARSES 
C THE HEADER AND READS THE NUMERICAL DATA INTO THE APPROPRIATE ARRAYS.
C EACH OF THESE REQUIRES A PRE-OPENED READ UNIT (IUNIT), AND
C REWINDS THE UNIT PRIOR TO RETURN. (SO MMINFON DOESN'T INTERFERE
C WITH MMREAD AND VICE VERSA).
C
C 18-OCT-96  KARIN A. REMINGTON, NIST ACMD (KARIN@CAM.NIST.GOV)
C 30-OCT-96  MINOR CHANGE TO CALLING SEQUENCES TO ACCOMMODATE NEW
C              INTEGER VALUE PARAMETER (IVAL)
C            + ASSOCIATED CHANGES TO CHECK FOR INTEGER VALUE DATA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      INTEGER Z_PARAM
      INTEGER N_PARAM
      INTEGER M_PARAM

      PARAMETER(N_PARAM = 263743)
      PARAMETER(M_PARAM = 263743)
      PARAMETER (Z_PARAM=1302464)

      INTEGER IUNIT,OUNIT
      CHARACTER REP*10
      CHARACTER FIELD*7
      CHARACTER SYMM*19
      CHARACTER IFILE*32,OFILE*32
      INTEGER IVAL(Z_PARAM)
      DOUBLE PRECISION RVAL(Z_PARAM)
      COMPLEX CVAL(Z_PARAM)


      INTEGER INDROW(Z_PARAM)
      INTEGER INDCOL(Z_PARAM)
      INTEGER IPNTR(N_PARAM+1)
      INTEGER JPNTR(N_PARAM+1)
      INTEGER IWA(6*N_PARAM+1)
      INTEGER MAXCLQ

C
      IF (IARGC() .GT. 0) THEN
         CALL GETARG(1,IFILE)
         IUNIT = 8
         OPEN(UNIT=IUNIT,FILE=IFILE)
      ELSE
         IUNIT = 5
      ENDIF
      PRINT *,'READING HEADER ONLY...'
      CALL MMINFO(IUNIT,REP,FIELD,SYMM,NROWS,NCOLS,NNZ)
      PRINT *,'  MATRIX IS TYPE: ',REP,' ',FIELD,' ',SYMM
      PRINT *,'  MATRIX SIZE: ',NROWS,' BY ',NCOLS,' WITH ',
     *                          NNZ,' NONZEROS.'
C
      PRINT *,'READING HEADER AND DATA...'
      CALL MMREAD(IUNIT,REP,FIELD,SYMM,NROWS,NCOLS,NNZ,Z_PARAM,
     *           INDROW,INDCOL,IVAL,RVAL,CVAL)
      PRINT *,'  MATRIX IS TYPE: ',REP,' ',FIELD,' ',SYMM
      PRINT *,'  MATRIX SIZE: ',NROWS,' BY ',NCOLS,' WITH ',
     *                          NNZ,' NONZEROS.'
C
      IF( REP .EQ. 'ARRAY' ) THEN
        PRINT *,'  DENSE ARRAY.'
        PRINT *,'  FIRST TWO ENTRIES:'
        IF ( FIELD .EQ. 'INTEGER' ) THEN
           PRINT *,(IVAL(I),I=1,2)
        ELSEIF ( FIELD .EQ. 'REAL') THEN
           PRINT *,(RVAL(I),I=1,2)
        ELSEIF ( FIELD .EQ. 'COMPLEX' ) THEN
           PRINT *,(CVAL(I),I=1,2)
        ENDIF
      ELSE
        PRINT *,'  SPARSE (COORDINATE) ARRAY.'
        PRINT *,'  FIRST TWO ENTRIES:'
        IF ( FIELD .EQ. 'INTEGER' ) THEN
          PRINT *,INDROW(1),INDCOL(1),IVAL(1)
          PRINT *,INDROW(2),INDCOL(2),IVAL(2)
        ELSEIF ( FIELD .EQ. 'REAL') THEN
          PRINT *,INDROW(1),INDCOL(1),RVAL(1)
          PRINT *,INDROW(2),INDCOL(2),RVAL(2)
        ELSEIF ( FIELD .EQ. 'COMPLEX' ) THEN
          PRINT *,INDROW(1),INDCOL(1),CVAL(1)
          PRINT *,INDROW(2),INDCOL(2),CVAL(2)
        ELSEIF ( FIELD .EQ. 'PATTERN' ) THEN
          PRINT *,INDROW(1),INDCOL(1)
          PRINT *,INDROW(2),INDCOL(2)
        ENDIF
      ENDIF

      CALL SRTDAT(NCOLS,NNZ,INDROW,INDCOL,JPNTR,IWA)
      CALL SETR(NROWS,NCOLS,INDROW,JPNTR,INDCOL,IPNTR,IWA)

      CALL DEGR(NCOLS,INDROW,JPNTR,INDCOL,IPNTR,IWA(5*NCOLS+1)
     *     ,IWA(NCOLS+1))



c     5*N +1 NDEG 
c     4*N + 1 list


      CALL SLO(NCOLS,INDROW,JPNTR,INDCOL,IPNTR,IWA(5*NCOLS+1)
     $  ,IWA(4*NCOLS+1), MAXCLQ,IWA(1)
     $     ,IWA(NCOLS+1),IWA(2*NCOLS+1), IWA(3*NCOLS+1))
c     CALL IDO IN THE ALTERNATIVE PROGRAM.
      CALL IDO(NROWS,NCOLS,INDROW,JPNTR,INDCOL,IPNTR,IWA(5*NCOLS+1),
     *     IWA(4*NCOLS+1),MAXCLQ,IWA(1),IWA(NCOLS+1),IWA(2*NCOLS+1),
     *     IWA(3*NCOLS+1))

      DO 40 I = 1,NCOLS
         PRINT *,'LIST;',I,';', IWA(4*NCOLS+I)
 40   CONTINUE
c$$$      IF (IARGC() .GT. 1) THEN
c$$$         CALL GETARG(2,OFILE)
c$$$         OUNIT = 9
c$$$         OPEN(UNIT=OUNIT,FILE=OFILE)
c$$$         PRINT *,'WRITING HEADER AND DATA...'
c$$$         CALL MMWRITE(OUNIT,REP,FIELD,SYMM,NROWS,NCOLS,NNZ,
c$$$     *                 INDROW,INDCOL,IVAL,RVAL,CVAL)
c$$$         CLOSE(9)
c$$$         PRINT *,'DONE WRITING TO FILE: ',OFILE
c$$$      ENDIF
      STOP
      END

