/*****************************************************************************
**
**  FUNCTION1
**  ---------
**
**  Purpose : Compute the sume of squared deviations from the put-call parity
**
**  Globals : PARAM1  -- Parameters
**            Y       -- Call price vector - put price vector
**            X       -- Properly stacked strike price matrix
**            XC      -- Properly stacked matrix of ones.
******************************************************************************/

  PROC (1)=FUNCT1(PARAM1);

  @ Local variables @

    LOCAL PV,I,Y1,X1,D,RES,SSE;

  @ Parameter specification @

    PARAM1=PARAM1./SF1;

    PV=ZEROS(ROWS(PARAM1)+1,1);
    I=1;
      DO WHILE I<=ROWS(PARAM1);
      PV[I+1,1]=PV[I,1]+ABS(PARAM1[I,1]);
      I=I+1;
    ENDO;

    Y1=Y+XC*PV;
    X1=ONES(ROWS(X),1)~X;
    D=INVPD(X1'X1)*(X1'Y1);
    RES=Y1-X1*D ;
    SSE=RES'RES;

  RETP(-SSE);

ENDP;
