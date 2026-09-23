!***********************************************************************
      SUBROUTINE SPLINE(N,X,Y,ALFA,BETA,TYP,A,B,C,D)
!-----------------------------------------------------------------------
!     INPUT:
!
!     N     NUMBER OF POINTS
!     X     ARRAY X VECTOR
!     Y     ARRAY Y VECTOR
!     ALFA  BOUNDARY CONDITION IN X(1)
!     BETA        "       IN X(N)
!     TYP   =  0  NOT-A-KNOT SPLINE
!              1  ALFA, BETA 1. ABLEITUNGEN VORGEGEBEN
!              2    "    "   2.     "           "
!              3    "    "   3.     "           "
!
!     BEMERKUNG: MIT TYP = 2 UND ALFA = BETA = 0 ERHAELT MAN
!           EINEN NATUERLICHEN SPLINE
!
!     OUTPUT:
!
!     A, B, C, D     ARRAYS OF SPLINE COEFFICIENTS
!       S = A(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2+ D(I)*(X-X(I))**3
!
!     BEI ANWENDUNGSFEHLERN WIRD DAS PROGRAMM MIT ENTSPRECHENDER
!     FEHLERMELDUNG ABGEBROCHEN
!-----------------------------------------------------------------------
!
!
      IMPLICIT NONE

      INTEGER  N, TYP
      REAL*8   X(N), Y(N), ALFA, BETA, A(N), B(N), C(N), D(N)
      INTEGER  I, IERR
      REAL*8   H(N)

      IF((TYP.LT.0).OR.(TYP.GT.3)) THEN
         WRITE(*,*) 'ERROR IN ROUTINE SPLINE: FALSE TYP'
         STOP
      ENDIF

      IF (N.LT.3) THEN
         WRITE(*,*) 'ERROR IN ROUTINE  SPLINE: N < 3'
         STOP
      ENDIF


!     BERECHNE DIFFERENZ AUFEINENDERFOLGENDER X-WERTE UND
!     UNTERSUCHE MONOTONIE
!
      DO I = 1, N-1
         H(I) = X(I+1)- X(I)
         IF(H(I).LE.0.d0) THEN
            WRITE(*,*) 'NON MONOTONIC ABCISSAE IN SPLINE: X(I-1)>=X(I)'
            STOP
         ENDIF
      ENDDO
!
!     AUFSTELLEN DES GLEICHUNGSSYSTEMS
!
      !DIR$ NOVECTOR
      DO 20 I = 1, N-2
         A(I) = 3.d0 * ((Y(I+2)-Y(I+1)) / H(I+1) - (Y(I+1)-Y(I)) / H(I))
         B(I) = H(I)
         C(I) = H(I+1)
         D(I) = 2.d0 * (H(I) + H(I+1))
   20 CONTINUE
!
!     BERUECKSICHTIGEN DER RANDBEDINGUNGEN

!     NOT-A-KNOT

      IF(TYP.EQ.0) THEN
         A(1)   = A(1) * H(2) / (H(1) + H(2))
         A(N-2) = A(N-2) * H(N-2) / (H(N-1) + H(N-2))
         D(1)   = D(1) - H(1)
         D(N-2) = D(N-2) - H(N-1)
         C(1)   = C(1) - H(1)
         B(N-2) = B(N-2) - H(N-1)
      ENDIF

!     1. ABLEITUNG VORGEGEBEN

      IF(TYP.EQ.1) THEN
         A(1)   = A(1) - 1.5d0 * ((Y(2)-Y(1)) / H(1) - ALFA)
         A(N-2) = A(N-2) - 1.5d0 * (BETA - (Y(N)-Y(N-1)) / H(N-1))
         D(1)   = D(1) - 0.5d0 * H(1)
         D(N-2) = D(N-2) - 0.5d0 * H(N-1)
      ENDIF
!
!     2. ABLEITUNG VORGEGEBEN
!
      IF(TYP.EQ.2) THEN
         A(1)   = A(1) - 0.5d0 * ALFA * H(1)
         A(N-2) = A(N-2) - 0.5d0 * BETA * H(N-1)
      ENDIF

!     3. ABLEITUNG VORGEGEBEN
!
      IF(TYP.EQ.3 ) THEN
         A(1)   = A(1) + 0.5d0 * ALFA * H(1) * H(1)
         A(N-2) = A(N-2) - 0.5d0 * BETA * H(N-1)* H(N-1)
         D(1)   = D(1) + H(1)
         D(N-2) = D(N-2) + H(N-1)
      ENDIF

!     BERECHNUNG DER KOEFFIZIENTEN
!
      CALL DGTSL(N-2,B,D,C,A,IERR)
      IF(IERR.NE.0) THEN
         WRITE(*,21)
         STOP
      ENDIF

!     UEBERSCHREIBEN DES LOESUNGSVEKTORS

      CALL DCOPY(N-2,A,1,C(2),1)
!
!     IN ABHAENGIGKEIT VON DEN RANDBEDINGUNGEN WIRD DER 1. UND
!     DER LETZTE WERT VON C KORRIGIERT
!
      IF(TYP.EQ.0) THEN
         C(1) = C(2) + H(1) * (C(2)-C(3)) / H(2)
         C(N) = C(N-1) + H(N-1) * (C(N-1)-C(N-2)) / H(N-2)
      ENDIF

      IF(TYP.EQ.1) THEN
         C(1) = 1.5d0*((Y(2)-Y(1)) / H(1) - ALFA) / H(1) - 0.5d0 * C(2)
         C(N) = -1.5d0*((Y(N)-Y(N-1)) / H(N-1)-BETA) / H(N-1) - 0.5d0*C(N-1)
      ENDIF

      IF(TYP.EQ.2) THEN
         C(1) = 0.5d0 * ALFA
         C(N) = 0.5d0 * BETA
      ENDIF

      IF(TYP.EQ.3) THEN
         C(1) = C(2) - 0.5d0 * ALFA * H(1)
         C(N) = C(N-1) + 0.5d0 * BETA * H(N-1)
      ENDIF

      CALL DCOPY(N,Y,1,A,1)

      DO I = 1, N-1
         B(I) = (A(I+1)-A(I)) / H(I) - H(I) * (C(I+1)+2.0d0 * C(I)) / 3.0d0
         D(I) = (C(I+1)-C(I)) / (3.0d0 * H(I))
      END DO

      B(N) = (3.0d0 * D(N-1) * H(N-1) + 2.0d0 * C(N-1)) * H(N-1) + B(N-1)

      RETURN

   21 FORMAT(1X,'ERROR IN SGTSL: MATRIX SINGULAR')
      END


!************************************************************************
      REAL FUNCTION SPWERT(N,XWERT,A,B,C,D,X,ABLTG)
!-----------------------------------------------------------------------
!     INPUT:
!
!     N           NUMBER OF GRID POINTS
!     XWERT       STELLE AN DER FUNKTIONSWERTE BERECHNET WERDEN
!    A, B, C, D  ARRAYS DER SPLINEKOEFFIZIENTEN (AUS SPLINE)
!     X           ARRAY DER KNOTENPUNKTE
!
!     OUTPUT:
!
!     SPWERT   FUNKTIONSWERT AN DER STELLE XWERT
!     ABLTG(I) I=1 : FIRST DERIVATIVE, ETC.
!-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER  N
      REAL*8   XWERT, A(N), B(N), C(N), D(N), X(N), ABLTG(3), XX
      INTEGER  I, K, M

!     SUCHE PASSENDES INTERVALL (BINAERE SUCHE)

      I = 1
      K = N

   10 M = (I+K) / 2

      IF(M.NE.I) THEN
         IF(XWERT.LT.X(M)) THEN
            K = M
         ELSE
            I = M
         ENDIF
         GOTO 10
      ENDIF

      XX = XWERT - X(I)

      ABLTG(1) = (3.0d0 * D(I) * XX + 2.0d0 * C(I)) * XX + B(I)
      ABLTG(2) = 6.0d0 * D(I) * XX + 2.0d0 * C(I)
      ABLTG(3) = 6.0d0 * D(I)

      SPWERT = ((D(I)*XX + C(I))*XX + B(I))*XX + A(I)

      RETURN
      END

!DECK DGTSL                                                             CAS02750
!** FROM NETLIB, TUE AUG 28 08:28:34 EDT 1990 ***
!** COPIED FROM SGTSL AND RENAMED             ***
!
      SUBROUTINE DGTSL(N,C,D,E,B,INFO)
      IMPLICIT NONE
      INTEGER N,INFO
      REAL*8 C(*),D(*),E(*),B(*)

!     SGTSL GIVEN A GENERAL TRIDIAGONAL MATRIX AND A RIGHT HAND
!     SIDE WILL FIND THE SOLUTION.
!
!     ON ENTRY
!
!        N       INTEGER
!                IS THE ORDER OF THE TRIDIAGONAL MATRIX.
!
!
!        C       REAL(N)
!                IS THE SUBDIAGONAL OF THE TRIDIAGONAL MATRIX.
!                C(2) THROUGH C(N) SHOULD CONTAIN THE SUBDIAGONAL.
!                ON OUTPUT C IS DESTROYED.
!
!        D       REAL(N)
!                IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.
!                ON OUTPUT D IS DESTROYED.
!
!        E       REAL(N)
!                IS THE SUPERDIAGONAL OF THE TRIDIAGONAL MATRIX.
!                E(1) THROUGH E(N-1) SHOULD CONTAIN THE SUPERDIAGONAL.
!                ON OUTPUT E IS DESTROYED.
!
!        B       REAL(N)
!                IS THE RIGHT HAND SIDE VECTOR.
!
!     ON RETURN
!
!        B       IS THE SOLUTION VECTOR.
!
!        INFO    INTEGER
!                = 0 NORMAL VALUE.
!                = K IF THE K-TH ELEMENT OF THE DIAGONAL BECOMES
!                    EXACTLY ZERO.  THE SUBROUTINE RETURNS WHEN
!                    THIS IS DETECTED.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
!
!     NO EXTERNALS
!     FORTRAN ABS
!
!     INTERNAL VARIABLES
!
      INTEGER K,KB,KP1,NM1,NM2
      REAL*8 T
!     BEGIN BLOCK PERMITTING ...EXITS TO 100

         INFO = 0
         C(1) = D(1)
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 40
            D(1) = E(1)
            E(1) = 0.0D0
            E(N) = 0.0D0

            DO 30 K = 1, NM1
               KP1 = K + 1

!              FIND THE LARGEST OF THE TWO ROWS

               IF (ABS(C(KP1)) .LT. ABS(C(K))) GO TO 10

!                 INTERCHANGE ROW

                  T = C(KP1)
                  C(KP1) = C(K)
                  C(K) = T
                  T = D(KP1)
                  D(KP1) = D(K)
                  D(K) = T
                  T = E(KP1)
                  E(KP1) = E(K)
                  E(K) = T
                  T = B(KP1)
                  B(KP1) = B(K)
                  B(K) = T
   10          CONTINUE

!              ZERO ELEMENTS

               IF (C(K) .NE. 0.0D0) GO TO 20
                  INFO = K
!     ............EXIT
                  GO TO 100
   20          CONTINUE
               T = -C(KP1)/C(K)
               C(KP1) = D(KP1) + T*D(K)
               D(KP1) = E(KP1) + T*E(K)
               E(KP1) = 0.0D0
               B(KP1) = B(KP1) + T*B(K)
   30       CONTINUE
   40    CONTINUE
         IF (C(N) .NE. 0.0D0) GO TO 50
            INFO = N
         GO TO 90
   50    CONTINUE

!           BACK SOLVE

            NM2 = N - 2
            B(N) = B(N)/C(N)
            IF (N .EQ. 1) GO TO 80
               B(NM1) = (B(NM1) - D(NM1)*B(N))/C(NM1)
               IF (NM2 .LT. 1) GO TO 70
               DO 60 KB = 1, NM2
                  K = NM2 - KB + 1
                  B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2))/C(K)
   60          CONTINUE
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE

      RETURN
      END

