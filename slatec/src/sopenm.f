*DECK SOPENM
      SUBROUTINE SOPENM (IPAGE, LPAGE)
C***BEGIN PROLOGUE  SOPENM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SPLP
C***LIBRARY   SLATEC
C***TYPE      ALL (SOPENM-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     1. OPEN UNIT NUMBER IPAGEF AS A RANDOM ACCESS FILE.
C
C     2. THE RECORD LENGTH IS CONSTANT=LPG.
C
C***SEE ALSO  SPLP
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   811215  DATE WRITTEN
C   890605  Corrected references to XERRWV.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C***END PROLOGUE  SOPENM
      CHARACTER*8 XERN1
C
C***FIRST EXECUTABLE STATEMENT  SOPENM
      IPAGEF=IPAGE
      LPG   =LPAGE
      OPEN(UNIT=IPAGEF,IOSTAT=IOS,ERR=100,STATUS='UNKNOWN',
     *ACCESS='DIRECT',FORM='UNFORMATTED',RECL=LPG)
      RETURN
C
 100  WRITE (XERN1, '(I8)') IOS
      CALL XERMSG ('SLATEC', 'SOPENM',
     *   'IN SPLP, OPEN HAS ERROR FLAG = ' // XERN1, 100, 1)
      RETURN
      END
