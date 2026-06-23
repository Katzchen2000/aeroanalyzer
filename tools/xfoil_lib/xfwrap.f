C=====================================================================72
C     xfwrap.f  --  minimal C-callable wrapper around the XFOIL 6.99
C     analysis core.  Drives the COMMON-block state in the same call
C     sequence the interactive OPER menu uses, so we get XFOIL-exact
C     viscous polars in-process with zero disk I/O.
C
C     Entry points (iso_c_binding, lowercase C names):
C        xfl_init()                         one-time defaults
C        xfl_set_airfoil(x,y,n)             load + panel a buffer airfoil
C        xfl_set_cond(re,ncrit,mach,xt,xb)  viscous operating conditions
C        xfl_solve(adeg, cl,cd,cm,conv)     one viscous point at alpha
C
C     NOT thread-safe (XFOIL uses global COMMON); serialize calls.
C=====================================================================72

      SUBROUTINE XFL_INIT() BIND(C, NAME="xfl_init")
      INCLUDE 'XFOIL.INC'
      CALL INIT
      LVISC  = .TRUE.
      LVCONV = .FALSE.
      RETYP  = 1
      MATYP  = 1
      MINF1  = 0.0
      ITMAX  = 100
      RETURN
      END

      SUBROUTINE XFL_SET_AIRFOIL(XIN, YIN, NIN)
     &           BIND(C, NAME="xfl_set_airfoil")
      USE ISO_C_BINDING
      INCLUDE 'XFOIL.INC'
      INTEGER(C_INT), VALUE :: NIN
      REAL(C_DOUBLE) :: XIN(NIN), YIN(NIN)
      INTEGER I
C---- copy coordinates into the buffer airfoil
      NB = NIN
      DO I = 1, NB
        XB(I) = XIN(I)
        YB(I) = YIN(I)
      ENDDO
C---- spline + geometry parameters (mirrors LOAD), then generate paneling
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB, W1,
     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKB,CAMBRB )
      CALL PANGEN(.FALSE.)
C---- a fresh geometry invalidates ALL cached solution state
      LGAMU  = .FALSE.
      LQAIJ  = .FALSE.
      LADIJ  = .FALSE.
      LWDIJ  = .FALSE.
      LIPAN  = .FALSE.
      LWAKE  = .FALSE.
      LBLINI = .FALSE.
      LVCONV = .FALSE.
      RETURN
      END

      SUBROUTINE XFL_SET_COND(RE, NCRIT, MACH, XTRT, XTRB)
     &           BIND(C, NAME="xfl_set_cond")
      USE ISO_C_BINDING
      INCLUDE 'XFOIL.INC'
      REAL(C_DOUBLE), VALUE :: RE, NCRIT, MACH, XTRT, XTRB
      LVISC = .TRUE.
      REINF1 = RE
      MINF1  = MACH
      ACRIT  = NCRIT
      XSTRIP(1) = XTRT
      XSTRIP(2) = XTRB
      CALL MRCL(1.0,MINF_CL,REINF_CL)
      CALL COMSET
C---- conditions changed: drop the warm BL so the next solve re-initializes
C     from scratch (a stale BL from a different Re/Ncrit can diverge to NaN)
      LBLINI = .FALSE.
      LVCONV = .FALSE.
      RETURN
      END

C---- force the next solve to cold-start its boundary layer (call between
C     unrelated operating points to avoid stale-warm-start divergence)
      SUBROUTINE XFL_RESET_BL() BIND(C, NAME="xfl_reset_bl")
      INCLUDE 'XFOIL.INC'
      LBLINI = .FALSE.
      LVCONV = .FALSE.
      RETURN
      END

      SUBROUTINE XFL_SOLVE(ADEGIN, CLO, CDO, CMO, ICONV)
     &           BIND(C, NAME="xfl_solve")
      USE ISO_C_BINDING
      INCLUDE 'XFOIL.INC'
      REAL(C_DOUBLE), VALUE :: ADEGIN
      REAL(C_DOUBLE) :: CLO, CDO, CMO
      INTEGER(C_INT) :: ICONV
      LALFA = .TRUE.
      ADEG  = ADEGIN
      ALFA  = DTOR*ADEG
      QINF  = 1.0
      CALL SPECAL
      LWAKE  = .FALSE.
      LVCONV = .FALSE.
      IF(LVISC) CALL VISCAL(ITMAX)
      CLO = CL
      CDO = CD
      CMO = CM
      IF(LVCONV) THEN
        ICONV = 1
      ELSE
        ICONV = 0
      ENDIF
      RETURN
      END
