C=====================================================================72
C     xfstubs.f -- no-op replacements for the plot / GUI / menu / polar /
C     inverse-design symbols that the XFOIL 6.99 analysis core references
C     (through PROGRAM XFOIL and the OPER/GDES/MDES/QDES menus) but never
C     actually reaches when driven headlessly through xfwrap.f.
C
C     Compiling these as empty subroutines lets us link a libxfoil.a that
C     contains ONLY the analysis path (panel + boundary layer + e^n), with
C     no dependency on the Xplot11 plot library or the interactive menus.
C
C     This set was verified sufficient: the analysis source set
C       aread naca sort spline userio xbl xblsys xfoil xgeom xoper
C       xpanel xsolve xutils
C     links cleanly against it with no undefined references. If a future
C     source addition pulls in a new plot/menu symbol the linker will name
C     it; add a matching empty stub here.
C=====================================================================72
      SUBROUTINE abcopy
      RETURN
      END
      SUBROUTINE airlim
      RETURN
      END
      SUBROUTINE annot
      RETURN
      END
      SUBROUTINE apcopy
      RETURN
      END
      SUBROUTINE axisadj
      RETURN
      END
      SUBROUTINE blplot
      RETURN
      END
      SUBROUTINE clrzoom
      RETURN
      END
      SUBROUTINE colorspectrumhues
      RETURN
      END
      SUBROUTINE cpaxes
      RETURN
      END
      SUBROUTINE cpvec
      RETURN
      END
      SUBROUTINE cpx
      RETURN
      END
      SUBROUTINE dash
      RETURN
      END
      SUBROUTINE dplot
      RETURN
      END
      SUBROUTINE gdes
      RETURN
      END
      SUBROUTINE getcolor
      RETURN
      END
      SUBROUTINE getxyf
      RETURN
      END
      SUBROUTINE mdes
      RETURN
      END
      SUBROUTINE newcolor
      RETURN
      END
      SUBROUTINE newcolorname
      RETURN
      END
      SUBROUTINE newfactor
      RETURN
      END
      SUBROUTINE newpen
      RETURN
      END
      SUBROUTINE oplset
      RETURN
      END
      SUBROUTINE panplt
      RETURN
      END
      SUBROUTINE plchar
      RETURN
      END
      SUBROUTINE plclose
      RETURN
      END
      SUBROUTINE plend
      RETURN
      END
      SUBROUTINE plflush
      RETURN
      END
      SUBROUTINE plgrid
      RETURN
      END
      SUBROUTINE plinitialize
      RETURN
      END
      SUBROUTINE plnumb
      RETURN
      END
      SUBROUTINE plot
      RETURN
      END
      SUBROUTINE plradd
      RETURN
      END
      SUBROUTINE plrcop
      RETURN
      END
      SUBROUTINE plrini
      RETURN
      END
      SUBROUTINE plrset
      RETURN
      END
      SUBROUTINE plrsrt
      RETURN
      END
      SUBROUTINE plrsum
      RETURN
      END
      SUBROUTINE plsymb
      RETURN
      END
      SUBROUTINE pltini
      RETURN
      END
      SUBROUTINE plxadd
      RETURN
      END
      SUBROUTINE plxini
      RETURN
      END
      SUBROUTINE polaxi
      RETURN
      END
      SUBROUTINE polplt
      RETURN
      END
      SUBROUTINE polread
      RETURN
      END
      SUBROUTINE polref
      RETURN
      END
      SUBROUTINE polwrit
      RETURN
      END
      SUBROUTINE ppaplt
      RETURN
      END
      SUBROUTINE prfcop
      RETURN
      END
      SUBROUTINE prfsum
      RETURN
      END
      SUBROUTINE qdes
      RETURN
      END
      SUBROUTINE replot
      RETURN
      END
      SUBROUTINE seqlab
      RETURN
      END
      SUBROUTINE seqplt
      RETURN
      END
      SUBROUTINE usetzoom
      RETURN
      END
      SUBROUTINE xaxis
      RETURN
      END
      SUBROUTINE xyline
      RETURN
      END
      SUBROUTINE yaxis
      RETURN
      END
