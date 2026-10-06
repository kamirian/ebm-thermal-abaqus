C ======================================================================
C  DFLUX - moving Gaussian electron-beam heat source (Abaqus/Standard)
C
C  Body heat flux deposited in the powder layer by an electron beam
C  that scans along +x on the top surface (z = Z0), on the symmetry
C  plane y = 0 of a half model.
C
C  Units: mm, s, mW (mm-tonne-s system). FLUX(1) is in mW/mm^3.
C
C  Surface intensity (Gaussian of size PHI):
C     HS = 2*U*IB/(PI*PHI^2) * exp(-2*((x-x0)^2 + (y-y0)^2)/PHI^2)
C  Depth profile over the penetration depth S (integrates to S):
C     IZ = (-2.25*(d/S)^2 + 1.5*(d/S) + 0.75)/0.75,   d = Z0 - z
C  Absorbed body flux in the powder layer (0 <= d <= S):
C     FLUX(1) = ETA*HS*IZ/S,   and 0 below it.
C ======================================================================
      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,
     1                 JLTYP,TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION FLUX(2), TIME(2), COORDS(3)
      CHARACTER*80 SNAME
C
C     PI     pi
C     PHI    Gaussian beam size (Phi_E in An et al., 2021), mm
C     ETA    absorption efficiency
C     U      accelerating voltage, V
C     IB     beam current, mA (U*IB is the beam power in mW)
C     V      scan speed, mm/s
C     S      penetration depth of the beam, mm
C     Z0     z of the top surface, mm
C     X00    start of the beam centre (PHI before the part), mm
      DOUBLE PRECISION PI, PHI, ETA, U, IB, V, S, Z0, X00
      PARAMETER (PI = 3.1415926D0, PHI = 0.55D0, ETA = 0.9D0,
     1           U = 60.0D3, IB = 6.7D0, V = 632.6D0,
     2           S = 0.062D0, Z0 = 10.07D0, X00 = -PHI)
C
      DOUBLE PRECISION X0, Y0, D, HS, IZ
C
C     Beam centre at the current total time
      X0 = X00 + V*TIME(2)
      Y0 = 0.0D0
C
C     Depth below the top surface
      D = Z0 - COORDS(3)
C
      FLUX(1) = 0.0D0
      FLUX(2) = 0.0D0
      IF (D .GE. 0.0D0 .AND. D .LE. S) THEN
         HS = 2.0D0*U*IB/(PI*PHI**2)
     1        *EXP(-2.0D0*((COORDS(1)-X0)**2 + (COORDS(2)-Y0)**2)
     2        /PHI**2)
         IZ = (-2.25D0*(D/S)**2 + 1.5D0*(D/S) + 0.75D0)/0.75D0
         FLUX(1) = ETA*HS*IZ/S
      END IF
C
      RETURN
      END
