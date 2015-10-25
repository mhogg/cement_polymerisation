
      module sharedata

C     Define parameters
      DOUBLE PRECISION, PARAMETER :: w_init=1.0D-6
      INTEGER, PARAMETER :: elm_min=589, elm_max=798

C     Define and set initial values of shared variables
      DOUBLE PRECISION :: w(elm_min:elm_max,4), wpi(elm_min:elm_max,4)
      DOUBLE PRECISION :: wmax_peak(elm_min:elm_max,4)
      INTEGER :: nnoel, nnpt, kkinc, whichstep(elm_min:elm_max,4)

      CONTAINS

        subroutine initialise()

          w=w_init
          wmax_peak=w_init
          whichstep=0

        end subroutine initialise

      end module sharedata

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      subroutine USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,TIME,DTIME,
     1      CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,
     2      NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)

      use sharedata

      include 'ABA_PARAM.INC'

      character*80 CMNAME,ORNAME
      character*3 FLGRAY(15)
      dimension FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),T(3,3),TIME(2)
      dimension ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C ----------------------------------------------------------------------

C     User variables
      INTEGER :: kounter=0

C ----------------------------------------------------------------------

C     NOTES:
C     ------
C     (1) One of the material properties must be a function of a field
C         variable otherwise this subroutine is not called. This can be
C         either density, conductivity or specific heat

C ----------------------------------------------------------------------

C     Call subroutine initialise in USDFLD as subroutine USDFLD is called
C     before subroutine HETVAL
      if (kounter/=1) then
        CALL initialise()
        kounter=1
      end if

      if (CMNAME(1:6)=='CEMENT') then

C       Subroutine HETVAL does not have NOEL, NPT or KINC passed in so
C       make them available via module sharedata
        nnoel = NOEL
        nnpt  = NPT
        kkinc = KINC

      end if

C ----------------------------------------------------------------------

      return
      end

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE HETVAL(CMNAME,TEMP,TIME,DTIME,STATEV,FLUX,PREDEF,DPRED)

      use sharedata

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION TEMP(2),STATEV(*),PREDEF(*),TIME(2),FLUX(2),DPRED(*)
     
C ----------------------------------------------------------------------

C     User variables
      INTEGER :: ref=1
      REAL, PARAMETER :: delT=1.0
      DOUBLE PRECISION :: w_dum,wpi_dum,w1,w2,wpeak
      REAL :: Qout, dwdt,temperature,dt
      REAL :: Qout1,dwdt1,T1,Qout2,dwdt2,T2

C ----------------------------------------------------------------------

      if (kkinc>=1) then

C       At the start of each increment set wpi to the value of w at the 
C       previous increment
        if (whichstep(nnoel,nnpt)/=kkinc) then
          wpi(nnoel,nnpt)=w(nnoel,nnpt)
          whichstep(nnoel,nnpt)=kkinc
        end if

C       Calculate rate of change of polyermization
        wpi_dum=wpi(nnoel,nnpt)
        wpeak=wmax_peak(nnoel,nnpt)
        dt=DTIME
        temperature=TEMP(1)
        CALL poly(Qout,dwdt,w_dum,wpi_dum,dt,temperature,ref,wpeak)

C       Capture output from subroutine polymerization
        w(nnoel,nnpt) = w_dum
        wmax_peak(nnoel,nnpt)=wpeak

C       Set FLUX(1) to the heat generated by the polymerization process
        FLUX(1)=Qout

C       Calculate FLUX(2) - the rate of change of the flux wrt temperature
C       at this point - to force solver to iterate and help convergence
        T1=TEMP(1)-delT
        T2=TEMP(1)+delT
        CALL poly(Qout1,dwdt1,w1,wpi,dt,T1,ref,wpeak)
        CALL poly(Qout2,dwdt2,w2,wpi,dt,T2,ref,wpeak)
        FLUX(2)=(Qout2-Qout1)/(T2-T1)

        STATEV(1) = w(nnoel,nnpt)
        STATEV(2) = dwdt
        STATEV(3) = FLUX(1)
        STATEV(4) = FLUX(2)

      end if

C ----------------------------------------------------------------------

      return
      end

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      subroutine poly(Qout,dwdt,w,wpi,dt,temperature,ref,wpeak)

      implicit none

C     This subroutine is used to calculate the rate of change of the fraction
C     of cement polymerization that has occurred. From this the heat 
C     generated by this exothermic reaction can be calculated. 

C     All equations and variable values are taken from:
C     Quarini et al, Numerical predictions of the thermal behaviour and resultant
C     effects of grouting cements while setting prosthetic components in bone,
C     J Engineering in Medicine, 2006, 220:625-634

      INTEGER:: ref, p, iteration
      DOUBLE PRECISION, PARAMETER :: TOL=1.0D-12
      DOUBLE PRECISION :: w, wpi, wpeak, w_in, resid
      REAL :: Qout, dwdt, wmax, temperature, dt
      REAL:: R, m, n, K, Ko, Eo, H, density_cement, Tg_deg, Tg
      LOGICAL :: flag

C     Set common variables
      density_cement = 1.1e-9       ! tonne/mm3 
      R = 8.314e+3                  ! mJ/(mol.K)

C     Set polymerization variables     
      select case (ref)
      case(1)
C       Brozacchiello et al, 1998 [7] - Listing 1
        Ko = 1.2088e+04             ! 1/s
        Eo = 3.3240e+07             ! mJ/mol
        m = 0.98                    ! non-dimensional
        n = 1.14                    ! non-dimensional
        p = 0                       ! non-dimensional
        H = 1.0e+11                 ! mJ/tonne
      case(2)
C       Brozacchiello et al, 1998 [7] - Listing 2
        Ko = 9.4840e+03             ! 1/s
        Eo = 3.1270e+07             ! mJ/mol
        m = 0.98                    ! non-dimensional
        n = 1.064                   ! non-dimensional
        p = 0                       ! non-dimensional
        H = 1.23e+11                ! mJ/tonne
      case(3)
C       Li et al, 2003 [13]
        Ko = 5.206e+09              ! 1/s
        Eo = 6.7891e+07             ! mJ/mol
        m = 1.05                    ! non-dimensional
        n = 1.10                    ! non-dimensional
        p = 0                       ! non-dimensional
        H = 1.0e+11                 ! mJ/tonne    
      end select

C     Set polymerization variable wmax
      select case (ref)
      case(1,2)
        Tg_deg=105.0                ! deg C
        Tg=273.0+Tg_deg             ! deg K
        if (temperature <= Tg) then
          wmax=(temperature/Tg)
        else
          wmax=1.0D0
        end if
        wpeak=MAX(wmax,wpeak)
        wmax=wpeak
      case(3)
        wmax=1.0
      end select

C     Calculate function K
      K=Ko*exp(-Eo/(R*temperature))

C     Update polymerization fraction using a fully implicit formulation with
C     a backwards difference wrt time for the rate of change of w, dwdt
      if (wpi>=wmax) then
        w=wmax
      else
        iteration=0
        flag=.FALSE.
        w_in=wpi
        do while(flag==.FALSE.)
          w = wpi + dt*K*(w_in**m)*((wmax-w_in)**n)*(1.0D0-w_in)**p
          w = MIN(w,wmax)
          resid=abs(w-w_in)
          if (resid<TOL .OR. iteration>10) then
            flag=.TRUE.
          else
            w_in=w
            iteration = iteration+1
            if (iteration>5) then
              write(*,*) iteration, wpi, w, wmax, resid
            end if
          end if
        end do
      end if

C     Calculate rate of change of polymerization fraction
      dwdt = (w-wpi)/dt

C     Calculate heat generated by polymerization process at this point
      Qout = density_cement*H*dwdt 

      return
      end
