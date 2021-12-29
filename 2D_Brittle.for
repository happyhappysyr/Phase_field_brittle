c     Author: Mojtaba Abdolkhani 
c     Email: mojtababdolkhani@gmail.com
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
c
      INCLUDE 'ABA_PARAM.INC'
c
      CHARACTER*20 CMNAME
c
c      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
c     1 STRAN(NTENS), DSTRAN(NTENS),PROPS(NPROPS), DROT(3, 3),TIME(2)
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),time(2),
     2 predef(1),dpred(1),props(nprops),coords(3),drot(3,3),dfgrd0(3,3),
     3 dfgrd1(3,3),jstep(4)
c
      PARAMETER(zero=0.D0, one=1.D0, two=2.D0, three=3.D0)
c      
      parameter (zeta = 1.d0)
c
      DIMENSION one_hydro(ntens,ntens), dstran_trace(ntens),
     1 one_dev(ntens,ntens), STRESS_MINUS(NTENS), 
     2 DDSDDE_MINUS(NTENS,NTENS)
c      
      REAL*8 EMOD,ENU,EBULK3,EBULK,EG,EG2,EG3,ELAM  
C ======================================================================
c        STATEV array explained:
c         STATEV(1) - H 
c
c     The energy decomposition is triggered via material PROPS(6) in input.
c     Two options are available:   ELASTOPLASTIC_ISO(41) or ELASTOPLASTIC_SD(40)  
c     
c             ELASTOPLASTIC_ISO - no deformation energy decomposition
c             ELASTOPLASTIC_SD  - spherical-deviatoric energy decomposition
C ======================================================================     
      phi=temp+dtemp    
C ======================================================================       
      !Material properties from input file
      EMOD=PROPS(1)
      ENU=PROPS(2)
	  xl=props(3) ! Phase field length scale
      Gc=props(4) ! Toughness
	  jeltype1 = PROPS(5) ! AT-2(30) OR TH(31)
	  type2 = PROPS(6) ! SD(40) OR ISO(41)
      xK = 0.00000001d0
c      
      !Degradation function
      degfnc = zero
      degfnc = (one - xK)*(one - phi)**2 + xK
	  dg=-2.d0*(1.d0-phi)
      ddg=2.d0
c      
      !Computing lame constants
      EBULK3=EMOD/(one-two*ENU)
      EBULK = EBULK3/three
      EG2=EMOD/(one+ENU)
      EG=EG2/two
      EG3=three*EG
      ELAM=(EBULK3-EG2)/three
c      
C ======================================================================  
      !Building material matrix
      DO 20 K1=1,NTENS
        DO 11 K2=1,NTENS
           DDSDDE(K2,K1)= zero
 11     CONTINUE
 20   CONTINUE
C
      DO 40 K1=1,NDI
        DO 30 K2=1,NDI
           DDSDDE(K2,K1)=ELAM
 30     CONTINUE
        DDSDDE(K1,K1)=EG2+ELAM
 40   CONTINUE
      DO 50 K1=NDI+1,NTENS
        DDSDDE(K1,K1)=EG
 50   CONTINUE
C ======================================================================
c      
      STRAN = STRAN + DSTRAN
      STRESS = zero
c
      !Stress predictor
      DO K1=1, NTENS
         DO K2=1, NTENS
             STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*(STRAN(K1))
         END DO
      END DO
c      
	  if (type2.eq.40)then      
          ! DECOMPOSITION - VOLUMETRIC DEVIATORIC
          ONE_hydro = zero
          A_trace = zero
          dstran_trace = zero
c          
          do i = 1, ndi
              do j = 1, ndi
                  ONE_hydro(i,j) = one/three
              end do
          end do
c          
          ONE_Dev = zero
          do i = 1, ntens
              do j = 1, ntens
                  if (i .eq. j) ONE_dev (i,j) = ONE
              end do
          end do
c      
          ONE_dev = ONE_dev - ONE_hydro      
          dstran_trace = matmul(ONE_hydro, STRAN)
c         
          if (dstran_trace(1) .le. zero) then
              A_trace = one 
c      
              DDSDDE_MINUS = A_trace*matmul(ONE_hydro, DDSDDE)
              DDSDDE = DDSDDE - DDSDDE_MINUS
c      
              STRESS_MINUS = A_trace*matmul(ONE_hydro, STRESS)
              STRESS = STRESS - STRESS_MINUS
          else
              STRESS_MINUS = zero
              DDSDDE_MINUS = zero
          end if
	  elseif (type2.eq.41)then  
          STRESS_MINUS = zero
          DDSDDE_MINUS = zero
          A_trace = zero
	  else
          write(6,*)""
          write(6,*)"WARNING: CHECK TYPE2 PROP INPUT!"
	  end if
C ======================================================================         
      DDSDDE = degfnc*DDSDDE
      STRESS = degfnc*STRESS
c         
      DDSDDE = DDSDDE + DDSDDE_MINUS
      STRESS = STRESS + STRESS_MINUS
C ======================================================================         
      SSE = zero      !elastic deformation energy array
c      
      SSE = 0.5d0*ELAM*(one-A_trace)*(stran(1) + stran(2) + stran(3))**2
     1 + EG*(stran(1)**2 + stran(2)**2 + stran(3)**2) +
     2 0.5d0*EG*(stran(4)**2)
c	
c   ************************************************************************************************************   
c      
	  !SSE - Elastic energy
	  psiE_old = statev(1)                  
	  psiE_new = SSE 
c      
	  !SPD - Plastic energy
	  psiE_SPD = zero
c                  
	  psiE_SSE = zero
	  psiE = zero
c                 
c	  ! AT-2 or TH model
	  if (jeltype1.eq.30)then      !AT-2 model
		  psiE_crit = zero
	  elseif (jeltype1.eq.31)then  !TH model
		  psiE_crit = 0.2652d0*(Gc/xL)
	  else
		  write(6,*)"***ERROR - check element types"
	  end if
c                 
	  ! Elastic energy history field
	  if (psiE_new.gt.psiE_old) then
		  psiE_SSE = psiE_new
	  else
		  psiE_SSE = psiE_old
	  end if
c                  
	  if ((psiE_SSE + psiE_SPD).gt.psiE_crit) then
		  psiE = psiE_SSE + psiE_SPD - psiE_crit
		  psiE = zeta*psiE
	  else
		  psiE = zero
	  end if
c                  
	  statev(1) = psiE_SSE
c ***************************************************************************************				  
      w=phi**2
      dw=2.d0*phi
      ddw=2.d0
      cw=0.5d0
c ***************************************************************************************
	  if (jeltype1.eq.30)then      !AT-2 model
		  rpl=((psiE/(Gc/(two*xl)))*((1.d0-phi)/(xl**2)))-(phi/(xl**2))
          drpldt=-((1.d0+(psiE/(Gc/(two*xl))))/(xl**2))
	  elseif (jeltype1.eq.31)then  !TH model
		  rpl=((psiE/psiE_crit)*((1.d0-phi)/(xl**2)))-(phi/(xl**2))
          drpldt=-((1.d0+(psiE/psiE_crit))/(xl**2))
	  else
		  write(6,*)"***ERROR - check element types"
	  end if
c 
      RETURN
      END