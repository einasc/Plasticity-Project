       SUBROUTINE VUMAT(
! READ ONLY - DO NOT MODIFY
     . NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL, STEPTIME,
     . TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH, PROPS, DENSITY,
     . STRAININC, RELSPININC, TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     . STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD, TEMPNEW,
     . STRETCHNEW, DEFGRADNEW, FIELDNEW,
! WRITE ONLY - DO NOT READ
     . STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
!-----------------------------------------------------------------------
!     ABAQUS implicit variable declaration included in VABA_PARAM.INC
!     states the following:
!     a to h are real variables
!     o to z are real variables
!     i to n are integer variables
!-----------------------------------------------------------------------
      INCLUDE 'VABA_PARAM.INC'
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     ABAQUS variables 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK,*),
     . CHARLENGTH(*), STRAININC(NBLOCK,NDIR+NSHR), RELSPININC(*),
     . TEMPOLD(*), STRETCHOLD(*), DEFGRADOLD(*), FIELDOLD(*),
     . STRESSOLD(NBLOCK,NDIR+NSHR), STATEOLD(NBLOCK,NSTATEV),
     . ENERINTERNOLD(NBLOCK),  ENERINELASOLD(NBLOCK), TEMPNEW(*),
     . STRETCHNEW(*), DEFGRADNEW(*), FIELDNEW(*),
     . STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV),
     . ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)
C
      CHARACTER*80 CMNAME
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Internal UMAT variables 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DIMENSION STRESS(NDIR+NSHR) ! Stress tensor inside UMAT
      DIMENSION DFDS(NDIR+NSHR)   ! Derivative of the yield function
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Internal UMAT variables could be declared as follow
!     but it is not required due to the VABA_PARAM.INC file
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!   Material parameters
!      REAL YOUNG       ! Young's modulus 
!      REAL POISS       ! Poisson's ratio
!      REAL C11,C12,C44 ! Elasticity matrix C components
!      REAL SIGMA0      ! Initial yield stress
!      REAL ET          ! Tangent modulus
!	   REAL ALPHA		! Pressure dependency parameter
!   Internal variables
!      REAL P           ! Equivalent plastic strain
!   Plasticity variables
!      REAL DLAMBDA     ! Plastic multiplier
!      REAL DDLAMBDA    ! Increment in plastic multiplier
!   Yield function variables
!      REAL F           ! Yield function
!      REAL PHI         ! Equivalent stress
!      REAL SIGMAY      ! Yield stress
!   Computational variables
!      REAL RESNOR      ! Convergence criterion
!      REAL TOL         ! Tolerance for the RMAP algorithm
!      INTEGER MXITER   ! Maximum number of iteration for the RMAP
!      INTEGER ITER     ! Number of iteration for the RMAP
!-----------------------------------------------------------------------
!     Read parameters from ABAQUS material card
!-----------------------------------------------------------------------
      YOUNG  = props(1)
      POISS  = props(2)
      SIGMA0 = props(3)
      ET     = props(4)
      TOL    = props(5)
      MXITER = props(6)
	  ALPHA  = props(7)
!-----------------------------------------------------------------------
!     Compute elasticity matrix
!-----------------------------------------------------------------------
      C11    = YOUNG*(1.0-POISS)/((1.0+POISS)*(1.0-2.0*POISS))
      C12    = POISS*C11/(1.0-POISS)
      C44    = 0.5*YOUNG/(1.0+POISS)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Loop over integration points
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DO i=1,NBLOCK
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!        If time = 0 then pure elastic computation
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
         IF(TOTALTIME.eq.0.0)THEN
            STRESSNEW(i,1) = STRESSOLD(i,1)+C11*STRAININC(i,1)
     .                       +C12*STRAININC(i,2)+C12*STRAININC(i,3)
            STRESSNEW(i,2) = STRESSOLD(i,2)+C12*STRAININC(i,1)
     .                       +C11*STRAININC(i,2)+C12*STRAININC(i,3)
            STRESSNEW(i,3) = STRESSOLD(i,3)+C12*STRAININC(i,1)
     .                       +C12*STRAININC(i,2)+C11*STRAININC(i,3)
            STRESSNEW(i,4) = STRESSOLD(i,4)+C44*STRAININC(i,4)*2.0
            STRESSNEW(i,5) = STRESSOLD(i,5)+C44*STRAININC(i,5)*2.0
            STRESSNEW(i,6) = STRESSOLD(i,6)+C44*STRAININC(i,6)*2.0
         ELSE
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!        Elastic-predictor-corrector scheme
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------------------------------------
!           Elastic prediction
!-----------------------------------------------------------------------
            STRESS(1) = STRESSOLD(i,1)+C11*STRAININC(i,1)
     .                                +C12*STRAININC(i,2)
     .                                +C12*STRAININC(i,3)
            STRESS(2) = STRESSOLD(i,2)+C12*STRAININC(i,1)
     .                                +C11*STRAININC(i,2)
     .                                +C12*STRAININC(i,3)
            STRESS(3) = STRESSOLD(i,3)+C12*STRAININC(i,1)
     .                                +C12*STRAININC(i,2)
     .                                +C11*STRAININC(i,3)
            STRESS(4) = STRESSOLD(i,4)+C44*STRAININC(i,4)*2.0
            STRESS(5) = STRESSOLD(i,5)+C44*STRAININC(i,5)*2.0
            STRESS(6) = STRESSOLD(i,6)+C44*STRAININC(i,6)*2.0   
!-----------------------------------------------------------------------
!           Von Mises stress
!-----------------------------------------------------------------------
            VMSTRESS   = sqrt(STRESS(1)*STRESS(1)
     .                      +STRESS(2)*STRESS(2)
     .                      +STRESS(3)*STRESS(3)
     .                      -STRESS(1)*STRESS(2)
     .                      -STRESS(2)*STRESS(3)
     .                      -STRESS(3)*STRESS(1)
     .                 +3.0*(STRESS(4)*STRESS(4)
     .                      +STRESS(5)*STRESS(5)
     .                      +STRESS(6)*STRESS(6)))
!-----------------------------------------------------------------------
!           Hydrostatic stress
!-----------------------------------------------------------------------
            HYDSTRESS   = (STRESS(1)+STRESS(2)+STRESS(3))/3.0
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!           Equivalent stress
!-----------------------------------------------------------------------
            PHI   = sqrt(1.0/(1.0+(ALPHA/3.0)**2)
     .					   *(VMSTRESS**2+ALPHA**2*HYDSTRESS**2)) 
!-----------------------------------------------------------------------
!           Equivalent plastic strain from previous time step
!-----------------------------------------------------------------------
            P        = STATEOLD(i,1)
!-----------------------------------------------------------------------
!           Update yield stress
!-----------------------------------------------------------------------
            SIGMAY   = SIGMA0+ET*P
!-----------------------------------------------------------------------
!           Compute yield function
!-----------------------------------------------------------------------
            F        = PHI-SIGMAY
!-----------------------------------------------------------------------
!           Initialize the plastic multiplier
!-----------------------------------------------------------------------
            DLAMBDA  = 0.0
!-----------------------------------------------------------------------
!           Find the derivative of the yield function
!-----------------------------------------------------------------------
            IF(PHI.eq.0)THEN
			   DENOM   = 1.0
            ELSE
			   DENOM   = PHI
            ENDIF        
c			
!-----------------------------------------------------------------------
!           Compute the derivative of yield function
!-----------------------------------------------------------------------
            DFDS(1) = 3.0/(2.0*DENOM*(1.0+(ALPHA/3)**2))
     .				*(STRESS(1)+(2*(ALPHA/3)**2-1)*HYDSTRESS)
            DFDS(2) = 3.0/(2.0*DENOM*(1.0+(ALPHA/3)**2))
     .				*(STRESS(2)+(2*(ALPHA/3)**2-1)*HYDSTRESS) 
            DFDS(3) = 3.0/(2.0*DENOM*(1.0+(ALPHA/3)**2))
     .				*(STRESS(3)+(2*(ALPHA/3)**2-1)*HYDSTRESS) 
            DFDS(4) = 3.0/(2.0*DENOM*(1.0+(ALPHA/3)**2))*(STRESS(4))
            DFDS(5) = 3.0/(2.0*DENOM*(1.0+(ALPHA/3)**2))*(STRESS(5))
            DFDS(6) = 3.0/(2.0*DENOM*(1.0+(ALPHA/3)**2))*(STRESS(6))	
!			write(*,*) 'DFDS1',DFDS(1)
!			write(*,*) 'DFDS2',DFDS(2)
!			write(*,*) 'DFDS4',DFDS(4)
!			write(*,*) 'DFDS5',DFDS(5)
!			write(*,*) 'DFDS6',DFDS(6)
!-----------------------------------------------------------------------
!           Helping parameter H
!-----------------------------------------------------------------------
            H = DFDS(1)*C11*DFDS(1)+DFDS(2)*C11*DFDS(2)
     .			+DFDS(3)*C11*DFDS(3)
     .			+2.0*DFDS(1)*C12*DFDS(2)+2.0*DFDS(1)*C12*DFDS(3)
     .			+2.0*DFDS(2)*C12*DFDS(3)
     .          +DFDS(4)*C44*DFDS(4)+DFDS(5)*C44*DFDS(5)
     .			+DFDS(6)*C44*DFDS(6)
!	            write(*,*) 'H',H
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!           Check for plasticity
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------------------------------------
!           Plastic flow
!-----------------------------------------------------------------------
            IF(F.GT.0.0)THEN ! Plastic flow
               DO ITER=1,MXITER
!-----------------------------------------------------------------------
!                 Compute increment in plastic multiplier
!-----------------------------------------------------------------------
                  DDLAMBDA = F/(H-ET)
!				  write(*,*) 'F', F
!				  write(*,*) 'H', H
!				  write(*,*) 'ET', ET
!-----------------------------------------------------------------------
!                 Update plastic multiplier
!-----------------------------------------------------------------------
                  DLAMBDA  = DLAMBDA+DDLAMBDA
!-----------------------------------------------------------------------
!                 Update equivalent plastic strain
!-----------------------------------------------------------------------
                  P      = P+DDLAMBDA
!-----------------------------------------------------------------------
!                    Increment stress tensor
!-----------------------------------------------------------------------
                     STRESS(1)= STRESS(1)-DDLAMBDA*(C11*DFDS(1)
     .                                                 +C12*DFDS(2)
     .                                                 +C12*DFDS(3))
                     STRESS(2)= STRESS(2)-DDLAMBDA*(C12*DFDS(1)
     .                                                 +C11*DFDS(2)
     .                                                 +C12*DFDS(3))
                     STRESS(3)= STRESS(3)-DDLAMBDA*(C12*DFDS(1)
     .                                                 +C12*DFDS(2)
     .                                                 +C11*DFDS(3))
                     STRESS(4)= STRESS(4)-DDLAMBDA*C44*DFDS(4)
                     STRESS(5)= STRESS(5)-DDLAMBDA*C44*DFDS(5)
                     STRESS(6)= STRESS(6)-DDLAMBDA*C44*DFDS(6) 			
!-----------------------------------------------------------------------
!           Calculate intermediate Von Mises stress
!-----------------------------------------------------------------------
            VMSTRESS   = sqrt(STRESS(1)*STRESS(1)
     .                      +STRESS(2)*STRESS(2)
     .                      +STRESS(3)*STRESS(3)
     .                      -STRESS(1)*STRESS(2)
     .                      -STRESS(2)*STRESS(3)
     .                      -STRESS(3)*STRESS(1)
     .                 +3.0*(STRESS(4)*STRESS(4)
     .                      +STRESS(5)*STRESS(5)
     .                      +STRESS(6)*STRESS(6)))
!-----------------------------------------------------------------------
!           Calculate intermediate Hydrostatic stress
!-----------------------------------------------------------------------
            HYDSTRESS   = (STRESS(1)+STRESS(2)+STRESS(3))/3.0
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!           Equivalent stress
!-----------------------------------------------------------------------
            PHI   = sqrt(1.0/(1.0+(ALPHA/3.0)**2)
     .					   *(VMSTRESS**2+ALPHA**2*HYDSTRESS**2)) 
!-----------------------------------------------------------------------					 
!-----------------------------------------------------------------------
!                 Update Intermediate Yield stress
!-----------------------------------------------------------------------
                  SIGMAY = SIGMA0+ET*P
!-----------------------------------------------------------------------
!                 Update Intermediate Yield function
!-----------------------------------------------------------------------
                  F      = PHI-SIGMAY
!-----------------------------------------------------------------------
!                 Compute convergence criterion
!-----------------------------------------------------------------------
                  RESNOR = ABS(F/SIGMAY)
!-----------------------------------------------------------------------
!                 Check for convergence
!-----------------------------------------------------------------------
                  IF(RESNOR.LE.TOL)THEN ! RMAP has converged
!-----------------------------------------------------------------------
!                    Update the stress tensor
!-----------------------------------------------------------------------
                     STRESSNEW(i,1)= STRESS(1)
                     STRESSNEW(i,2)= STRESS(2)
                     STRESSNEW(i,3)= STRESS(3)
                     STRESSNEW(i,4)= STRESS(4)
                     STRESSNEW(i,5)= STRESS(5)
                     STRESSNEW(i,6)= STRESS(6) 
!-----------------------------------------------------------------------
!                    Update the history variables
!-----------------------------------------------------------------------
                     STATENEW(i,1) = P
                     STATENEW(i,2) = PHI
                     STATENEW(i,3) = F
                     STATENEW(i,4) = SIGMAY
                     STATENEW(i,5) = DGAMA
                     STATENEW(i,6) = ITER 
                     GOTO 90
                  ELSE ! RMAP has not converged yet
                     IF(ITER.eq.MXITER)THEN
                        write(*,*) 'RMAP has not converged'
                        write(*,*) 'Integration point',i
                        write(*,*) 'Convergence',RESNOR
                        write(*,*) 'dlambda',DLAMBDA,P
                        STOP
                     ENDIF
                  ENDIF          
               ENDDO
!-----------------------------------------------------------------------
!           Elastic point
!-----------------------------------------------------------------------
            ELSE
!-----------------------------------------------------------------------
!              Update the stress tensor
!-----------------------------------------------------------------------
               STRESSNEW(i,1) = STRESS(1)
               STRESSNEW(i,2) = STRESS(2)
               STRESSNEW(i,3) = STRESS(3)
               STRESSNEW(i,4) = STRESS(4)
               STRESSNEW(i,5) = STRESS(5)
               STRESSNEW(i,6) = STRESS(6)
!-----------------------------------------------------------------------
!              Update the history variables
!-----------------------------------------------------------------------
               STATENEW(i,1) = P
               STATENEW(i,2) = PHI
               STATENEW(i,3) = F
               STATENEW(i,4) = SIGMAY
               STATENEW(i,5) = 0.0
               STATENEW(i,6) = 0
            ENDIF
!-----------------------------------------------------------------------
!        End of loop over integration points
!-----------------------------------------------------------------------
         ENDIF           
  90     CONTINUE  
      ENDDO
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      RETURN
      END