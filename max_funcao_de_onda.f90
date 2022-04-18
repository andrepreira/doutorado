PROGRAM matriz1e
 
!  IMPLICIT none 
  INTEGER :: contE,k,j
  INTEGER :: seed, amostra, linha, coluna, tamanho_matriz
  DOUBLE PRECISION :: largura_desordem, largura_desordem_hopping
  INTEGER, PARAMETER :: Nmax = 3000, MMAX=10000000, MMAX1=10000000
  DOUBLE PRECISION :: WORK(MMAX)
  DOUBLE PRECISION:: PAR, max = 0.
  INTEGER:: IWORK(MMAX1),LWORK,LIWORK,INFO
  DOUBLE PRECISION :: Ep(Nmax),T(Nmax),H(Nmax, Nmax) = 0.0d0
  CHARACTER*60 C1,C2,C3,C4,C5,C6,C7,C8
  CHARACTER*1 COMPZ,UPLO, JOBZ
  ! COMPZ="V"
  JOBZ="V"
   UPLO="U"
! gfortran max_funcao_de_onda.f90 -O2 -o teste.exe -L/usr/local/lib -llapack -lblas

! Dimensão da matriz

! OPEN(unit=100, file='dados.txt', status='old', action='read')
! read(100,*)tamanho_matriz,largura_desordem,largura_desordem_hopping,seed

   READ(*,*)  tamanho_matriz, largura_desordem, largura_desordem_hopping, seed

!    PRINT *, " Digite o numero de Sitios:"
!    READ *, tamanho_matriz

! ! Digite a largura da desordem de hopping
!     PRINT *, "Digite a largura da desordem de potencial."
!     READ *, largura_desordem
! ! Digite a largura da desordem de hopping
!     PRINT *, "Digite a largura da desordem do hopping."
!     READ *, largura_desordem_hopping
        
! !Valor da seed
!    PRINT *, " Digite o valor da seed: "
!    READ *, seed

!#########################################################################!        
!############# Registre a energia potencial em um arquivo histórico ######!
              
          C2="L_desordem="
          C3=""
          C4=""
          C5=".dat"
          C6=""
	  C7="AMTRA="
	  
!	  C1="UoVDENERL" C2 = "WX"
! CALL SUBARCH(10,C1,tamanho_matriz,C2,largura_desordem,C7,seed,C5)
!	  C1="UoVDpartiL" C2 = "WX"
! CALL SUBARCH(20,C1,tamanho_matriz,C2,largura_desordem,C7,seed,C5)
!	  C1="UoVDENERL"
!   CALL SUBARCH(10,C1,tamanho_matriz,C2,largura_desordem,C7,seed,C5)
!	  C1="UoVDpartiL"
!   CALL SUBARCH(20,C1,tamanho_matriz,C2,largura_desordem,C7,seed,C5)
  !         C1="PART_VS_N="
  !  CALL SUBARCH(90,C1,tamanho_matriz,C2,largura_desordem,C7,seed,C5)
    !        C1="PART_VS_E_N="
    ! CALL SUBARCH(20,C1,tamanho_matriz,C2,largura_desordem,C7,seed,C5)
          C1="MAX_WAVEFUNC_VS_E_N="
   CALL SUBARCH(80,C1,tamanho_matriz,C2,largura_desordem,C7,seed,C5)

!#########################################################################!        
!############# Registre a energia potencial em um arquivo histórico ######!


!#########################################################################!
!#################### DISTRIBUIÇÃO ORDENADA DE POTENCIAL #################!	
!       DO linha=1, tamanho_matriz
!          Ep(linha)= 0.0d0       
!       END DO
!#########################################################################!
!#################### DISTRIBUIÇÃO ORDENADA DE POTENCIAL #################!

!Condições iniciais 
   amostra = -seed
!#########################################################################!
!################# DISTRIBUIÇÃO DESORDENADA DE POTENCIAL #################!

!Gerador de números aleatórios

    Do  linha =1, tamanho_matriz
        !  Ep(linha) = 0.0d0
         Ep(linha) =(ran1(amostra)-0.5d0)*largura_desordem
         T(linha) =(ran1(amostra)-0.5d0)*largura_desordem_hopping
!         PRINT *, "Desordem da Energia Potencial"
!         PRINT *, Ep(linha)
!         PRINT *,
    END DO
  
!#########################################################################!
!################# DISTRIBUIÇÃO DESORDENADA DE POTENCIAL #################!


!#########################################################################!
!############ DISTRIBUIÇÃO DESORDENADA DE HOPPPING UNIFORME ############# !

!#########################################################################!
!######################### HAMILTONIANO #################################!

     Do coluna = 1, tamanho_matriz
      Do linha = 1, tamanho_matriz
          ! H( linha , coluna ) =  0.0d0
          IF( linha == coluna) THEN
            H( linha, coluna ) = Ep(linha)
          ELSE
              IF ( linha .EQ. coluna + 1)THEN
                H( linha, coluna ) = T(linha)
                H( coluna, linha ) =H( linha, coluna )
              END IF
          END IF
       END DO  
     END DO
    !condicao de cadeia fechada
     H(1,tamanho_matriz)=T(tamanho_matriz)
     H(tamanho_matriz,1)= H(1,tamanho_matriz)
     
!Imprima a matriz uma linha de cada vez
 
    ! Do linha = 1, tamanho_matriz
    !    PRINT *, ( H( linha, coluna ), coluna = 1, tamanho_matriz)	
    ! END DO



!#########################################################################!
!######################### HAMILTONIANO #################################!


 10966 continue
!#########################################################################!
!####################### ROTINA DE DIAGONALIZAÇÃO ########################!
  
 	LWORK=10+6*tamanho_matriz+2*tamanho_matriz**2
        LIWORK=10+5*tamanho_matriz

! Chama as rotinas de diagonalização Ep= diagonal principal, 
! tamanho_matriz= tamanho da cadeia, Nmax= tamanho maximo da matriz,
! H = matriz a ser diagonalizada
        INFO = 0

        CALL DSYEVD(JOBZ,UPLO,tamanho_matriz,H,Nmax,Ep,WORK,LWORK, &
                                                      IWORK,LIWORK,INFO)
        ! PRINT * , INFO
  
        IF (INFO.ne.0) then 
          go to 10966 
        END IF
          
 
!#########################################################################!
!####################### ROTINA DE DIAGONALIZAÇÃO ########################!
    ! Do linha = 1, tamanho_matriz
    !    PRINT *, ( H( linha, coluna ), coluna = 1, tamanho_matriz)	
    ! END DO

!#########################################################################!
!################## Escreve os auto-estados de energia da Matriz #########!
       DO linha = 1, tamanho_matriz
         WRITE(10,*) Ep(linha)
         if (Ep(linha).ne. 0.0d0) contE=contE+1
       END DO
!#########################################################################!
!################### Calculo da função participação ######################!



       do j=1,contE
         PAR=0.0d0
         do k=1,tamanho_matriz
           PAR=PAR+H(k,j)**4.0d0
         enddo
	        write(20,*)Ep(j),1.0d0/PAR
          WRITE(90,*)1.0d0/PAR
       enddo


!#########################################################################!
!################### Calculo da função participação ######################!


!#########################################################################!
!################### Calculo da função de onda máxima ao quadrado ######################!

  ! print *, " "    

  ! Do j = 1, contE
  !   PRINT *, ( H( j, k ), k = 1, tamanho_matriz)
  ! END DO

  ! PRINT *, "O numero de autoestados não nulos é: "
  ! PRINT *, contE     

  max = 0.0d0
  do j=1,contE
    max = H(1,j)**2
    do k=1,tamanho_matriz
      if (max < H(k,j)**2) THEN
        max = H(k, j)**2
      END IF
    enddo
    WRITE(80,*)Ep(j),max
  enddo


!#########################################################################!
!################### Calculo da função de onda máxima ao quadrado ######################!

END PROGRAM matriz1e
          
!Subrotina


!#########################################################################!
	!C1,tamanho_matriz,C2,largura_desordem,C7,seed,C5)
      SUBROUTINE SUBARCH(NOUT,C1,i1,C2,r2,C7,i7,C5)
      CHARACTER*300 FOUT,FDUMMY
      CHARACTER*42 C1,C2,C3,C4,C5,C6,C7
      integer i1,i2,i3,i4,i5,i6,i7
      DOUBLE PRECISION r1,r2,r3,r4,r5,r6,r7,r8
     
      FDUMMY  = ' '
      FOUT    = ' ' !C1,N,C2,WX,C7,SEED3,C5)
      WRITE(FDUMMY,'(A,I6,A,F7.2,A,I6,A)')C1,i1,C2,r2,C7,i7,C5
!     REMOVE THE BLANKS FROM NAME
      IPOS = 0
          DO 41 I=1,LEN(FDUMMY)
              IF(FDUMMY(I:I).NE.' ')THEN
                  IPOS = IPOS+1
                  FOUT(IPOS:IPOS) = FDUMMY(I:I)
              ENDIF
  41       CONTINUE  
      OPEN(NOUT,FILE=FOUT,STATUS='REPLACE')
      RETURN
      END


!#########################################################################!
          
         
         FUNCTION ran1(idum)
          INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
          REAL ran1,AM,EPS,RNMX
          PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836)
          PARAMETER (NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

          INTEGER j,k,iv(NTAB),iy
          SAVE iv,iy
          DATA iv /NTAB*0/, iy /0/
          if (idum.le.0.or.iy.eq.0) then
          idum=max(-idum,1)
          do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11              continue
          iy=iv(1)
          endif
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          j=1+iy/NDIV
          iy=iv(j)
          iv(j)=idum
          ran1=min(AM*iy,RNMX)
          return
          END
