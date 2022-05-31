PROGRAM matriz1e
 
!  IMPLICIT none 
  INTEGER :: contE,k,j,g
  INTEGER :: seed, amostra, linha, coluna, tamanho_matriz, seed1, fator_amostra, n_amostras
  REAL :: golden_ratio
  DOUBLE PRECISION, PARAMETER :: PI=4.D0*DATAN(1.D0)
  INTEGER, PARAMETER :: Nmax = 6000, MMAX=100000000, MMAX1=100000000
  !COLOCAR UMA VALIDAÇÃO PARA AVISAR QUANDO MMAX E MMAX1 FOR ATINGIDO !!
  DOUBLE PRECISION :: WORK(MMAX)
  DOUBLE PRECISION:: PAR,PAR_MED, max = 0.0d0, segundo_max=0.0d0, max_wave_med = 0.0d0
  INTEGER:: IWORK(MMAX1),LWORK,LIWORK,INFO
  DOUBLE PRECISION :: Ep(Nmax),T(Nmax),H(Nmax, Nmax) = 0.0d0
  CHARACTER*60 C1,C2,C3,C4,C5,C6,C7,C8
  CHARACTER*1 COMPZ,UPLO, JOBZ
  ! COMPZ="V"
  JOBZ="V"
  UPLO="U"
  golden_ratio=(sqrt(5.)-1)/2
! gfortran participacao_max_funcao_media.f90 -O2 -o teste.exe -L/usr/local/lib -llapack -lblas
tamanho_matriz = 4000
n_amostras = 1

! Dimensão da matriz
PRINT *, " Digite o valor do fator multiplicador da amostra"
READ *, fator_amostra

C2="n_amostras="
C3="v_zero="
C4=""
C5=".dat"
C6=""
C7="E_PAR_MAX_WAVE_FUC_MEDIO"

C1="ARQUIVO_MEDIAS_N="
!C1,tamanho_matriz,C2,n_amostras,C7,C5
CALL SUBARCH_ARQUIVO_MEDIAS(100,C1,tamanho_matriz,C2,n_amostras,C7,C5)

Do  seed1 =1, n_amostras
    seed = seed1*fator_amostra

    !#########################################################################!        
    !############# Registre a energia potencial em um arquivo histórico ######!
                  
      C2="v_zero="
      C3="L_desordem_hopping="
      C4=""
      C5=".dat"
      C6=""
      C7="AMTRA="
    
      ! C1="WAVEFUNC_VS_E_N="
      ! CALL SUBARCH(80,C1,tamanho_matriz,C2,v_zero,C7,seed,C5)

    !#########################################################################!        
    !############# Registre a energia potencial em um arquivo histórico ######!


    !#########################################################################!
    !#################### DISTRIBUIÇÃO ORDENADA DE POTENCIAL #################!	
    !       DO linha=1, tamanho_matriz
    !          Ep(linha)= 0.0d0       
    !       END DO
    !#########################################################################!
    !#################### DISTRIBUIÇÃO ORDENADA DE POTENCIAL #################!
      ! Condições iniciais 
      amostra = -seed

    !#########################################################################!
    !################# DISTRIBUIÇÃO DESORDENADA DE POTENCIAL #################!

    !Gerador de números aleatórios
      DO g =10,40
              v_zero = g/10.0d0
                  Do  i =1, tamanho_matriz
                      !  Ep(i) = 0.0d0
                      Ep(i) =v_zero*cos(2*PI*golden_ratio*i)
                      T(i) =1.
              !         PRINT *, "Desordem da Energia Potencial"
              !         PRINT *, Ep(i)
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
                  ! PRINT *, ' '



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
                      ! WRITE(10,*) Ep(linha)
                      if (Ep(linha).ne. 0.0d0) contE=contE+1
                    END DO
                    ! print *,  ' '
                    ! print *,  contE
                    ! print *,  ' '
              !#########################################################################!
              !################### Calculo da função participação ######################!



                    ! do j=1,contE
                    !   PAR=0.0d0
                    !   do k=1,tamanho_matriz
                    !     PAR=PAR+H(k,j)**4.0d0
                    !   enddo
                    !     write(20,*)Ep(j),1.0d0/PAR
                    !     WRITE(90,*)1.0d0/PAR
                    ! enddo


              !#########################################################################!
              !################### Calculo da função participação ######################!
              
              !#########################################################################!
              !################### Calculo da função de onda ao quadrado ######################!
            

                    ! psi = 0.0d0
                    ! do j=1,contE
                    !       do k=1,tamanho_matriz
                    !               psi=H(k,j)**2.0d0
                    !               WRITE(80,*)Ep(j),psi
                    !       enddo
                    ! enddo
              !#########################################################################!
              !################### Calculo da função de onda ao quadrado ######################!

              !#########################################################################!
              !################### Calculo da função de onda máxima ao quadrado ######################!

                ! print *, " "    

                ! Do j = 1, contE
                !   PRINT *, ( H( j, k ), k = 1, tamanho_matriz)
                ! END DO

                ! PRINT *, "O numero de autoestados não nulos é: "
                ! PRINT *, contE     

                max = 0.0d0
                segundo_max = 0.0d0
                PAR_MED = 0.0d0
                max_wave_med = 0.0d0
                do j=1,contE
                      max = H(1,j)**2
                      PAR=0.0d0
                      do k=1,tamanho_matriz
                              PAR=PAR+H(k,j)**4.0d0
                              if (max < H(k,j)**2) THEN
                                segundo_max = max
                                max = H(k, j)**2
                              ELSE 
                                  if (segundo_max < H(k,j)**2) segundo_max = H(k, j)**2
                              END IF

                      enddo
                      PAR_MED = PAR_MED + 1.0d0/PAR
                      max_wave_med = max_wave_med + max
                    enddo
                    WRITE(100,*)PAR_MED/contE,max_wave_med/contE, v_zero
                contE = 0
                PAR_MED = 0.0d0
                max_wave_med = 0.0d0
                ! WRITE(100,*) ' '
                Do coluna = 1, tamanho_matriz
                  Do linha = 1, tamanho_matriz
                      H( linha , coluna ) =  0.0d0
                  END DO  
                END DO

                Do linha = 1, tamanho_matriz
                  Ep(linha) = 0.0d0
                  T(linha) = 0.0d0
                END DO 
          ENDDO
          !#########################################################################!
          !################### Calculo da função de onda máxima ao quadrado ######################!
    END DO
END PROGRAM matriz1e
          
!Subrotina


!#########################################################################!
	!(80,C1,tamanho_matriz,C2,v_zero,C7,seed,C5)
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

!#########################################################################!
	!(100,C1,tamanho_matriz,C2,n_amostras,C7,C5
      SUBROUTINE SUBARCH_ARQUIVO_MEDIAS(NOUT,C1,i1,C2,i2,C7,C5)
      CHARACTER*300 FOUT,FDUMMY
      CHARACTER*42 C1,C2,C3,C4,C5,C6,C7
      integer i1,i2,i3,i4,i5,i6,i7
      REAL r1,r2,r3,r4,r5,r6,r7,r8
     
      FDUMMY  = ' '
      FOUT    = ' ' !!C1,tamanho_matriz,C2,n_amostras,C7,C5
      WRITE(FDUMMY,'(A,I6,A,I6,A,A)')C1,i1,C2,i2,C7,C5
!     REMOVE THE BLANKS FROM NAME
      IPOS = 0
          DO 41 I=1,LEN(FDUMMY)
              IF(FDUMMY(I:I).NE.' ')THEN
                  IPOS = IPOS+1
                  FOUT(IPOS:IPOS) = FDUMMY(I:I)
              ENDIF
  41       CONTINUE  
      OPEN(NOUT,FILE=FOUT,STATUS='REPLACE', ACTION='WRITE', ACCESS ='SEQUENTIAL' )
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
