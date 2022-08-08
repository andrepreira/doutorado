
MODULE procedures

contains
    !#########################################################################!
        !(100,C1,tamanho_matriz,C3,v_zero,C4,,C7,C5)
    SUBROUTINE SUBARCH_HISTOGRAMA(NOUT,C1,i1,C3,r3,C4,i4,C7,C5)
        CHARACTER*300 FOUT,FDUMMY
        CHARACTER*42 C1,C2,C3,C4,C5,C6,C7
        integer i1,i2,i3,i4,i5,i6,i7
        REAL*8 r1,r2,r3,r4,r5,r6,r7,r8
    
        FDUMMY  = ' '
        FOUT    = ' ' !!C1,tamanho_matriz,C3,v_zero,C7,C5
        WRITE(FDUMMY,'(A,I6,A,F7.2,A,I6,A,A)')C1,i1,C3,r3,C4,i4,C7,C5
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


    SUBROUTINE sort(a)
        REAL(KIND=8), INTENT(in out), DIMENSION(:) :: a
        REAL(KIND=8) :: temp
        INTEGER :: i, j
        LOGICAL :: swapped
    
        DO j = SIZE(a)-1, 1, -1
        swapped = .FALSE.
        DO i = 1, j
            IF (a(i) > a(i+1)) THEN
            temp = a(i)
            a(i) = a(i+1)
            a(i+1) = temp
            swapped = .TRUE.
            END IF
        END DO
        IF (.NOT. swapped) EXIT
        END DO
    END SUBROUTINE sort   

END MODULE

PROGRAM histograma

    ! IMPLICIT none
    USE procedures
    INTEGER :: i, j, k, io, nlines, tamanho_matriz, qtd_pontos_pular
    CHARACTER (len=255)  :: nome_do_arquivo
    REAL(KIND=8) :: x_line,y_line,v_zero
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: histo_shape, histo_calda
    ! gfortran histograma.f90 -O2 -o histo.exe

    ! SUBARCH
!---------------------------CRIANDO ARQUIVOS DE SAIDA----------------------------------------------
    CHARACTER*60 C1,C2,C3,C4,C5,C6,C7,C8
    
    ! tamanho_matriz = 1000
    ! v_zero = 1.0
    qtd_pontos_pular = 50

    PRINT *, "Digite o nome do arquivo"
    READ (*,*) nome_do_arquivo
    OPEN (1,file=nome_do_arquivo)
    
    print *, nome_do_arquivo

    ! tamanho_matriz
    PRINT *, " tamanho_matriz = "
    READ *, tamanho_matriz

    ! v_zero
    PRINT *, " v_zero = "
    READ *, v_zero

    ! ! qtd_pontos_pular
    ! PRINT *, " qtd_pontos_pular = "
    ! READ *, qtd_pontos_pular
    
    C2="n_amostras="
    C3="v_zero="
    C4=""
    C5=".dat"
    C6=""
    C7="_SHAPE"
    
    C1="HISTOGRAMA_N="
    !C1,tamanho_matriz,C3,v_zero,C7,C5
    CALL SUBARCH_HISTOGRAMA(100,C1,tamanho_matriz,C3,v_zero,C4,qtd_pontos_pular,C7,C5)
    C7="_CALDA"
    CALL SUBARCH_HISTOGRAMA(200,C1,tamanho_matriz,C3,v_zero,C4,qtd_pontos_pular,C7,C5)
    
!------------------------CRIANDO ARQUIVOS DE SAIDA-------------------------------------------------
    
    nlines = 0 
    x_line=0
    y_line=0
    DO
        READ(1,*,iostat=io)
        IF (io/=0) EXIT
        nlines = nlines + 1
    END DO

    REWIND(1)


    ALLOCATE(histo_shape(nlines))
    ALLOCATE(histo_calda(nlines))

    DO j =1,nlines
        READ(1,*,iostat=io)x_line,y_line
        histo_shape(j) = LOG(y_line)
        histo_calda(j) = y_line
    END DO

    ! print *," "

    ! DO i = 1,nlines, 40
    !     print *, i, histo_shape(i), histo_calda(i)
    ! END DO

    !! algoritmo Bubble sort

    CALL sort(histo_shape)
    CALL sort(histo_calda)

    ! print *," "

    ! DO i = 1,nlines, 40
    !     print *, i, histo_shape(i), histo_calda(i)
    ! END DO

    !ler uma quantidade menor de pontos e escrever arquivo
    DO i = 1,nlines,qtd_pontos_pular
        WRITE(100,*) i,histo_shape(i)
        WRITE(200,*) i,histo_calda(i)
    END DO

    DEALLOCATE(histo_shape)
    DEALLOCATE(histo_calda)
    CLOSE (1)


END PROGRAM histograma

