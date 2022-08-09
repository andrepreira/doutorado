PROGRAM test_do
 
    !  IMPLICIT none 
      INTEGER :: contE,k,j,g

      ! gfortran teste_do.f90 -O2 -o teste.exe

      DO j=2000, 16000, 2000
            PRINT * , j
      END DO


END PROGRAM test_do