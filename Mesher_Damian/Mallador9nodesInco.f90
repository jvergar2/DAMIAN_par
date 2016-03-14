!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  PROGRAMA PARA GENERAR ARCHIVO DE ENTRADA PARA EL PROGRAMA DE ELEMENTOS    C
!C                   FINITOS DESARROLLADO POR JUAN CARLOS.                    C
!C                                                                            C
!C EL PROGRAMA SE ALIMENTA POR EL ARCHIVO *.msh GENERADO POR EL PROGRAMA GMSH C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C
PROGRAM MALLAS
!
!
USE OMP_LIB
IMPLICIT NONE
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    DECLARACIÓN DE VARIABLES
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!
!VARIABLES GENERALES.
INTEGER I, J, K, L, M, N, O, P, Q, F
INTEGER II, JJ
INTEGER NP, ND, NEL, NMAT, PROBLEM
INTEGER NGF, NGLXY, NGLX, NGLY
INTEGER NX, NY, NLF, NTC
INTEGER NG  !NÚMERO DE PUNTOS DE GAUSS.
INTEGER PARA    !PARA SABER SI SE IMPRIME O NO INFORMACIÓN PARA EL PARAVIEW.
INTEGER CONTA1, CONTA2, CONTA3, CONTA4, CONTA5, CONTA6
INTEGER NP1
REAL*8 JUAN
REAL*8 NUHS, BETAHS, RHOHS                  !PROPIEDADES DEL HALF-SPACE
REAL*8, ALLOCATABLE, DIMENSION(:):: X, Y    !COORDINATES POINTS.
INTEGER, ALLOCATABLE, DIMENSION(:):: IDN    !LABEL POINTS.
INTEGER, ALLOCATABLE, DIMENSION(:):: TP     !ELEMENT TYPE 2: TRIANGLE, 3: CUADRILATERAL
INTEGER, ALLOCATABLE, DIMENSION(:):: TPB    !AUXILIAR
INTEGER, ALLOCATABLE, DIMENSION(:):: TL     !NÚMERO DE LINEAS A LAS CUALES SE LES APLICA UN GRUPO FISICO.
INTEGER, ALLOCATABLE, DIMENSION(:):: TLRX   !NÚMERO DE LINEAS CON RESTRICCIONES EN X.
INTEGER, ALLOCATABLE, DIMENSION(:):: TLRY   !NÚMERO DE LINEAS CON RESTRICCIONES EN Y.
INTEGER, ALLOCATABLE, DIMENSION(:):: TLF    !NÚMERO DE LINEAS CARGADAS.
INTEGER, ALLOCATABLE, DIMENSION(:):: BOUNDARY, NODEIZ, NODEDER, NODECEN
INTEGER, ALLOCATABLE, DIMENSION(:):: PROP   !VECTOR PARA DETERMINAR LAS PROPIEDADES DE LOS ELEMENTOS.
INTEGER, ALLOCATABLE, DIMENSION(:,:):: MIE  !CONECTIVIDADES.
INTEGER, ALLOCATABLE, DIMENSION(:):: UX, UY !GRADOS DE LIBERTAD DESPLAZAMIENTOS.
INTEGER, ALLOCATABLE, DIMENSION(:,:):: NODEX, NODEY !GRADOS DE LIBERTAD RESTRINGIDOS EN CADA DIRECCIÓN.
INTEGER, ALLOCATABLE, DIMENSION(:,:):: NF   !NODOS CARGADOS
INTEGER, ALLOCATABLE, DIMENSION(:,:):: TPLOAD   !TIPO DE CARGA SOBRE EL NODO. HACE REFERENCIA AL NÚMERO DE CARGA Y SE COMPLEMENTA CON EL VALOR Y LA DIRECCIONALIDAD.
REAL*8, ALLOCATABLE, DIMENSION(:,:):: MATERIAL  !VECTOR CON LA INFORMACIÓN DE LOS MATERIALES A UTILIZAR. POISSON Y MÓDULO DE ELASTICIDAD.
REAL*8, ALLOCATABLE, DIMENSION(:):: BETA    !VELOCIDAD ONDA DE CORTE.
REAL*8, ALLOCATABLE, DIMENSION(:):: NU      !COEFICIENTE DE POISSON.
REAL*8, ALLOCATABLE, DIMENSION(:):: RHO     !DENSIDAD.
INTEGER NT      !NÚMERO DE TIEMPOS PARA EL PULSO DE RICKER.
REAL*8 FC       !FRECUENCIA CARACTERISTICA DEL PULSO DE RICKER
REAL*8 C        !CENTRO DEL PULSO DE RICKER.
REAL*8 TF       !TIEMPO TOTAL DE LA SEÑAL
INTEGER SUM
INTEGER, ALLOCATABLE, DIMENSION(:):: IDPR
INTEGER, ALLOCATABLE, DIMENSION(:):: U
REAL*8 X1, X2, X3, X4, Y1, Y2, Y3, Y4
REAL*8 LONG, ALT
!
!
INTEGER, ALLOCATABLE, DIMENSION(:,:):: EXPMESH
INTEGER, ALLOCATABLE, DIMENSION(:):: HMEPN

INTEGER, ALLOCATABLE, DIMENSION(:,:):: INCO1, INCO
INTEGER, ALLOCATABLE, DIMENSION(:):: POSREP

CHARACTER(15):: NCOMPLE
INTEGER NTHREADS
INTEGER POINTS_PER_THREAD
INTEGER THREAD_NUM
INTEGER ISTART, IEND
INTEGER NPL

INTEGER PARAL                                   !VARIABLE DONDE LE DIGO EL NÚMERO DE HILOS
INTEGER REP


REAL*8 elapsed_time
INTEGER tclock1, tclock2, clock_rate
call system_clock(tclock1)  ! start wall timer
!
!
OPEN(10,FILE='Entrada/canonquad.msh')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(1I10)') NP        !NUMBER OF NODES.
!
ALLOCATE(X(NP), Y(NP), IDN(NP))
!
DO I=1, NP
!
    READ(10,*) IDN(I), X(I), Y(I), JUAN
!
END DO
!
READ(10,*)
READ(10,*)
!
READ(10,'(1I10)') ND        !DATA TO READ.
!
!
!ALLOCATE(TPB(ND))
!
NEL=0
!TPB=1000
!
!ACÁ VOY A CONTAR LOS ELEMENTOS DE LA FRONTERA SOBRE LOS CUALES SE VA A PONER UN ELEMENTO
!ABSORBENTE. DE PASO SALE EL NÚMERO DE ELEMENTOS DEL PROBLEMA.
!
DO I=1, ND
!
    READ(10,*) J, L
!
    IF (L>8) THEN
!
        NEL=NEL+1
!
    END IF
!
END DO
!
CLOSE(10)
!
!
OPEN (10, FILE='Coord.txt')
!
DO I=1, NP
    WRITE(10, '(2F12.8)') X(I), Y(I)
END DO
!
CLOSE(10)
!
X=0.0D0
Y=0.0D0
!
OPEN (10, FILE='Coord.txt')
!
DO I=1, NP
    READ(10, '(2F12.8)') X(I), Y(I)
END DO
!
X=X/10.0D0
Y=Y/10.0D0
!
CLOSE(10, STATUS='DELETE')
!
LONG=MAXVAL(X)
ALT=MAXVAL(Y)
!
NPL=0
DO I=1, NP
    IF (Y(I)==ALT) THEN
        NPL=NPL+1
    END IF
END DO
!
!
!NEL:NUMERO DE ELEMENTOS.
!
OPEN(10,FILE='Entrada/canonquad.msh')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
!
DO I=1, NP
!
    READ(10,*)
!
END DO
!
READ(10,*)
READ(10,*)
READ(10,*)
!
ALLOCATE(BOUNDARY(ND-NEL), NODEIZ(ND-NEL), NODECEN(ND-NEL), &
         NODEDER(ND-NEL), TL(ND-NEL))
!
NODEIZ=0
NODEDER=0
NODECEN=0
!
NGLX=0
NGLY=0
!
!ACÁ YA IDENTIFICO POR COMPLETO LOS ELEMENTOS SOBRE LOS CUALES SE VAN A PONER LAS
!FRONTERAS ABSORBENTES.
!
DO I=1, ND-NEL
!
    READ(10,*) J, K, L, BOUNDARY(I), TL(I), NODEIZ(I), NODEDER(I), NODECEN(I)   
                            !LOS VECTORES NODEIZ Y NODEDER ME DICEN LOS NODOS CON 
                            !FRONTERAS ABSORBENTES.
!
END DO
!
!
DEALLOCATE(BOUNDARY)
!
!
!ACÁ VOY A ENCONTRAR LA MATRIZ INDICADORA DE ECUACIÓN Y LAS PROPIEDADES DE CADA
!ELEMENTO DEL DOMINIO
!
ALLOCATE(TP(NEL), MIE(9,NEL), PROP(NEL))
!
TP=0
MIE=0
PROP=0
!
DO I=1, NEL
!
    READ(10,*) J, TP(I), K, PROP(I), L, MIE(1,I), MIE(2,I), MIE(3,I), MIE(4,I), &
                                        MIE(5,I), MIE(6,I), MIE(7,I), MIE(8,I), &
                                        MIE(9,I)
    PROP(I)=PROP(I)/1000
!
END DO
!
!
NMAT=1
!
!
!NMAT: NÚMERO DE MATERIALES A UTILIZAR.
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
ALLOCATE(INCO1(NEL,3))
!
INCO1(:,1)=3
INCO1(:,2)=1
!
DO I=1, ND-NEL
!
    NX=NODEIZ(I)
    NY=NODEDER(I)
!
    DO J=1, NEL
!
        IF (NX==MIE(1,J) .AND. NY==MIE(2,J)) THEN
!
            INCO1(J,1)=5
            INCO1(J,2)=2
            INCO1(J,3)=3
            EXIT
!
        END IF
!
!
        IF (NX==MIE(2,J) .AND. NY==MIE(3,J)) THEN
!
            INCO1(J,1)=5
            INCO1(J,2)=3
            INCO1(J,3)=4
            EXIT
!
        END IF
!
!
        IF (NX==MIE(3,J) .AND. NY==MIE(4,J)) THEN
!
            INCO1(J,1)=5
            INCO1(J,2)=4
            INCO1(J,3)=1
            EXIT
!
        END IF
!
!
        IF (NX==MIE(4,J) .AND. NY==MIE(1,J)) THEN
!
            INCO1(J,1)=5
            INCO1(J,2)=5
            INCO1(J,3)=2
            EXIT
!
        END IF
!
!
    END DO
!
END DO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!
TP=TP+1     !AJUSTA EL TIPO DE ELEMENTO AL QUE NECESITA EL PROGRAMA DE FEM DE NOSOTROS
            !PUES PARA EL TRIANGLE:3 NODES Y CUADRILATERAL: 4 NODES
!
!
OPEN(10,FILE='../canon9.inp')
!
WRITE(10, '(1A20)') 'simple square mesh'
WRITE(10, '(2I16, 1I6, 1F6.2, 7I6)') NP, ND, 9, 25.00D0, 50001, 2, 9, 18, 5, 1, 0
!
DO I=1, NP
!
    WRITE(10, '(1I16, 3I6, 2F12.6)') I, 2, 0, 0, X(I), Y(I)
!
END DO
!
!
WRITE(10, '(3I4, 5F10.2, 1I4)') 1, 5, 0, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0
WRITE(10, '(3I4, 5F10.2, 1I4)') 2, 5, 1, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 3
WRITE(10, '(3I4, 5F10.2, 1I4)') 3, 5, 1, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 4
WRITE(10, '(3I4, 5F10.2, 1I4)') 4, 5, 1, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 1
WRITE(10, '(3I4, 5F10.2, 1I4)') 5, 5, 1, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 2
WRITE(10, '(3I4, 5F10.2, 1I4)') 6, 5, 1, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 5
WRITE(10, '(3I4, 5F10.2, 1I4)') 7, 5, 1, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 6
WRITE(10, '(3I4, 5F10.2, 1I4)') 8, 5, 1, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 7
WRITE(10, '(3I4, 5F10.2, 1I4)') 9, 5, 1, 2.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 8
!
!1: TIPO DE ONDA, 1: P-WAVE, 2: S-WAVE
!2: TIEMPO TOTAL DE LA SEÑAL
!3: TIEMPO CENTRAL DEL PUSLO
!4: FRECUENCIA CARACTERISTICA DEL PULSO
!5: AMPLITUD DE LA SEÑAL
!6: ÁNGULO DE INCIDENCIA DE LA ONDA
!
WRITE(10, '(1I6, 5F8.4)') 1, 25.0D0, 10.15D0, 8.0D0, 1.0D0, 0.0D0
!
!
DO I=1, NEL
!
    IF (MIE(1,I)==11) THEN
        WRITE(10,'(1I16, 4I8, 9I16)') I, INCO1(I,1), 18, 8, 9, MIE(1,I), MIE(2,I), &
                           MIE(3,I), MIE(4,I), MIE(5,I), MIE(6,I), &
                           MIE(7,I), MIE(8,I), MIE(9,I)
                           
        WRITE(*,'(1I16, 4I8, 9I16)')I, INCO1(I,1), 18, 8, 9, MIE(1,I), MIE(2,I), &
                           MIE(3,I), MIE(4,I), MIE(5,I), MIE(6,I), &
                           MIE(7,I), MIE(8,I), MIE(9,I)
    ELSE
        IF (MIE(3,I)==8) THEN
            WRITE(10,'(1I16, 4I8, 9I16)') I, INCO1(I,1), 18, 6, 9, MIE(1,I), &
                   MIE(2,I), MIE(3,I), MIE(4,I), MIE(5,I), MIE(6,I), &
                   MIE(7,I), MIE(8,I), MIE(9,I)
                   
            WRITE(*,'(1I16, 4I8, 9I16)')I, INCO1(I,1), 18, 8, 9,MIE(1,I),MIE(2,I),&
                           MIE(3,I), MIE(4,I), MIE(5,I), MIE(6,I), &
                           MIE(7,I), MIE(8,I), MIE(9,I)
        ELSE
            WRITE(10,'(1I16, 4I8, 9I16)') I, INCO1(I,1), 18, INCO1(I,2), 9, &
                   MIE(1,I), MIE(2,I), &
                   MIE(3,I), MIE(4,I), MIE(5,I), MIE(6,I), &
                   MIE(7,I), MIE(8,I), MIE(9,I)
        END IF
    END IF
END DO
!
DO I=1, ND-NEL
!
    WRITE(10,'(1I16, 4I8, 3I16)') I+NEL, 7, 6, 1, 3, &
                                    NODEDER(I), NODECEN(I), NODEIZ(I)
!
END DO
!
! DO I=1, NP
!     IF (Y(I)==ALT) THEN
!         WRITE(10, '(2I6)') I, 2
!     END IF
! END DO
!
!
call system_clock(tclock2, clock_rate)
elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
print 11, elapsed_time
11 format("Elapsed time = ",f12.4, " seconds")

END PROGRAM MALLAS
























