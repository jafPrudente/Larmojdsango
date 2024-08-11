FLAGS = -O2 -ffree-form -w 

OBJS = arrays.o saverRT.o saverT.o evolveDirac.o evolveEscalar.o evolveProca.o fuenteDirac.o fuenteEscalar.o \
		fuenteProca.o initialDirac.o initialEscalar.o initialProca.o main.o vars.o metricEscalar.o \
		metricDirac.o metricProca.o perturbacionEscalar.o perturbacionDirac.o

exe : folder link

folder : 
	@ mkdir -p data

link: $(OBJS)
	@ echo 'Compilando (...)'
	@ gfortran $(FLAGS) -o estrella $(OBJS)
	@ echo 'Limpiando archivos extra (...)'
	@ rm *.o
	@ echo 'Moviendo el ejecutable a la carpeta 'data' (...)'
	@ mv estrella data
	@ echo 'Copiando los parámetros de entrada a la carpeta 'data' (...)'
	@ cp input.par data
	@ echo 'Todo está listo para ejecutar UwUr'

clean :
	@ echo 'Limpiando el espacio de trabajo (...)'
	@ rm -r data
	@ echo 'Todo quedó como antes de compilar UwUr'

arrays.o : arrays.f90
	@ gfortran $(FLAGS) -c arrays.f90

saverRT.o : saverRT.f90
	@ gfortran $(FLAGS) -c saverRT.f90

saverT.o : saverT.f90
	@ gfortran $(FLAGS) -c saverT.f90

evolveDirac.o : evolveDirac.f90
	@ gfortran $(FLAGS) -c evolveDirac.f90

evolveEscalar.o : evolveEscalar.f90
	@ gfortran $(FLAGS) -c evolveEscalar.f90

evolveProca.o : evolveProca.f90
	@ gfortran $(FLAGS) -c evolveProca.f90

fuenteDirac.o : fuenteDirac.f90
	@ gfortran $(FLAGS) -c fuenteDirac.f90

fuenteEscalar.o : fuenteEscalar.f90
	@ gfortran $(FLAGS) -c fuenteEscalar.f90

fuenteProca.o : fuenteProca.f90
	@ gfortran $(FLAGS) -c fuenteProca.f90

initialDirac.o : initialDirac.f90
	@ gfortran $(FLAGS) -c initialDirac.f90

initialEscalar.o : initialEscalar.f90
	@ gfortran $(FLAGS) -c initialEscalar.f90

initialProca.o : initialProca.f90
	@ gfortran $(FLAGS) -c initialProca.f90

metricEscalar.o : metricEscalar.f90
	@ gfortran $(FLAGS) -c metricEscalar.f90

metricDirac.o : metricDirac.f90
	@ gfortran $(FLAGS) -c metricDirac.f90

metricProca.o : metricProca.f90
	@ gfortran $(FLAGS) -c metricProca.f90

perturbacionEscalar.o : perturbacionEscalar.f90
	@ gfortran $(FLAGS) -c perturbacionEscalar.f90

perturbacionDirac.o : perturbacionDirac.f90
	@ gfortran $(FLAGS) -c perturbacionDirac.f90

main.o : main.f90
	@ gfortran $(FLAGS) -c main.f90

vars.o : vars.f90
	@ gfortran $(FLAGS) -c vars.f90