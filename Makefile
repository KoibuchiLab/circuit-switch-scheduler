#
# Makefile	Feb-26-2019 huyao@nii.ac.jp
#
# GNU make
#
# v1.0

#compiler
CC = g++

#object files
objects = circuit-switch-scheduler.o 

all : css

css : circuit-switch-scheduler.o 
	${CC} circuit-switch-scheduler.o -pthread -o css.out

clean : 
	-rm ${objects}

