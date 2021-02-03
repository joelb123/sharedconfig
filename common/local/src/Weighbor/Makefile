#
# must have ansi compiler if not
# gcc then substitute something else
#
#
CC=cc

#
# Set flags to be given to the complier. Currently set to optimize
#

CFLAGS=-O3 -ffast-math

##########################################################################
##########################################################################
####
#### End of user editable options
####
#### There should be no need to edit anything past here
####
##########################################################################


WOBJS = tree.o matrix.o io.o build.o noise.o calcq.o calcb.o calcz.o calerf.o 

weighbor: weighbor.o $(WOBJS)
	$(CC) $(CFLAGS) -o weighbor weighbor.o $(WOBJS) -lm


clean: 
	rm -rf *.o core *~ *.dvi *.log *.toc *.bak *.ps Doc/*~

spotless:	
	rm -rf *.o core *~ *.dvi *.log *.toc *.bak *.ps
	rm -rf weighbor weighbor.out Doc/*~
	rm -rf config.* configure* weighbor.tar

# DO NOT DELETE

build.o: tree.h matrix.h weighbor.h io.h
calcb.o: tree.h matrix.h weighbor.h
calcq.o: tree.h matrix.h weighbor.h
calcz.o: tree.h matrix.h weighbor.h
io.o: weighbor.h matrix.h tree.h
matrix.o: matrix.h
noise.o: weighbor.h matrix.h tree.h
tree.o: tree.h matrix.h
weighbor.o: tree.h matrix.h io.h weighbor.h
