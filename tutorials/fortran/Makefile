
SRC=testFilter1d.f
OBJS=$(SRC:.f=.o)

goTestFilter1d: $(OBJS)
	$(FC) $(OBJS) -o @

.o:.f
	$(FC) $(FOPT) -c < @
