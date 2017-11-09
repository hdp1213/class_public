CC = icc
FC = ifort #-nofor_main
AR = ar rvs
CCFLAG = -O3 -march=native# -DCAMB -DHYRECPATH=\"$(PWD)/\" 
FCFLAG = -O3 -march=native
LDFLAG = -O3 -march=native
DBUG = -DDBUG# -DPBH_DEBUG

%.o: %.c 
	$(CC) $(CCFLAG) $(DBUG) -c $*.c -o $*.o

vpath %.f dierckx

%.o: %.f
	$(FC) $(FCFLAG) $(DBUG) -c $< -o $*.o

HYREC_SRC = hyrectools.o helium.o hydrogen.o history.o energy_injection.o hyrec_21cm.o bispev.o fpbisp.o fpbspl.o pbh.o
HYREC_EXE = hyrec.o 

default: libhyrec.a

clean: 
	rm -rf *.o

hyrec: $(HYREC_SRC) $(HYREC_EXE) 
	$(CC) $(LDFLAG) $(DBUG) -o hyrec $^ -lm


libhyrec.a: $(HYREC_SRC)
	$(AR) $@ $(HYREC_SRC)
