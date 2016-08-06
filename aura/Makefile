# Makefile for Aura

F90 = ifort

FLAGS = -align all -ipo

OBJS = mod_phys.f90 mod_infodefs.f90 mod_info.f90 mod_sync.f90 mod_ic.f90 mod_brem.f90 mod_pp.f90 mod_emis.f90 mod_io.f90 mod_user.f90 aura.f90

default: aura

aura:	$(OBJS)
	$(F90) $(FLAGS) $(OBJS) -o aura

clean:
	rm *.mod
	rm aura
