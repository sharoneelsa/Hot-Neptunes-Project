include $(MESA_DIR)/utils/makefile_header

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_JOB_DIR = $(MESA_DIR)/star/job

RUN_STAR_EXTRA_OBJS =
RUN_STAR_EXTRA_OBJS += custom_eos.o
# RUN_STAR_EXTRA_OBJS += custom_kap.o   # uncomment if you need it

#OBJS = $(RUN_STAR_EXTRA_OBJS) $(STAR_WORK_OBJS) run_star_extras.o run_star.o run.o
OBJS = $(RUN_STAR_EXTRA_OBJS) $(STAR_WORK_OBJS) run.o

WORK_COMPILE = \
   $(LOADER) $(FCbasic) $(FCopenmp) $(FCchecks) $(FCopt) $(FCdebug) $(FCfree) $(FCstandard) \
   -I$(MESA_INCLUDE_DIR) -I$(WORK_SRC_DIR) -I$(MESA_DIR)/star/make -c

STAR = star

all: $(STAR)

$(STAR): $(OBJS)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(STAR) $(OBJS) $(LOAD_MESA_STAR)

%.o: $(WORK_SRC_DIR)/%.f90
	$(WORK_COMPILE) $<

%.o: $(STAR_JOB_DIR)/%.f90
	$(WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod *.smod *.so *.a $(WORK_DIR)/$(STAR)
