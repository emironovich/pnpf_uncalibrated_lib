###########################################################################
## Makefile generated for MATLAB file/project 'p35p_solver'. 
## 
## Makefile     : p35p_solver_rtw.mk
## Generated on : Fri Oct 04 01:44:20 2019
## MATLAB Coder version: 4.3 (R2019b)
## 
## Build Info:
## 
## Final product: ./p35p_solver.a
## Product type : static-library
## 
###########################################################################

###########################################################################
## MACROS
###########################################################################

# Macro Descriptions:
# PRODUCT_NAME            Name of the system to build
# MAKEFILE                Name of this makefile
# MODELLIB                Static library target

PRODUCT_NAME              = p35p_solver
MAKEFILE                  = p35p_solver_rtw.mk
MATLAB_ROOT               = /home/elizaveta/Matlab
MATLAB_BIN                = /home/elizaveta/Matlab/bin
MATLAB_ARCH_BIN           = $(MATLAB_BIN)/glnxa64
MASTER_ANCHOR_DIR         = 
START_DIR                 = /home/elizaveta/WORK/p3.5p/P3.5P/codegen/lib/p35p_solver
TGT_FCN_LIB               = ISO_C++
SOLVER_OBJ                = 
CLASSIC_INTERFACE         = 0
MODEL_HAS_DYNAMICALLY_LOADED_SFCNS = 
RELATIVE_PATH_TO_ANCHOR   = .
C_STANDARD_OPTS           = -fwrapv -ansi -pedantic -Wno-long-long
CPP_STANDARD_OPTS         = -fwrapv -std=c++03 -pedantic -Wno-long-long
MODELLIB                  = p35p_solver.a

###########################################################################
## TOOLCHAIN SPECIFICATIONS
###########################################################################

# Toolchain Name:          GNU gcc/g++ | gmake (64-bit Linux)
# Supported Version(s):    
# ToolchainInfo Version:   2019b
# Specification Revision:  1.0
# 
#-------------------------------------------
# Macros assumed to be defined elsewhere
#-------------------------------------------

# C_STANDARD_OPTS
# CPP_STANDARD_OPTS

#-----------
# MACROS
#-----------

WARN_FLAGS         = -Wall -W -Wwrite-strings -Winline -Wstrict-prototypes -Wnested-externs -Wpointer-arith -Wcast-align
WARN_FLAGS_MAX     = $(WARN_FLAGS) -Wcast-qual -Wshadow
CPP_WARN_FLAGS     = -Wall -W -Wwrite-strings -Winline -Wpointer-arith -Wcast-align
CPP_WARN_FLAGS_MAX = $(CPP_WARN_FLAGS) -Wcast-qual -Wshadow

TOOLCHAIN_SRCS = 
TOOLCHAIN_INCS = 
TOOLCHAIN_LIBS = 

#------------------------
# BUILD TOOL COMMANDS
#------------------------

# C Compiler: GNU C Compiler
CC = gcc

# Linker: GNU Linker
LD = g++

# C++ Compiler: GNU C++ Compiler
CPP = g++

# C++ Linker: GNU C++ Linker
CPP_LD = g++

# Archiver: GNU Archiver
AR = ar

# MEX Tool: MEX Tool
MEX_PATH = $(MATLAB_ARCH_BIN)
MEX = "$(MEX_PATH)/mex"

# Download: Download
DOWNLOAD =

# Execute: Execute
EXECUTE = $(PRODUCT)

# Builder: GMAKE Utility
MAKE_PATH = %MATLAB%/bin/glnxa64
MAKE = "$(MAKE_PATH)/gmake"


#-------------------------
# Directives/Utilities
#-------------------------

CDEBUG              = -g
C_OUTPUT_FLAG       = -o
LDDEBUG             = -g
OUTPUT_FLAG         = -o
CPPDEBUG            = -g
CPP_OUTPUT_FLAG     = -o
CPPLDDEBUG          = -g
OUTPUT_FLAG         = -o
ARDEBUG             =
STATICLIB_OUTPUT_FLAG =
MEX_DEBUG           = -g
RM                  = @rm -f
ECHO                = @echo
MV                  = @mv
RUN                 =

#----------------------------------------
# "Faster Builds" Build Configuration
#----------------------------------------

ARFLAGS              = ruvs
CFLAGS               = -c $(C_STANDARD_OPTS) -fPIC \
                       -O0
CPPFLAGS             = -c $(CPP_STANDARD_OPTS) -fPIC \
                       -O0
CPP_LDFLAGS          = -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)"
CPP_SHAREDLIB_LDFLAGS  = -shared -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)" -Wl,--no-undefined
DOWNLOAD_FLAGS       =
EXECUTE_FLAGS        =
LDFLAGS              = -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)"
MEX_CPPFLAGS         =
MEX_CPPLDFLAGS       =
MEX_CFLAGS           =
MEX_LDFLAGS          =
MAKE_FLAGS           = -f $(MAKEFILE)
SHAREDLIB_LDFLAGS    = -shared -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)" -Wl,--no-undefined



###########################################################################
## OUTPUT INFO
###########################################################################

PRODUCT = ./p35p_solver.a
PRODUCT_TYPE = "static-library"
BUILD_TYPE = "Static Library"

###########################################################################
## INCLUDE PATHS
###########################################################################

INCLUDES_BUILDINFO = -I$(START_DIR) -I/home/elizaveta/WORK/p3.5p/P3.5P -I$(MATLAB_ROOT)/extern/include -I$(MATLAB_ROOT)/simulink/include -I$(MATLAB_ROOT)/rtw/c/src -I$(MATLAB_ROOT)/rtw/c/src/ext_mode/common -I$(MATLAB_ROOT)/rtw/c/ert

INCLUDES = $(INCLUDES_BUILDINFO)

###########################################################################
## DEFINES
###########################################################################

DEFINES_CUSTOM = 
DEFINES_STANDARD = -DMODEL=p35p_solver -DHAVESTDIO -DUSE_RTMODEL -DUNIX

DEFINES = $(DEFINES_CUSTOM) $(DEFINES_STANDARD)

###########################################################################
## SOURCE FILES
###########################################################################

SRCS = $(START_DIR)/rt_nonfinite.cpp $(START_DIR)/rtGetNaN.cpp $(START_DIR)/rtGetInf.cpp $(START_DIR)/p35p_solver_rtwutil.cpp $(START_DIR)/p35p_solver_data.cpp $(START_DIR)/p35p_solver_initialize.cpp $(START_DIR)/p35p_solver_terminate.cpp $(START_DIR)/p35p_solver.cpp $(START_DIR)/quadruple_constraint.cpp $(START_DIR)/equations_for_groebner.cpp $(START_DIR)/mult_poly22.cpp $(START_DIR)/mult_poly42.cpp $(START_DIR)/mult_for_groebner.cpp $(START_DIR)/mldivide.cpp $(START_DIR)/eig.cpp $(START_DIR)/schur.cpp $(START_DIR)/xnrm2.cpp $(START_DIR)/sqrt.cpp $(START_DIR)/xzlarf.cpp $(START_DIR)/xgerc.cpp $(START_DIR)/xdhseqr.cpp $(START_DIR)/xdlanv2.cpp $(START_DIR)/xrot.cpp $(START_DIR)/xzggev.cpp $(START_DIR)/xzggbal.cpp $(START_DIR)/xzgghrd.cpp $(START_DIR)/xzlartg.cpp $(START_DIR)/xzhgeqz.cpp $(START_DIR)/xztgevc.cpp $(START_DIR)/xzggbak.cpp $(START_DIR)/find_f.cpp $(START_DIR)/qr.cpp $(START_DIR)/cat.cpp

ALL_SRCS = $(SRCS)

###########################################################################
## OBJECTS
###########################################################################

OBJS = rt_nonfinite.o rtGetNaN.o rtGetInf.o p35p_solver_rtwutil.o p35p_solver_data.o p35p_solver_initialize.o p35p_solver_terminate.o p35p_solver.o quadruple_constraint.o equations_for_groebner.o mult_poly22.o mult_poly42.o mult_for_groebner.o mldivide.o eig.o schur.o xnrm2.o sqrt.o xzlarf.o xgerc.o xdhseqr.o xdlanv2.o xrot.o xzggev.o xzggbal.o xzgghrd.o xzlartg.o xzhgeqz.o xztgevc.o xzggbak.o find_f.o qr.o cat.o

ALL_OBJS = $(OBJS)

###########################################################################
## PREBUILT OBJECT FILES
###########################################################################

PREBUILT_OBJS = 

###########################################################################
## LIBRARIES
###########################################################################

LIBS = 

###########################################################################
## SYSTEM LIBRARIES
###########################################################################

SYSTEM_LIBS =  -lm -lstdc++

###########################################################################
## ADDITIONAL TOOLCHAIN FLAGS
###########################################################################

#---------------
# C Compiler
#---------------

CFLAGS_BASIC = $(DEFINES) $(INCLUDES)

CFLAGS += $(CFLAGS_BASIC)

#-----------------
# C++ Compiler
#-----------------

CPPFLAGS_BASIC = $(DEFINES) $(INCLUDES)

CPPFLAGS += $(CPPFLAGS_BASIC)

###########################################################################
## INLINED COMMANDS
###########################################################################

###########################################################################
## PHONY TARGETS
###########################################################################

.PHONY : all build clean info prebuild download execute


all : build
	@echo "### Successfully generated all binary outputs."


build : prebuild $(PRODUCT)


prebuild : 


download : $(PRODUCT)


execute : download


###########################################################################
## FINAL TARGET
###########################################################################

#---------------------------------
# Create a static library         
#---------------------------------

$(PRODUCT) : $(OBJS) $(PREBUILT_OBJS)
	@echo "### Creating static library "$(PRODUCT)" ..."
	$(AR) $(ARFLAGS)  $(PRODUCT) $(OBJS)
	@echo "### Created: $(PRODUCT)"


###########################################################################
## INTERMEDIATE TARGETS
###########################################################################

#---------------------
# SOURCE-TO-OBJECT
#---------------------

%.o : %.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : %.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : /home/elizaveta/WORK/p3.5p/P3.5P/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : /home/elizaveta/WORK/p3.5p/P3.5P/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(MATLAB_ROOT)/rtw/c/src/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : $(MATLAB_ROOT)/rtw/c/src/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


rt_nonfinite.o : $(START_DIR)/rt_nonfinite.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


rtGetNaN.o : $(START_DIR)/rtGetNaN.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


rtGetInf.o : $(START_DIR)/rtGetInf.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


p35p_solver_rtwutil.o : $(START_DIR)/p35p_solver_rtwutil.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


p35p_solver_data.o : $(START_DIR)/p35p_solver_data.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


p35p_solver_initialize.o : $(START_DIR)/p35p_solver_initialize.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


p35p_solver_terminate.o : $(START_DIR)/p35p_solver_terminate.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


p35p_solver.o : $(START_DIR)/p35p_solver.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


quadruple_constraint.o : $(START_DIR)/quadruple_constraint.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


equations_for_groebner.o : $(START_DIR)/equations_for_groebner.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


mult_poly22.o : $(START_DIR)/mult_poly22.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


mult_poly42.o : $(START_DIR)/mult_poly42.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


mult_for_groebner.o : $(START_DIR)/mult_for_groebner.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


mldivide.o : $(START_DIR)/mldivide.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


eig.o : $(START_DIR)/eig.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


schur.o : $(START_DIR)/schur.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xnrm2.o : $(START_DIR)/xnrm2.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


sqrt.o : $(START_DIR)/sqrt.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xzlarf.o : $(START_DIR)/xzlarf.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xgerc.o : $(START_DIR)/xgerc.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xdhseqr.o : $(START_DIR)/xdhseqr.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xdlanv2.o : $(START_DIR)/xdlanv2.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xrot.o : $(START_DIR)/xrot.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xzggev.o : $(START_DIR)/xzggev.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xzggbal.o : $(START_DIR)/xzggbal.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xzgghrd.o : $(START_DIR)/xzgghrd.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xzlartg.o : $(START_DIR)/xzlartg.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xzhgeqz.o : $(START_DIR)/xzhgeqz.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xztgevc.o : $(START_DIR)/xztgevc.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


xzggbak.o : $(START_DIR)/xzggbak.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


find_f.o : $(START_DIR)/find_f.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


qr.o : $(START_DIR)/qr.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


cat.o : $(START_DIR)/cat.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


###########################################################################
## DEPENDENCIES
###########################################################################

$(ALL_OBJS) : rtw_proj.tmw $(MAKEFILE)


###########################################################################
## MISCELLANEOUS TARGETS
###########################################################################

info : 
	@echo "### PRODUCT = $(PRODUCT)"
	@echo "### PRODUCT_TYPE = $(PRODUCT_TYPE)"
	@echo "### BUILD_TYPE = $(BUILD_TYPE)"
	@echo "### INCLUDES = $(INCLUDES)"
	@echo "### DEFINES = $(DEFINES)"
	@echo "### ALL_SRCS = $(ALL_SRCS)"
	@echo "### ALL_OBJS = $(ALL_OBJS)"
	@echo "### LIBS = $(LIBS)"
	@echo "### MODELREF_LIBS = $(MODELREF_LIBS)"
	@echo "### SYSTEM_LIBS = $(SYSTEM_LIBS)"
	@echo "### TOOLCHAIN_LIBS = $(TOOLCHAIN_LIBS)"
	@echo "### CFLAGS = $(CFLAGS)"
	@echo "### LDFLAGS = $(LDFLAGS)"
	@echo "### SHAREDLIB_LDFLAGS = $(SHAREDLIB_LDFLAGS)"
	@echo "### CPPFLAGS = $(CPPFLAGS)"
	@echo "### CPP_LDFLAGS = $(CPP_LDFLAGS)"
	@echo "### CPP_SHAREDLIB_LDFLAGS = $(CPP_SHAREDLIB_LDFLAGS)"
	@echo "### ARFLAGS = $(ARFLAGS)"
	@echo "### MEX_CFLAGS = $(MEX_CFLAGS)"
	@echo "### MEX_CPPFLAGS = $(MEX_CPPFLAGS)"
	@echo "### MEX_LDFLAGS = $(MEX_LDFLAGS)"
	@echo "### MEX_CPPLDFLAGS = $(MEX_CPPLDFLAGS)"
	@echo "### DOWNLOAD_FLAGS = $(DOWNLOAD_FLAGS)"
	@echo "### EXECUTE_FLAGS = $(EXECUTE_FLAGS)"
	@echo "### MAKE_FLAGS = $(MAKE_FLAGS)"


clean : 
	$(ECHO) "### Deleting all derived files..."
	$(RM) $(PRODUCT)
	$(RM) $(ALL_OBJS)
	$(ECHO) "### Deleted all derived files."


