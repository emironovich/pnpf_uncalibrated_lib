###########################################################################
## Makefile generated for MATLAB file/project 'p35p_solver'. 
## 
## Makefile     : p35p_solver_rtw.mk
## Generated on : Mon Sep 30 04:08:11 2019
## MATLAB Coder version: 4.2 (R2019a)
## 
## Build Info:
## 
## Final product: ./p35p_solver_double.lib
## Product type : static-library
## 
###########################################################################

###########################################################################
## MACROS
###########################################################################

# Macro Descriptions:
# PRODUCT_NAME            Name of the system to build
# MAKEFILE                Name of this makefile
# COMPUTER                Computer type. See the MATLAB "computer" command.
# COMPILER_COMMAND_FILE   Compiler command listing model reference header paths
# CMD_FILE                Command file

PRODUCT_NAME              = p35p_solver
MAKEFILE                  = p35p_solver_rtw.mk
COMPUTER                  = PCWIN64
MATLAB_ROOT               = C:/PROGRA~1/MATLAB/R2019a
MATLAB_BIN                = C:/PROGRA~1/MATLAB/R2019a/bin
MATLAB_ARCH_BIN           = $(MATLAB_BIN)/win64
MASTER_ANCHOR_DIR         = 
START_DIR                 = C:/Users/Елизавета/UNI/WORK/p3.5p/P3.5P/codegen/lib/p35p_solver
ARCH                      = win64
TGT_FCN_LIB               = ISO_C++
RELATIVE_PATH_TO_ANCHOR   = .
COMPILER_COMMAND_FILE     = p35p_solver_rtw_comp.rsp
CMD_FILE                  = p35p_solver_rtw.rsp
C_STANDARD_OPTS           = -ansi -pedantic -Wno-long-long -fwrapv
CPP_STANDARD_OPTS         = -std=c++98 -pedantic -Wno-long-long -fwrapv

###########################################################################
## TOOLCHAIN SPECIFICATIONS
###########################################################################

# Toolchain Name:          MinGW64 | gmake (64-bit Windows)
# Supported Version(s):    6.x
# ToolchainInfo Version:   R2019a
# Specification Revision:  1.0
# 
#-------------------------------------------
# Macros assumed to be defined elsewhere
#-------------------------------------------

# C_STANDARD_OPTS
# CPP_STANDARD_OPTS
# MINGW_ROOT
# MINGW_C_STANDARD_OPTS

#-----------
# MACROS
#-----------

WARN_FLAGS            = -Wall -W -Wwrite-strings -Winline -Wstrict-prototypes -Wnested-externs -Wpointer-arith -Wcast-align
WARN_FLAGS_MAX        = $(WARN_FLAGS) -Wcast-qual -Wshadow
CPP_WARN_FLAGS        = -Wall -W -Wwrite-strings -Winline -Wpointer-arith -Wcast-align
CPP_WARN_FLAGS_MAX    = $(CPP_WARN_FLAGS) -Wcast-qual -Wshadow
MEX_OPTS_FILE         = $(MATLAB_ROOT)/bin/win64/mexopts/mingw64.xml
MEX_CPP_OPTS_FILE     = $(MATLAB_ROOT)/bin/win64/mexopts/mingw64_g++.xml
MW_EXTERNLIB_DIR      = $(MATLAB_ROOT)/extern/lib/win64/mingw64
SHELL                 = %SystemRoot%/system32/cmd.exe

TOOLCHAIN_SRCS = 
TOOLCHAIN_INCS = 
TOOLCHAIN_LIBS = -lws2_32

#------------------------
# BUILD TOOL COMMANDS
#------------------------

# C Compiler: GNU C Compiler
CC_PATH = $(MINGW_ROOT)
CC = "$(CC_PATH)/gcc"

# Linker: GNU Linker
LD_PATH = $(MINGW_ROOT)
LD = "$(LD_PATH)/g++"

# C++ Compiler: GNU C++ Compiler
CPP_PATH = $(MINGW_ROOT)
CPP = "$(CPP_PATH)/g++"

# C++ Linker: GNU C++ Linker
CPP_LD_PATH = $(MINGW_ROOT)
CPP_LD = "$(CPP_LD_PATH)/g++"

# Archiver: GNU Archiver
AR_PATH = $(MINGW_ROOT)
AR = "$(AR_PATH)/ar"

# MEX Tool: MEX Tool
MEX_PATH = $(MATLAB_BIN)/win64
MEX = "$(MEX_PATH)/mex"

# Download: Download
DOWNLOAD =

# Execute: Execute
EXECUTE = $(PRODUCT)

# Builder: GMAKE Utility
MAKE_PATH = %MATLAB%/bin/win64
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
RM                  = @del
ECHO                = @echo
MV                  = @move
RUN                 =

#----------------------------------------
# "Faster Builds" Build Configuration
#----------------------------------------

ARFLAGS              = ruvs
CFLAGS               = -c $(MINGW_C_STANDARD_OPTS) -m64 \
                       -O0
CPPFLAGS             = -c $(CPP_STANDARD_OPTS) -m64 \
                       -O0
CPP_LDFLAGS          = -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)" -static -m64
CPP_SHAREDLIB_LDFLAGS  = -shared -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)" -Wl,--no-undefined \
                         -Wl,--out-implib,$(basename $(PRODUCT)).lib
DOWNLOAD_FLAGS       =
EXECUTE_FLAGS        =
LDFLAGS              = -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)" -static -m64
MEX_CPPFLAGS         = -R2018a -MATLAB_ARCH=$(ARCH) $(INCLUDES) \
                         \
                       CXXOPTIMFLAGS="$(MINGW_C_STANDARD_OPTS)  \
                       -O0 \
                        $(DEFINES)" \
                         \
                       -silent
MEX_CPPLDFLAGS       = LDFLAGS=='$$LDFLAGS'
MEX_CFLAGS           = -R2018a -MATLAB_ARCH=$(ARCH) $(INCLUDES) \
                         \
                       COPTIMFLAGS="$(MINGW_C_STANDARD_OPTS)  \
                       -O0 \
                        $(DEFINES)" \
                         \
                       -silent
MEX_LDFLAGS          = LDFLAGS=='$$LDFLAGS'
MAKE_FLAGS           = -f $(MAKEFILE)
SHAREDLIB_LDFLAGS    = -shared -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)" -Wl,--no-undefined \
                       -Wl,--out-implib,$(basename $(PRODUCT)).lib



###########################################################################
## OUTPUT INFO
###########################################################################

PRODUCT = ./p35p_solver_double.lib
PRODUCT_TYPE = "static-library"
BUILD_TYPE = "Static Library"

###########################################################################
## INCLUDE PATHS
###########################################################################

INCLUDES_BUILDINFO = 

INCLUDES = $(INCLUDES_BUILDINFO)

###########################################################################
## DEFINES
###########################################################################

DEFINES_ = -D__USE_MINGW_ANSI_STDIO=1 -DMODEL=p35p_solver_double -DHAVESTDIO -DUSE_RTMODEL
DEFINES_CUSTOM = 
DEFINES_STANDARD = -DMODEL=p35p_solver_double -DHAVESTDIO -DUSE_RTMODEL

DEFINES = $(DEFINES_) $(DEFINES_CUSTOM) $(DEFINES_STANDARD)

###########################################################################
## SOURCE FILES
###########################################################################

SRCS = $(START_DIR)/p35p_solver_rtwutil.cpp $(START_DIR)/p35p_solver_initialize.cpp $(START_DIR)/p35p_solver_terminate.cpp $(START_DIR)/p35p_solver.cpp $(START_DIR)/init_F.cpp $(START_DIR)/quadruple_constraint.cpp $(START_DIR)/equations_for_groebner.cpp $(START_DIR)/mult_poly22.cpp $(START_DIR)/mult_poly42.cpp $(START_DIR)/mult_for_groebner.cpp $(START_DIR)/mldivide.cpp $(START_DIR)/xgetrf.cpp $(START_DIR)/xger.cpp $(START_DIR)/make_mult_matrix.cpp $(START_DIR)/eig.cpp $(START_DIR)/anyNonFinite.cpp $(START_DIR)/isfinite.cpp $(START_DIR)/ishermitian.cpp $(START_DIR)/schur.cpp $(START_DIR)/triu.cpp $(START_DIR)/xgehrd.cpp $(START_DIR)/xzlarfg.cpp $(START_DIR)/xnrm2.cpp $(START_DIR)/sqrt.cpp $(START_DIR)/xzlarf.cpp $(START_DIR)/xgemv.cpp $(START_DIR)/xgerc.cpp $(START_DIR)/xungorghr.cpp $(START_DIR)/xhseqr.cpp $(START_DIR)/xdhseqr.cpp $(START_DIR)/xdlanv2.cpp $(START_DIR)/rsf2csf.cpp $(START_DIR)/xzgeev.cpp $(START_DIR)/xzggev.cpp $(START_DIR)/xzlangeM.cpp $(START_DIR)/xzlascl.cpp $(START_DIR)/xzggbal.cpp $(START_DIR)/xzgghrd.cpp $(START_DIR)/xzlartg.cpp $(START_DIR)/xzhgeqz.cpp $(START_DIR)/xztgevc.cpp $(START_DIR)/xzggbak.cpp $(START_DIR)/find_f.cpp $(START_DIR)/qr.cpp $(START_DIR)/cat.cpp $(START_DIR)/rt_nonfinite.cpp $(START_DIR)/rtGetNaN.cpp $(START_DIR)/rtGetInf.cpp

ALL_SRCS = $(SRCS)

###########################################################################
## OBJECTS
###########################################################################

OBJS = p35p_solver_rtwutil.obj p35p_solver_initialize.obj p35p_solver_terminate.obj p35p_solver.obj init_F.obj quadruple_constraint.obj equations_for_groebner.obj mult_poly22.obj mult_poly42.obj mult_for_groebner.obj mldivide.obj xgetrf.obj xger.obj make_mult_matrix.obj eig.obj anyNonFinite.obj isfinite.obj ishermitian.obj schur.obj triu.obj xgehrd.obj xzlarfg.obj xnrm2.obj sqrt.obj xzlarf.obj xgemv.obj xgerc.obj xungorghr.obj xhseqr.obj xdhseqr.obj xdlanv2.obj rsf2csf.obj xzgeev.obj xzggev.obj xzlangeM.obj xzlascl.obj xzggbal.obj xzgghrd.obj xzlartg.obj xzhgeqz.obj xztgevc.obj xzggbak.obj find_f.obj qr.obj cat.obj rt_nonfinite.obj rtGetNaN.obj rtGetInf.obj

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

SYSTEM_LIBS = 

###########################################################################
## ADDITIONAL TOOLCHAIN FLAGS
###########################################################################

#---------------
# C Compiler
#---------------

CFLAGS_BASIC = $(DEFINES) $(INCLUDES) @$(COMPILER_COMMAND_FILE)

CFLAGS += $(CFLAGS_BASIC)

#-----------------
# C++ Compiler
#-----------------

CPPFLAGS_BASIC = $(DEFINES) $(INCLUDES) @$(COMPILER_COMMAND_FILE)

CPPFLAGS += $(CPPFLAGS_BASIC)

#---------------------
# MEX C++ Compiler
#---------------------

MEX_CPP_Compiler_BASIC =  @$(COMPILER_COMMAND_FILE)

MEX_CPPFLAGS += $(MEX_CPP_Compiler_BASIC)

#-----------------
# MEX Compiler
#-----------------

MEX_Compiler_BASIC =  @$(COMPILER_COMMAND_FILE)

MEX_CFLAGS += $(MEX_Compiler_BASIC)

###########################################################################
## INLINED COMMANDS
###########################################################################


ifdef SIM_TARGET_BUILD
MINGW_C_STANDARD_OPTS = $(filter-out -ansi, $(C_STANDARD_OPTS))
else
MINGW_C_STANDARD_OPTS = $(C_STANDARD_OPTS)
endif


###########################################################################
## PHONY TARGETS
###########################################################################

.PHONY : all build clean info prebuild download execute


all : build
	@echo "### Successfully generated all binary outputs."


build : prebuild $(PRODUCT)


prebuild : 


download : build


execute : download


###########################################################################
## FINAL TARGET
###########################################################################

#---------------------------------
# Create a static library         
#---------------------------------

$(PRODUCT) : $(OBJS) $(PREBUILT_OBJS)
	@echo "### Creating static library "$(PRODUCT)" ..."
	$(AR) $(ARFLAGS)  $(PRODUCT) @$(CMD_FILE)
	@echo "### Created: $(PRODUCT)"


###########################################################################
## INTERMEDIATE TARGETS
###########################################################################

#---------------------
# SOURCE-TO-OBJECT
#---------------------

%.obj : %.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.obj : %.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.obj : $(RELATIVE_PATH_TO_ANCHOR)/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.obj : $(RELATIVE_PATH_TO_ANCHOR)/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.obj : C:/Users/Елизавета/UNI/WORK/p3.5p/P3.5P/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.obj : C:/Users/Елизавета/UNI/WORK/p3.5p/P3.5P/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.obj : $(START_DIR)/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.obj : $(START_DIR)/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.obj : $(MATLAB_ROOT)/rtw/c/src/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.obj : $(MATLAB_ROOT)/rtw/c/src/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


###########################################################################
## DEPENDENCIES
###########################################################################

$(ALL_OBJS) : rtw_proj.tmw $(COMPILER_COMMAND_FILE) $(MAKEFILE)


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
	$(RM) $(subst /,\,$(PRODUCT))
	$(RM) $(subst /,\,$(ALL_OBJS))
	$(ECHO) "### Deleted all derived files."


