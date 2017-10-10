
SRCDIR = src
BUILDDIR = obj

SRCNAMES := $(shell find $(SOURCEDIR) -name '*.cpp' -type f -exec basename {} \;)

HNAMES := $(shell find $(SOURCEDIR) -name '*.h' -type f -exec basename {} \;)

OBJECTS := $(addprefix $(BUILDDIR)/, $(SRCNAMES:%.cpp=%.o))
SRCS := $(addprefix $(SRCDIR)/, $(SRCNAMES))

LIKDIR1=/usr/local/likwid/
LIKDIR2=/home/soft/likwid/
LIBS =  -DLIKWID_PERFMON -I$(LIKDIR2)include -L$(LIKDIR2)lib -llikwid -lm

# warnings and flags
WARN = -Wall
WNO = -Wno-comment  -Wno-sign-compare
FLAGS = -O3 $(WARN) $(WNO)

# Executable filename
bin = invmat

#compiler
compiler = g++ -std=c++11

all: pre obj_dir list_srcnames $(bin)

pre:


doc:
	doxygen doxyconfig

set_debug:
	$(eval FLAGS = -O0 -g $(WARN))

debug: set_debug all

rebuild: clean all

buildclean: all clean

obj_dir:
	mkdir -p $(BUILDDIR)

list_srcnames:
	@echo
	@echo "Found source files to compile:"
	@echo $(SRCNAMES)
	@echo "Found header files:"
	@echo $(HNAMES)
	@echo

# Tool invocations
$(bin): $(OBJECTS)
	@echo 'Building target: $@'
	@echo 'Invoking Linker'
	$(compiler) -o "$(bin)" $(OBJECTS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking Compiler'
	$(compiler) $(FLAGS) $(LIBS) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo

clean:
	rm -rf  ./$(BUILDDIR)/*.o  ./$(BUILDDIR)/*.d $(bin)

cleanAll: clean
	rm -rf  $(bin)