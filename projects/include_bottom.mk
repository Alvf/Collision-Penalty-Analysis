# object files
STRIPPED = $(subst ../../, , $(SOURCES))
OBJS = $(patsubst %.cpp, ../build/%.o , $(STRIPPED))
MAIN_OBJ = $(patsubst %.cpp, ./build/%.o , $(MAIN))
EXECUTABLE = $(patsubst %.cpp, ../../bin/%, $(MAIN))

DIRS = ext \
			 ext/glvu \
			 ext/catch2 \
			 ext/solvePoly \
			 src \
			 src/Collision \
			 src/util \
			 src/Timestepper \
			 src/Timestepper/Volume \
			 src/Timestepper/Shell \
			 src/Timestepper/Strand \
			 src/Geometry \
			 src/Damping \
			 src/Damping/Volume \
			 src/Hyperelastic/ \
			 src/Hyperelastic/Strand \
			 src/Hyperelastic/Volume \
			 src/Hyperelastic/Shell

# how to make the main target
$(EXECUTABLE): $(OBJS) $(MAIN_OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: release
release: $(EXECUTABLE)

# how to compile each file
.SUFFIXES:
../build/%.o:
	$(CC) $(CFLAGS) -o $@ $<

.SUFFIXES:
./build/%.o:
	$(CC) $(CFLAGS) -o $@ $<
$(MAIN_OBJ):
	$(CC) $(CFLAGS) -o $@ $(MAIN)

# cleaning up
.PHONY: clean
clean:
	-find ./build/ -name '*.o' -delete
	-find ../build/ -name '*.o' -delete
	-rm -f $(EXECUTABLE)

# dependencies are automatically generated
.PHONY: depend
depend:
	-mkdir ../build
	-mkdir build
	-rm -f build/depend
	-mkdir $(foreach dir,$(DIRS), ../build/$(dir))
	$(foreach srcfile,$(SOURCES),$(CC) $(INCLUDES) -MM $(srcfile) -MT $(patsubst %.cpp, ../build/%.o, $(subst ../../,,$(srcfile))) >> build/depend;)
	$(CC) $(INCLUDES) -MM $(MAIN) -MT $(MAIN_OBJ) >> build/depend;

-include build/depend
