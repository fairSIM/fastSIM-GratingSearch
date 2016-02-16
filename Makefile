
.PHONY: all clean

# Options for the java compiler
JFLAGS = -g -Xlint:unchecked -extdirs "./external"
JC = javac6
JAR = jar

# remove command to clean up
RM = rm -rf

# java -> class rules
.SUFFIXES: .java .class
.java.class:
	$(JC) $(JFLAGS) $*.java


# just build the project
all:	$(wildcard de/bio_photonics/*/*.java)
	$(JC) $(JFLAGS) $(wildcard de/bio_photonics/*/*.java)


# create a jar-file
jar: all
	$(JAR) -cvf SLM_GratingSearch_$(shell date +%Y%m%d-%H%M).jar \
	de/*/*/*.class


clean :
	$(RM) ./de/*/*/*.class SLM_GratingSearch_*.jar
