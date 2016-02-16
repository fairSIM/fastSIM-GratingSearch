
.PHONY: all clean

# Options for the java compiler
JFLAGS = -g -d ./ -Xlint:unchecked -extdirs "./external"
JC = javac6
JAR = jar

# remove command to clean up
RM = rm -rf

# java -> class rules
.SUFFIXES: .java .class
.java.class:
	$(JC) $(JFLAGS) $*.java


# just build the project
all:	$(wildcard *.java)
	$(JC) $(JFLAGS) $(wildcard *.java)


# create a jar-file
jar: all
	$(JAR) -cvf SLM_GratingSearch_$(shell date +%Y%m%d-%H%M).jar \
	de/*/*/*.class


clean :
	$(RM) ./de SLM_GratingSearch_*.jar
