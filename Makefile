
.PHONY: all clean

# Options for the java compiler
JFLAGS = -g -Xlint:unchecked -extdirs "./external" -d ./
JFLAGS+= -target 1.6 -source 1.6 -bootclasspath ./external/rt-1.6.jar
JC = javac
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
jar:	all git-version
	$(JAR) -cvfm SLM_GratingSearch_$(shell head -c 10 git-version.txt).jar \
	Manifest.txt plugins.config de/*/*/*.class

git-version :
	git rev-parse HEAD > ./git-version.txt  ; \
	git tag --contains >> ./git-version.txt ; \
	echo "n/a" >> ./git-version.txt
	 	

clean :
	$(RM) ./de/*/*/*.class SLM_GratingSearch_*.jar
