CC=g++
CFLAGS= -g -c -Wall -I. -I/usr/X11/include
LDFLAGS=
SOURCES=elastic_membrane.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=elastic_membrane
INCLUDE=
LIBS = -L /opt/X11/lib/ -lm
#LIBS= -L /opt/X11/lib/ -lglut -lGL -lGLU -lX11 -lXmu -lXi -lm
all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $(INCLUDE) $@

clean:
	@rm $(EXECUTABLE) $(OBJECTS)
cleanall: 
	@rm -f $(EXECUTABLE) $(OBJECTS) *.txt
