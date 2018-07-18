CC = cc
TARGET = tool
OBJECTS = main.o
#IDIR = -I/usr/local/include
#LDIR = -L/usr/local/lib
CFLAGS = -std=c99 $(IDIR)
LINK = -lpng -lm

all : $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDIR) $(LINK)

clean :
	rm $(TARGET) $(OBJECTS)
