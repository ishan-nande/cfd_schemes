TARGET = example_4_1_HVWM

BINDIR = bin
CFLAGS = gcc -std=c11 -g -Wall 
CC = cc
LDFLAGS = -l lis -l m

$(TARGET): $(TARGET)
	$(CFLAGS) $(TARGET).c -o ./$(BINDIR)/$(TARGET) $(LDFLAGS)

#remove executbale file 
rm_exe:
	rm $(TARGET)
clean:	
	rm *.o
