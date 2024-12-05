CC ?= gcc
APP = cbl
CLIBS = -lgsl -lm
CFLAGS = -Wall -I/usr/include/gsl
OBJS = cbl.o main.o

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

$(APP): $(OBJS)
	$(CC) -o $@ $^ $(CLIBS)
	./$(APP)

clean:
	rm -f $(APP) *.o *~
