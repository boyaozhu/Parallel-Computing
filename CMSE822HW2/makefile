
CC := gcc
TAUCC := tau_cc.sh

#CFLAGS := -std=gnu99 -Wall -Wextra -Wpedantic -Wunused-parameter -O3
CFLAGS := -std=gnu99 -O3 #-qopenmp
SRCDIR := src
CFILES := $(filter-out ./$(SRCDIR)/main.c,$(wildcard ./$(SRCDIR)/*.c))
OFILES := $(CFILES:.c=.o)
TAUOFILES := $(CFILES:.c=.tau.o)

LIB := -lm


multiplication: $(OFILES)
	$(CC) -o multiplication.x ./$(SRCDIR)/main.c $(OFILES) $(CFLAGS) $(LIB)

$(SRCDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c -o $@ $< $(CFLAGS) $(LIB)


multiplication-tau: $(TAUOFILES)
	$(TAUCC) -o multiplication-tau.x ./$(SRCDIR)/main.c $(TAUOFILES) $(CFLAGS) $(LIB)

$(SRCDIR)/%.tau.o: $(SRCDIR)/%.c
	$(TAUCC) -c -o $@ $< $(CFLAGS) $(LIB)


clean:
	-rm  *.x $(OFILES) $(TAUOFILES) main.o
