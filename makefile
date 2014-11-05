CC=gcc
CFLAGS= -Wall

temp_dist2D: temp_dist2D.c
	$(CC) $(CFLAGS) -o temp_dist2D temp_dist2D.c
