CFLAGS = -O3
CC = g++

p2v: p2v_utils.o p2v_input_reader.o p2v.o
	$(CC) $(CFLAGS) -o p2v p2v_utils.o p2v_input_reader.o p2v.o

p2v_input_reader.o: p2v_input_reader.cpp
	$(CC) $(CFLAGS) -c p2v_input_reader.cpp

p2v_utils.o: p2v_utils.cpp
	$(CC) $(CFLAGS) -c p2v_utils.cpp

p2v.o: p2v.cpp
	$(CC) $(CFLAGS) -c p2v.cpp

clean:
	rm -f core *.o 
