CFLAGS = -O3 -c
LDFLAGS= 
CC = g++

#SRC = p2v.cpp p2v_utils.cpp p2v_input_reader.cpp p2v_parameters.cpp p2v_chromosome.cpp
#OBJ = $(SRC:.cpp = .o)

#all: $(SRC) p2v

#p2v: $(OBJ)
#	$(CC) $(LDFLAGS) -o p2v $(OBJ)

#.cpp.o:
#	$(CC) $(CFLAGS) $< -o $@

p2v: p2v_utils.o p2v_input_reader.o p2v.o p2v_parameters.o p2v_chromosome.o p2v_genome.o p2v_svdata.o
	$(CC) $(LDFLAGS) -o p2vbin p2v_utils.o p2v_input_reader.o p2v.o p2v_parameters.o p2v_chromosome.o p2v_genome.o p2v_svdata.o

p2v_input_reader.o: p2v_input_reader.cpp
	$(CC) $(CFLAGS) -c p2v_input_reader.cpp

p2v_utils.o: p2v_utils.cpp
	$(CC) $(CFLAGS) -c p2v_utils.cpp

p2v_parameters.o: p2v_parameters.cpp
	$(CC) $(CFLAGS) -c p2v_parameters.cpp

p2v_chromosome.o: p2v_chromosome.cpp
	$(CC) $(CFLAGS) -c p2v_chromosome.cpp

p2v_genome.o: p2v_genome.cpp
	$(CC) $(CFLAGS) -c p2v_genome.cpp

p2v_svdata.o: p2v_svdata.cpp
	$(CC) $(CFLAGS) -c p2v_svdata.cpp

#p2v_genotype.o: p2v_genotype.cpp
#	$(CC) $(CFLAGS) -c p2v_genotype.cpp

p2v.o: p2v.cpp
	$(CC) $(CFLAGS) -c p2v.cpp

clean:
	rm -f core *.o 

