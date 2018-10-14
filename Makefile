TARGET=./debug/GeneSimilarity

OBJS = ./obj/main.o \
		./obj/lfc.o \
		./obj/lfc_generate.o \
		./obj/gene.o \
		./obj/term.o \
		./obj/term_index.o

FLAGS = -Wall -pg -g -std=c++11

vpath %cpp ./src/
vpath %.h  ./include/


$(TARGET) : $(OBJS)
	g++ $(FLAGS) -o $(TARGET) $(OBJS)

./obj/main.o : main.cpp anno.h  defs.h  sim_gene.h  sim_lfc.h  sim_term.h  term.h
	gcc -c $(FLAGS) ./src/main.cpp -o ./obj/main.o

./obj/lfc.o : lfc.cpp sim_lfc.h
	gcc -c $(FLAGS) ./src/lfc.cpp -o ./obj/lfc.o

./obj/gene.o : gene.cpp sim_gene.h
	gcc -c $(FLAGS) ./src/gene.cpp -o ./obj/gene.o

./obj/term.o : term.cpp sim_term.h
	gcc -c $(FLAGS) ./src/term.cpp -o ./obj/term.o

./obj/term_index.o : term_index.cpp sim_term.h
	gcc -c $(FLAGS) ./src/term_index.cpp -o ./obj/term_index.o

./obj/lfc_generate.o : lfc_gene_generate.cpp sim_lfc.h
	gcc -c $(FLAGS) ./src/lfc_gene_generate.cpp -o ./obj/lfc_generate.o

.PHONY:run
.PHONY:clean
.PHONY:line
.PHONY:mem
line:
	@wc -l ./src/*.cpp
	@wc -l ./include/*.h

run:
	@./debug/GeneSimilarity

clean:
	@-rm $(TARGET)
	@-rm $(OBJS)
	@-rm ./debug/.fuse*
	@echo "finished clean"

mem:
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all $(TARGET)
