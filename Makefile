TARGET=./debug/GeneSimilarity

OBJS = ./obj/main.o \
		./obj/lfc.o \
		./obj/gene.o \
		./obj/term.o

FLAGS = -Wall -pg -g

vpath %cpp ./src/
vpath %.h  ./include/


$(TARGET) : $(OBJS)
	g++ $(FLAGS) -o $(TARGET) $(OBJS)

./obj/main.o : main.cpp sim.h
	gcc -c $(FLAGS) ./src/main.cpp -o ./obj/main.o

./obj/lfc.o : lfc.cpp sim.h
	gcc -c $(FLAGS) ./src/lfc.cpp -o ./obj/lfc.o

./obj/gene.o : gene.cpp sim.h
	gcc -c $(FLAGS) ./src/gene.cpp -o ./obj/gene.o

./obj/term.o : term.cpp sim.h
	gcc -c $(FLAGS) ./src/term.cpp -o ./obj/term.o

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
