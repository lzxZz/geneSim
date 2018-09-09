TARGET=./debug/GeneSimilarity

OBJECTS= ./obj/data.o \
		./obj/evaluator.o \
		./obj/calculator.o \
		./obj/getter.o \
		./obj/main.o \
		./obj/log.o

#./obj/shared.o \

#vpath %.o ./obj
vpath %.cpp ./src
vpath %.h ./include
FLAGS= -Wall -I./include -g

$(TARGET) : $(OBJECTS)
	g++ $(FLAGS) -o $(TARGET) $(OBJECTS)

./obj/getter.o : getter.cpp getter.h
	gcc -c $(FLAGS) ./src/getter.cpp -o ./obj/getter.o

./obj/log.o : log.cpp log.h
	gcc -c $(FLAGS) ./src/log.cpp -o ./obj/log.o

./obj/data.o : data.cpp data.h term.h edge.h anno.h log.h
	gcc -c $(FLAGS) ./src/data.cpp -o ./obj/data.o

# ./obj/shared.o : shared.cpp shared.h
# 	gcc -c $(FLAGS) ./src/shared.cpp -o ./obj/shared.o

./obj/evaluator.o : evaluator.cpp calc.h
	gcc -c $(FLAGS) ./src/evaluator.cpp -o ./obj/evaluator.o

./obj/calculator.o : calculator.cpp calc.h
	gcc -c $(FLAGS) ./src/calculator.cpp -o ./obj/calculator.o

./obj/main.o : main.cpp term.h
	gcc -c $(FLAGS) ./src/main.cpp -o ./obj/main.o


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
	@-rm $(OBJECTS)
	@echo "finished clean"

mem:
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all $(TARGET)

