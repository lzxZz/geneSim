TARGET=./debug/GeneSimilarity

OBJECTS= ./obj/main.o \
		./obj/data.o \
		./obj/shared.o \
		./obj/evaluator.o \
		./obj/calculator.o

#vpath %.o ./obj
vpath %.cpp ./src
vpath %.h ./include
FLAGS= -Wall -I./include -g

$(TARGET) : $(OBJECTS)
	g++ $(FLAGS) -o $(TARGET) $(OBJECTS)


./obj/main.o : main.cpp term.h
	gcc -c $(FLAGS) ./src/main.cpp -o ./obj/main.o

./obj/data.o : data.cpp data.h term.h edge.h anno.h
	gcc -c $(FLAGS) ./src/data.cpp -o ./obj/data.o

./obj/shared.o : shared.cpp 
	gcc -c $(FLAGS) ./src/shared.cpp -o ./obj/shared.o

./obj/evaluator.o : evaluator.cpp calc.h
	gcc -c $(FLAGS) ./src/evaluator.cpp -o ./obj/evaluator.o

./obj/calculator.o : calculator.cpp calc.h
	gcc -c $(FLAGS) ./src/calculator.cpp -o ./obj/calculator.o




.PHONY:run
.PHONY:clean
.PHONY:line

line:
	@wc -l ./src/*.cpp
	@wc -l ./include/*.h

run:
	@./debug/GeneSimilarity

clean:
	@-rm $(TARGET)
	@-rm $(OBJECTS)
	@echo "finished clean"