TARGET=./debug/GeneSimilarity

OBJECTS= ./obj/main.o ./obj/data.o
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





.PHONY:run
.PHONY:clean




run:
	@./debug/GeneSimilarity

clean:
	@rm $(TARGET)
	@rm $(OBJECTS)
	@echo "finished clean"