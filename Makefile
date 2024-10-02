# имя программы:
TARGET = BOX_FLAT_LIP

# исходники для сборки:
SOURCES = \
       BOX_FLAT_LIP.f90

OBJECTS=$(SOURCES:%.f90=%.o)

TYPE=master

# про флаги:
# https://www.opennet.ru/docs/RUS/cpp/cpp-10.html

# простая сборка:
all: $(TARGET)

$(OBJECTS): $(SOURCES)

$(TARGET): $(OBJECTS)
	gfortran -O3 -o $(TARGET) -cpp $(SOURCES) 
clean:
	$(RM) $(TARGET) *.mod
	
.PHONY: all clean