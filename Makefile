SRC_DIR := ./src
OBJ_DIR := ./src
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
ROOT_DIR := .
CFLAGS = -I $(ROOT_DIR)/cudd/cudd -I $(ROOT_DIR)/cudd/util
LFLAGS = $(ROOT_DIR)/cudd/cudd/.libs/libcudd.a $(ROOT_DIR)/abc/libabc.a -lm -ldl -lreadline -ltinfo -lpthread -lboost_program_options

VanQiRA: $(OBJ_FILES)
	g++ -o $@ $^ $(LFLAGS) 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	g++ -c -o $@ $< $(CFLAGS)

clean:
	rm $(OBJ_FILES)
