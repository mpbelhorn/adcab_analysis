EXEC = Fit3d
SRC = Fit3D.cpp main.cpp
OBJECTS = $(patsubst %.cpp,%.o,$(SRC))

SHARE_SRC = DileptonEvents.cpp EventSelectors.cpp
SHARE_DIR = /home/matt/research/belle/adcab/analysis/DileptonEvents/includes
SHARE_SRCS = $(patsubst %.cpp,$(SHARE_DIR)/%.cpp,$(SHARE_SRC))
SHARE_OBJECTS = $(patsubst %.cpp,%.o,$(SHARE_SRCS))

CC = g++

CFLAGS = -pthread -m64 -I/usr/include/root
LDFLAGS = -L/usr/lib64/root -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -pthread -lm -ldl -rdynamic

$(EXEC): $(OBJECTS) $(SHARE_OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -I$(SHARE_DIR) -c $< -o $@

clean: 
	rm -rf *.o a.out 