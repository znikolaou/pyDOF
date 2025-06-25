
BUILD_PATH=$(PWD)
SRC=$(PYDOF_CSRC)/pydof.cpp $(PYDOF_CSRC)/pydofFilter.cpp 
SRC_TEST1=test1.cpp 
SRC_TEST2=test2.cpp 


OBJS_TEST1=$(SRC_TEST1:%.cpp=%.o)
OBJS_TEST2=$(SRC_TEST2:%.cpp=%.o)

OBJS1=$(patsubst $(PYDOF_CSRC)/%,$(BUILD_PATH)/%,$(SRC:%.cpp=%.o)) $(OBJS_TEST1)
OBJS2=$(patsubst $(PYDOF_CSRC)/%,$(BUILD_PATH)/%,$(SRC:%.cpp=%.o)) $(OBJS_TEST2)

goTest1: $(OBJS1)
	$(CXX) $(CFLAGS) $(OBJS1) -o $@

goTest2: $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) -o $@


%.o:%.cpp 
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<  

$(BUILD_PATH)/%.o:$(PYDOF_CSRC)/%.cpp $(PYDOF_CSRC)/%.h
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< 

clean:
	@rm -f *.o  a.out* go*

path:
	@echo $(BUILD_PATH)

printObjs:
	@echo $(OBJS)


