## How to compile and generate synthetic data



## Linux or any other environment that simulates the linux terminal experience (such as WSL)

1) Change the number "64" on line 49 in main.cpp to the number of threads that you will use
2) Make sure that you have a g++ compiler that supports c++23 or later (GCC 13.2.0 was used for this project)
3) Navigate to the "Simulator" folder and simply type "make", it will create a binary called "main"
4) Run the synthetic data generation by executing "./main"



## Windows or any other environment

Not directly supported. Ask your favorite LLM for help on how to translate the contents of "Makefile" into a successful compilation on your system.