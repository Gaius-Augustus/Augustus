################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Compute_UTRs.cpp \
../src/Coord_Transform.cpp \
../src/Genomic_Data.cpp \
../src/Splice_Sites.cpp \
../src/Supporting_Methods.cpp \
../src/Test.cpp \
../src/UTRs.cpp \
../src/main.cpp 

OBJS += \
./src/Compute_UTRs.o \
./src/Coord_Transform.o \
./src/Genomic_Data.o \
./src/Splice_Sites.o \
./src/Supporting_Methods.o \
./src/Test.o \
./src/UTRs.o \
./src/main.o 

CPP_DEPS += \
./src/Compute_UTRs.d \
./src/Coord_Transform.d \
./src/Genomic_Data.d \
./src/Splice_Sites.d \
./src/Supporting_Methods.d \
./src/Test.d \
./src/UTRs.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/boost -O0 -g3 -pedantic -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


