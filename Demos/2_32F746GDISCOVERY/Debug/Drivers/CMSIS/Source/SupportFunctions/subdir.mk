################################################################################
# Automatically-generated file. Do not edit!
# Toolchain: GNU Tools for STM32 (11.3.rel1)
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Drivers/CMSIS/Source/SupportFunctions/arm_copy_f32.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_copy_q15.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_copy_q31.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_copy_q7.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_fill_f32.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_fill_q15.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_fill_q31.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_fill_q7.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q15.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q31.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q7.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_float.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q31.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q7.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_float.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q15.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q7.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_float.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q15.c \
../Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q31.c 

OBJS += \
./Drivers/CMSIS/Source/SupportFunctions/arm_copy_f32.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q15.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q31.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q7.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_fill_f32.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q15.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q31.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q7.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q15.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q31.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q7.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_float.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q31.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q7.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_float.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q15.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q7.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_float.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q15.o \
./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q31.o 

C_DEPS += \
./Drivers/CMSIS/Source/SupportFunctions/arm_copy_f32.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q15.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q31.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q7.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_fill_f32.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q15.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q31.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q7.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q15.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q31.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q7.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_float.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q31.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q7.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_float.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q15.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q7.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_float.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q15.d \
./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q31.d 


# Each subdirectory must supply rules for building sources it contributes
Drivers/CMSIS/Source/SupportFunctions/%.o Drivers/CMSIS/Source/SupportFunctions/%.su Drivers/CMSIS/Source/SupportFunctions/%.cyclo: ../Drivers/CMSIS/Source/SupportFunctions/%.c Drivers/CMSIS/Source/SupportFunctions/subdir.mk
	arm-none-eabi-gcc "$<" -mcpu=cortex-m7 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32F746xx -DARM_MATH_CM7 -c -I../Core/Inc -I"C:/Users/ozand/STM32CubeIDE/dip_workspace/2_32F746GDISCOVERY/Drivers/CMSIS/Include" -I../Drivers/STM32F7xx_HAL_Driver/Inc -I../Drivers/STM32F7xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32F7xx/Include -I../Drivers/CMSIS/Include -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" --specs=nano.specs -mfpu=fpv5-sp-d16 -mfloat-abi=hard -mthumb -o "$@"

clean: clean-Drivers-2f-CMSIS-2f-Source-2f-SupportFunctions

clean-Drivers-2f-CMSIS-2f-Source-2f-SupportFunctions:
	-$(RM) ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_f32.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_f32.d ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_f32.o ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_f32.su ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q15.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q15.d ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q15.o ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q15.su ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q31.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q31.d ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q31.o ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q31.su ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q7.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q7.d ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q7.o ./Drivers/CMSIS/Source/SupportFunctions/arm_copy_q7.su ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_f32.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_f32.d ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_f32.o ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_f32.su ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q15.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q15.d ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q15.o ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q15.su ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q31.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q31.d ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q31.o ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q31.su ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q7.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q7.d ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q7.o ./Drivers/CMSIS/Source/SupportFunctions/arm_fill_q7.su ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q15.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q15.d ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q15.o ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q15.su ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q31.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q31.d ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q31.o ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q31.su ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q7.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q7.d ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q7.o ./Drivers/CMSIS/Source/SupportFunctions/arm_float_to_q7.su ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_float.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_float.d ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_float.o ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_float.su ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q31.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q31.d ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q31.o ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q31.su ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q7.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q7.d ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q7.o ./Drivers/CMSIS/Source/SupportFunctions/arm_q15_to_q7.su ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_float.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_float.d ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_float.o ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_float.su ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q15.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q15.d ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q15.o ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q15.su ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q7.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q7.d ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q7.o ./Drivers/CMSIS/Source/SupportFunctions/arm_q31_to_q7.su ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_float.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_float.d ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_float.o ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_float.su ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q15.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q15.d ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q15.o ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q15.su ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q31.cyclo ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q31.d ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q31.o ./Drivers/CMSIS/Source/SupportFunctions/arm_q7_to_q31.su

.PHONY: clean-Drivers-2f-CMSIS-2f-Source-2f-SupportFunctions

