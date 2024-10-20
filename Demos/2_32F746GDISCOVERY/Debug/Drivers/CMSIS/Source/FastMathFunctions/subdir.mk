################################################################################
# Automatically-generated file. Do not edit!
# Toolchain: GNU Tools for STM32 (11.3.rel1)
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Drivers/CMSIS/Source/FastMathFunctions/arm_cos_f32.c \
../Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q15.c \
../Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q31.c \
../Drivers/CMSIS/Source/FastMathFunctions/arm_sin_f32.c \
../Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q15.c \
../Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q31.c \
../Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q15.c \
../Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q31.c 

OBJS += \
./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_f32.o \
./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q15.o \
./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q31.o \
./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_f32.o \
./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q15.o \
./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q31.o \
./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q15.o \
./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q31.o 

C_DEPS += \
./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_f32.d \
./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q15.d \
./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q31.d \
./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_f32.d \
./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q15.d \
./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q31.d \
./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q15.d \
./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q31.d 


# Each subdirectory must supply rules for building sources it contributes
Drivers/CMSIS/Source/FastMathFunctions/%.o Drivers/CMSIS/Source/FastMathFunctions/%.su Drivers/CMSIS/Source/FastMathFunctions/%.cyclo: ../Drivers/CMSIS/Source/FastMathFunctions/%.c Drivers/CMSIS/Source/FastMathFunctions/subdir.mk
	arm-none-eabi-gcc "$<" -mcpu=cortex-m7 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32F746xx -DARM_MATH_CM7 -c -I../Core/Inc -I"C:/Users/ozand/STM32CubeIDE/dip_workspace/2_32F746GDISCOVERY/Drivers/CMSIS/Include" -I../Drivers/STM32F7xx_HAL_Driver/Inc -I../Drivers/STM32F7xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32F7xx/Include -I../Drivers/CMSIS/Include -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" --specs=nano.specs -mfpu=fpv5-sp-d16 -mfloat-abi=hard -mthumb -o "$@"

clean: clean-Drivers-2f-CMSIS-2f-Source-2f-FastMathFunctions

clean-Drivers-2f-CMSIS-2f-Source-2f-FastMathFunctions:
	-$(RM) ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_f32.cyclo ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_f32.d ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_f32.o ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_f32.su ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q15.cyclo ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q15.d ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q15.o ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q15.su ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q31.cyclo ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q31.d ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q31.o ./Drivers/CMSIS/Source/FastMathFunctions/arm_cos_q31.su ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_f32.cyclo ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_f32.d ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_f32.o ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_f32.su ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q15.cyclo ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q15.d ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q15.o ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q15.su ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q31.cyclo ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q31.d ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q31.o ./Drivers/CMSIS/Source/FastMathFunctions/arm_sin_q31.su ./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q15.cyclo ./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q15.d ./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q15.o ./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q15.su ./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q31.cyclo ./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q31.d ./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q31.o ./Drivers/CMSIS/Source/FastMathFunctions/arm_sqrt_q31.su

.PHONY: clean-Drivers-2f-CMSIS-2f-Source-2f-FastMathFunctions

