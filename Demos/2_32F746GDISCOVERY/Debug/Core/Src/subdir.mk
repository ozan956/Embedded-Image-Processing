################################################################################
# Automatically-generated file. Do not edit!
# Toolchain: GNU Tools for STM32 (11.3.rel1)
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Core/Src/DIP.c \
../Core/Src/DIP_LTDC.c \
../Core/Src/DIP_SERIAL.c \
../Core/Src/DIP_UTIL.c \
../Core/Src/main.c \
../Core/Src/mandrill.c \
../Core/Src/stm32f7xx_hal_msp.c \
../Core/Src/stm32f7xx_it.c \
../Core/Src/syscalls.c \
../Core/Src/sysmem.c \
../Core/Src/system_stm32f7xx.c \
../Core/Src/t1_480x270.c \
../Core/Src/t2_480x270.c \
../Core/Src/t3_480x270.c 

OBJS += \
./Core/Src/DIP.o \
./Core/Src/DIP_LTDC.o \
./Core/Src/DIP_SERIAL.o \
./Core/Src/DIP_UTIL.o \
./Core/Src/main.o \
./Core/Src/mandrill.o \
./Core/Src/stm32f7xx_hal_msp.o \
./Core/Src/stm32f7xx_it.o \
./Core/Src/syscalls.o \
./Core/Src/sysmem.o \
./Core/Src/system_stm32f7xx.o \
./Core/Src/t1_480x270.o \
./Core/Src/t2_480x270.o \
./Core/Src/t3_480x270.o 

C_DEPS += \
./Core/Src/DIP.d \
./Core/Src/DIP_LTDC.d \
./Core/Src/DIP_SERIAL.d \
./Core/Src/DIP_UTIL.d \
./Core/Src/main.d \
./Core/Src/mandrill.d \
./Core/Src/stm32f7xx_hal_msp.d \
./Core/Src/stm32f7xx_it.d \
./Core/Src/syscalls.d \
./Core/Src/sysmem.d \
./Core/Src/system_stm32f7xx.d \
./Core/Src/t1_480x270.d \
./Core/Src/t2_480x270.d \
./Core/Src/t3_480x270.d 


# Each subdirectory must supply rules for building sources it contributes
Core/Src/%.o Core/Src/%.su Core/Src/%.cyclo: ../Core/Src/%.c Core/Src/subdir.mk
	arm-none-eabi-gcc "$<" -mcpu=cortex-m7 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32F746xx -DARM_MATH_CM7 -c -I../Core/Inc -I"C:/Users/ozand/STM32CubeIDE/dip_workspace/2_32F746GDISCOVERY/Drivers/CMSIS/Include" -I../Drivers/STM32F7xx_HAL_Driver/Inc -I../Drivers/STM32F7xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32F7xx/Include -I../Drivers/CMSIS/Include -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" --specs=nano.specs -mfpu=fpv5-sp-d16 -mfloat-abi=hard -mthumb -o "$@"

clean: clean-Core-2f-Src

clean-Core-2f-Src:
	-$(RM) ./Core/Src/DIP.cyclo ./Core/Src/DIP.d ./Core/Src/DIP.o ./Core/Src/DIP.su ./Core/Src/DIP_LTDC.cyclo ./Core/Src/DIP_LTDC.d ./Core/Src/DIP_LTDC.o ./Core/Src/DIP_LTDC.su ./Core/Src/DIP_SERIAL.cyclo ./Core/Src/DIP_SERIAL.d ./Core/Src/DIP_SERIAL.o ./Core/Src/DIP_SERIAL.su ./Core/Src/DIP_UTIL.cyclo ./Core/Src/DIP_UTIL.d ./Core/Src/DIP_UTIL.o ./Core/Src/DIP_UTIL.su ./Core/Src/main.cyclo ./Core/Src/main.d ./Core/Src/main.o ./Core/Src/main.su ./Core/Src/mandrill.cyclo ./Core/Src/mandrill.d ./Core/Src/mandrill.o ./Core/Src/mandrill.su ./Core/Src/stm32f7xx_hal_msp.cyclo ./Core/Src/stm32f7xx_hal_msp.d ./Core/Src/stm32f7xx_hal_msp.o ./Core/Src/stm32f7xx_hal_msp.su ./Core/Src/stm32f7xx_it.cyclo ./Core/Src/stm32f7xx_it.d ./Core/Src/stm32f7xx_it.o ./Core/Src/stm32f7xx_it.su ./Core/Src/syscalls.cyclo ./Core/Src/syscalls.d ./Core/Src/syscalls.o ./Core/Src/syscalls.su ./Core/Src/sysmem.cyclo ./Core/Src/sysmem.d ./Core/Src/sysmem.o ./Core/Src/sysmem.su ./Core/Src/system_stm32f7xx.cyclo ./Core/Src/system_stm32f7xx.d ./Core/Src/system_stm32f7xx.o ./Core/Src/system_stm32f7xx.su ./Core/Src/t1_480x270.cyclo ./Core/Src/t1_480x270.d ./Core/Src/t1_480x270.o ./Core/Src/t1_480x270.su ./Core/Src/t2_480x270.cyclo ./Core/Src/t2_480x270.d ./Core/Src/t2_480x270.o ./Core/Src/t2_480x270.su ./Core/Src/t3_480x270.cyclo ./Core/Src/t3_480x270.d ./Core/Src/t3_480x270.o ./Core/Src/t3_480x270.su

.PHONY: clean-Core-2f-Src

