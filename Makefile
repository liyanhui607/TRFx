# 编译器设置
CC      = gcc
NVCC    = nvcc
CFLAGS  = -g -Wall -O2 -Wc++-compat -Wno-unused-function
#CFLAGS  =  -g -O0 -fno-omit-frame-pointer 
LDFLAGS = -lstdc++ -lm -lz -lpthread

# 增强版GPU检测（同时检测驱动和工具链）
CUDA_PATH          ?= /usr/local/cuda
GPU_DEVICE_PRESENT := $(shell which nvidia-smi >/dev/null 2>&1 && echo 1 || echo 0)
CUDA_INSTALLED     := $(shell [ -d "$(CUDA_PATH)/bin" ] && echo 1 || echo 0)
HAVE_GPU           := $(shell if [ $(GPU_DEVICE_PRESENT) -eq 1 ] && [ $(CUDA_INSTALLED) -eq 1 ]; then echo 1; else echo 0; fi)

# 条件编译设置
ifeq ($(HAVE_GPU),1)
    OBJS        = kthread.o bseq.o map.o trfx.o GetTopPeriods_cuda.o
    INCLUDES    = -I$(CUDA_PATH)/include
    LDFLAGS    += -L$(CUDA_PATH)/lib64 -lcuda -lcudart
    NVCC_FLAGS  = -arch=sm_75  # 根据实际GPU架构调整
    CFLAGS     += -DHAVE_GPU   # 关键：传递宏定义给源代码
else
    OBJS        = kthread.o bseq.o map.o trfx.o
    INCLUDES    = 
endif

PROG       = trfx
PROG_EXTRA = sdust trfx-lite

.SUFFIXES: .c .o .cu

# 编译规则
.c.o:
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

ifeq ($(HAVE_GPU),1)
.cu.o:
	$(NVCC) -c $(NVCC_FLAGS) $(INCLUDES) $< -o $@
endif

# 主目标
all: $(PROG)

# 构建主程序
trfx: trfx.o libtrfx.a
	$(CC) $(CFLAGS) $< -o $@ -L. -ltrfx $(LDFLAGS)

# 构建静态库
libtrfx.a: $(OBJS)
	$(AR) -csrD $@ $(OBJS)

# 依赖关系
trfx.o: trfx.c trfx.h bseq.h
bseq.o: bseq.h kseq.h
map.o: bseq.h trfx.h trfrun.h tr30dat.h
GetTopPeriods_cuda.o: GetTopPeriods_cuda.cu
ifeq ($(HAVE_GPU),1)
tr30dat.o: tr30dat.c tr30dat.h GetTopPeriods_cuda.o
else
tr30dat.o: tr30dat.c tr30dat.h
endif

# 辅助目标
trfx-lite: example.o libtrfx.a
	$(CC) $(CFLAGS) $< -o $@ -L. -ltrfx $(LDFLAGS)

clean:
	rm -fr *.o *.a $(PROG) $(PROG_EXTRA) *~ *.dSYM session* output

.PHONY: all clean