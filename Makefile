# TRFx Smart Cross-Platform Makefile
# Automatic GPU Detection with Fallback to CPU Version
# Using GCC only (no G++)

# Basic compiler settings - Use GCC for everything
CC      = gcc
CFLAGS  = -g -Wall -O2 -Wc++-compat -Wno-unused-function
LDFLAGS = -lm -lz -lpthread  # Linux不需要-lstdc++，macOS会在检测后添加

# Detect operating system
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Darwin)
    # macOS settings
    CFLAGS  += -mmacosx-version-min=10.9
    CUDA_PATH ?= /usr/local/cuda
    LDFLAGS += -lstdc++  # macOS需要显式链接C++标准库
else
    # Linux settings
    CUDA_PATH ?= /usr/local/cuda
endif

# ============================================================================
# Enhanced Universal GPU Detection Solution
# ============================================================================

# 1. Smart CUDA path detection with multiple fallbacks
ifdef CUDA_PATH
    # Use user-defined CUDA path first
    CUDA_PATH := $(CUDA_PATH)
else
    # Try multiple detection strategies
    ifneq ($(wildcard /usr/local/cuda/include/cuda.h),)
        CUDA_PATH := /usr/local/cuda
    else ifneq ($(wildcard /usr/lib/cuda/include/cuda.h),)
        CUDA_PATH := /usr/lib/cuda
    else ifneq ($(wildcard /usr/include/cuda.h),)
        CUDA_PATH := /usr
    else ifneq ($(wildcard /opt/cuda/include/cuda.h),)
        CUDA_PATH := /opt/cuda
    else
        # Strategy 2: Infer from nvcc path if available
        NVCC_PATH := $(shell which nvcc 2>/dev/null)
        ifneq ($(NVCC_PATH),)
            NVCC_REALPATH := $(shell readlink -f $(NVCC_PATH) 2>/dev/null || echo $(NVCC_PATH))
            INFERRED_CUDA_PATH := $(shell dirname "$(shell dirname "$(NVCC_REALPATH)")")
            ifneq ($(wildcard $(INFERRED_CUDA_PATH)/include/cuda.h),)
                CUDA_PATH := $(INFERRED_CUDA_PATH)
            else
                CUDA_PATH := /usr/local/cuda
            endif
        else
            CUDA_PATH := /usr/local/cuda
        endif
    endif
endif

# 2. Enhanced component detection
GPU_DRIVER_PRESENT := $(shell which nvidia-smi > /dev/null 2>&1 && nvidia-smi > /dev/null 2>&1 && echo 1 || echo 0)
NVCC_AVAILABLE := $(shell which nvcc > /dev/null 2>&1 && echo 1 || echo 0)

CUDA_HEADERS_PRESENT := $(shell \
    if [ -f "$(CUDA_PATH)/include/cuda.h" ] || \
       [ -f "$(CUDA_PATH)/include/cuda/cuda.h" ] || \
       [ -f "/usr/include/cuda.h" ] || \
       [ -f "/usr/local/cuda/include/cuda.h" ]; then \
        echo 1; \
    else \
        echo 0; \
    fi \
)

CUDA_LIBS_PRESENT := $(shell \
    if [ -f "$(CUDA_PATH)/lib64/libcudart.so" ] || \
       [ -f "$(CUDA_PATH)/lib/libcudart.so" ] || \
       [ -f "$(CUDA_PATH)/lib64/libcudart.dylib" ] || \
       [ -f "$(CUDA_PATH)/lib/libcudart.dylib" ] || \
       [ -f "/usr/lib/x86_64-linux-gnu/libcudart.so" ] || \
       [ -f "/usr/lib64/libcudart.so" ] || \
       [ -f "/usr/lib/libcudart.so" ]; then \
        echo 1; \
    else \
        echo 0; \
    fi \
)

# 3. Comprehensive evaluation - Always build at least CPU version
HAS_CUDA_DEVELOPMENT_ENV := $(shell \
    if [ $(NVCC_AVAILABLE) -eq 1 ] && [ $(CUDA_HEADERS_PRESENT) -eq 1 ] && [ $(CUDA_LIBS_PRESENT) -eq 1 ]; then \
        echo 1; \
    else \
        echo 0; \
    fi \
)

# Final decision: Try GPU if possible, but always fallback to CPU
ifeq ($(GPU_DRIVER_PRESENT)$(HAS_CUDA_DEVELOPMENT_ENV), 11)
    ENABLE_CUDA = 1
    # Smart GPU architecture detection
    ifneq ("$(shell which nvidia-smi 2>/dev/null)","")
        GPU_ARCH_RAW := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1)
        ifneq ($(GPU_ARCH_RAW),)
            GPU_ARCH_MAJOR := $(shell echo "$(GPU_ARCH_RAW)" | cut -d'.' -f1)
            GPU_ARCH_MINOR := $(shell echo "$(GPU_ARCH_RAW)" | cut -d'.' -f2)
            GPU_ARCH := $(GPU_ARCH_MAJOR)$(GPU_ARCH_MINOR)
            
            ifeq ($(shell expr $(GPU_ARCH_MAJOR) \>= 8), 1)
                NVCC_FLAGS = -arch=sm_80
            else ifeq ($(shell expr $(GPU_ARCH_MAJOR) \>= 7), 1)
                ifeq ($(GPU_ARCH_MINOR),5)
                    NVCC_FLAGS = -arch=sm_75
                else
                    NVCC_FLAGS = -arch=sm_70
                endif
            else ifeq ($(shell expr $(GPU_ARCH_MAJOR) \>= 6), 1)
                NVCC_FLAGS = -arch=sm_60
            else
                NVCC_FLAGS = -arch=sm_50
            endif
        else
            NVCC_FLAGS = -arch=sm_75
        endif
    else
        NVCC_FLAGS = -arch=sm_75
    endif
else
    ENABLE_CUDA = 0
endif

# 4. Enhanced diagnostic output with clear build intention
ifneq ($(MAKECMDGOALS),clean)
$(info ===========================================)
$(info Operating System: $(UNAME_S))
$(info 🔍 GPU Environment Detection Report)
$(info ===========================================)
$(info NVIDIA Driver: $(if $(filter 1,$(GPU_DRIVER_PRESENT)),✅,❌))
$(info CUDA Compiler: $(if $(filter 1,$(NVCC_AVAILABLE)),✅,❌))
$(info CUDA Headers: $(if $(filter 1,$(CUDA_HEADERS_PRESENT)),✅,❌))
$(info CUDA Libraries: $(if $(filter 1,$(CUDA_LIBS_PRESENT)),✅,❌))
$(info -------------------------------------------)

# Always build at least one version
ifneq ($(ENABLE_CUDA),1)
$(info 🚨 GPU support not available, building CPU version)
$(info 💡 To enable GPU support:)
$(info   • Install NVIDIA drivers and CUDA Toolkit)
$(info   • Or use 'make gpu' to try GPU version anyway)
else
$(info ✅ GPU support available, building GPU-accelerated version)
$(info 🎯 Detected GPU Architecture: $(NVCC_FLAGS))
$(info 💡 Use 'make cpu' to force CPU-only version)
endif
$(info ===========================================)
endif

# ============================================================================
# Build Configuration - Always produce executable
# ============================================================================

# Default to CPU version if GPU detection fails, but allow explicit GPU build
ifeq ($(ENABLE_CUDA),1)
    # GPU version configuration
    OBJS = kthread.o bseq.o map.o trfx.o GetTopPeriods_cuda.o
    INCLUDES = -I$(CUDA_PATH)/include
    CFLAGS += -DHAVE_GPU
    NVCC = nvcc
    
    # Set library paths based on platform
    ifeq ($(UNAME_S), Darwin)
        CUDA_LIB_PATH = $(CUDA_PATH)/lib
        LDFLAGS += -L$(CUDA_LIB_PATH) -lcuda -lcudart
        RPATH_FLAG = -Wl,-rpath,$(CUDA_LIB_PATH)
    else
        ifneq ($(wildcard $(CUDA_PATH)/lib64),)
            CUDA_LIB_PATH = $(CUDA_PATH)/lib64
        else
            CUDA_LIB_PATH = $(CUDA_PATH)/lib
        endif
        LDFLAGS += -L$(CUDA_LIB_PATH) -lcuda -lcudart
        RPATH_FLAG = -Wl,-rpath,$(CUDA_LIB_PATH)
    endif
else
    # CPU version configuration (fallback)
    OBJS = kthread.o bseq.o map.o trfx.o
    INCLUDES = 
    CFLAGS += -DNO_GPU
    RPATH_FLAG =
endif

PROG = trfx

# ============================================================================
# Build Rules - Ensure executable is always produced
# ============================================================================

# Default target - always build an executable
all: $(PROG)

# Compile C files
%.o: %.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

# Compile CUDA files (only when CUDA is enabled and file exists)
ifeq ($(ENABLE_CUDA),1)
ifneq ($(wildcard GetTopPeriods_cuda.cu),)
%.o: %.cu
	$(NVCC) -c $(NVCC_FLAGS) $(INCLUDES) $< -o $@
else
$(warning ⚠️  GetTopPeriods_cuda.cu not found, falling back to CPU version)
ENABLE_CUDA = 0
OBJS = kthread.o bseq.o map.o trfx.o
CFLAGS += -DNO_GPU
endif
endif

# Main target - build based on availability (USING GCC FOR LINKING)
$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(RPATH_FLAG)
	@echo "✅ Build completed: $(PROG) ($(if $(filter 1,$(ENABLE_CUDA)),GPU-accelerated,CPU) version) - Built with GCC"

# Build GPU version (explicitly specified)
gpu: 
	@if [ "$(GPU_DRIVER_PRESENT)" = "0" ]; then \
		echo "❌ Cannot build GPU version: NVIDIA driver not found"; \
		echo "💡 Building CPU version instead..."; \
		$(MAKE) cpu; \
	elif [ "$(NVCC_AVAILABLE)" = "0" ]; then \
		echo "❌ Cannot build GPU version: nvcc compiler not found"; \
		echo "💡 Building CPU version instead..."; \
		$(MAKE) cpu; \
	elif [ "$(CUDA_HEADERS_PRESENT)" = "0" ]; then \
		echo "❌ Cannot build GPU version: CUDA headers not found"; \
		echo "💡 Building CPU version instead..."; \
		$(MAKE) cpu; \
	elif [ "$(CUDA_LIBS_PRESENT)" = "0" ]; then \
		echo "❌ Cannot build GPU version: CUDA libraries not found"; \
		echo "💡 Building CPU version instead..."; \
		$(MAKE) cpu; \
	else \
		echo "🔧 Building GPU version..."; \
		$(MAKE) ENABLE_CUDA=1 OBJS="kthread.o bseq.o map.o trfx.o GetTopPeriods_cuda.o" \
			CFLAGS="$(CFLAGS) -DHAVE_GPU" INCLUDES="-I$(CUDA_PATH)/include" $(PROG); \
	fi

# Build CPU version (explicitly specified)
cpu:
	@echo "🔧 Building CPU version..."
	@$(MAKE) ENABLE_CUDA=0 OBJS="kthread.o bseq.o map.o trfx.o" \
		CFLAGS="$(filter-out -DHAVE_GPU,$(CFLAGS)) -DNO_GPU" INCLUDES="" $(PROG)

# Create static library
libtrfx.a: $(OBJS)
	$(AR) rcs $@ $(OBJS)
	@echo "✅ Static library created: $@"

# Clean
clean:
	rm -f *.o *.a $(PROG) *~ *.dSYM core.*
	@echo "🧹 All build files cleaned"

# Install
install: $(PROG)
	cp $(PROG) /usr/local/bin/
	@echo "📦 Installed to /usr/local/bin/$(PROG)"

# Detailed information
info:
	@echo "=== TRFx Build Configuration ==="
	@echo "🔧 Operating System: $(UNAME_S)"
	@echo "🔧 C Compiler: $(CC)"
	@echo "🔧 CUDA Support: $(if $(filter 1,$(ENABLE_CUDA)),✅ Enabled,❌ Disabled)"
	@if [ "$(ENABLE_CUDA)" = "1" ]; then \
		echo "🔧 CUDA Path: $(CUDA_PATH)"; \
		echo "🔧 GPU Architecture: $(NVCC_FLAGS)"; \
	fi
	@echo "🔧 Object Files: $(OBJS)"
	@echo "🔧 Compilation Flags: $(CFLAGS)"
	@echo "🔧 Link Flags: $(LDFLAGS)"
	@echo ""
	@echo "=== GPU Detection Details ==="
	@echo "🔍 NVIDIA Driver: $(if $(filter 1,$(GPU_DRIVER_PRESENT)),✅,❌)"
	@echo "🔍 CUDA Compiler: $(if $(filter 1,$(NVCC_AVAILABLE)),✅,❌)"
	@echo "🔍 CUDA Headers: $(if $(filter 1,$(CUDA_HEADERS_PRESENT)),✅,❌)"
	@echo "🔍 CUDA Libraries: $(if $(filter 1,$(CUDA_LIBS_PRESENT)),✅,❌)"

# Test
test: $(PROG)
	@echo "🧪 Testing $(PROG) ($(if $(filter 1,$(ENABLE_CUDA)),GPU,CPU) version)..."
	@echo ">test_seq" > test.fasta
	@echo "ACGTACGTACGTACGTACGTACGTACGTACGT" >> test.fasta
	@./$(PROG) test.fasta 2>&1 | head -5 || echo "⚠️  Program execution failed, please check compilation"
	@rm -f test.fasta
	@echo "✅ Test completed"

# Performance test
benchmark: $(PROG)
	@echo "⚡ Performance testing $(PROG) ($(if $(filter 1,$(ENABLE_CUDA)),GPU,CPU) version)..."
	@time ./$(PROG) test/input.fasta 2>&1 | tail -5 || echo "⚠️  Test file does not exist"

# GPU diagnostic command
gpu-info:
	@echo "=== GPU System Diagnosis ==="
	@echo "🔧 Check nvidia-smi:"
	@which nvidia-smi > /dev/null 2>&1 && nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader || echo "❌ nvidia-smi not available"
	@echo ""
	@echo "🔧 Check nvcc:"
	@which nvcc > /dev/null 2>&1 && nvcc --version || echo "❌ nvcc not available"
	@echo ""
	@echo "🔧 Check CUDA Path: $(CUDA_PATH)"
	@echo "🔧 Check CUDA Headers: $(if $(filter 1,$(CUDA_HEADERS_PRESENT)),✅Present,❌Missing)"
	@echo "🔧 Check CUDA Libraries: $(if $(filter 1,$(CUDA_LIBS_PRESENT)),✅Present,❌Missing)"

# Help
help:
	@echo "Available commands:"
	@echo "  make           - Auto-detect and build appropriate version"
	@echo "  make gpu       - Force build GPU version (if available)"
	@echo "  make cpu       - Force build CPU version"
	@echo "  make clean     - Clean build files"
	@echo "  make info      - Show detailed build configuration"
	@echo "  make gpu-info  - GPU system diagnostic information"
	@echo "  make test      - Run simple test"
	@echo "  make install   - Install to /usr/local/bin"
	@echo ""
	@echo "Current configuration:"
	@echo "  Operating System: $(UNAME_S)"
	@echo "  CUDA Available: $(if $(filter 1,$(ENABLE_CUDA)),✅Yes,❌No)"

.PHONY: all gpu cpu clean install info test benchmark help gpu-info
