CC ?= gcc
PKG_CONFIG ?= pkg-config
CSTD ?= -std=c17
OPT ?= -O3
NATIVE ?= 1
WARN ?= -Wall -Wextra -Wpedantic

ifeq ($(NATIVE),1)
ARCH := -march=native
endif

CFLAGS ?= $(OPT) $(ARCH) $(CSTD) $(WARN)
JSON_C_CFLAGS := $(shell $(PKG_CONFIG) --cflags json-c)
JSON_C_LIBS := $(shell $(PKG_CONFIG) --libs json-c)
LDLIBS_COMMON := -lgmp -lm

SRC := src
BUILD := build
COMMON := $(SRC)/ramanujan_c_common.c
HDR := $(SRC)/ramanujan_c_common.h
CERT := $(BUILD)/Ramanujan_10K4_Implicit_Certificate_C.json
CERT_TXT := $(BUILD)/Ramanujan_10K4_Implicit_Certificate_C.txt
PI1000 := $(BUILD)/pi1000.txt
PI1000_EXPLICIT := $(BUILD)/pi1000_explicit.txt
PI1000_TRACE := $(BUILD)/pi1000_newton_trace.txt
EXPECTED_SHA256 := e898fea26734a6d3af5396b9f4c60ae5dcc88fc40944d835911a9ee8a672ea1b

.PHONY: all clean core certify restore1000 explicit1000 trace1000 dimensionbench loworderbench d2basicbench smoke

all: $(BUILD)/ramanujan_core $(BUILD)/ramanujan_certifier $(BUILD)/ramanujan_restore $(BUILD)/ramanujan_explicit_pi

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/ramanujan_core: $(SRC)/ramanujan_nested_core.c $(COMMON) $(HDR) | $(BUILD)
	$(CC) $(CFLAGS) -I$(SRC) $(SRC)/ramanujan_nested_core.c $(COMMON) $(LDLIBS_COMMON) -o $@

$(BUILD)/ramanujan_certifier: $(SRC)/ramanujan_10k4_certifier.c $(COMMON) $(HDR) | $(BUILD)
	$(CC) $(CFLAGS) $(JSON_C_CFLAGS) -I$(SRC) $(SRC)/ramanujan_10k4_certifier.c $(COMMON) $(JSON_C_LIBS) $(LDLIBS_COMMON) -o $@

$(BUILD)/ramanujan_restore: $(SRC)/ramanujan_certificate_restore.c $(COMMON) $(HDR) | $(BUILD)
	$(CC) $(CFLAGS) $(JSON_C_CFLAGS) -I$(SRC) $(SRC)/ramanujan_certificate_restore.c $(COMMON) $(JSON_C_LIBS) $(LDLIBS_COMMON) -o $@

$(BUILD)/ramanujan_explicit_pi: $(SRC)/ramanujan_explicit_pi.c $(COMMON) $(HDR) | $(BUILD)
	$(CC) $(CFLAGS) -I$(SRC) $(SRC)/ramanujan_explicit_pi.c $(COMMON) $(LDLIBS_COMMON) -o $@

$(BUILD)/ramanujan_dimension_bench: validation/benchmark_modular_dimension_sweetspot.c $(COMMON) $(HDR) | $(BUILD)
	$(CC) $(CFLAGS) -I$(SRC) validation/benchmark_modular_dimension_sweetspot.c $(COMMON) $(LDLIBS_COMMON) -o $@

$(BUILD)/ramanujan_low_order_bench: validation/benchmark_low_order_explicit_vs_chud.c $(COMMON) $(HDR) | $(BUILD)
	$(CC) $(CFLAGS) -I$(SRC) validation/benchmark_low_order_explicit_vs_chud.c $(COMMON) $(LDLIBS_COMMON) -o $@

$(BUILD)/ramanujan_d2_basic_bench: validation/benchmark_original_d2_basic_vs_chud.c $(COMMON) $(HDR) | $(BUILD)
	$(CC) $(CFLAGS) -I$(SRC) validation/benchmark_original_d2_basic_vs_chud.c $(COMMON) $(LDLIBS_COMMON) -o $@

core: $(BUILD)/ramanujan_core
	$(BUILD)/ramanujan_core

certify: $(BUILD)/ramanujan_certifier
	$(BUILD)/ramanujan_certifier $(CERT) $(CERT_TXT)

restore1000: certify $(BUILD)/ramanujan_restore
	$(BUILD)/ramanujan_restore $(CERT) --digits 1000 --guard 140 --output $(PI1000)

explicit1000: $(BUILD)/ramanujan_explicit_pi
	$(BUILD)/ramanujan_explicit_pi --digits 1000 --guard 120 --output $(PI1000_EXPLICIT) --quiet

trace1000: certify $(BUILD)/ramanujan_restore
	$(BUILD)/ramanujan_restore $(CERT) --digits 1000 --guard 140 --output $(PI1000_TRACE) --trace-newton

dimensionbench: $(BUILD)/ramanujan_dimension_bench
	$(BUILD)/ramanujan_dimension_bench

loworderbench: $(BUILD)/ramanujan_low_order_bench
	$(BUILD)/ramanujan_low_order_bench 1000 10000 30000 100000 300000

d2basicbench: $(BUILD)/ramanujan_d2_basic_bench
	$(BUILD)/ramanujan_d2_basic_bench 1000 10000 30000 100000 200000

smoke: restore1000 explicit1000
	@restored=$$(sha256sum $(PI1000) | awk '{print $$1}'); \
	 explicit=$$(sha256sum $(PI1000_EXPLICIT) | awk '{print $$1}'); \
	 echo "expected: $(EXPECTED_SHA256)"; \
	 echo "restore:  $$restored"; \
	 echo "explicit: $$explicit"; \
	 test "$$restored" = "$(EXPECTED_SHA256)"; \
	 test "$$explicit" = "$(EXPECTED_SHA256)"; \
	 cmp $(PI1000) $(PI1000_EXPLICIT); \
	 echo "SMOKE TEST: PASS"

clean:
	rm -rf $(BUILD)
