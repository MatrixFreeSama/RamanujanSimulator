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
EXPECTED_SHA256 := e898fea26734a6d3af5396b9f4c60ae5dcc88fc40944d835911a9ee8a672ea1b

.PHONY: all clean core certify restore1000 smoke

all: $(BUILD)/ramanujan_core $(BUILD)/ramanujan_certifier $(BUILD)/ramanujan_restore

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/ramanujan_core: $(SRC)/ramanujan_nested_core.c $(COMMON) $(HDR) | $(BUILD)
	$(CC) $(CFLAGS) -I$(SRC) $(SRC)/ramanujan_nested_core.c $(COMMON) $(LDLIBS_COMMON) -o $@

$(BUILD)/ramanujan_certifier: $(SRC)/ramanujan_10k4_certifier.c $(COMMON) $(HDR) | $(BUILD)
	$(CC) $(CFLAGS) $(JSON_C_CFLAGS) -I$(SRC) $(SRC)/ramanujan_10k4_certifier.c $(COMMON) $(JSON_C_LIBS) $(LDLIBS_COMMON) -o $@

$(BUILD)/ramanujan_restore: $(SRC)/ramanujan_certificate_restore.c $(COMMON) $(HDR) | $(BUILD)
	$(CC) $(CFLAGS) $(JSON_C_CFLAGS) -I$(SRC) $(SRC)/ramanujan_certificate_restore.c $(COMMON) $(JSON_C_LIBS) $(LDLIBS_COMMON) -o $@

core: $(BUILD)/ramanujan_core
	$(BUILD)/ramanujan_core

certify: $(BUILD)/ramanujan_certifier
	$(BUILD)/ramanujan_certifier $(CERT) $(CERT_TXT)

restore1000: certify $(BUILD)/ramanujan_restore
	$(BUILD)/ramanujan_restore $(CERT) --digits 1000 --guard 140 --output $(PI1000)

smoke: restore1000
	@actual=$$(sha256sum $(PI1000) | awk '{print $$1}'); \
	 echo "expected: $(EXPECTED_SHA256)"; \
	 echo "actual:   $$actual"; \
	 test "$$actual" = "$(EXPECTED_SHA256)" && echo "SMOKE TEST: PASS"

clean:
	rm -rf $(BUILD)
