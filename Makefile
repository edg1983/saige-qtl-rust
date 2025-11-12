.PHONY: help build test clean docker-build docker-test docker-push all

# Variables
IMAGE_NAME := saige-qtl-rust
IMAGE_TAG := latest
REGISTRY := ghcr.io/edg1983
FULL_IMAGE := $(REGISTRY)/$(IMAGE_NAME):$(IMAGE_TAG)

help:
	@echo "SAIGE-QTL-Rust Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  build         - Build Rust binaries (release mode)"
	@echo "  test          - Run Rust tests"
	@echo "  clean         - Clean build artifacts"
	@echo "  docker-build  - Build Docker image"
	@echo "  docker-test   - Test Docker image"
	@echo "  docker-push   - Push Docker image to registry"
	@echo "  all           - Build and test everything"

build:
	cargo build --release

test:
	cargo test --release

clean:
	cargo clean
	rm -f *.sif

docker-build:
	docker build -t $(IMAGE_NAME):$(IMAGE_TAG) .

docker-test:
	@echo "Testing Docker image..."
	docker run --rm $(IMAGE_NAME):$(IMAGE_TAG) step1-fit-null --help
	docker run --rm $(IMAGE_NAME):$(IMAGE_TAG) step2-run-tests --help
	@echo "Docker image tests passed!"

docker-push: docker-build
	docker tag $(IMAGE_NAME):$(IMAGE_TAG) $(FULL_IMAGE)
	docker push $(FULL_IMAGE)

singularity:
	singularity build $(IMAGE_NAME).sif docker-daemon://$(IMAGE_NAME):$(IMAGE_TAG)

all: build test docker-build docker-test
	@echo "All builds and tests completed successfully!"
