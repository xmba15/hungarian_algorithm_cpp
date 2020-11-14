UTEST=OFF
BUILD_EXAMPLES=OFF
CMAKE_ARGS:=$(CMAKE_ARGS)

default:
	@mkdir -p build
	@cd build && cmake .. -DBUILD_TEST=$(UTEST) -DBUILD_EXAMPLES=$(BUILD_EXAMPLES) -DCMAKE_BUILD_TYPE=Release $(CMAKE_ARGS) && make

apps:
	@make default BUILD_EXAMPLES=ON

unittest:
	@make default UTEST=ON
	@cd build && make test

clean:
	@rm -rf build*
