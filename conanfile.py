from conans import ConanFile, CMake

class FlashMathConan(ConanFile):
    name = "flash_math"
    version = "0.1"
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"
    #default_options = ""
    exports_sources = "*"

    def imports(self):
        self.copy("*.dll", dst="bin", src="bin") # From bin to bin
        self.copy("*.dylib*", dst="bin", src="lib") # From lib to bin

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        #cmake.configure(source_dir="include")
        cmake.build()

    def package(self):
        self.copy("*.h", "include", "include")
        self.copy("*.lib", "lib", "", keep_path=False)
        self.copy("*.a", "lib", "", keep_path=False)

    def package_info(self):
        self.cpp_info.libs = ["flash_math"]
    
