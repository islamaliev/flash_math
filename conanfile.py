from conans import ConanFile, CMake

class FlashMathConan(ConanFile):
    name = "flash_math"
    version = "0.1"
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"
    #default_options = ""
    exports_sources = "*"

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        self.copy("*.h", "include", "src")
        self.copy("*.lib", "lib", "", keep_path=False)
        self.copy("*.a", "lib", "", keep_path=False)

