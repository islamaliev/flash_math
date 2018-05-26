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

    def package_info(self):
        libs = ['flash_math']

        self.cpp_info.libs = []
        for lib in libs:
            if self.settings.os == "Windows" and self.settings.build_type == "Debug":
                suffix = "d"
            else:
                suffix = ""
            self.cpp_info.libs += ["%s%s" % (lib, suffix)]

