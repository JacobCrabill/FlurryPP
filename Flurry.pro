TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += src/global.cpp \
           src/input.cpp \
		   src/ele.cpp \
		   src/polynomials.cpp \
		   src/operators.cpp \
    src/mesh.cpp \
    src/output.cpp
		   
HEADERS += include/global.hpp \
           include/input.hpp \
		   include/ele.hpp \
		   include/polynomials.hpp \
		   include/operators.hpp \
    include/mesh.hpp \
    include/output.hpp

DISTFILES += \
    README.md