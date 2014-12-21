TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += src/global.cpp \
           src/input.cpp \
		   src/ele.cpp \
		   src/polynomials.cpp \
		   src/operators.cpp \
    src/mesh.cpp \
    src/output.cpp \
    src/face.cpp \
    src/flux.cpp \
    geo.cpp \
    src/geo.cpp
		   
HEADERS += include/global.hpp \
           include/input.hpp \
		   include/ele.hpp \
		   include/polynomials.hpp \
		   include/operators.hpp \
    include/mesh.hpp \
    include/output.hpp \
    include/face.hpp \
    include/flux.hpp \
    geo.hpp

DISTFILES += \
    README.md