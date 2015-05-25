TEMPLATE = app
CONFIG += console debug
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++11

INCLUDEPATH += $$PWD/include

SOURCES += src/global.cpp \
    src/matrix.cpp \
    src/input.cpp \
    src/ele.cpp \
    src/polynomials.cpp \
    src/operators.cpp \
    src/geo.cpp \
    src/output.cpp \
    src/face.cpp \
    src/flux.cpp \
    src/flurry.cpp \
    src/solver.cpp \
    include/geo.inl \
    src/intFace.cpp \
    src/boundFace.cpp \
    src/mpiFace.cpp
		   
HEADERS += include/global.hpp \
    include/matrix.hpp \
    include/input.hpp \
    include/ele.hpp \
    include/polynomials.hpp \
    include/operators.hpp \
    include/geo.hpp \
    include/output.hpp \
    include/face.hpp \
    include/flux.hpp \
    include/flurry.hpp \
    include/solver.hpp \
    include/error.hpp \
    include/intFace.hpp \
    include/boundFace.hpp \
    include/mpiFace.hpp

DISTFILES += \
    README.md \

OTHER_FILES += \
    bin/input_supwall \
    tests/euler/cylinder/input_cyl \
    planning \
    makefile
