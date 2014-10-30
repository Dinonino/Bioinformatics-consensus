#-------------------------------------------------
#
# Project created by QtCreator 2014-10-30T20:27:54
#
#-------------------------------------------------

QT       -= gui

TARGET = MakeConsensus
TEMPLATE = lib

DEFINES += MAKECONSENSUS_LIBRARY

SOURCES += makeconsensus.cpp

HEADERS += makeconsensus.h\
        makeconsensus_global.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
