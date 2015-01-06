#-------------------------------------------------
#
# Project created by QtCreator 2014-12-11T22:14:30
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = MakeConsensus
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    read.cpp \
    unitig.cpp \
    realigner.cpp \
    nucleic_codes.cpp \
    consensus.cpp

HEADERS += \
    read.h \
    unitig.h \
    realigner.h \
    nucleic_codes.h \
    consensus.h
