#PARTIE A MODIFIER : Liste des fichiers .cpp (et uniquement les .cpp) à compiler
SOURCES=projet.cpp test.cpp
#FIN DE LA PARTIE A MODIFIER
#Nom du compilateur
CXX=g++
# Arguments pour l'etape d'édition de liens : c++11 standart et tous les avertissements
LDFLAGS=-std=c++11 -Wall -Wextra
# Arguments pour l'etape de compilation : c++11 standart, tous les avertissements(et même plus) sauf "ignored-attributes" (présence d'un bug dans GCC version 6et plus)
CPPFLAGS=-std=c++11 -Wall -Wextra -Wno-ignored-attributes
# Librairies : none
LDLIBS=
# Liste des fichiers objets (*.o), générée automagiquement
OBJETS=$(subst .cpp,.o,$(SOURCES))

all: build

build: $(OBJETS)
	$(CXX) $(LDFLAGS) -o run $(OBJETS) $(LDLIBS)

depend: .depend
.depend: $(SOURCES)
	rm -f ./.depend
	$(CXX) $(CPPFLAGS)-MM $^>>./.depend;

clean:
	rm $(OBJETS)

dist-clean:clean
	rm -f *~ .depend

#include .depend
