#
#  Makefile for Linux 
#
#     make all   (construction de l'executable)
#     make clean (effacement des fichiers objets et de l'executable)
#
#  A adapter en fonction des ordinateurs/environnements 
#  Compilateur, edition de liens, 
#
#
CC       = gcc  
LD       = gcc
CFLAGS   = -O3 -Dgraphic -Wall -g
LFLAGS   = -Wall -O3 -g
LIBS     = -lglfw -lm -lGL -lGLU
#
PROG     = myFem
LISTEOBJ = \
  tsunami.o main.o homework.o
# ATTENTION... aucun caractere apres le caractere de continuation "\"
#
# compilation
#
.c.o :
	$(CC) -c  $(CFLAGS) -o $@ $<
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# dependances
#
all        : $(PROG)
homework.o : homework.c tsunami.h
tsunami.o  : tsunami.c tsunami.h
main.o     : main.c tsunami.h

#
# edition de lien
#
$(PROG) : $(LISTEOBJ)
	$(LD) -o $(PROG) $(LFLAGS) $(LISTEOBJ) $(LIBS)
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# effacement des fichiers intermediaires
#
clean :
	rm -vf $(PROG) $(LISTEOBJ) core a.out
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# ATTENTION... il faut une ligne vide a la fin du fichier.


