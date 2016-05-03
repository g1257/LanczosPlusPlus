CPPFLAGS =  -I../src/Engine -I../src/Models/FeBasedSc  -I../src/Models/Heisenberg  
CPPFLAGS += -I../src/Models/HubbardOneOrbital  -I../src/Models/Immm  -I../src/Models/TjMultiOrb
CPPFLAGS += -I../../PsimagLite/src -I../../PsimagLite/src/Geometry 
CXX = g++ 

all: manual.pdf

manual.tex: manual.ptex 
	find ../src -iname "*.h" -or -iname "*.cpp" |\
              perl ../../PsimagLite/scripts/doc.pl manual.ptex

manual.pdf: manual.tex 
	pdflatex manual.tex
	bibtex manual.aux
	pdflatex manual.tex
	pdflatex manual.tex

clean:
	rm -f manual.tex
	rm -f manual.pdf


