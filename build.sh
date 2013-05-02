# build images
cd images; ./build.sh; cd ..
# build tex
pdflatex main.tex
bibtex   main.aux
pdflatex main.tex
pdflatex main.tex
