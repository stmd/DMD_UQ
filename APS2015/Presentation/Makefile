filename=APS2015

all: 
	pdflatex --shell-escape ${filename}
	bibtex ${filename}||true
	pdflatex --shell-escape ${filename}
	pdflatex --shell-escape ${filename}

read:
	evince ${filename}.pdf &

clean:
	rm -f ${filename}.{ps,pdf,log,aux,out,dvi,bbl,blg,glo,pdfsync}
