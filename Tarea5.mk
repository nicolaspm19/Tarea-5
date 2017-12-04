
CC=g++ -std=c++11 -O2  
NAME=CurvaRotacion

all:
	$(CC) $(NAME).c -o $(NAME)
	./$(NAME)
	python Plots.py
	pdflatex -interaction=nonstopmode "Results_hw5".tex
clean:
	rm *.aux *.log *.pdf *png