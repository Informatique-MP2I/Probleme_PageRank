SRC := Nom_PageRank

EXE := pagerank
CFLAGS := -Wall -O2

.PHONY: clean

$(EXE): $(SRC).c
	gcc $(CFLAGS) -o $@ $< -lm

clean:
	rm $(EXE) *~ 
