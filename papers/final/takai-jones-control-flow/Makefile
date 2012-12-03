# This is here just to remember how to create an image from the `dot'
# file. Graphviz and `dot' are required.

.PHONY: all clean

all: diagram.png

diagram.png: diagram.dot
	dot -Tpng diagram.dot -o diagram.png

clean:
	rm diagram.png
