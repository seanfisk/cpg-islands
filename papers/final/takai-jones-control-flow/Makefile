# This is here just to remember how to create an image from the `dot'
# file. Graphviz and `dot' are required.

.PHONY: all clean

all: diagram.png diagram-nogoto.png

diagram.png: diagram.dot
	dot -Tpng diagram.dot -o diagram.png

diagram-nogoto.png: diagram-nogoto.dot
	dot -Tpng diagram-nogoto.dot -o diagram-nogoto.png

clean:
	rm *.png
