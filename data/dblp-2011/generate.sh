# data file basename
basename=dblp-2011
# generate .obl .offsets files
java -cp "../../lib/*" it.unimi.dsi.webgraph.BVGraph -o -O -L $basename
java -cp "../../lib/*" it.unimi.dsi.webgraph.BVGraph -o -O -L $basename-t
# generate sym file
java -cp "../../lib/*" it.unimi.dsi.webgraph.Transform union $basename $basename-t $basename-sym

