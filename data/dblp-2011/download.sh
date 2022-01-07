basename=dblp-2011
for ext in .properties .graph .md5sums; do
    wget -c http://data.law.di.unimi.it/webdata/$basename/$basename$ext
    wget -c http://data.law.di.unimi.it/webdata/$basename/$basename-t$ext
done

