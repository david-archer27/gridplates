#!/bin/tcsh
rm -r upload.$1
mkdir upload.$1
sed "s/output_tag/$1/g" gridplates.html > upload.$1/index.html
cd outfiles.$1/animations
cp *mp4 ../../upload.$1
foreach f (*)
    echo $1 $f
	#cp *mp4 ../upload.$1/
	pwd
	echo cp $f/img.001.png ../../upload.$1/$f.0.png
	cp $f/img.001.png ../../upload.$1/$f.0.png
end
cd ../charts
cp * ../../upload.$1



