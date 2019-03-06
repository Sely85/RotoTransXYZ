RotoTransXYZ
============

Welcome to RotoTransXYZ code!
This code will rotate (around an axis) or translate your xyz configuration file.


BUILD (Linux)
-------------
g++ -lm RotoTransXYZ.cpp -o RotoTransXYZ


USAGE
-----
./RotoTransXYZ <file.xyz> <options>
where options are 
     translate: -t x y z 
     OR 
     rotate: -r axis degrees


EXAMPLES
--------
Rotation of a benzene molecule around z axis of 30°  
./RotoTransXYZ benzene.xyz -r z 30

Translation vector (10, 5, 15) applied to a benzene molecule
./RotoTransXYZ benzene.xyz -t 10 5 15
