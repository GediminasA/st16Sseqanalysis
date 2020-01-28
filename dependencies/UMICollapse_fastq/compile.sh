rm -r umicollapse
rm test.jar 
javac -cp lib/htsjdk-2.19.0.jar:lib/snappy-java-1.1.7.3.jar  -sourcepath src -d .  -verbose src/umicollapse/**/*.java
jar cfmv  test.jar Manifest.txt umicollapse/*
