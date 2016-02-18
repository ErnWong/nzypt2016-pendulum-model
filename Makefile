model.exe: model.cpp
	g++ -std=gnu++11 -O2 -o $@ $^ -lm -Wall -Wextra

program.exe: main.cpp
	g++ -std=gnu++11 -O2 -o $@ $^ -lm -Wall -Wextra

program.mypers.exe: main.mypers.cpp
	g++ -std=gnu++11 -O2 -o $@ $^ -lm -Wall -Wextra

tocsv.exe: tocsv.cpp
	g++ -std=gnu++11 -O2 -o $@ $^ -lm -Wall -Wextra

testargs.exe: testargs.cpp
	g++ -std=gnu++11 -O2 -o $@ $^ -lm -Wall -Wextra
