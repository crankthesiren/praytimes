# praytimes library port from praytimes.org to C

Port of praytimes.org C++ code to C; lightweight offline use  for embedded system like Raspberry Pi


Original source: http://praytimes.org/wiki/Code

# How to compile

Run the following in your raspberry pi terminal:

gcc prayertimes.c -lm -o prayertimes 

# How to use compiled executable

Run executable with longitude and lattitude of location to get prayer times: 

./prayertimes -l 44 -n 44

# Additional parameters

Please run program without parameters, it should print out all available options
