
*** To compile:

 * Linux / GNU:

 Type "make heatflow".

 * Other systems:

 Compile and link *.cpp.


*** To run:

 Run 'heatflow' and redirect one of the config files in,
 e.g. in unix: 
 
   ./heatflow < config/input1d

 or in Windows:

   heatflow < config\input1d

 See sample input files for input file format.


*** Output:

 Output (T and effective k field) will be in CSV format, outputted 
 in to working directory. Any outputted temperature field (T.nnn.csv) can
 be renamed to Tini.csv and used as an input for initial temperature field (see
 sample config files for the option.) CSV field order is x,y,z,value (number
 of coordinate fields = number of dimensions in the problem).


