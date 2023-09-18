To run code please execute :

make clean
make lib && make main.exe && (cd bin && echo && echo && ./main.exe ; echo && echo && cd ..)

where main represents which main(main_wavelets,main_artificial_signal,main_nmreader) you wish to execute, for example: main_wavelets

In /bin we can find:
-object files and executables of the program;

In /Data.AMS we can find:
-some data of the AMS experience;

In /lib we can find:
-the library being used to run this code;

In /main we can find:
-the 3 main programs with some examples for an artificial signal (main_artificial_signal),NMReader(main_nmreader) and wavelets(main_wavelets). 

 In /src we can find:
-the class NMReader(both .C and .h files) being implemented to read neutron monitor data online and some methods which can study the data(correlations and moving averages);
-Also the class wavelets(both .C and .h files) which has some methods for studying periodicity of a given time-series. 

  

