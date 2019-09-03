echo %date% >>log_corn.txt
set month=%date:~4,2%
echo month=%month%  >>log_corn.txt
set day=%date:~7,2%
echo day=%day%  >>log_corn.txt
set year=%date:~10,4%
echo year=%year%  >>log_corn.txt

cd C:\projects\IFPRI\FSP\MATLAB\outputfolder

set dateoftheweek=%date:~0,3%
echo dateoftheweek=%dateoftheweek%  

REM CREATE INPUT FILE NAME 

set input_file_name=corn%month%%day%%year%.txt

echo inputfilename=%input_file_name% >>log_corn.txt
REM COPY THE INPUT FILE FROM THE DROPBOX, AND CHANGE THE FILE NAME INTO {commodity}.txt

copy "C:\Dropbox (IFPRI)\gauss\input\%input_file_name%" "C:\projects\IFPRI\FSP\MATLAB\input\corn.txt" /y 

REM RUN THE GAUSS PROGRAM 
"C:\Program Files\MATLAB\R2019a\bin\matlab.exe" -nosplash -noFigureWindows -r "try; run('C:\Users\SOONHOKIM\Documents\MATLAB\summer_2019\sbk_2Lags_corn.m'); catch; end; quit" 


timeout /t 3600 /nobreak

set filename="output_corn_%year%_%month%_%day%.txt"

echo %filename%  >>log_corn.txt

REM CHANGE THE LOCATION OF FOLDER 
cd \
cd C:\projects\IFPRI\FSP\MATLAB\outputfolder

REM COPY FILES 

rename output_corn.txt %filename%  


copy C:\projects\IFPRI\FSP\MATLAB\outputfolder\%filename% "C:\Dropbox (IFPRI)\MATLAB\output\%filename%" /y 
copy C:\projects\IFPRI\FSP\MATLAB\outputfolder\corn.png "C:\Dropbox (IFPRI)\MATLAB\output\corn.png" /y

REM GO BACK TO THE BATCH FOLDER  
cd \
cd C:\projects\IFPRI\FSP\MATLAB\batch

echo done_in_%year%_%month%_%day%  >>log_corn.txt

