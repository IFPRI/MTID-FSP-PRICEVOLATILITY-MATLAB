echo %date% >>log_wheatkcbt.txt
set month=%date:~4,2%
echo month=%month%  >>log_wheatkcbt.txt
set day=%date:~7,2%
echo day=%day%  >>log_wheatkcbt.txt
set year=%date:~10,4%
echo year=%year%  >>log_wheatkcbt.txt

cd C:\projects\IFPRI\FSP\MATLAB\outputfolder

set dateoftheweek=%date:~0,3%
echo dateoftheweek=%dateoftheweek%  

REM CREATE INPUT FILE NAME 

set input_file_name=wheatkcbt%month%%day%%year%.txt

echo inputfilename=%input_file_name% >>log_wheatkcbt.txt
REM COPY THE INPUT FILE FROM THE DROPBOX, AND CHANGE THE FILE NAME INTO {commodity}.txt

copy "C:\Dropbox (IFPRI)\gauss\input\%input_file_name%" "C:\projects\IFPRI\FSP\MATLAB\input\wheatkcbt.txt" /y 

REM RUN THE GAUSS PROGRAM 
"C:\Program Files\MATLAB\R2019a\bin\matlab.exe" -nosplash -noFigureWindows -r "try; run('C:\Users\SOONHOKIM\Documents\MATLAB\summer_2019\sbk_2Lags_wheatkcbt.m'); catch; end; quit" 


timeout /t 3600 /nobreak

set filename="output_wheatkcbt_%year%_%month%_%day%.txt"

echo %filename%  >>log_wheatkcbt.txt

REM CHANGE THE LOCATION OF FOLDER 
cd \
cd C:\projects\IFPRI\FSP\MATLAB\outputfolder

REM COPY FILES 

rename output_wheatkcbt.txt %filename%  


copy C:\projects\IFPRI\FSP\MATLAB\outputfolder\%filename% "C:\Dropbox (IFPRI)\MATLAB\output\%filename%" /y 
copy C:\projects\IFPRI\FSP\MATLAB\outputfolder\wheatkcbt.png "C:\Dropbox (IFPRI)\MATLAB\output\wheatkcbt.png" /y

REM GO BACK TO THE BATCH FOLDER  
cd \
cd C:\projects\IFPRI\FSP\MATLAB\batch

echo done_in_%year%_%month%_%day%  >>log_wheatkcbt.txt

