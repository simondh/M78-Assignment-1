/*
    Swansea University
    CS M78, High Performance Computing in C and C++
    Assignment 1
    2 dimensional temperature simulation
    This is variant of tempSim, that repeatedly re-runs with increasing iterations, and records the time taken
    to file tempSim.csv
    Simon Hewitt 806068
    November 2014
 
 temp_dist2D.c
 
 Speed test version, takes several command line parameters:
 
 temp_dist2D Cx Cy nTimes datafile [outputFileDir||""] [finalNTimes]
 
 where
 Cx, Cy are trye temperature siumulation coefficients
 nTimes is the number of timesliec simulations to execute
 datafile is the CSV file with the initial conditions
 outputDirFile is thedirectory to hold the output files. If blank (""), no files are written. NOTE 1 file is written per iteration,
 I recommend this is OFF for large iteration timing runs
 testruns :- temp_dist2D will repeat the cycle from nTimes to finalNTimes increasing x10 each time
 eg 1000 1000000 will generate simulation runs of 1000 10,000 100,000 and 1,000,000. The run times are output for each cycle
 
 Examples:
 Simple run, no output files:  temp_dist2D 0.1 0.1 initialdata.csv 1000
 As above but write data files for each iteration:  temp_dist2D 0.1 0.1 initialdata.csv 1000 ./data
 Large iteration timing run: temp_dist2D 0.1 0.1 initialdata.csv 1000 "" 1000000
 
 
*/

// Includes
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <time.h>



// Global constants 
// The matrix dimensions, 100 X 100
#define XDIM 100
#define YDIM 100
#define FIRST_RUN 10
#define LAST_RUN 1000001


// Error codes returned by tempSim.
enum _tempSim_error
{
    E_SUCCESS = 0,
    E_WRONG_PARAMS = -1,    // Incorrect number of parameters supplied
    E_Cx_ERROR = -2,        // Cx must be integer >= 1,
    E_Cy_ERROR = -3,        // dito Cy
    E_TIMESTEP_ERROR = -4,  // Timesteps must be an integer >= 1
    E_FILENOTFOUNT = -5,    // Data file not found or cannot be opened
    E_DATAERROR = -6,       // Data file contains invalid data
    E_INTERNALERROR = -7,   // internal logic error, should never happen...
    E_FILEOPENFAIL = -8,    // cannot open performance log file for writing
    E_PARAMERROR = -9       // something wrong with command line parameters
};

// type to provide in your API
typedef enum _tempSim_error  error_t;


// Typedefs etc
// a single copy of this struct will  be created in main to make globally significant data easily available
typedef struct {
    double Cx;  /* temperature constant, X axis */
    double Cy;  /* ditto, Y axis */
    unsigned long timeSteps;  /* How many time steps to simulate, */
				/* set on command line invocation of the program */
    unsigned long finalTimeSteps;
    char outputFileName[250]; // will be set to 'output_data_yymmdd-hhmmss-tttttt.csv' where ttt is the timeslice number
    } CommonDataStruct;



// function forward definitions

error_t checkParams (int argc, char *argv[], CommonDataStruct *commonData, FILE **ipFile);
error_t loadFromFile (double temperatureArray [XDIM][YDIM], FILE *ipFile);
void    calculateNextStep (double arrayNow[XDIM][YDIM], double arrayNext[XDIM][YDIM], CommonDataStruct *params);
void    writeData (unsigned long iter,  double arrayNow[XDIM][YDIM], CommonDataStruct *params);
void    setFileName (char *prefix, char* postfix, char *filenameBuffer);
void    recordClock (int command, unsigned long tag);



/*========== main ============*/
/*  returns 0 on success, -ve numbers indicate fatal errors, see enum _tempSim_error (above) and program documentation
    invoke as tempSim Cx, Cy, nTime, dataFile
    where Cx and Cy are the temperature coefficient constants, positive non zero real numbers (double)
	nTime is the number of time slices to simulate, < 65000
	datafile is te CSV file containing the starting data, as a 100X100 numeric array
    example tempSim 0.1 0.1 1000 initial_data.csv

    tempSimSpeed changes - timeslices is ignored (but still needed to avoid changing that routine, any integer)
    And the run is repeated for FIRST_RUN to LAST_RUN in decades, i.e. X10 each time
    (so make sure the end is a power of 10 times the start)
    The performance figures are written to tempSim.csv
 */

int main (int argc, char* argv[])
{
    FILE *ipFile;
    CommonDataStruct commonData;
    error_t localError;

    // the time arrays for the simulation.
    double temperatureArray_1 [XDIM][YDIM];
    double temperatureArray_2 [XDIM][YDIM];
    double *tArray_1_Ptr = (double *)temperatureArray_1;
    double *tArray_2_Ptr = (double *)temperatureArray_2;
    double *tempArrayPtr;
    // To build the next time slice, we only need the current time slice, hence two arrays
    // To avoid copying the arrays , the POINTER tArray_i_ptr is always the current time slice, and _2 the next one,
    // and the pointers are switched after each iteration.

    unsigned long timeCount;
    unsigned long ticker;
    int runCount = 1;
    

    // First check the parameters have been correctly supplied
    // Also opens the file, or reports error if it cannot be opened.
    if ((localError = checkParams(argc, argv, &commonData, &ipFile)))
        {
        fprintf (stderr, "Fatal error incorrect programme parameters\n");
        fprintf (stderr, "tempSim Cx Cy nTime dataFileName [outputFileDir | \"\"] finalNtime\n");
        exit(localError);
        }

    // At this stage, got the program parameters and the data file is open ready to read
    if (commonData.timeSteps < commonData.finalTimeSteps) {
        printf ("tempSim performance running with parameters Cx: %f Cy: %f time slices from: %lu to: %lu\ndata file : %s\n",
                commonData.Cx, commonData.Cy, commonData.timeSteps, commonData.finalTimeSteps, argv[4]);
    } else {
        printf ("tempSim running with parameters Cx: %f Cy: %f time slices : %lu \ndata file : %s\n",
            commonData.Cx, commonData.Cy, commonData.timeSteps, argv[4]);
    }

        if ((localError = loadFromFile (temperatureArray_1, ipFile)))
        {
            fprintf (stderr, "Cannot load data\n"); // detailed error output printed in the function
            exit(localError);
        }
    
    // Got the data in temperatureArray_1, ready to start the time simulation
    // first copy array 1 to array 2, so the edge values are set at zero
    memcpy(tArray_2_Ptr, tArray_1_Ptr, sizeof(double) * XDIM * YDIM);
    recordClock (0,0);

    
    while (commonData.timeSteps <= commonData.finalTimeSteps+1)
    // this loop runs jsut once UNLESS asked to do a performance run, where final param specifies maximum number of steps
    {
        printf ("Starting run: %d for %lu time slices\n", runCount, commonData.timeSteps);
        recordClock (1, commonData.timeSteps);
        
        ticker = commonData.timeSteps / 20; // ticker messages, so we can see something is happening. set to 20 per run.
        
        for (timeCount = 0; timeCount < commonData.timeSteps; timeCount++) {
            calculateNextStep ((double(*)[YDIM])tArray_1_Ptr, (double(*)[YDIM])tArray_2_Ptr, &commonData);
            
            // and swap A and B (1 and 2) for the next round
            // this makes sure array ptr 1 is always for time t, array ptr 2 is overwritten with new t+1 values
            tempArrayPtr = tArray_2_Ptr;
            tArray_2_Ptr = tArray_1_Ptr;
            tArray_1_Ptr = tempArrayPtr;
            if (!(timeCount % ticker)) printf ("Time slice %06lu\n", timeCount); // so we can see if we are still alive...
            // not done for speedtest// writeData (timeCount, (double(*)[YDIM])tArray_1_Ptr, &commonData);
            writeData(timeCount, (double(*)[YDIM])tArray_1_Ptr, &commonData);

        }
        recordClock (2, commonData.timeSteps);

        commonData.timeSteps *= 10; // try again 10 times more iterations!
        runCount++;
    }

	return (E_SUCCESS);
} //=== main ===


void recordClock (int command, unsigned long tag)
{
    // command can be:  0 -> initialise (open file etc)
    //                  1 -> start clock
    //                  2 -> end clock, record result
    // done this way to hide the implementation (hey lets use C++ instead?)
    
    static time_t clockStart, clockEnd;
    static FILE *dataFile;
    static int setup = 0;
    char fname[160];
    
    if ((setup == 0) && (command != 0)) {
        fprintf (stderr, "recordClock function %d called before initialisation\n\n", command);
        exit (E_INTERNALERROR);
    }
    
    switch (command) {
        case 0:
            setFileName ("run_data", ".csv", fname);
            dataFile = fopen (fname, "w");
            if (dataFile == NULL) {
                fprintf (stderr, "Cannot open output file %s\n", fname);
                exit (E_FILEOPENFAIL);
            }
            printf ("Writing run performance data to: %s\n", fname);
            setup = 1;
            break;
            
        case 1:
            clockStart = clock();
            break;
            
        case 2:
            clockEnd = clock();
            fprintf (dataFile, "Run, %lu, took, %f, seconds\n", tag, (double)(clockEnd - clockStart) / CLOCKS_PER_SEC);
            break;
            
        default:
            exit(E_INTERNALERROR); // no reason it should reach here
            break;
    }
    
    
}

 void    writeData (unsigned long iter, double arrayNow[XDIM][YDIM], CommonDataStruct *params)
{
    // write the timeslice to a CSV file, numbers in the form "1234.000, "
    // Data files are written only if parameter 5 is specified, which must be a valid existing directory ("." works of course)
    // The directory must exist. the string must end with a '/', eg /Users/simonhewitt/temp/HPCdata/
    // in bash export temp_dist_data="/Users/simonhewitt/temp/HPCdata/"
    
    int row, col;
    static int errorCount = 0;
    static int first = 1;
    static int justSayNo = 0; // set if output file dir is empty, indicating no write wanted.
    FILE *opFile;
    char fname[160];
    
    if (justSayNo || ( strlen(params->outputFileName) == 0)) {
        justSayNo = 1;  // this just avoids calling strlen every time.
        return; // no output directory specified (run parameter #5) so do not write data files.
    }
    
    sprintf (fname, "%s-%06lu.csv", (char *)params->outputFileName, iter);
    opFile = fopen(fname, "w");
    if (first) {
        printf ("Writing one data file for each timeslice, 1st is:\n%s\n", fname);
        first = 0;
    }
    
    if (opFile == NULL)
    {
        if (errorCount < 3) // don't write potentially thousands of error messages
            fprintf (stderr, "Cannot create timeslice file %s\n", fname);
        errorCount++;
        return;
    }
    
    for (row = 0; row < XDIM; row++) {
        for (col = 0; col < YDIM; col++) {
            fprintf (opFile, "%10.3f, ", arrayNow[row][col]);
        }
        fprintf (opFile, "\n");
    }
    
    fclose (opFile);
    
}// === writeData ===


void    calculateNextStep (double arrayNow[XDIM][YDIM], double arrayNext[XDIM][YDIM], CommonDataStruct *params)
{
    // formula is (from assignment paper)
    // Ux,y,t = Ux,y,t−1 + Cx × (Ux+1,y,t−1 + Ux−1,y,t−1 − 2 × Ux,y,t−1) + Cy × (Ux,y+1,t−1 + Ux,y−1,t−1 − 2 × Ux,y,t−1)
    // x is the row, y the col. t-1 is arrayNow, t is arrayNext
    // the edges start at zero and stay at zero, which makes the calculation much easier
    // we will process [1..98][1..98], knowing the these cells ALWAYS have left, right, up, down cells
    // i.e. no need for special processing for edge cells.
    // Ensure the arrayNext has edges set to zero before starting
    
    int row, col;
    double newValue;
    double Cx = params->Cx;
    double Cy = params->Cy;
    
    
    for (row = 1; row < 99; row++)
        for (col = 1; col < 99; col++) {
            newValue = arrayNow[row][col]
                    + Cx * (arrayNow[row][col+1] + arrayNow[row][col-1] - (2 * arrayNow[row][col]))
                    + Cy * (arrayNow[row+1][col] + arrayNow[row-1][col] - (2 * arrayNow[row][col]));
            arrayNext[row][col] = newValue;
        }
    
} // === calculateNextStep ===


error_t loadFromFile (double temperatureArray [XDIM][YDIM], FILE *ipFile)
{
    // input data file has been opened, with ipFile
    // read YDIM data lines each of XDIM data elements, which are double numbers (100 X 100) and store in timeArray
    // the file must be simple CSV in the form ddd.dd, ddd.dd  -  no quotes are handled. Can start with whitespace ,
    // and any amount of whitespace between digits.
    
#define BUFFSIZE 16000
    // The largest number in the test data is 6002500.000, then <,  >, 14 chars each, make it 16 to be safe, X 100 = 16000
    
    int rows = 0, cols = 0;
    char readBuff[BUFFSIZE];
    double t;
    char *dataPtr, *endPtr;
    
    while (fgets(readBuff, BUFFSIZE, ipFile)) // for each row of data read from file
    {
        if (rows >=  100) {
            fprintf (stderr, "Too many data rows ( : %d) \n", rows);
            return (E_DATAERROR);
        }
        
        cols = 0;
        dataPtr = readBuff;
        
        while (*dataPtr)  // for each data element in the row
        {
            t = strtod(dataPtr, &endPtr);
            temperatureArray[rows][cols] = t;
            cols++;
            if (cols > 100) {
                fprintf (stderr, "Too many data elements in row : %d\n", rows);
                return (E_DATAERROR);
            }
            
            dataPtr = endPtr;
            while ((*dataPtr) && !isdigit (*dataPtr)) dataPtr++;  // stops if \0 end of line or is any digit
        }
        
        if (cols < 100){
            fprintf (stderr, "Too few data elements (: %d) in row : %d\n", cols, rows);
            return (E_DATAERROR);
            }
        rows++;

    }
    
    if (rows != 100) {
        fprintf (stderr, "Too many or too few data rows ( : %d) \n", rows);
        return (E_DATAERROR);
    }
    
    return (E_SUCCESS);
}  // ===  loadFromFile ===




errno_t checkParams (int argc, char *argv[], CommonDataStruct *commonData, FILE **ipFile)
{
    // see notes at top of this file
    // checking for params:
    //  temp_dist2D Cx Cy nTimes datafile [outputFileDir||""] [finalNTimes]
    
    double i;
    unsigned long l;
    time_t rawtime;
    struct tm *timeStruct;
    char buffer[1600];
    char dirBuffer[160];
    
    if ((argc < 5) || (argc >> 7)) // remember parameter[0] is the name of the program
    {
        fprintf (stderr, "Incorrect number of parameters\n");
        return (E_WRONG_PARAMS);
    }
    
    i = atof (argv[1]);
    if (i == 0)
    {
        fprintf (stderr, "Cx is zero or not a number\n");
        return (E_Cx_ERROR);
    }
    commonData->Cx = i;
    
    i = atof (argv[2]);
    if (i == 0)
    {
        fprintf (stderr, "Cy is zero or not a number\n");
        return (E_Cy_ERROR);
    }
    commonData->Cy = i;
    
    i = atol (argv[3]);
    if (i == 0)
    {
        fprintf (stderr, "nTime is zero or not a number\n");
        return (E_TIMESTEP_ERROR);
    }
    commonData->timeSteps = i;
    
    *ipFile = fopen (argv[4], "r");
    if (*ipFile == NULL)
    {
        fprintf (stderr, "Cannot open data file : %s\n", argv[4]);
        perror ("temp_sim2D");
        return (E_FILENOTFOUNT);
    }
    
    if ((argc >= 6) && (strlen(argv[5]) != 0))
    {
        // next, set up the output filename template, based on clock time
        // output_dir/ must be param 5, must be a valid existing directory
        // here we set up output_dir/output_data_yymmdd-hhmmss, the time-slice writer must add the time slice number and filetype -tttttt.csv
        time( &rawtime );
        timeStruct = localtime( &rawtime );
        strcpy (dirBuffer, (char *)argv[5]);
        if (dirBuffer[strlen(dirBuffer)-1] == '/')
        {
            dirBuffer[strlen(dirBuffer)-1] = '\0'; //strip off folder '/', (then add it below), so its DEFINATELY there
        }
        
        sprintf(buffer, "%s/output_data_%02d%02d%02d-%02d%02d%02d", dirBuffer, timeStruct->tm_year%100, timeStruct->tm_mon + 1, timeStruct->tm_mday,
                timeStruct->tm_hour, timeStruct->tm_min, timeStruct->tm_sec);
        strcpy((char *)commonData->outputFileName, buffer);
    } else {
        commonData->outputFileName[0] = '\0';
    }
    
    commonData->finalTimeSteps = commonData->timeSteps;
    if (argc == 7)
    {
        // Yes we are going for a repeated speed run
        l = atol(argv[6]);
        if (l  <= commonData->timeSteps)
        {
            fprintf (stderr, "finalNTimes (%lu) must be greater than nTimes (%lu)\n", l , commonData->timeSteps);
            return (E_PARAMERROR);
        }
        commonData->finalTimeSteps = l;
    }
    
    // All parameters OK, ready to go...
    
    return (E_SUCCESS);
    
} // === checkParams ===


            
void setFileName (char *prefix, char* postfix, char *filenameBuffer)
{
    // create a filename with prefix, postfix and yymmdd-hhmmss in between eg
    // setFileName ("MyData", ".csv", buffer) sets buffer to MyData_141101-213300.csv
    time_t rawtime;
    struct tm *timeStruct;
    char buffer[50];
    
    time( &rawtime );
    timeStruct = localtime( &rawtime );
    sprintf(buffer, "%s_%02d%02d%02d-%02d%02d%02d%s", prefix, timeStruct->tm_year%100,
            timeStruct->tm_mon + 1, timeStruct->tm_mday,
            timeStruct->tm_hour, timeStruct->tm_min, timeStruct->tm_sec, postfix);
    strcpy(filenameBuffer, buffer);
    
}
                              
                                  
                              

