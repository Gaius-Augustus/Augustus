/*	Functions that prints in seconds the elapsed time to 
	run a task
	
	Created: 22-September-2011	
	Last modified: 14-November-2011
*/

#include <time.h>
#include <iomanip>
#include <iostream> 
#include <stdio.h>

using namespace std;  

int printElapsedTime(int tEnd, int tStart)
{
  int tElapsed;
  tElapsed = difftime(tEnd, tStart);     
  
  printf("Elapsed time: %.3d seconds.\n", tElapsed); 
  // cout << "Elapsed time: " << setw(3) << tElapsed << " seconds.\n"; 

  return tElapsed;
}
