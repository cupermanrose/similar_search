#ifndef _MYTIME_H
#define _MYTIME_H

#include <iostream>
#include <time.h>
#include <algorithm>
#include <string>

clock_t Start_Time, End_Time;

void init_time() {
	srand((int)time(0));
	Start_Time = clock();
}

void out_time(std::string time_description) {
	End_Time = clock();
	fout << time_description << " Total Time : " << (double)(End_Time - Start_Time) / CLOCKS_PER_SEC << "s" << std::endl;
	fout << std::endl;
}

#endif
