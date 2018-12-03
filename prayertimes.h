/*-------------------------- In the name of God ----------------------------*\

    libprayertimes 1.0
    Islamic prayer times calculator library
    Based on PrayTimes 1.1 JavaScript library

----------------------------- Copyright Block --------------------------------

Copyright (C) 2007-2010 PrayTimes.org

Developed By: Mohammad Ebrahim Mohammadi Panah <ebrahim at mohammadi dot ir>
Based on a JavaScript Code By: Hamid Zarrabi-Zadeh
Ported from C++ By: Igor Pylypiv

License: GNU LGPL v3.0

TERMS OF USE:
    Permission is granted to use this code, with or
    without modification, in any website or application
    provided that credit is given to the original work
    with a link back to PrayTimes.org.

This program is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY.

PLEASE DO NOT REMOVE THIS COPYRIGHT BLOCK.

------------------------------------------------------------------------------

User's Manual:
http://praytimes.org/manual

Calculating Formulas:
http://praytimes.org/calculation

Code Repository:
http://code.ebrahim.ir/prayertimes/

\*--------------------------------------------------------------------------*/

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

enum
{
	VERSION_MAJOR = 1,
	VERSION_MINOR = 0,
};

// Calculation Methods
enum CalculationMethod
{
	Jafari, 	// Ithna Ashari
	Karachi,	// University of Islamic Sciences, Karachi
	ISNA,   	// Islamic Society of North America (ISNA)
	MWL,    	// Muslim World League (MWL)
	Makkah, 	// Umm al-Qura, Makkah
	Egypt,  	// Egyptian General Authority of Survey
	Custom, 	// Custom Setting

	CalculationMethodsCount
};

// Juristic Methods
enum JuristicMethod
{
	Shafii,    // Shafii (standard)
	Hanafi,    // Hanafi
};

// Adjusting Methods for Higher Latitudes
enum AdjustingMethod
{
	None,      	// No adjustment
	MidNight,  	// middle of night
	OneSeventh,	// 1/7th of night
	AngleBased,	// angle/60th of night
};

// Time IDs
enum TimeID
{
	Fajr,
	Sunrise,
	Dhuhr,
	Asr,
	Sunset,
	Maghrib,
	Isha,

	TimesCount
};



/* ---------------------- Misc Functions ----------------------- */

/* range reduce hours to 0..23 */
static double fix_hour(double a)
{
	a = a - 24.0 * floor(a / 24.0);
	a = a < 0.0 ? a + 24.0 : a;
	return a;
}

/* compute the difference between two times  */
static double time_diff(double time1, double time2)
{
	return fix_hour(time2 - time1);
}

/* ---------------------- Julian Date Functions ----------------------- */

/* calculate julian date from a calendar date */
double get_julian_date(int year, int month, int day)
{
	if (month <= 2)
	{
		year -= 1;
		month += 12;
	}

	double a = floor(year / 100.0);
	double b = 2 - a + floor(a / 4.0);

	return floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + b - 1524.5;
}

/* convert a calendar date to julian date (second method) */
double calc_julian_date(int year, int month, int day)
{
	double j1970 = 2440588.0;
	struct tm date = { 0 };
	date.tm_year = year - 1900;
	date.tm_mon = month - 1;
	date.tm_mday = day;
	date.tm_isdst = -1;		// determine it yourself from system
	time_t ms = mktime(&date);		// seconds since midnight Jan 1, 1970
	double days = floor(ms / (double) (60 * 60 * 24));
	return j1970 + days - 0.5;
}

/* ---------------------- Trigonometric Functions ----------------------- */

/* degree to radian */
static double deg2rad(double d)
{
	return d * M_PI / 180.0;
}

/* radian to degree */
static double rad2deg(double r)
{
	return r * 180.0 / M_PI;
}

/* degree sin */
static double dsin(double d)
{
	return sin(deg2rad(d));
}

/* degree cos */
static double dcos(double d)
{
	return cos(deg2rad(d));
}

/* degree tan */
static double dtan(double d)
{
	return tan(deg2rad(d));
}

/* degree arcsin */
static double darcsin(double x)
{
	return rad2deg(asin(x));
}

/* degree arccos */
static double darccos(double x)
{
	return rad2deg(acos(x));
}

/* degree arctan */
//static double darctan(double x)
//{
//	return rad2deg(atan(x));
//}

/* degree arctan2 */
static double darctan2(double y, double x)
{
	return rad2deg(atan2(y, x));
}

/* degree arccot */
static double darccot(double x)
{
	return rad2deg(atan(1.0 / x));
}

/* range reduce angle in degrees. */
static double fix_angle(double a)
{
	a = a - 360.0 * floor(a / 360.0);
	a = a < 0.0 ? a + 360.0 : a;
	return a;
}


/* get hours and minutes parts of a float time */
void get_float_time_parts(double time, int *hours, int *minutes)
{
	time = fix_hour(time + 0.5 / 60);		// add 0.5 minutes to round
	*hours = floor(time);
	*minutes = floor((time - *hours) * 60);
}

/* convert float hours to 24h format */
const char* float_time_to_time24(double time, char *buffer)
{
	if (isnan(time))
		return "";
	int hours, minutes;
	get_float_time_parts(time, &hours, &minutes);
	sprintf(buffer, "%2.2d:%2.2d", hours, minutes);
	return buffer;
}

/* convert float hours to 12h format */
const char* float_time_to_time12(double time, bool no_suffix, char *buffer)
{
	if (isnan(time))
		return "";
	int hours, minutes;
	get_float_time_parts(time, &hours, &minutes);
	const char* suffix = hours >= 12 ? "PM" : "AM";
	hours = (hours + 12 - 1) % 12 + 1;
	sprintf(buffer, "%d:%2.2d %s", hours, minutes, (no_suffix ? "" : suffix));
	return buffer;
}

/* convert float hours to 12h format with no suffix */
const char* float_time_to_time12ns(double time, char *buffer)
{
	return float_time_to_time12(time, true, buffer);
}

/* ---------------------- Time-Zone Functions ----------------------- */

/* compute local time-zone for a specific date */
double get_effective_timezone(time_t local_time)
{
	struct tm* tmp = localtime(&local_time);
	tmp->tm_isdst = 0;
	time_t local = mktime(tmp);
	tmp = gmtime(&local_time);
	tmp->tm_isdst = 0;
	time_t gmt = mktime(tmp);
	return (local - gmt) / 3600.0;
}

/* compute local time-zone for a specific date */
double get_effective_timezone_(int year, int month, int day)
{
	struct tm date = { 0 };
	date.tm_year = year - 1900;
	date.tm_mon = month - 1;
	date.tm_mday = day;
	date.tm_isdst = -1;		// determine it yourself from system
	time_t local = mktime(&date);		// seconds since midnight Jan 1, 1970
	return get_effective_timezone(local);
}


/* ------------------- Calc Method Parameters -------------------- */

struct MethodConfig
{
	double fajr_angle;
	bool   maghrib_is_minutes;
	double maghrib_value;		// angle or minutes
	bool   isha_is_minutes;
	double isha_value;		// angle or minutes
} method_params[CalculationMethodsCount] =
{
	{16.0, false, 4.0, false, 14.0},	// Jafari
	{18.0, true,  0.0, false, 18.0},	// Karachi
	{15.0, true,  0.0, false, 15.0},	// ISNA
	{18.0, true,  0.0, false, 17.0},	// MWL
	{19.0, true,  0.0, true,  90.0},	// Makkah
	{19.5, true,  0.0, false, 17.5},	// Egypt
	{18.0, true,  0.0, false, 17.0},	// Custom
};

enum CalculationMethod calc_method = Jafari;		// caculation method
enum JuristicMethod asr_juristic = Shafii;		// Juristic method for Asr
enum AdjustingMethod adjust_high_lats = MidNight;	// adjusting method for higher latitudes
double dhuhr_minutes = 0;			// minutes after mid-day for Dhuhr

double latitude;
double longitude;
double gtimezone;
double julian_date;

/* --------------------- Technical Settings -------------------- */

static const int NUM_ITERATIONS = 1;		// number of iterations needed to compute times

/* ---------------------- Calculation Functions ----------------------- */

/* References: */
/* http://www.ummah.net/astronomy/saltime   */
/* http://aa.usno.navy.mil/faq/docs/SunApprox.html */

struct sunpos {
	double declination;
	double equation; 
};

/* compute declination angle of sun and equation of time */
static struct sunpos sun_position(double jd)
{
	double d = jd - 2451545.0;
	double g = fix_angle(357.529 + 0.98560028 * d);
	double q = fix_angle(280.459 + 0.98564736 * d);
	double l = fix_angle(q + 1.915 * dsin(g) + 0.020 * dsin(2 * g));

	// double r = 1.00014 - 0.01671 * dcos(g) - 0.00014 * dcos(2 * g);
	double e = 23.439 - 0.00000036 * d;

	double dd = darcsin(dsin(e) * dsin(l));
	double ra = darctan2(dcos(e) * dsin(l), dcos(l)) / 15.0;
	ra = fix_hour(ra);
	double eq_t = q / 15.0 - ra;

	struct sunpos t;
	t.declination = dd;
	t.equation = eq_t;
	return t;
}

/* compute equation of time */
static double equation_of_time(double jd)
{
	return sun_position(jd).equation;
}

/* compute declination angle of sun */
static double sun_declination(double jd)
{
	return sun_position(jd).declination;
}

/* compute mid-day (Dhuhr, Zawal) time */
static double compute_mid_day(double _t)
{
	double t = equation_of_time(julian_date + _t);
	double z = fix_hour(12 - t);
	return z;
}

/* compute time for a given angle G */
static double compute_time(double g, double t)
{
	double d = sun_declination(julian_date + t);
	double z = compute_mid_day(t);
	double v = 1.0 / 15.0 * darccos((-dsin(g) - dsin(d) * dsin(latitude)) / (dcos(d) * dcos(latitude)));
	return z + (g > 90.0 ? - v :  v);
}

/* compute the time of Asr */
static double compute_asr(int step, double t)  // Shafii: step=1, Hanafi: step=2
{
	double d = sun_declination(julian_date + t);
	double g = -darccot(step + dtan(fabs(latitude - d)));
	return compute_time(g, t);
}

/* ----------------- Compute Prayer Times ----------------------- */

// array parameters must be at least of size TimesCount

/* convert hours to day portions  */
static void day_portion(double times[])
{
	for (int i = 0; i < TimesCount; ++i)
		times[i] /= 24.0;
}

/* compute prayer times at given julian date */
static void compute_times(double times[])
{
	day_portion(times);

	times[Fajr]    = compute_time(180.0 - method_params[calc_method].fajr_angle, times[Fajr]);
	times[Sunrise] = compute_time(180.0 - 0.833, times[Sunrise]);
	times[Dhuhr]   = compute_mid_day(times[Dhuhr]);
	times[Asr]     = compute_asr(1 + asr_juristic, times[Asr]);
	times[Sunset]  = compute_time(0.833, times[Sunset]);
	times[Maghrib] = compute_time(method_params[calc_method].maghrib_value, times[Maghrib]);
	times[Isha]    = compute_time(method_params[calc_method].isha_value, times[Isha]);
}

/* the night portion used for adjusting times in higher latitudes */
static double night_portion(double angle)
{
	switch (adjust_high_lats)
	{
		case AngleBased:
			return angle / 60.0;
		case MidNight:
			return 1.0 / 2.0;
		case OneSeventh:
			return 1.0 / 7.0;
		default:
			// Just to return something!
			// In original library nothing was returned
			// Maybe I should throw an exception
			// It must be impossible to reach here
			return 0;
	}
}


/* adjust Fajr, Isha and Maghrib for locations in higher latitudes */
static void adjust_high_lat_times(double times[])
{
	double night_time = time_diff(times[Sunset], times[Sunrise]);		// sunset to sunrise

	// Adjust Fajr
	double fajr_diff = night_portion(method_params[calc_method].fajr_angle) * night_time;
	if (isnan(times[Fajr]) || time_diff(times[Fajr], times[Sunrise]) > fajr_diff)
		times[Fajr] = times[Sunrise] - fajr_diff;

	// Adjust Isha
	double isha_angle = method_params[calc_method].isha_is_minutes ? 18.0 : method_params[calc_method].isha_value;
	double isha_diff = night_portion(isha_angle) * night_time;
	if (isnan(times[Isha]) || time_diff(times[Sunset], times[Isha]) > isha_diff)
		times[Isha] = times[Sunset] + isha_diff;

	// Adjust Maghrib
	double maghrib_angle = method_params[calc_method].maghrib_is_minutes ? 4.0 : method_params[calc_method].maghrib_value;
	double maghrib_diff = night_portion(maghrib_angle) * night_time;
	if (isnan(times[Maghrib]) || time_diff(times[Sunset], times[Maghrib]) > maghrib_diff)
		times[Maghrib] = times[Sunset] + maghrib_diff;
}



/* adjust times in a prayer time array */
static void adjust_times(double times[])
{
	for (int i = 0; i < TimesCount; ++i)
		times[i] += gtimezone - longitude / 15.0;
	times[Dhuhr] += dhuhr_minutes / 60.0;		// Dhuhr
	if (method_params[calc_method].maghrib_is_minutes)		// Maghrib
		times[Maghrib] = times[Sunset] + method_params[calc_method].maghrib_value / 60.0;
	if (method_params[calc_method].isha_is_minutes)		// Isha
		times[Isha] = times[Maghrib] + method_params[calc_method].isha_value / 60.0;

	if (adjust_high_lats != None)
		adjust_high_lat_times(times);
}


/* compute prayer times at given julian date */
static void compute_day_times(double times[])
{
	double default_times[] = { 5, 6, 12, 13, 18, 18, 18 };		// default times
	for (int i = 0; i < TimesCount; ++i)
		times[i] = default_times[i];

	for (int i = 0; i < NUM_ITERATIONS; ++i)
		compute_times(times);

	adjust_times(times);
}

//--------------

/* return prayer times for a given date */
void get_prayer_times_(int year, int month, int day, double _latitude, double _longitude, double _timezone, double times[])
{
	latitude = _latitude;
	longitude = _longitude;
	gtimezone = _timezone;
	julian_date = get_julian_date(year, month, day) - longitude / (double) (15 * 24);
	compute_day_times(times);
}

/* return prayer times for a given date */
void get_prayer_times(time_t date, double latitude, double longitude, double timezone, double times[])
{
	struct tm* t = localtime(&date);
	get_prayer_times_(1900 + t->tm_year, t->tm_mon + 1, t->tm_mday, latitude, longitude, timezone, times);
}

/* set the calculation method  */
void set_calc_method(enum CalculationMethod method_id)
{
	calc_method = method_id;
}

/* set the juristic method for Asr */
void set_asr_method(enum JuristicMethod method_id)
{
	asr_juristic = method_id;
}

/* set adjusting method for higher latitudes */
void set_high_lats_adjust_method(enum AdjustingMethod method_id)
{
	adjust_high_lats = method_id;
}

/* set the angle for calculating Fajr */
void set_fajr_angle(double angle)
{
	method_params[Custom].fajr_angle = angle;
	calc_method = Custom;
}

/* set the angle for calculating Maghrib */
void set_maghrib_angle(double angle)
{
	method_params[Custom].maghrib_is_minutes = false;
	method_params[Custom].maghrib_value = angle;
	calc_method = Custom;
}

/* set the angle for calculating Isha */
void set_isha_angle(double angle)
{
	method_params[Custom].isha_is_minutes = false;
	method_params[Custom].isha_value = angle;
	calc_method = Custom;
}

/* set the minutes after mid-day for calculating Dhuhr */
void set_dhuhr_minutes(double minutes)
{
	dhuhr_minutes = minutes;
}

/* set the minutes after Sunset for calculating Maghrib */
void set_maghrib_minutes(double minutes)
{
	method_params[Custom].maghrib_is_minutes = true;
	method_params[Custom].maghrib_value = minutes;
	calc_method = Custom;
}

/* set the minutes after Maghrib for calculating Isha */
void set_isha_minutes(double minutes)
{
	method_params[Custom].isha_is_minutes = true;
	method_params[Custom].isha_value = minutes;
	calc_method = Custom;
}
