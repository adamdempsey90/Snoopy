/*	
	This file is part of the Snoopy code.

    Snoopy code is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snoopy code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
*/

double time_shift(double t);
void remap(double complex qi[]);
void kvolve(const double tremap);
#ifdef WITH_ELLIPTICAL_VORTEX
void kvolvellipse(const double t); //AJB
void knumvolvellipse(const double t, const double dt, const int kflag); //AJB
void f ( double t, double y[], double yp[] );
#endif
