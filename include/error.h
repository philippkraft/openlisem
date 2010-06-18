/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*
 * error.h, global errorstring for exception throwing
 *
 *  Created on: Feb 28, 2010
 *      Author: jetten
 */

#ifndef ERROR_H_
#define ERROR_H_

extern QString ErrorString;

#define Error(s,s1,s2)  ErrorString=QString("%1 %2 %2").arg(s).arg(s1).arg(s2)

#endif /* ERROR_H_ */
