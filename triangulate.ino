/*
    Author: MaHir & ToniNii
    See also additional text later

    This script is used for wifi-positioning:
    First wifiAPs are scanned and the distances to
    those calculated. The coordinates of APs are known
    in advance. After three APs are been found,
    the position of scanner is calculated based on
    triangulate postioning method. Due to poor 
    quality of scan accurascy, error handling is made
    with "smallest enclosing circle" method.
	
	Code is written for ESP8266 module (Wemos D1 mini pro).
*/
#include <math.h>
#include <random>
#include "SmallestEnclosingCircle.h"
#include "ESP8266WiFi.h"

#define ARRAY_SIZE(someArray) (sizeof(someArray) / sizeof(someArray[0]))

// SSIDs and positions for all WiFi beacons:
int BeaconCount = 3;		// Amound of beacons
int BeaconInfoCount = 3;	// Amount of data variables for single beacon
struct Beacon {
	const char SSID;
	int x;
	int y;
	int r;
	} Beacons [BeaconCount];

Beacons[0].SSID = "SSID1";
Beacons[0].x = 10;
Beacons[0].y = 10;

Beacons[1].SSID = "SSID2";
Beacons[1].x = 25;
Beacons[1].y = 35;

Beacons[2].SSID = "SSID3";
Beacons[2].x = 32;
Beacons[2].y = 5;

//float Ax=10;
//float Ay=10;

//const char BeaconB = "SSID2";
//float Bx=25;
//float By=35;

//const char BeaconC = "SSID3";
//float Cx=32;
//float Cy=5;

//float Ar =16;
//float Br =16;
//float Cr =16; //Define these from rssi!

float error=0;
int nroOFintersections=0;

// TODO: dynamic intersection array
// now only for three towers (6 intersections)
Point intersectionsArray[6];

void setup() {
  Serial.begin(115200);
  
  // Set WiFi to station mode and disconnect from an AP if it was previously connected
  WiFi.mode(WIFI_STA);
  WiFi.disconnect();
  delay(100);
  
  Serial.println("Setup done");
}

void loop() {
  // STEP1: scan the wifi APs
  // System scans all networks and if three or more networks are found
  // it continues and adds the radius for all beacons.
  Serial.println("scan start");
  // WiFi.scanNetworks will return the number of networks found
  int n = WiFi.scanNetworks();
  int rssi=0;
  double to_exp=0;
  int = BcnCount = 0;
  Serial.println("scan done");
  if (n == 0) {
    Serial.println("no networks found");
  } 
  if (n < 3){Serial.println("Not enough beacons found, can not triangulate!")}
  else {
    Serial.print(n);
    Serial.println(" networks found");
    for (int i = 0; i < n; ++i) {
		for (int j = 0; j < BeaconCount; ++j){
			if (WiFi.SSID(i)==beacons[j].SSID){
				Serial.print(WiFi.SSID(i)+ " ,strength: "+ WiFi.RSSI(i) +"dBm");
				rssi=-WiFi.RSSI(i);
				to_exp=(27.55-(20*log10(2412))+rssi)/20;
				distance=100*pow(10,to_exp);
				beacons[j].r = distance;			//Insert radius of beacon to list of beacon info
				BcnCount++;
				Serial.print(" ,distance: ");
				Serial.print((float)(distance),1);
				Serial.println("cm");
			}
		}
    }
  }
  if (BcnCount < 3){Serial.println("Not enough beacons found, can not triangulate!");}
  
  Serial.print("scan done");
  Serial.print(BcnCount);
  Serial.print("towers found, radius for beacons added")

  // STEP2: calculate intersections
  // System calculates intersections of three
  // wifi APs based on previously saved array radius
  // TODO: add proper array for different tower positions
  // and radius
  Serial.println("starting to calculate intersections");
  int luku=-1;
  while (luku==-1){ // loops until funtion return something
    luku=intersections(Ax,Ay,Ar+error,Bx,By,Br+error); 
    // If circles are separate, funtion return 0
    // and error value is been increased
    if (luku==0){
      error=error+0.5;
    }
  }
  
  luku=-1;
  while (luku==-1){ // loops until funtion return something
    luku=intersections(Ax,Ay,Ar,Cx,Cy,Cr); 
    // If circles are separate, funtion return 0
    // and error value is been increased
    if (luku==0){
      error=error+0.5;
    }
  }
  luku=-1;
  while (luku==-1){ // loops until funtion return something
    luku=intersections(Bx,By,Br,Cx,Cy,Cr);  
    // If circles are separate, funtion return 0
    // and error value is been increased
    if (luku==0){
      error=error+0.5;
    }
  }
  Serial.println("intersections calculated"); 

  // DEBUG: prints all intersection coordinates
  long Count=ARRAY_SIZE(intersectionsArray);
  for (int i = 0; i < Count; i++) 
  {
    Point p = intersectionsArray[i];
    Serial.print("x: ");
    Serial.print(p.x);
    Serial.print(" y: ");
    Serial.println(p.y); 
  }

  // STEP3: find smallest enclosing circle
  // Goes through all combinations of possible 
  // intersectinons (which contain one from 
  // each pair of circles)
  Serial.println("smallest enclosing circle started to calculate");
  Circle circula; // the final result  
  Circle circulaTemp; // the temp result
  
  // TODO: effiecient way of making all combinations 
  // mentioned above (now just harcoded combo with three towers)  
  Point Parray[3]={intersectionsArray[0],intersectionsArray[2],intersectionsArray[4]};
  circula=makeSmallestEnclosingCircle(Parray,ARRAY_SIZE(Parray));
  // DEBUG: print result of first combo
  Serial.print("x: ");
  Serial.print( circula.c.x);
  Serial.print( " y: ");
  Serial.print(circula.c.y);
  Serial.print( " r: ");
  Serial.println(circula.r);
  
  Parray[0]=intersectionsArray[0];
  Parray[1]=intersectionsArray[2];
  Parray[2]=intersectionsArray[5];
  circulaTemp=makeSmallestEnclosingCircle(Parray, ARRAY_SIZE(Parray));
  if (circula.r>circulaTemp.r){
    circula=circulaTemp;
  }
  // DEBUG: print result of first combo
  Serial.print("x: ");
  Serial.print( circulaTemp.c.x);
  Serial.print( " y: ");
  Serial.print(circulaTemp.c.y);
  Serial.print( " r: ");
  Serial.println(circulaTemp.r);
  
  Parray[0]=intersectionsArray[0];
  Parray[1]=intersectionsArray[3];
  Parray[2]=intersectionsArray[4];
  circulaTemp=makeSmallestEnclosingCircle(Parray, ARRAY_SIZE(Parray));
  if (circula.r>circulaTemp.r){
    circula=circulaTemp;
  }
  // DEBUG: print result of first combo
  Serial.print("x: ");
  Serial.print( circulaTemp.c.x);
  Serial.print( " y: ");
  Serial.print(circulaTemp.c.y);
  Serial.print( " r: ");
  Serial.println(circulaTemp.r); 

  Parray[0]=intersectionsArray[0];
  Parray[1]=intersectionsArray[3];
  Parray[2]=intersectionsArray[5];
  circulaTemp=makeSmallestEnclosingCircle(Parray, ARRAY_SIZE(Parray));
  if (circula.r>circulaTemp.r){
    circula=circulaTemp;
  }
  // DEBUG: print result of first combo
  Serial.print("x: ");
  Serial.print( circulaTemp.c.x);
  Serial.print( " y: ");
  Serial.print(circulaTemp.c.y);
  Serial.print( " r: ");
  Serial.println(circulaTemp.r); 

  Parray[0]=intersectionsArray[1];
  Parray[1]=intersectionsArray[2];
  Parray[2]=intersectionsArray[4];
  circulaTemp=makeSmallestEnclosingCircle(Parray, ARRAY_SIZE(Parray));
  if (circula.r>circulaTemp.r){
    circula=circulaTemp;
  }
  // DEBUG: print result of first combo
  Serial.print("x: ");
  Serial.print( circulaTemp.c.x);
  Serial.print( " y: ");
  Serial.print(circulaTemp.c.y);
  Serial.print( " r: ");
  Serial.println(circulaTemp.r);

  Parray[0]=intersectionsArray[1];
  Parray[1]=intersectionsArray[2];
  Parray[2]=intersectionsArray[5];
  circulaTemp=makeSmallestEnclosingCircle(Parray, ARRAY_SIZE(Parray));
  if (circula.r>circulaTemp.r){
    circula=circulaTemp;
  }
  // DEBUG: print result of first combo
  Serial.print("x: ");
  Serial.print( circulaTemp.c.x);
  Serial.print( " y: ");
  Serial.print(circulaTemp.c.y);
  Serial.print( " r: ");
  Serial.println(circulaTemp.r);

  Parray[0]=intersectionsArray[1];
  Parray[1]=intersectionsArray[3];
  Parray[2]=intersectionsArray[4];
  circulaTemp=makeSmallestEnclosingCircle(Parray, ARRAY_SIZE(Parray));
  if (circula.r>circulaTemp.r){
    circula=circulaTemp;
  }
  // DEBUG: print result of first combo
  Serial.print("x: ");
  Serial.print( circulaTemp.c.x);
  Serial.print( " y: ");
  Serial.print(circulaTemp.c.y);
  Serial.print( " r: ");
  Serial.println(circulaTemp.r);

  Parray[0]=intersectionsArray[1];
  Parray[1]=intersectionsArray[3];
  Parray[2]=intersectionsArray[5];
  circulaTemp=makeSmallestEnclosingCircle(Parray, ARRAY_SIZE(Parray));
  if (circula.r>circulaTemp.r){
    circula=circulaTemp;
  }
  // DEBUG: print result of first combo
  Serial.print("x: ");
  Serial.print( circulaTemp.c.x);
  Serial.print( " y: ");
  Serial.print(circulaTemp.c.y);
  Serial.print( " r: ");
  Serial.println(circulaTemp.r);

  Serial.print("smallest enclosing circle found: x: ");
  Serial.print( circula.c.x);
  Serial.print( " y: ");
  Serial.print(circula.c.y);
  Serial.print( " r: ");
  Serial.println(circula.r);
  
}

// returns 1, if two intersections exist (saves those to intersectionsArray)
// returns 0, if circles are separate
// returns 2, if one circle is contained within the other
int intersections(float x0, float y0, float r0, float x1, float y1, float r1)
{
  double d;
  double a;
  double h;
  double x;
  double y;

   // distance between the center of the circles
   d=sqrt(pow((x1-x0),2)+pow((y1-y0),2));
   // circles are separate
   if ( d > r0 + r1 )
   {
      return 0;     
   }
   // one circle is contained within the other
   else if ( d < abs(r0 - r1) )
   {
      return 2;    
   }
   else
   { 
    a=(pow(r0,2)-pow(r1,2)+pow(d,2))/(2*d);
    h=sqrt(pow(r0,2)-pow(a,2));

    // coordinates of the first intersection 
    x=x0+a*(x1-x0)/d+h*(y1-y0)/d;
    y=y0+a*(y1-y0)/d-h*(x1-x0)/d;
    Point newPoint;
    newPoint.x=x;
    newPoint.y=y;
    intersectionsArray[nroOFintersections]=newPoint;
    nroOFintersections=nroOFintersections+1;

    // coordinates of the second intersection 
    x=x0+a*(x1-x0)/d-h*(y1-y0)/d;
    y=y0+a*(y1-y0)/d+h*(x1-x0)/d;
    newPoint.x=x;
    newPoint.y=y;
    intersectionsArray[nroOFintersections]=newPoint;
    nroOFintersections=nroOFintersections+1;
    return 1;
   }
}


/* This is modified version of Project Nayuki. See the text below 
/* 
 * Smallest enclosing circle - Library (C++)
 * 
 * Copyright (c) 2017 Project Nayuki
 * https://www.nayuki.io/page/smallest-enclosing-circle
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program (see COPYING.txt and COPYING.LESSER.txt).
 * If not, see <http://www.gnu.org/licenses/>.
 */

/*---- Members of struct Point ----*/

Point Point::subtract(const Point p) const {
  return Point{x - p.x, y - p.y};
}

double Point::distance(const Point p) const {
  return sqrt(pow((x - p.x),2)+pow(( y - p.y),2));
}

double Point::cross(const Point p) const {
  return x * p.y - y * p.x;
}


/*---- Members of struct Circle ----*/

const Circle Circle::INVALID{Point{0, 0}, -1};
const double Circle::MULTIPLICATIVE_EPSILON = 1 + 1e-14;


bool Circle::contains(const Point p) const{
  return c.distance(p) <= r * MULTIPLICATIVE_EPSILON;
}


bool Circle::contains(const Point pointArray[]) const {
  for (int i=0;i<sizeof(pointArray)/sizeof(Point);i++) {
    if (!contains(pointArray[i]))
      return false;
  }
  return true;
}


/*---- Smallest enclosing circle algorithm ----*/

//static Circle makeSmallestEnclosingCircleOnePoint(const Point pointArray[], int endi, const Point p);
//static Circle makeSmallestEnclosingCircleTwoPoints(const Point pointArray[], int endi, const Point p, const Point q);


// Initially: No boundary points known
Circle makeSmallestEnclosingCircle(const Point pointArray[], int countti) { 
  
  /*DEBUG
  for (int i = 0; i < countti; i++) 
  {
    Point p = pointArray[i];
    Serial.print("x: ");
    Serial.print(p.x);
    Serial.print(" y: ");
    Serial.println(p.y); 
  }
  */
  
  // Clone list to preserve the caller's data, randomize order
  Point shuffled[countti];
  memcpy(shuffled,pointArray,sizeof(shuffled));
  for (int i=0; i < countti; i++) 
  {
    int n = random(0, countti);  // Integer from 0 to questionCount-1
    Point temp = shuffled[n];
    shuffled[n] =  shuffled[i];
    shuffled[i] = temp;
  }
  
  // Progressively add points to circle or recompute circle
  Circle c(Circle::INVALID);
  for (int i = 0; i < countti; i++) {
    const Point p = shuffled[i];
    if (c.r < 0 || !c.contains(p))
      c = makeSmallestEnclosingCircleOnePoint(shuffled, i + 1, p);
  }
  return c;
}


// One boundary point known
Circle makeSmallestEnclosingCircleOnePoint(const Point pointArray[], int endi, const Point p) {
  Circle c{p, 0};
  for (int i = 0; i < endi; i++) {
    const Point q = pointArray[i];
    if (!c.contains(q)) {
      if (c.r == 0)
        //Serial.println("c.r==0");
        c = makeDiameter(p, q);
      else
        c = makeSmallestEnclosingCircleTwoPoints(pointArray, i + 1, p, q);
    }
  }
  return c;
}


// Two boundary points known
Circle makeSmallestEnclosingCircleTwoPoints(const Point pointArray[], int endi, const Point p, const Point q) {
  Circle circ = makeDiameter(p, q);
  Circle left = Circle::INVALID;
  Circle right = Circle::INVALID;
  
  // For each point not in the two-point circle
  Point pq = q.subtract(p);
  for (int i = 0; i < endi; i++) {
    const Point r = pointArray[i];
    if (circ.contains(r))
      continue;
    
    // Form a circumcircle and classify it on left or right side
    double cross = pq.cross(r.subtract(p));
    Circle c = makeCircumcircle(p, q, r);
    if (c.r < 0)
      continue;
    else if (cross > 0 && (left.r < 0 || pq.cross(c.c.subtract(p)) > pq.cross(left.c.subtract(p))))
      left = c;
    else if (cross < 0 && (right.r < 0 || pq.cross(c.c.subtract(p)) < pq.cross(right.c.subtract(p))))
      right = c;
  }
  
  // Select which circle to return
  if (left.r < 0 && right.r < 0)
    return circ;
  else if (left.r < 0)
    return right;
  else if (right.r < 0)
    return left;
  else
    return left.r <= right.r ? left : right;
}


Circle makeDiameter(const Point a, const Point b) {
  Point c{(a.x + b.x) / 2, (a.y + b.y) / 2};
  return Circle{c, max(c.distance(a), c.distance(b))};
}


Circle makeCircumcircle(const Point a, const Point b, const Point c) {
  // Mathematical algorithm from Wikipedia: Circumscribed circle
  double ox = (min(min(a.x, b.x), c.x) + max(min(a.x, b.x), c.x)) / 2;
  double oy = (min(min(a.y, b.y), c.y) + max(min(a.y, b.y), c.y)) / 2;
  double ax = a.x - ox, ay = a.y - oy;
  double bx = b.x - ox, by = b.y - oy;
  double cx = c.x - ox, cy = c.y - oy;
  double d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2;
  if (d == 0)
    return Circle::INVALID;
  double x = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d;
  double y = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d;
  Point p{ox + x, oy + y};
  double r = max(max(p.distance(a), p.distance(b)), p.distance(c));
  return Circle{p, r};
}

