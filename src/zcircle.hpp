#ifndef ZCIRCLE_HPP_
#define ZCIRCLE_HPP_

#include "dielectric.hpp"

  // Local functions used by Zcircle()
  void Set_to_zero(double *x,double *y) {
    *x=0.0; *y=0.0;
  }

  void Auxilliary(int m, int *i,int *j) { 
    *i=*i+1; *j=*j-1;
    if (*i>m) *i=*i-m;
    if (*j<1) *j=*j+m;
  }

  void Interpolation(int m, int i,int j,double *X,double *Y,
	double *U,double *V,double *xx,double *yy) {
    j=i-1;
    if (j<1) j=j+m;
    //Regula Falsi interpolation for the zero
    *xx=(X[i]*U[j]+X[j]*U[i])/(U[i]+U[j]);
    *yy=(Y[i]*V[j]+Y[j]*V[i])/(V[i]+V[j]);
  } // End local functions

/************************************************
*        Complex root search subroutine         *
* --------------------------------------------- *
* This routine searches for the complex roots   *
* of an analytical function by encircling the   *
* zero and estimating where it is. The circle   *
* is subsequently tightened by a factor e, and  *
* a new estimate made.                          *
* The inputs to the subroutine are;             *
*     w      initial radius of search circle    *
*     x0,y0  the initial guess                  *
*     e      factor by which circle is reduced  *
*     n      maximum number of iterations       *
*     m      evaluation points per quadrant     *
* The routine returns Z=X+IY (X,Y), and the     *
* number of iterations, k.                      *
************************************************/
void Zcircle(double w, double e, double x0, double y0, int n, int m, 
				double &xx, double &yy, dielectric &epsilon)  {
  //Labels: 50,100,200,250,300,310,320,330,340,350,400,450
  //        460,470,500,510,520,521
  int i,j,m3,m4;
  double X[80],Y[80],U[80],V[80];
  double a,b1,b2,PI,x1,x2,xm1,xm2;
  double x3,y3,x4,y1,y2,y4,uu,vv;
  int k;
  double yy0;
  m=m*4; k=1;
  PI = 4*atan(1);
  a=2*PI/m;
e50:for (i=1; i<m+1; i++) {
    X[i]=w*cos(a*(i-1))+x0;
    Y[i]=w*sin(a*(i-1))+yy0;
	}
  // Determine the corresponding U[i] and V[i]
  for (i=1; i<m+1; i++) {
    xx=X[i] ; yy=Y[i];
    //Eval(xx,yy,&uu,&vv);
	epsilon.determinant(xx,yy,&uu,&vv);
    U[i]=uu;
    V[i]=vv;
  }
  // Find the position at which uu changes sign
  // in the counterclockwise direction
  i=1; uu=U[i];
e100: Auxilliary(m,&i,&j);
  if (uu*U[i]<0) goto e200;
  // Guard against infinite loop
  if (i==1) Set_to_zero(&xx,&yy);
  goto e100;
  // Transition found
e200: xm1=i;
  // Search for the other transition, starting
  // on the other side of the circle
  i=(int) xm1+m/2;
  if (i>m) i=i-m;
  j=i;
  uu=U[i];
  // Flip directions alternately
e250: Auxilliary(m,&i,&j);
  if (uu*U[i]<0) goto e300;
  if (uu*U[j]<0) goto e310;
  // Guard against infinite loop
  if (i==xm1+(xm2/2)) Set_to_zero(&xx,&yy);
  if (j==xm1+(xm2/2)) Set_to_zero(&xx,&yy);
  goto e250;
  // Transition found
e300: m3=i;
  goto e320;
  // Transition found
e310: if (j==m) j=0;
  m3=j+1;
  // xm1 and m3 have been determined, now for xm2 and m4
  // now for the vv transitions                         
e320: i=(int) xm1+m/4;
  if (i>m) i=i-m;
  j=i; vv=V[i];
e330: Auxilliary(m,&i,&j);
  if (vv*V[i]<0) goto e340;
  if (vv*V[j]<0) goto e350;
  // Guard against infinite loop
  if (i==xm1+m/4) Set_to_zero(&xx,&yy);
  if (j==xm1+m/4) Set_to_zero(&xx,&yy);
  goto e330;
  // xm2 has been found
e340: xm2=i;
  goto e400;
e350: if (j==m) j=0;
  xm2=j+1;
  // xm2 has been found, now for m4
e400: i=(int) xm2+m/2;
  if (i>m) i=i-m;
  j=i; vv=V[i];
e450: Auxilliary(m,&i,&j);
  if (uu*V[i]<0) goto e460;
  if (vv*V[j]<0) goto e470;
  // Guard against infinite loop
  if (i==xm2+m/2) Set_to_zero(&xx,&yy);
  if (j==xm2+m/2) Set_to_zero(&xx,&yy);
  goto e450;
e460: m4=i;
  goto e500;
e470: if (j==m) j=0;
  m4=j+1;
  // All the intersections have been determine
  // Interpolate to find the four (x,y) coordinates
e500: i=(int) xm1;
  Interpolation(m,i,j,X,Y,U,V,&xx,&yy);
  x1=xx ; y1=yy ; i=(int) xm2;
  Interpolation(m,i,j,X,Y,U,V,&xx,&yy);
  x2=xx ; y2=yy ; i=m3;
  Interpolation(m,i,j,X,Y,U,V,&xx,&yy);
  x3=xx ; y3=yy ; i=m4;
  Interpolation(m,i,j,X,Y,U,V,&xx,&yy);
  x4=xx ; y4=yy;
  // Calculate the intersection of the lines
  // Guard against a divide by zero
  if (x1!=x3) goto e510;
  xx=x1; yy=(y1+y3)/2.0;
  goto e520;
e510: xm1=(y3-y1)/(x3-x1);
  if (x2!=x4) goto e520;
  xm2=1e8;
  goto e521;
e520: xm2=(y2-y4)/(x2-x4);
e521: b1=y1-xm1*x1;
  b2=y2-xm2*x2;
  xx=-(b1-b2)/(xm1-xm2);
  yy=(xm1*b2+xm2*b1)/(xm1+xm2);
  // is another iteration in order ?
  if (k==n) return;
  x0=xx ; yy0=yy ; k=k+1 ; w=w*e;
  goto e50;
} // Zcircle

#endif
