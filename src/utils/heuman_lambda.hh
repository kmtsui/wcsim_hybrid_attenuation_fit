#ifndef HEUMAN_LAMBDA_HH
#define HEUMAN_LAMBDA_HH

# include <cmath>
# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <algorithm>

constexpr double PI = 3.14159265358979323846;

// TOMS577, a C++ code which evaluates Carlson's elliptic integral functions RC, RD, RF and RJ. 
// Copied from https://people.sc.fsu.edu/~jburkardt/cpp_src/toms577/toms577.html

//****************************************************************************80

double rc ( double x, double y, double errtol, int &ierr )

//****************************************************************************80
//
//  Purpose:
//
//    RC computes the elementary integral RC(X,Y).
//
//   Discussion:
//
//     This function computes the elementary integral
//
//       RC(X,Y) = Integral ( 0 <= T < oo )
//
//                   -1/2     -1
//         (1/2)(T+X)    (T+Y)  DT,
//
//     where X is nonnegative and Y is positive.  The duplication
//     theorem is iterated until the variables are nearly equal,
//     and the function is then expanded in Taylor series to fifth
//     order.
//
//     Logarithmic, inverse circular, and inverse hyperbolic
//     functions can be expressed in terms of RC.
//
//     Check by addition theorem:
//
//       RC(X,X+Z) + RC(Y,Y+Z) = RC(0,Z),
//       where X, Y, and Z are positive and X * Y = Z * Z.
//
//   Licensing:
//
//     This code is distributed under the GNU LGPL license.
//
//   Modified:
//
//     02 June 2018
//
//   Author:
//
//     Original FORTRAN77 version by Bille Carlson, Elaine Notis.
//     This C++ version by John Burkardt.
//
//   Reference:
//
//     Bille Carlson,
//     Computing Elliptic Integrals by Duplication,
//     Numerische Mathematik,
//     Volume 33, 1979, pages 1-16.
//
//     Bille Carlson, Elaine Notis,
//     Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
//     ACM Transactions on Mathematical Software,
//     Volume 7, Number 3, pages 398-403, September 1981.
//
//   Parameters:
//
//     Input, double X, Y, the arguments in the integral.
//
//     Input, double ERRTOL, the error tolerance.
//     Relative error due to truncation is less than
//       16 * ERRTOL ^ 6 / (1 - 2 * ERRTOL).
//     Sample choices:
//       ERRTOL   Relative truncation error less than
//       1.D-3    2.D-17
//       3.D-3    2.D-14
//       1.D-2    2.D-11
//       3.D-2    2.D-8
//       1.D-1    2.D-5
//
//     Output, int &IERR, the error flag.
//     0, no error occurred.
//     1, abnormal termination.
//
{
  double c1;
  double c2;
  double lamda;
  const double lolim = 3.0E-78;
  double mu;
  double s;
  double sn;
  const double uplim = 1.0E+75;
  double value;
  double xn;
  double yn;
//
//   LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
//   LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
//   UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.
//
  if ( 
    x < 0.0 || 
    y <= 0.0 || 
    ( x + y ) < lolim || 
    uplim < x || 
    uplim < y )
  {
    std::cout << "\n";
    std::cout << "RC - Error!\n";
    std::cout << "  Invalid input arguments.\n";
    std::cout << "  X = " << x << "\n";
    std::cout << "  Y = " << y << "\n";
    std::cout << "\n";
    ierr = 1;
    value = 0.0;
    return value;
  }

  ierr = 0;
  xn = x;
  yn = y;

  while ( true )
  {
    mu = ( xn + yn + yn ) / 3.0;
    sn = ( yn + mu ) / mu - 2.0;

    if ( fabs ( sn ) < errtol )
    {
      c1 = 1.0 / 7.0;
      c2 = 9.0 / 22.0;
      s = sn * sn * ( 0.3 + sn * ( c1 + sn * ( 0.375 + sn * c2 ) ) );
      value = ( 1.0 + s ) / sqrt ( mu );
      return value;
    }

    lamda = 2.0 * sqrt ( xn ) * sqrt ( yn ) + yn;
    xn = ( xn + lamda ) * 0.25;
    yn = ( yn + lamda ) * 0.25;
  }

}
//****************************************************************************80

double rd ( double x, double y, double z, double errtol, int &ierr )

//****************************************************************************80
//
//  Purpose:
//
//    RD computes an incomplete elliptic integral of the second kind, RD(X,Y,Z).
//
//   Discussion:
//
//     This function computes an incomplete elliptic integral of the second kind.
//
//     RD(X,Y,Z) = Integral ( 0 <= T < oo )
//
//                     -1/2     -1/2     -3/2
//           (3/2)(T+X)    (T+Y)    (T+Z)    DT,
//
//     where X and Y are nonnegative, X + Y is positive, and Z is positive.
//
//     If X or Y is zero, the integral is complete.
//
//     The duplication theorem is iterated until the variables are
//     nearly equal, and the function is then expanded in Taylor
//     series to fifth order.
//
//     Check:
//
//       RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y) = 3 / sqrt ( X * Y * Z ),
//       where X, Y, and Z are positive.
//
//   Licensing:
//
//     This code is distributed under the GNU LGPL license.
//
//   Modified:
//
//     02 June 2018
//
//   Author:
//
//     Original FORTRAN77 version by Bille Carlson, Elaine Notis.
//     This C++ version by John Burkardt.
//
//   Reference:
//
//     Bille Carlson,
//     Computing Elliptic Integrals by Duplication,
//     Numerische Mathematik,
//     Volume 33, 1979, pages 1-16.
//
//     Bille Carlson, Elaine Notis,
//     Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
//     ACM Transactions on Mathematical Software,
//     Volume 7, Number 3, pages 398-403, September 1981.
//
//   Parameters:
//
//     Input, double X, Y, Z, the arguments in the integral.
//
//     Input, double ERRTOL, the error tolerance.
//     The relative error due to truncation is less than
//       3 * ERRTOL ^ 6 / (1-ERRTOL) ^ 3/2.
//     Sample choices:
//       ERRTOL   Relative truncation error less than
//       1.D-3    4.D-18
//       3.D-3    3.D-15
//       1.D-2    4.D-12
//       3.D-2    3.D-9
//       1.D-1    4.D-6
//
//     Output, int &IERR, the error flag.
//     0, no error occurred.
//     1, abnormal termination.
//
{
  double c1;
  double c2;
  double c3;
  double c4;
  double ea;
  double eb;
  double ec;
  double ed;
  double ef;
  double epslon;
  double lamda;
  const double lolim = 6.0E-51;
  double mu;
  double power4;
  double sigma;
  double s1;
  double s2;
  const double uplim = 1.0E+48;
  double value;
  double xn;
  double xndev;
  double xnroot;
  double yn;
  double yndev;
  double ynroot;
  double zn;
  double zndev;
  double znroot;
//
//   LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
//   LOLIM IS NOT LESS THAN 2 / (MACHINE MAXIMUM) ^ (2/3).
//   UPLIM IS NOT GREATER THAN (0.1 * ERRTOL / MACHINE
//   MINIMUM) ^ (2/3), WHERE ERRTOL IS DESCRIBED BELOW.
//   IN THE FOLLOWING TABLE IT IS ASSUMED THAT ERRTOL WILL
//   NEVER BE CHOSEN SMALLER THAN 1.D-5.
//
  if ( 
    x < 0.0 || 
    y < 0.0 || 
    x + y < lolim || 
    z < lolim || 
    uplim < x || 
    uplim < y || 
    uplim < z )
  {
    std::cout << "\n";
    std::cout << "RD - Error!\n";
    std::cout << "  Invalid input arguments.\n";
    std::cout << "  X = " << x << "\n";
    std::cout << "  Y = " << y << "\n";
    std::cout << "  Z = " << z << "\n";
    std::cout << "\n";
    ierr = 1;
    value = 0.0;
    return value;
  }

  ierr = 0;
  xn = x;
  yn = y;
  zn = z;
  sigma = 0.0;
  power4 = 1.0;

  while ( true )
  {
    mu = ( xn + yn + 3.0 * zn ) * 0.2;
    xndev = ( mu - xn ) / mu;
    yndev = ( mu - yn ) / mu;
    zndev = ( mu - zn ) / mu;
    epslon = std::max( fabs ( xndev ), 
      std::max( fabs ( yndev ), fabs ( zndev ) ) );

    if ( epslon < errtol )
    {
      c1 = 3.0 / 14.0;
      c2 = 1.0 / 6.0;
      c3 = 9.0 / 22.0;
      c4 = 3.0 / 26.0;
      ea = xndev * yndev;
      eb = zndev * zndev;
      ec = ea - eb;
      ed = ea - 6.0 * eb;
      ef = ed + ec + ec;
      s1 = ed * ( - c1 + 0.25 * c3 * ed - 1.5 * c4 * zndev * ef );
      s2 = zndev  * ( c2 * ef + zndev * ( - c3 * ec + zndev * c4 * ea ) );
      value = 3.0 * sigma  + power4 * ( 1.0 + s1 + s2 ) / ( mu * sqrt ( mu ) );
      return value;
    }

    xnroot = sqrt ( xn );
    ynroot = sqrt ( yn );
    znroot = sqrt ( zn );
    lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot;
    sigma = sigma + power4 / ( znroot * ( zn + lamda ) );
    power4 = power4 * 0.25;
    xn = ( xn + lamda ) * 0.25;
    yn = ( yn + lamda ) * 0.25;
    zn = ( zn + lamda ) * 0.25;
  }

}
//****************************************************************************80

double rf ( double x, double y, double z, double errtol, int &ierr )

//****************************************************************************80
//
//  Purpose:
//
//    RF computes an incomplete elliptic integral of the first kind, RF(X,Y,Z).
//
//   Discussion:
//
//     This function computes the incomplete elliptic integral of the first kind.
//
//     RF(X,Y,Z) = Integral ( 0 <= T < oo )
//
//                     -1/2     -1/2     -1/2
//           (1/2)(T+X)    (T+Y)    (T+Z)    DT,
//
//     where X, Y, and Z are nonnegative and at most one of them is zero.
//
//     If X or Y or Z is zero, the integral is complete.
//
//     The duplication theorem is iterated until the variables are
//     nearly equal, and the function is then expanded in Taylor
//     series to fifth order.
//
//     Check by addition theorem:
//
//       RF(X,X+Z,X+W) + RF(Y,Y+Z,Y+W) = RF(0,Z,W),
//       where X, Y, Z, W are positive and X * Y = Z * W.
//
//   Licensing:
//
//     This code is distributed under the GNU LGPL license.
//
//   Modified:
//
//     02 June 2018
//
//   Author:
//
//     Original FORTRAN77 version by Bille Carlson, Elaine Notis.
//     This C++ version by John Burkardt.
//
//   Reference:
//
//     Bille Carlson,
//     Computing Elliptic Integrals by Duplication,
//     Numerische Mathematik,
//     Volume 33, 1979, pages 1-16.
//
//     Bille Carlson, Elaine Notis,
//     Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
//     ACM Transactions on Mathematical Software,
//     Volume 7, Number 3, pages 398-403, September 1981.
//
//   Parameters:
//
//     Input, double X, Y, Z, the arguments in the integral.
//
//     Input, double ERRTOL, the error tolerance.
//     Relative error due to truncation is less than
//       ERRTOL ^ 6 / (4 * (1 - ERRTOL)).
//     Sample choices:
//       ERRTOL   Relative truncation error less than
//       1.D-3    3.D-19
//       3.D-3    2.D-16
//       1.D-2    3.D-13
//       3.D-2    2.D-10
//       1.D-1    3.D-7
//
//     Output, int &IERR, the error flag.
//     0, no error occurred.
//     1, abnormal termination.
//
{
  double c1;
  double c2;
  double c3;
  double e2;
  double e3;
  double epslon;
  double lamda;
  const double lolim = 3.0E-78;
  double mu;
  double s;
  const double uplim = 1.0E+75;
  double value;
  double xn;
  double xndev;
  double xnroot;
  double yn;
  double yndev;
  double ynroot;
  double zn;
  double zndev;
  double znroot;
//
//   LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
//   LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
//   UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.
//
  if ( 
    x < 0.0 || 
    y < 0.0 || 
    z < 0.0 || 
    x + y < lolim || 
    x + z < lolim || 
    y + z < lolim || 
    uplim <= x || 
    uplim <= y || 
    uplim <= z )
  {
    std::cout << "\n";
    std::cout << "RF - Error!\n";
    std::cout << "  Invalid input arguments.\n";
    std::cout << "  X = " << x << "\n";
    std::cout << "  Y = " << y << "\n";
    std::cout << "  Z = " << z << "\n";
    std::cout << "\n";
    ierr = 1;
    value = 0.0;
    return value;
  }

  ierr = 0;
  xn = x;
  yn = y;
  zn = z;

  while ( true )
  {
    mu = ( xn + yn + zn ) / 3.0;
    xndev = 2.0 - ( mu + xn ) / mu;
    yndev = 2.0 - ( mu + yn ) / mu;
    zndev = 2.0 - ( mu + zn ) / mu;
    epslon = std::max( fabs ( xndev ), 
      std::max( fabs ( yndev ), fabs ( zndev ) ) );

    if ( epslon < errtol )
    {
      c1 = 1.0 / 24.0;
      c2 = 3.0 / 44.0;
      c3 = 1.0 / 14.0;
      e2 = xndev * yndev - zndev * zndev;
      e3 = xndev * yndev * zndev;
      s = 1.0 + ( c1 * e2 - 0.1 - c2 * e3 ) * e2 + c3 * e3;
      value = s / sqrt ( mu );
      return value;
    }

    xnroot = sqrt ( xn );
    ynroot = sqrt ( yn );
    znroot = sqrt ( zn );
    lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot;
    xn = ( xn + lamda ) * 0.25;
    yn = ( yn + lamda ) * 0.25;
    zn = ( zn + lamda ) * 0.25;
  }

}
//****************************************************************************80

double rj ( double x, double y, double z, double p, double errtol, int &ierr )

//****************************************************************************80
//
//  Purpose:
//
//    RJ computes an incomplete elliptic integral of the third kind, RJ(X,Y,Z,P).
//
//   Discussion:
//
//     This function computes an incomplete elliptic integral of the third kind.
//
//     RJ(X,Y,Z,P) = Integral ( 0 <= T < oo )
//
//                   -1/2     -1/2     -1/2     -1
//         (3/2)(T+X)    (T+Y)    (T+Z)    (T+P)  DT,
//
//     where X, Y, and Z are nonnegative, at most one of them is
//     zero, and P is positive.
//
//     If X or Y or Z is zero, then the integral is complete.
//
//     The duplication theorem is iterated until the variables are nearly equal,
//     and the function is then expanded in Taylor series to fifth order.
//
//     Check by addition theorem:
//
//       RJ(X,X+Z,X+W,X+P)
//       + RJ(Y,Y+Z,Y+W,Y+P) + (A-B) * RJ(A,B,B,A) + 3 / sqrt ( A)
//       = RJ(0,Z,W,P), where X,Y,Z,W,P are positive and X * Y
//       = Z * W,  A = P * P * (X+Y+Z+W),  B = P * (P+X) * (P+Y),
//       and B - A = P * (P-Z) * (P-W).
//
//     The sum of the third and fourth terms on the left side is 3 * RC(A,B).
//
//   Licensing:
//
//     This code is distributed under the GNU LGPL license.
//
//   Modified:
//
//     02 June 2018
//
//   Author:
//
//     Original FORTRAN77 version by Bille Carlson, Elaine Notis.
//     This C++ version by John Burkardt.
//
//   Reference:
//
//     Bille Carlson,
//     Computing Elliptic Integrals by Duplication,
//     Numerische Mathematik,
//     Volume 33, 1979, pages 1-16.
//
//     Bille Carlson, Elaine Notis,
//     Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
//     ACM Transactions on Mathematical Software,
//     Volume 7, Number 3, pages 398-403, September 1981.
//
//   Parameters:
//
//     Input, double X, Y, Z, P, the arguments in the integral.
//
//     Input, double ERRTOL, the error tolerance.
//     Relative error due to truncation of the series for rj
//     is less than 3 * ERRTOL ^ 6 / (1 - ERRTOL) ^ 3/2.
//     An error tolerance (ETOLRC) will be passed to the subroutine
//     for RC to make the truncation error for RC less than for RJ.
//     Sample choices:
//       ERRTOL   Relative truncation error less than
//       1.D-3    4.D-18
//       3.D-3    3.D-15
//       1.D-2    4.D-12
//       3.D-2    3.D-9
//       1.D-1    4.D-6
//
//     Output, int &IERR, the error flag.
//     0, no error occurred.
//     1, abnormal termination.
//
{
  double alfa;
  double beta;
  double c1;
  double c2;
  double c3;
  double c4;
  double ea;
  double eb;
  double ec;
  double e2;
  double e3;
  double epslon;
  double etolrc;
  double lamda;
  const double lolim = 2.0E-26;
  double mu;
  double pn;
  double pndev;
  double power4;
  double sigma;
  double s1;
  double s2;
  double s3;
  const double uplim = 3.0E+24;
  double value;
  double xn;
  double xndev;
  double xnroot;
  double yn;
  double yndev;
  double ynroot;
  double zn;
  double zndev;
  double znroot;
//
//   LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
//   LOLIM IS NOT LESS THAN THE CUBE ROOT OF THE VALUE
//   OF LOLIM USED IN THE SUBROUTINE FOR RC.
//   UPLIM IS NOT GREATER THAN 0.3 TIMES THE CUBE ROOT OF
//   THE VALUE OF UPLIM USED IN THE SUBROUTINE FOR RC.
//
  if ( 
    x < 0.0 || 
    y < 0.0 || 
    z < 0.0 || 
    x + y < lolim || 
    x + z < lolim || 
    y + z < lolim || 
    p < lolim || 
    uplim < x || 
    uplim < y || 
    uplim < z || 
    uplim < p )
  {
    std::cout << "\n";
    std::cout << "RJ - Error!\n";
    std::cout << "  Invalid input arguments.\n";
    std::cout << "  X = " << x << "\n";
    std::cout << "  Y = " << y << "\n";
    std::cout << "  Z = " << z << "\n";
    std::cout << "  P = " << p << "\n";
    std::cout << "\n";
    ierr = 1;
    value = 0.0;
    return value;
  }

  ierr = 0;
  xn = x;
  yn = y;
  zn = z;
  pn = p;
  sigma = 0.0;
  power4 = 1.0;
  etolrc = 0.5 * errtol;

  while ( true )
  {
    mu = ( xn + yn + zn + pn + pn ) * 0.2;
    xndev = ( mu - xn ) / mu;
    yndev = ( mu - yn ) / mu;
    zndev = ( mu - zn ) / mu;
    pndev = ( mu - pn ) / mu;
    epslon = std::max( fabs ( xndev ), 
      std::max( fabs ( yndev ), 
      std::max( fabs ( zndev ), fabs ( pndev ) ) ) );

    if ( epslon < errtol )
    {
      c1 = 3.0 / 14.0;
      c2 = 1.0 / 3.0;
      c3 = 3.0 / 22.0;
      c4 = 3.0 / 26.0;
      ea = xndev * ( yndev + zndev ) + yndev * zndev;
      eb = xndev * yndev * zndev;
      ec = pndev * pndev;
      e2 = ea - 3.0 * ec;
      e3 = eb + 2.0 * pndev * ( ea - ec );
      s1 = 1.0 + e2 * ( - c1 + 0.75 * c3 * e2 - 1.5 * c4 * e3 );
      s2 = eb * ( 0.5 * c2 + pndev * ( - c3 - c3 + pndev * c4 ) );
      s3 = pndev * ea * ( c2 - pndev * c3 ) - c2 * pndev * ec;
      value = 3.0 * sigma + power4 * ( s1 + s2 + s3 ) / ( mu * sqrt ( mu ) );
      return value;
    }

    xnroot = sqrt ( xn );
    ynroot = sqrt ( yn );
    znroot = sqrt ( zn );
    lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot;
    alfa = pn * ( xnroot + ynroot + znroot ) + xnroot * ynroot * znroot;
    alfa = alfa * alfa;
    beta = pn * ( pn + lamda ) * ( pn + lamda );
    sigma = sigma + power4 * rc ( alfa, beta, etolrc, ierr );

    if ( ierr != 0 )
    {
      value = 0.0;
      return value;
    }

    power4 = power4 * 0.25;
    xn = ( xn + lamda ) * 0.25;
    yn = ( yn + lamda ) * 0.25;
    zn = ( zn + lamda ) * 0.25;
    pn = ( pn + lamda ) * 0.25;
  }

}

//****************************************************************************

// Evaluates the Heuman Lambda Function
// Only valid for fabs(phi)<Pi/2
// Related to Carlson's integrals as explained in https://www.boost.org/doc/libs/1_62_0/libs/math/doc/html/math_toolkit/ellint/heuman_lambda.html
double heuman_lambda(double phi, double k2)
{
    double delta2 = 1-(1-k2)*sin(phi)*sin(phi);
    double delta = sqrt(1-(1-k2)*sin(phi)*sin(phi));

    double errtol = 1.e-1;
    int ierr;
    double rfval = rf(0,1-k2,1,errtol,ierr);
    double rjval = rj(0,1-k2,1,1-k2/delta2,errtol,ierr);

    return 2./PI*(1.-k2)*sin(phi)*cos(phi)/delta*(rfval+k2/3./delta2*rjval); 
}

#endif // HEUMAN_LAMBDA_HH