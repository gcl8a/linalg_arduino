/*
 * 2018.07.05: "Hard-coded" L2 norm, since it's so common
 * 2014.12.20: Converted to Arduino for library. All lowers now equal zero.
 * 2011-06-21: operator >> now indexes using operator [], and does not read the vector size
 * [this last part is a hack that needs to be undone]
 * TBD: Average() and StDev() functions
 * 2011-05-23: added statistical calculations Max(), Min(), and Range() [nb, Max and Min don't use fabs]
 * 2011-05-01: added Zero() function that I thought I already had
 * 2011-04-28: added Last() to vector
10/10/04: added Norm funcitons. CHANGED CalcMax to return infinity norm
09/15/04: added Average function
09/15/04: added Correl function, which emulates Excel CORREL (Covariance/<PI>Variances)
09/09/04: added CalcMax function [uses fabs]
06/16/03: removed XVectorError in favor of XError*/
//01/09/03: changed NO_VAILDATEs to VALIDATE
//11/11/02: fixed some ANSI incompatibilities
//1/22/02: made -1 default upper for null vector
//to avoid error when copying null vector
//12/04/01: added transpose

/*My very own vector template. Incorporates vector addition, etc.
Allows for vectors that start at something other than zero
Takes care of the memory management

Note that at the end, I define a few frequently used classes:

vector --> vector of doubles
fvector --> vector of floats
///////ivector --> vector of ints
*/

#ifndef __LINALG_VECTOR_H
#define __LINALG_VECTOR_H

//using namespace std;
/*
template < class T >
T MIN( T t1, T t2 ) { return t1 < t2 ? t1 : t2; }

template < class T >
T MAX( T t1, T t2 )
   { return t1 > t2 ? t1 : t2; }
*/

//see matrix.h


/*
Define a vector that runs from indices _lower->_upper, inclusive.
_tData is a pointer to the block allocated in memory, but
_tIndex is offset so that it can be accessed properly
*/


template < class T >
class TVector
   {
private:
   int _upper;
   T * _tData;
   T * _tIndex;

public:
   TVector( void ) : _upper( -1 ) { _tData = _tIndex = 0; }
   TVector( int length );
   TVector( T * data, int length );

   //copy constructor and assignment routines
   TVector( const TVector < T > & );

   TVector < T > & operator = ( const TVector & );
   TVector < T > & operator = ( const T& );

   //destructor
   ~TVector( void )
      {
      if ( _tData ) delete[] _tData;
      _tData = 0;
      }

   //access to an entry is by reference to allow assignment as an lvalue
   T & operator[] ( int i ) const
      {
#ifdef __VALIDATE__
      if ( i < 0 || i > _upper )
         throw XError( "Array out-of-bounds" );
#endif
      return _tIndex[i];
      }

   unsigned int Length( void )const { return _upper + 1; }
   T Last(void) const {return (*this)[Length()-1];}
   void Zero(void)
   {
	   if(_tData)
	   {
		   for ( int i = 0; i <= _upper; i++ ) _tIndex[i] = 0;
	   }
   }

   //checks to see if two vectors are compatible for comparison, addition, whatever
   int IsCompatible( const TVector < T > & t ) const
      { return ( _upper == t._upper ); }

   //comparison operators
   int operator == ( const TVector < T > & ) const;
   //overload an operator for "nearly" equals
   int operator %= ( const TVector < T > & ) const;

   TVector < T > operator + ( const TVector < T > & ) const;
   TVector < T > operator - ( const TVector < T > & ) const;
   TVector < T > & operator += ( const TVector < T > & );
   TVector < T > & operator -= ( const TVector < T > & );

   TVector < T > operator * ( const T & ) const;
   TVector < T > operator / ( const T & ) const;
   TVector < T > & operator *= ( const T & );
   TVector < T > & operator /= ( const T & );

   TVector < T > DotMultiply( const TVector < T > & ) const;
   TVector < T > DotDivide( const TVector < T > & ) const;

   //TVector < T > operator / ( const TMatrix < T > & tM ) const;
   //{TVector x(Length()); tM.Solve(x, *this); return x;}

   T Dot( const TVector < T > & ) const;
   T Average( void )const;
   T Correl( const TVector < T > & ) const;

   T Range (void) const {return FindMax()-FindMin();}
   T FindMin( void )const;
   int IndexOfMax( void )const;
   T FindMax( void )const;

   T AbsMax( void )const;
   T CalcNorm( int )const;

   T CalcL1Norm( void )const { return CalcNorm( 1 ); } //temp -- make faster later
       T CalcL2Norm( void ) const;// { return CalcNorm( 2 ); }

   //ostream& Write (ostream&) const;
   //istream& Read (istream&);
 
   };

#define dvector TVector<double>
#define fvector TVector<float>
//#define ivector TVector<int>

template < class T >
TVector < T >::TVector( int length )
: _upper( length - 1 )
   {
   _tData = new T[length];
   _tIndex = _tData;
   for ( int i = 0; i <= _upper; i++ ) _tIndex[i] = 0;
   }

//create a vector from an existing array, data.
//memory management of data is taken over by this object

template < class T >
TVector < T >::TVector( T * data, int length )
: _upper( length - 1 )
   {
   _tData = 0;
   _tIndex = data;
   }

//creates a copy

template < class T >
TVector < T >::TVector( const TVector < T > & t )
: _upper( t._upper )
   {
   _tData = new T[_upper + 1];
   _tIndex = _tData;

   for ( int i = 0; i <= _upper; i++ ) _tIndex[i] = t._tIndex[i];
   }

//assignment operator

template < class T >
TVector < T > & TVector < T >::operator = ( const TVector < T > & tV )
   {
   _upper = tV._upper;

   //if we have memory already allocated, delete it
   if ( _tData ) delete[] _tData;

   //now start anew
   _tData = new T[_upper + 1];
   _tIndex = _tData;

   for ( int i = 0; i <= _upper; i++ ) _tIndex[i] = tV._tIndex[i];
   return * this;
   }

template < class T >
TVector < T > & TVector < T >::operator = ( const T & t )
{
	for ( int i = 0; i <= _upper; i++ )
			_tData[i] = t;

	return * this;
}


template < class T >
int TVector < T >::operator == ( const TVector < T > & tV ) const
   {
#ifdef __VALIDATE__
   if ( !IsCompatible( tV ) ) return 0;
#endif

   for ( int i = 0; i <= _upper; i++ )
      if ( ( * this ) [i] != tV[i] ) return 0;
   return 1;
   }

#define ALLOWABLE_ERROR 1E-4

//checks for "nearly" equal, within factor of ALLOWABLE_ERROR
//scaled by the infinity norm of this vector

template < class T >
TVector < T > TVector < T >::operator + ( const TVector < T > & tV ) const
   {
   TVector < T > tSum( * this );
   return tSum += tV;
   }

template < class T >
TVector < T > TVector < T >::operator - ( const TVector < T > & tV ) const
   {
#ifdef __VALIDATE__
   if ( !IsCompatible( tV ) )
       throw XError( "Invalid vector addition!" );
#endif

   TVector < T > tDiff( Length() );
   for ( int i = 0; i <= _upper; i++ ) tDiff[i] = _tIndex[i] - tV[i];
   return tDiff;
   }

template < class T >
TVector < T > & TVector < T >::operator += ( const TVector < T > & tV )
   {
#ifdef __VALIDATE__
   if ( !IsCompatible( tV ) )
       throw XError( "Incompatible vector addition!" );
#endif

   for ( int i = 0; i <= _upper; i++ ) _tIndex[i] += tV[i];
   return * this;
   }

template < class T >
TVector < T > & TVector < T >::operator -= ( const TVector < T > & tV )
   {
#ifdef __VALIDATE__
   if ( !IsCompatible( tV ) )
       throw XError( "Invalid vector subtraction!" );
#endif

   for ( int i = 0; i <= _upper; i++ ) _tIndex[i] -= tV[i];
   return * this;
   }

template < class T >
TVector < T > TVector < T >::operator * ( const T & t ) const
   {
   TVector < T > tProduct( * this );
   return tProduct *= t;
   }

template < class T >
TVector < T > TVector < T >::operator / ( const T & t ) const
   {
   TVector < T > tProduct( * this );
   return tProduct /= t;
   }

template < class T >
TVector < T > & TVector < T >::operator *= ( const T & t )
   {
   for ( int i = 0; i <= _upper; i++ ) _tIndex[i] *= t;
   return * this;
   }

template < class T >
TVector < T > & TVector < T >::operator /= ( const T & t )
   {
   for ( int i = 0; i <= _upper; i++ ) _tIndex[i] /= t;
   return * this;
   }

template < class T >
TVector < T > TVector < T >::DotMultiply( const TVector < T > & tV ) const
   {
   TVector < T > tProduct( * this );
   for ( int i = 0; i <= _upper; i++ ) tProduct[i] *= tV[i];
   return tProduct;
   }

template < class T >
TVector < T > TVector < T >::DotDivide( const TVector < T > & tV ) const
   {
   TVector < T > tProduct( * this );
   for ( int i = 0; i <= _upper; i++ ) tProduct[i] /= tV[i];
   return tProduct;
   }

template < class T >
T TVector < T >::Dot( const TVector < T > & tV ) const
   {
#ifdef __VALIDATE__
   if ( !IsCompatible( tV ) )
       throw XError( "Invalid vector dot product!" );
#endif

   T tDot = 0;
   for ( int i = 0; i <= _upper; i++ ) tDot += _tIndex[i] * tV[i];
   return tDot;
   }

template < class T >
T TVector < T >::Average( void )const
   {
   T avg = 0;
   for ( int i = 0; i <= _upper; i++ )
      { avg += _tIndex[i]; }

   return avg / Length();
   }

template < class T >
T TVector < T >::Correl( const TVector < T > & tV ) const
   {

   T cov = 0;
   T var1 = 0;
   T var2 = 0;

   T avg1 = Average();
   T avg2 = tV.Average();


   for ( int i = 0; i <= _upper; i++ )
      {
      cov += ( _tIndex[i] - avg1 ) * ( tV[i] - avg2 );
      var1 += ( _tIndex[i] - avg1 ) * ( _tIndex[i] - avg1 );
      var2 += ( tV[i] - avg2 ) * ( tV[i] - avg2 );
      }

   return cov / sqrt( var1 * var2 );
   }

template < class T >
T TVector < T >::FindMin( void ) const
   {
   T tMin = _tIndex[0];
   for ( int i = 1; i <= _upper; i++ )
      {
      tMin = _tIndex[i] < tMin ? _tIndex[i] : tMin;
      }
   return tMin;
   }

template < class T >
int TVector < T >::IndexOfMax( void )const
   {
   int max = 0;
   for ( int i = 1; i <= _upper; i++ )
      {
      max = _tIndex[i] > _tIndex[max] ? i : max;
      }
   return max;
   }

template < class T >
T TVector < T >::FindMax( void )const
   {
   T tMax = 0;
   for ( int i = 0; i <= _upper; i++ )
      {
      tMax = _tIndex[i] > tMax ? _tIndex[i] : tMax;
      }
   return tMax;
   }

template < class T >
T TVector < T >::AbsMax( void )const
   {
   T tMax = 0;
   for ( int i = 0; i <= _upper; i++ )
      {
      T fabsi = fabs( _tIndex[i] );
      tMax = fabsi > tMax ? fabsi : tMax;
      }
   return tMax;
   }

template < class T >
T TVector < T >::CalcNorm( int order ) const
{
    T tNorm = 0;
    if ( !order ) //zero is surrogate for infinity norm
        return AbsMax();
    
    for ( int i = 0; i <= _upper; i++ ) tNorm += pow( fabs( _tIndex[i] ), order );
    return pow( tNorm, 1.0 / ( double )order );
}

//"hard-code" L2 norm, since it's so common
template < class T > T TVector < T >::CalcL2Norm( void ) const
{
    T tNorm = 0;
    for ( int i = 0; i <= _upper; i++ ) tNorm += _tIndex[i] * _tIndex[i];
    return sqrt(tNorm);
}

#endif
